#!/usr/bin/env python3
"""Replace the mesh shadings of a one-page vector figure by a raster image of the painted field.

    python3 shading_to_image.py in.pdf out.pdf [--dpi 600] [--jpeg Q | --flate]

The headless IST Vector Postprocessor paints a coloured surface as one free-form Gouraud triangle mesh
(PDF ShadingType 4) with a vertex per mesh node (49 554 for the upper face of the penny crack, figs/pdf/penny_opening.pdf);
PDF viewers rasterise such a mesh at every redraw, which takes seconds. This tool renders the shadings alone
(pdftocairo, white background) at <dpi> dots per inch of the figure's own coordinates, crops the raster to the painted
area and draws it as an image XObject just before the form that held the shadings; the colour bar, the triad and the
text stay vector. The sh operators and the /Shading resources are removed. The image is Flate-coded with the PNG Paeth
predictor (lossless) unless --jpeg Q is given; the default picks JPEG only when the lossless stream exceeds 1.5 MB.
"""
import argparse
import io
import os
import re
import subprocess
import sys
import tempfile
import zlib

import numpy as np
from PIL import Image
from pypdf import PdfReader, PdfWriter
from pypdf.generic import (ArrayObject, BooleanObject, DictionaryObject, FloatObject, NameObject, NumberObject,
                           StreamObject)

PAINT = {b"f", b"F", b"f*", b"B", b"B*", b"b", b"b*", b"S", b"s", b"n"}
PATH = {b"m", b"l", b"c", b"v", b"y", b"h", b"re"}
KEEP = {b"q", b"Q", b"cm", b"gs", b"sh"}
TOKEN = re.compile(rb"/[^\s/\[\]<>(){}]*|\[|\]|<<|>>|\([^)]*\)|<[0-9A-Fa-f\s]*>|[^\s\[\]<>(){}/]+")


def shading_only(content):
    """Keep q/Q/cm/gs/sh and clipping paths; drop every painting operator and all text."""
    out, path, clip, intext = [], [], False, False
    for tok in TOKEN.findall(content):
        if tok == b"BT":
            intext = True
            continue
        if tok == b"ET":
            intext = False
            continue
        if intext:
            continue
        if tok in PATH:
            path.append(tok)
            continue
        if tok in (b"W", b"W*"):
            clip = tok
            continue
        if tok in PAINT:
            if clip:
                out.extend(path + [clip, b"n"])
            path, clip = [], False
            continue
        if tok in KEEP:
            out.extend(path + [tok])
            path = []
            continue
        if tok[:1] in b"/-.0123456789[]<(":
            path.append(tok)
            continue
        path = []  # other operators (colour, text state, Do): dropped with their operands
    return b" ".join(out) + b"\n"


def strip_shading(content):
    return re.sub(rb"/[^\s]+\s+sh\b", b"", content)


def forms_with_shading(page):
    res = page.get("/Resources") or {}
    out = []
    for name, ref in (res.get("/XObject") or {}).items():
        xo = ref.get_object()
        if xo.get("/Subtype") == "/Form" and (xo.get("/Resources") or {}).get("/Shading"):
            out.append((name, xo))
    return out


def new_stream(template, data, compress=True):
    s = StreamObject()
    for k, v in template.items():
        if k not in ("/Length", "/Filter", "/DecodeParms"):
            s[NameObject(k)] = v
    s.set_data(data)
    return s.flate_encode() if compress else s


def paeth_flate(rgb):
    """PNG Paeth-filtered, zlib-compressed scanlines of an (H, W, 3) uint8 array."""
    h, w, _ = rgb.shape
    x = rgb.astype(np.int16)
    a = np.zeros_like(x); a[:, 1:] = x[:, :-1]
    b = np.zeros_like(x); b[1:] = x[:-1]
    c = np.zeros_like(x); c[1:, 1:] = x[:-1, :-1]
    p = a + b - c
    pa, pb, pc = np.abs(p - a), np.abs(p - b), np.abs(p - c)
    pred = np.where((pa <= pb) & (pa <= pc), a, np.where(pb <= pc, b, c))
    filt = ((x - pred) & 0xFF).astype(np.uint8).reshape(h, w * 3)
    rows = np.concatenate([np.full((h, 1), 4, np.uint8), filt], axis=1)
    return zlib.compress(rows.tobytes(), 9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--jpeg", type=int, default=None, help="JPEG quality (lossy); default lossless Flate")
    ap.add_argument("--flate", action="store_true", help="force lossless Flate whatever its size")
    ap.add_argument("--max-flate", type=float, default=1.5, help="MB above which JPEG 92 is used instead")
    a = ap.parse_args()

    reader = PdfReader(a.src)
    if len(reader.pages) != 1:
        sys.exit("one-page figures only")
    page = reader.pages[0]
    forms = forms_with_shading(page)
    if not forms:
        sys.exit("no form XObject with a /Shading resource in %s" % a.src)
    mb = page.mediabox
    pw, ph = float(mb.width), float(mb.height)

    # 1. the shadings alone, rendered on white
    w0 = PdfWriter(clone_from=reader)
    p0 = w0.pages[0]
    xo = p0["/Resources"]["/XObject"]
    nverts = 0
    for name, xob in forms:
        src_form = xo[name].get_object()
        for sname, sref in src_form["/Resources"]["/Shading"].items():
            sh = sref.get_object()
            if "/BitsPerFlag" in sh:
                bpv = (int(sh["/BitsPerFlag"]) + 2 * int(sh["/BitsPerCoordinate"]) + 3 * int(sh["/BitsPerComponent"])) // 8
                nverts += len(sh.get_data()) // bpv
        xo[NameObject(name)] = w0._add_object(new_stream(src_form, shading_only(src_form.get_data())))
    with tempfile.TemporaryDirectory() as td:
        tmp = os.path.join(td, "sh.pdf")
        with open(tmp, "wb") as f:
            w0.write(f)
        subprocess.run(["pdftocairo", "-png", "-r", str(a.dpi), "-singlefile", tmp, os.path.join(td, "sh")], check=True)
        im = np.asarray(Image.open(os.path.join(td, "sh.png")).convert("RGB"))
    ink = np.any(im < 250, axis=2)
    ys, xs = np.where(ink)
    if len(xs) == 0:
        sys.exit("the shadings painted nothing")
    m = 2
    x0, x1 = max(xs.min() - m, 0), min(xs.max() + m + 1, im.shape[1])
    y0, y1 = max(ys.min() - m, 0), min(ys.max() + m + 1, im.shape[0])
    crop = np.ascontiguousarray(im[y0:y1, x0:x1])
    H, W = crop.shape[:2]
    s = 72.0 / a.dpi
    bx, by = x0 * s, ph - y1 * s           # lower-left corner in page points
    bw, bh = W * s, H * s
    print("shading vertices %d; raster %dx%d px at %d dpi; placed at (%.2f, %.2f) size %.2f x %.2f pt"
          % (nverts, W, H, a.dpi, bx, by, bw, bh))

    # 2. encode the image
    flate = paeth_flate(crop)
    use_jpeg = a.jpeg is not None or (not a.flate and len(flate) > a.max_flate * 1e6)
    img = StreamObject()
    img.update({NameObject("/Type"): NameObject("/XObject"), NameObject("/Subtype"): NameObject("/Image"),
                NameObject("/Width"): NumberObject(W), NameObject("/Height"): NumberObject(H),
                NameObject("/ColorSpace"): NameObject("/DeviceRGB"), NameObject("/BitsPerComponent"): NumberObject(8),
                NameObject("/Interpolate"): BooleanObject(True)})
    if use_jpeg:
        q = a.jpeg or 92
        buf = io.BytesIO()
        Image.fromarray(crop).save(buf, format="JPEG", quality=q, subsampling=0, optimize=True)
        img._data = buf.getvalue()
        img[NameObject("/Filter")] = NameObject("/DCTDecode")
        print("image: JPEG quality %d, %.0f kB (lossless Flate would be %.0f kB)" % (q, len(img._data) / 1e3, len(flate) / 1e3))
    else:
        img._data = flate
        img[NameObject("/Filter")] = NameObject("/FlateDecode")
        img[NameObject("/DecodeParms")] = DictionaryObject({
            NameObject("/Predictor"): NumberObject(15), NameObject("/Colors"): NumberObject(3),
            NameObject("/BitsPerComponent"): NumberObject(8), NameObject("/Columns"): NumberObject(W)})
        print("image: lossless Flate (Paeth), %.0f kB" % (len(flate) / 1e3))

    # 3. the figure without the shadings, the image drawn first
    w = PdfWriter(clone_from=reader)
    p = w.pages[0]
    xo = p["/Resources"]["/XObject"]
    for name, _ in forms:
        form = xo[name].get_object()
        data = strip_shading(form.get_data())
        form_res = DictionaryObject({NameObject(k): v for k, v in form["/Resources"].items() if k != "/Shading"})
        ns = new_stream(form, data)
        ns[NameObject("/Resources")] = form_res
        xo[NameObject(name)] = w._add_object(ns)
    iname = "/ImField"
    xo[NameObject(iname)] = w._add_object(img)
    draw = ("q %.4f 0 0 %.4f %.4f %.4f cm %s Do Q\n" % (bw, bh, bx, by, iname)).encode()
    old = p.get_contents().get_data()
    p.replace_contents(new_stream({}, draw + old))
    w.add_metadata({"/Producer": "shading_to_image.py (pypdf)"})
    with open(a.dst, "wb") as f:
        w.write(f)
    print("wrote %s (%.0f kB)" % (a.dst, os.path.getsize(a.dst) / 1e3))


if __name__ == "__main__":
    main()
