#!/usr/bin/env python3
"""Colour sequence of the growing corner crack (Step2.Bench corner), from the light EnSight Gold case.

    python3 fig_corner_color.py <ens dir> [steps=0,2,4,6,8,11] [out dir] [scale=0] [view=x,y,z]

<ens dir>: corner_fronts.csv (fronts per step) and crack/crack.geoNNNN, crack.openNNNN, crack.dispNNNN (ASCII Gold,
parts crack_face (tria3), front (bar2)). One PNG per step in a fixed view, the upper face of the crack coloured by the
normal opening divided by 8 (1 - nu^2) sigma a / (pi E), a the mean radius of that step's front (as in the opening
figures of the manuscript), the earlier fronts in black, the block edges in gray; plus the contact sheet
corner_sequence.png (3 x 2). A painter's algorithm over the triangles (no GL, no matplotlib). The colour is clipped at the prescribed (smooth) front of the step
and the pixels behind it that the cut mesh leaves uncovered (the front is not fitted to the mesh) take the nearest colour: the region
behind the front is the union of the quadrilaterals between successive fronts and of the fan from the corner to the first one.
"""
import csv
import math
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

NU = 0.3


def turbo(t):
    """Google's turbo colormap (polynomial approximation), t in [0, 1] -> (n, 3) uint8."""
    t = np.clip(np.asarray(t, float), 0, 1)
    r = 0.13572138 + t * (4.61539260 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))))
    g = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))))
    b = 0.10667330 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))))
    return (np.clip(np.stack([r, g, b], axis=-1), 0, 1) * 255).astype(np.uint8)


def read_geo(path):
    """Parts of an ASCII Gold geometry: {name: (xyz (n, 3), etype, conn (m, k) 0-based)}."""
    with open(path) as f:
        lines = f.read().split("\n")
    i = 4  # description (2), node id, element id
    parts = {}
    while i < len(lines) and lines[i].strip() == "part":
        pid = int(lines[i + 1])
        name = lines[i + 2].strip()
        assert lines[i + 3].strip() == "coordinates"
        n = int(lines[i + 4])
        vals = np.array(lines[i + 5:i + 5 + 3 * n], dtype=float)
        xyz = vals.reshape(3, n).T
        i += 5 + 3 * n
        etype = lines[i].strip()
        m = int(lines[i + 1])
        k = {"tria3": 3, "bar2": 2, "tetra4": 4}[etype]
        conn = np.array([lines[i + 2 + q].split() for q in range(m)], dtype=np.int64) - 1 if m else np.zeros((0, k), np.int64)
        i += 2 + m
        parts[name] = (pid, xyz, etype, conn)
    return parts


def read_var(path, sizes, ncomp):
    """Node variable per part id: {pid: (n, ncomp)}."""
    with open(path) as f:
        lines = f.read().split("\n")
    i = 1
    out = {}
    while i < len(lines) and lines[i].strip() == "part":
        pid = int(lines[i + 1])
        assert lines[i + 2].strip() == "coordinates"
        n = sizes[pid]
        vals = np.array(lines[i + 3:i + 3 + ncomp * n], dtype=float)
        out[pid] = vals.reshape(ncomp, n).T
        i += 3 + ncomp * n
    return out


def view_basis(eye):
    w = np.array(eye, float)
    w /= np.linalg.norm(w)
    up = np.array([0.0, 0.0, 1.0])
    r = np.cross(up, w)
    r /= np.linalg.norm(r)
    v = np.cross(w, r)
    return r, v, w


def font(size):
    for p in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def main():
    ens = sys.argv[1]
    steps = [int(s) for s in (sys.argv[2] if len(sys.argv) > 2 else "0,2,4,6,8,11").split(",")]
    out = sys.argv[3] if len(sys.argv) > 3 else ens
    scale = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    eye = [float(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [0.55, 0.62, 0.55]
    fit = sys.argv[6] if len(sys.argv) > 6 else "all"
    os.makedirs(out, exist_ok=True)
    r, v, w = view_basis(eye)

    fronts = {}
    fcsv = os.path.join(ens, "corner_fronts.csv")
    if not os.path.exists(fcsv):
        fcsv = os.path.join(os.path.dirname(os.path.abspath(ens)), "corner_fronts.csv")
    for row in csv.DictReader(open(fcsv)):
        fronts.setdefault(int(row["step"]), []).append((float(row["x"]), float(row["y"]), float(row["z"])))
    fronts = {k: np.array(val) for k, val in fronts.items()}

    def region_polys(upto):
        """Triangles (3D) of the region behind front `upto`: fan from the corner to front 0 and strips between successive fronts."""
        tris = []
        f0 = fronts[0]
        c0 = np.array([0.0, 0.0, 0.0])
        for j in range(len(f0) - 1):
            tris.append((c0, f0[j], f0[j + 1]))
        for s2 in range(upto):
            a, b = fronts[s2], fronts[s2 + 1]
            m = min(len(a), len(b))
            for j in range(m - 1):
                tris.append((a[j], a[j + 1], b[j]))
                tris.append((a[j + 1], b[j + 1], b[j]))
        return tris

    frames = []
    for st in steps:
        tag = "%04d" % st
        geo = read_geo(os.path.join(ens, "crack", "crack.geo" + tag))
        pid, xyz, _, conn = geo["crack_face"]
        sizes = {p[0]: len(p[1]) for p in geo.values()}
        op = read_var(os.path.join(ens, "crack", "crack.open" + tag), sizes, 1)[pid][:, 0]
        if scale > 0:
            u = read_var(os.path.join(ens, "crack", "crack.disp" + tag), sizes, 3)[pid]
            xyz = xyz + scale * u
        onface = (np.all(np.abs(xyz[conn][:, :, 0]) < 1e-9, axis=1) | np.all(np.abs(xyz[conn][:, :, 1]) < 1e-9, axis=1) |
                  np.all(np.abs(xyz[conn][:, :, 0] - 1) < 1e-9, axis=1) | np.all(np.abs(xyz[conn][:, :, 1] - 1) < 1e-9, axis=1))
        conn = conn[~onface]
        f = fronts[st]
        seg = np.linalg.norm(np.diff(f, axis=0), axis=1).sum()
        a = seg / (0.5 * math.pi)
        ref = 8 * (1 - NU ** 2) * a / math.pi  # sigma = E = 1
        frames.append(dict(step=st, xyz=xyz, conn=conn, val=np.maximum(op, 0) / ref, a=a, front=f))
    allv = np.concatenate([fr["val"] for fr in frames])
    vmax = float(np.percentile(allv, 99.5))
    print("opening / reference: 99.5 percentile %.3f, max %.3f" % (vmax, allv.max()))

    # fixed screen window from all frames, the fronts and the corner region of the block
    pts = np.concatenate([fr["xyz"] for fr in frames] + [fr["front"] for fr in frames])
    sx, sy = pts @ r, pts @ v
    pad = 0.04 * max(np.ptp(sx), np.ptp(sy))
    x0, x1, y0, y1 = sx.min() - pad, sx.max() + pad, sy.min() - pad, sy.max() + pad
    W = 1100
    H = int(W * (y1 - y0) / (x1 - x0))
    ss = 2
    fnt, fnt_s = font(34), font(26)
    zlo, zhi = float(pts[:, 2].min()) - 0.12, float(pts[:, 2].max()) + 0.12
    box = [(x, y, z) for z in (zlo, zhi) for y in (0, 1) for x in (0, 1)]
    edges = [(0, 1), (2, 3), (0, 2), (1, 3), (4, 5), (6, 7), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7)]

    win = [x0, x1, y0, y1]

    def to_px(p):
        p = np.atleast_2d(p)
        return np.stack([(p @ r - win[0]) / (win[1] - win[0]) * W * ss, (win[3] - p @ v) / (win[3] - win[2]) * H * ss], axis=1)

    cb_w = 70
    images = []
    for fr in frames:
        if fit == "frame":
            q = np.concatenate([fr["xyz"], fr["front"]])
            qx, qy = q @ r, q @ v
            cx, cy = 0.5 * (qx.min() + qx.max()), 0.5 * (qy.min() + qy.max())
            hx, hy = 0.5 * np.ptp(qx) * 1.18, 0.5 * np.ptp(qy) * 1.18
            asp = (x1 - x0) / (y1 - y0)
            if hx / hy < asp:
                hx = hy * asp
            else:
                hy = hx / asp
            win[:] = [cx - hx, cx + hx, cy - hy, cy + hy]
        img = Image.new("RGB", (W * ss, H * ss), "white")
        d = ImageDraw.Draw(img)
        bp = to_px(np.array(box))
        for i, j in edges:
            d.line([tuple(bp[i]), tuple(bp[j])], fill=(190, 190, 190), width=2 * ss)
        xyz, conn, val = fr["xyz"], fr["conn"], fr["val"]
        depth = (xyz @ w)[conn].mean(axis=1)
        order = np.argsort(depth)
        px = to_px(xyz)
        col = turbo(np.clip(val / vmax, 0, 1))
        tcol = col[conn].mean(axis=1).astype(np.uint8)
        cov = Image.new("L", img.size, 0)
        dc = ImageDraw.Draw(cov)
        for t in order:
            c = conn[t]
            poly = [tuple(px[c[0]]), tuple(px[c[1]]), tuple(px[c[2]])]
            d.polygon(poly, fill=tuple(int(q) for q in tcol[t]))
            dc.polygon(poly, fill=255)
        # region behind the smooth front: clip the colour to it and fill what the mesh leaves uncovered
        mask = Image.new("L", img.size, 0)
        dm = ImageDraw.Draw(mask)
        for tri in region_polys(fr["step"]):
            q = to_px(np.array(tri))
            dm.polygon([tuple(q[0]), tuple(q[1]), tuple(q[2])], fill=255)
        A = np.array(img).astype(np.float32)
        have = np.array(cov) > 0
        inside = np.array(mask) > 0
        have &= inside
        for it in range(40 * ss):
            miss = inside & ~have
            if not miss.any():
                break
            acc = np.zeros_like(A)
            cnt = np.zeros(have.shape, np.float32)
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    h2 = np.roll(np.roll(have, dy, 0), dx, 1)
                    c2 = np.roll(np.roll(A, dy, 0), dx, 1)
                    acc += c2 * h2[..., None]
                    cnt += h2
            fill = miss & (cnt > 0)
            A[fill] = acc[fill] / cnt[fill][:, None]
            have |= fill
        base = Image.new("RGB", img.size, "white")
        db = ImageDraw.Draw(base)
        for i, j in edges:
            db.line([tuple(bp[i]), tuple(bp[j])], fill=(190, 190, 190), width=2 * ss)
        outimg = np.where((inside)[..., None], A, np.array(base)).astype(np.uint8)
        img = Image.fromarray(outimg)
        d = ImageDraw.Draw(img)
        for st2, f in sorted(fronts.items()):
            if st2 <= fr["step"]:
                p = to_px(f)
                d.line([tuple(q) for q in p], fill=(0, 0, 0) if st2 < fr["step"] else (255, 255, 255), width=(2 if st2 < fr["step"] else 3) * ss)
                if st2 == fr["step"]:
                    d.line([tuple(q) for q in p], fill=(0, 0, 0), width=2 * ss)
        img = img.resize((W, H), Image.LANCZOS)
        d = ImageDraw.Draw(img)
        d.text((16, 10), "step %d    mean front radius a = %.3f" % (fr["step"], fr["a"]), fill=(20, 20, 20), font=fnt)
        # axes triad
        ox, oy = 70, H - 70
        for k, (name, e) in enumerate((("x", (1, 0, 0)), ("y", (0, 1, 0)), ("z", (0, 0, 1)))):
            dx, dy = float(np.dot(e, r)), float(-np.dot(e, v))
            n = math.hypot(dx, dy) or 1
            d.line([(ox, oy), (ox + 48 * dx / n * min(1, n * 1.6), oy + 48 * dy / n * min(1, n * 1.6))], fill=(60, 60, 60), width=3)
            d.text((ox + 56 * dx / n * min(1, n * 1.6) - 6, oy + 56 * dy / n * min(1, n * 1.6) - 12), name, fill=(60, 60, 60), font=fnt_s)
        # colour bar
        bx, by0, by1 = W - 120, 90, H - 90
        for yy in range(by0, by1):
            t = 1 - (yy - by0) / (by1 - by0)
            d.line([(bx, yy), (bx + cb_w // 2, yy)], fill=tuple(int(q) for q in turbo([t])[0]))
        d.rectangle([bx, by0, bx + cb_w // 2, by1], outline=(80, 80, 80))
        for k in range(5):
            t = k / 4
            yy = by1 - t * (by1 - by0)
            d.line([(bx + cb_w // 2, yy), (bx + cb_w // 2 + 8, yy)], fill=(80, 80, 80))
            d.text((bx + cb_w // 2 + 12, yy - 13), "%.2f" % (t * vmax), fill=(40, 40, 40), font=fnt_s)
        d.text((bx - 90, by0 - 44), "opening / ref.", fill=(40, 40, 40), font=fnt_s)
        path = os.path.join(out, "corner_seq_%02d%s.png" % (fr["step"], "_zoom" if fit == "frame" else ""))
        img.save(path)
        images.append(img)
        print("  wrote", path)
    cols = 3
    rows = (len(images) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * W, rows * H), "white")
    for k, im in enumerate(images):
        sheet.paste(im, ((k % cols) * W, (k // cols) * H))
    sheet.save(os.path.join(out, "corner_sequence%s.png" % ("_zoom" if fit == "frame" else "")))
    print("  wrote", os.path.join(out, "corner_sequence%s.png" % ("_zoom" if fit == "frame" else "")))


if __name__ == "__main__":
    main()
