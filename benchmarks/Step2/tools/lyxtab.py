"""Cell-level edits of native LyX tables (refuse-or-apply).
   cells(s, label) -> list of rows, each a list of dicts {start,end,text} where start:end spans the cell's Plain Layout body.
   retext(s, label, {(r,c): 'mini markup'}) -> new s.  Mini markup: \\, thin space, ~ nonbreaking space, $...$ formula, everything else literal."""
import re

THIN = "\\begin_inset space \\thinspace{}\n\\end_inset\n\n"
NBSP = "\\begin_inset space ~\n\\end_inset\n\n"

def render(text):
    """mini markup -> LyX lines for the body of a Plain Layout (without the layout begin/end lines)."""
    out = []; buf = ""
    i = 0
    def flush():
        nonlocal buf
        if buf != "":
            out.append(buf + "\n"); buf = ""
    while i < len(text):
        if text.startswith("\\,", i):
            flush(); out.append(THIN); i += 2
        elif text[i] == "~":
            flush(); out.append(NBSP); i += 1
        elif text[i] == "$":
            j = text.index("$", i + 1)
            flush(); out.append("\\begin_inset Formula $%s$\n\\end_inset\n\n" % text[i + 1:j]); i = j + 1
        else:
            buf += text[i]; i += 1
    flush()
    return "".join(out)

def table_span(s, label):
    k = s.index('name "%s"\n' % label)
    a = s.index("<lyxtabular", k); b = s.index("</lyxtabular>", a) + len("</lyxtabular>")
    return a, b

def cells(s, label):
    a, b = table_span(s, label)
    t = s[a:b]
    rows = []
    for rm in re.finditer(r"<row[^>]*>\n(.*?)</row>", t, re.S):
        row = []
        for cm in re.finditer(r"<cell[^>]*>\n\\begin_inset Text\n\n\\begin_layout Plain Layout\n(.*?)\\end_layout\n\n\\end_inset\n</cell>", rm.group(1), re.S):
            gs = a + rm.start(1) + cm.start(1); ge = a + rm.start(1) + cm.end(1)
            row.append({"start": gs, "end": ge, "text": cm.group(1)})
        rows.append(row)
    return rows

def retext(s, label, edits):
    rows = cells(s, label)
    for (r, c) in edits:
        assert r < len(rows) and c < len(rows[r]), ("no such cell", label, r, c)
    for (r, c) in sorted(edits, key=lambda rc: -rows[rc[0]][rc[1]]["start"]):
        cl = rows[r][c]
        m = re.match(r"(\n?\\size \w+\n)", cl["text"])
        s = s[:cl["start"]] + (m.group(1) if m else "") + render(edits[(r, c)]) + s[cl["end"]:]
    return s

def plain(cell_text):
    """readable form of a cell body: formulas as $..$, thin space as \\,"""
    t = re.sub(r"^\n?\\size \w+\n", "", cell_text).replace(THIN, "\\,").replace(NBSP, "~")
    t = re.sub(r"\\begin_inset Formula (\$.*?\$)\n\\end_inset\n\n", r"\1", t)
    return t.replace("\n", "")

def dump(s, label):
    for r, row in enumerate(cells(s, label)):
        print(r, [plain(c["text"]).strip() for c in row])

if __name__ == "__main__":
    import sys
    dump(open(sys.argv[1], encoding="utf-8", newline="").read(), sys.argv[2])


def _cell(text, align, attrs=""):
    return ('<cell alignment="%s" valignment="top"%s usebox="none">\n\\begin_inset Text\n\n\\begin_layout Plain Layout\n%s\\end_layout\n\n\\end_inset\n</cell>\n'
            % (align, attrs, render(text)))

def float_table(caption, label, header, rows, aligns, placement="document"):
    """Complete Float table markup in the house form: centering ERT, caption (mini markup), label, booktabs tabular.
       header: list of cell strings; rows: list of lists; aligns: list of 'left'/'right'/'center' per column."""
    ncol = len(header); nrow = len(rows) + 1
    t = '\\begin_inset Float table\nplacement %s\nalignment document\nwide false\nsideways false\nstatus open\n\n' % placement
    t += '\\begin_layout Plain Layout\n\\begin_inset ERT\nstatus collapsed\n\n\\begin_layout Plain Layout\n\n\\backslash\ncentering\n\\end_layout\n\n\\end_inset\n\n \n'
    t += '\\begin_inset Caption Standard\n\n\\begin_layout Plain Layout\n' + render(caption) + '\\end_layout\n\n\\end_inset\n\n\n'
    t += '\\begin_inset CommandInset label\nLatexCommand label\nname "%s"\n\n\\end_inset\n\n \n' % label
    t += '\\begin_inset Tabular\n<lyxtabular version="3" rows="%d" columns="%d">\n<features booktabs="true" tabularvalignment="middle">\n' % (nrow, ncol)
    for a in aligns: t += '<column alignment="%s" valignment="top">\n' % a
    for ri, row in enumerate([header] + rows):
        t += "<row>\n"
        attrs = ' topline="true"' if ri <= 1 else (' bottomline="true"' if ri == nrow - 1 else "")
        for ci, c in enumerate(row):
            t += _cell(c, aligns[ci], attrs)
        t += "</row>\n"
    t += '</lyxtabular>\n\n\\end_inset\n\n\n\\end_layout\n\n\\end_inset\n\n\n'
    return t
