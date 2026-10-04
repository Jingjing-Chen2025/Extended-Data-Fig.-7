# EXTENDED DATA FIG. 7 - merge the spatial-map panel onto the text layer.
#
# Clipping is lossless: the source page's content stream is wrapped in a clip path, so nothing is rasterised
# and the result stays fully editable vector art.
#
# ED7_PANEL names the page the panel is clipped out of, and the clip box below is for the laid-out artwork used
# in the manuscript (a 612 x 792 pt page carrying the panel in its upper half). Running R/04_abundance_proximity.R
# from scratch writes `04_spatial_maps.pdf` at 12 x 8 in instead, so the box has to be re-measured for that page.
#
# Usage: ED7_PANEL="<path to the page holding the spatial maps>" python python/compose.py
import os
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject

REPO  = os.environ.get("ED7_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT   = os.environ.get("ED7_OUT", os.path.join(REPO, "figures"))
BASE  = os.path.join(OUT, "ED7 base.pdf")                      # from R/06_assemble_base.R
PANEL = os.environ.get("ED7_PANEL", os.path.join(OUT, "04_spatial_maps.pdf"))
FINAL = os.path.join(OUT, "Extended Data Fig. 7.pdf")

BOX  = tuple(float(v) for v in os.environ.get("ED7_BOX", "100,478,550,736").split(","))
SLOT = (30, 40, 582, 378)

def clipped(path, box):
    w = PdfWriter(clone_from=path); pg = w.pages[0]
    x0, y0, x1, y1 = box
    s = DecodedStreamObject()
    s.set_data(b"q %f %f %f %f re W n\n" % (x0, y0, x1 - x0, y1 - y0) + pg.get_contents().get_data() + b"\nQ")
    pg.replace_contents(s)
    return pg

def place(target, pg, box, slot):
    x0, y0, x1, y1 = box; a, b, c, d = slot
    sw, sh = c - a, d - b
    s = min(sw / (x1 - x0), sh / (y1 - y0))
    target.merge_transformed_page(pg, (s, 0, 0, s,
                                       a + (sw - s * (x1 - x0)) / 2 - s * x0,
                                       b + (sh - s * (y1 - y0)) / 2 - s * y0))

w = PdfWriter(); w.append(PdfReader(BASE)); t = w.pages[0]
place(t, clipped(PANEL, BOX), BOX, SLOT)
with open(FINAL, "wb") as fh: w.write(fh)
print("written", FINAL)
