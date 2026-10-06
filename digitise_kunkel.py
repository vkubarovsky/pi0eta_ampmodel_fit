"""Digitise CLAS g12 pi0 photoproduction, Kunkel PRC 98 015207 figure 4.

Red filled circles = the g12 measurement.  Error bars are thin vertical red
lines; the markers are blobs, so an erosion separates them.  Axis calibration
is read from the tick marks, not guessed: decades at y = 333/423/514/604 in the
top row (90.33 px per decade) and |t| = 0 at the left frame with 75.5 px/GeV^2.
"""
import numpy as np
from PIL import Image
from scipy import ndimage

a = np.asarray(Image.open("pg6-6.png").convert("RGB")).astype(int)
R, G, B = a[:,:,0], a[:,:,1], a[:,:,2]
red = (R > 120) & (R - G > 60) & (R - B > 60)

COLS = [(898,1502),(1502,2111),(2111,2716)]
ROWS = [(289,649),(649,1011)]
# (panel label, W in GeV) in the paper's order (a)..(f)
PANELS = [("a",2.490),("b",2.635),("c",2.790),("d",2.940),("e",3.080),("f",3.170)]
YREF = {0:(333,3.0), 1:(785,2.0)}      # row -> (pixel, log10 sigma)
PXDEC = 90.33

out = {}
for k,(lab,W) in enumerate(PANELS):
    r, c = k//3, k%3
    x0,x1 = COLS[c]; y0,y1 = ROWS[r]
    m = np.zeros_like(red); m[y0+2:y1-2, x0+2:x1-2] = red[y0+2:y1-2, x0+2:x1-2]
    core = ndimage.binary_erosion(m, np.ones((5,5)))
    lbl, n = ndimage.label(core)
    yref, lref = YREF[r]; sx = (x1-x0)/8.0
    pts = []
    for i in range(1, n+1):
        ys, xs = np.where(lbl == i)
        if len(ys) < 8: continue
        t = (xs.mean() - x0)/sx
        s = 10**(lref + (yref - ys.mean())/PXDEC)
        pts.append((t, s))
    pts.sort()
    out[lab] = (W, pts)
    print(f"panel ({lab})  W = {W:.3f} GeV   {len(pts)} markers")

with open("kunkel_g12_pi0_photoproduction.txt","w") as f:
    f.write("# CLAS g12 gamma p -> pi0 p, Kunkel et al. PRC 98 015207 (2018), figure 4\n")
    f.write("# DIGITISED from the published figure at 400 dpi: red filled circles,\n")
    f.write("# axis calibration from the tick marks (90.33 px/decade, 75.5 px/GeV^2).\n")
    f.write("# Accuracy of the digitisation is a few per cent on a log scale; the\n")
    f.write("# published uncertainties (statistical only in that figure) are NOT captured.\n")
    f.write("#  W [GeV]   |t| [GeV^2]   dsigma/dt [nb/GeV^2]\n")
    for lab,(W,pts) in out.items():
        for t,s in pts:
            f.write(f"  {W:8.3f} {t:12.4f} {s:16.4g}\n")
print("\n-> kunkel_g12_pi0_photoproduction.txt")
for lab,(W,pts) in out.items():
    sel=[(t,s) for t,s in pts if t<1.6]
    if sel: print(f"  ({lab}) W={W:.2f}: " + "  ".join(f"|t|={t:.2f}->{s:.0f}" for t,s in sel[:6]))
