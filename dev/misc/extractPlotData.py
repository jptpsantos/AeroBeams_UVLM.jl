"""
Extract data from a plot image (raster, e.g. a figure of a paper) into CSV.

Works on line plots with markers and/or plain curves. Each series is identified by its color
(optionally restricted to a rectangle in data coordinates, which allows telling apart series that
share a color, e.g. solid vs dashed lines of the same color).

Steps:
  1. The plot box (axes frame) is detected automatically as the outermost long horizontal/vertical lines
     of the image (axes lines or gridlines), and mapped to the axes limits given by --xlim/--ylim (the values at the box edges).
     Alternatively, give two reference pixels per axis with --xcal/--ycal.
  2. The series colors are detected automatically (or given with --series).
  3. Markers are located by a "cross test": a pixel is a marker center if the series color extends at
     least --marker-size pixels both vertically and horizontally through it. This rejects the lines
     (thin and mostly shallow or steep) and accepts squares, circles, diamonds, '+' and so on.
     With --curve, the curve itself is also extracted (mean y of the series pixels at each pixel column).
  4. A check image with the detected points drawn over the original is written, to be inspected.

Usage (run from anywhere; the image is only read, never modified):
  python3 -I extractPlotData.py figure.png --xlim 0 10 --ylim 500 4000
  python3 -I extractPlotData.py figure.png --xlim 0 10 --ylim 500 4000 \
      --series "S10 TE=#0072BD:,,2500,4000" "S10 LE=#0072BD:,,,2500" "S20 TE=#D95319:,,2500,4000"

The --series syntax is NAME=COLOR[:XMIN,XMAX,YMIN,YMAX] with limits in data units (blank = open).
Run first without --series to see the detected colors and marker positions, then split by regions.

Options worth knowing:
  --marker-size N   side length in pixels (default 9) of the cross test; must be larger than the line
                    thickness (otherwise lines are detected) and smaller than the marker
  --tol T           RGB distance for a pixel to belong to a series color (default 60)
  --curve           also extract the curves (written to <out>_curve.csv)
  --box L R T B     plot box in pixels, instead of automatic detection
  --xcal PX1 V1 PX2 V2 / --ycal PY1 V1 PY2 V2   pixel calibration, instead of the box + limits
  --xlog / --ylog   logarithmic axes
  --legend X0 Y0 X1 Y1  legend rectangle in pixels to ignore (default: detected automatically)
  --out PREFIX      output prefix (default: outputs/<image name>)

Limitations: markers hidden below other markers (or cut by them) may be missed or biased (compare with
the check image); points in the legend are ignored only if the legend has a dark frame.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi


def load_rgb(path):
    im = Image.open(path)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im)
    return np.asarray(im.convert("RGB")).astype(float)


def luma(rgb):
    return rgb @ np.array([0.299, 0.587, 0.114])


def longest_runs(mask, axis):
    """Longest run of True along `axis` for each line of the other axis."""
    m = mask if axis == 1 else mask.T
    out = np.zeros(m.shape[0], int)
    for i, row in enumerate(m):
        if not row.any():
            continue
        d = np.diff(np.concatenate(([0], row.astype(int), [0])))
        s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        out[i] = (e - s).max()
    return out


def group_centers(idx):
    """Centers of runs of consecutive indices."""
    if len(idx) == 0:
        return []
    groups = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    return [g.mean() for g in groups]


def detect_box(rgb):
    # Axes lines and gridlines at the limits (possibly light gray) are the longest non-white lines
    ink = luma(rgb) < 245
    h, w = ink.shape
    rowrun = longest_runs(ink, 1)
    colrun = longest_runs(ink, 0)
    rows = group_centers(np.flatnonzero(rowrun > 0.5 * w))
    cols = group_centers(np.flatnonzero(colrun > 0.5 * h))
    if len(rows) < 2 or len(cols) < 2:
        sys.exit("Could not detect the plot box; give it with --box L R T B or use --xcal/--ycal")
    return cols[0], cols[-1], rows[0], rows[-1]


def detect_legend(rgb, box):
    """Pixel rectangles of dark frames inside the plot box, other than the box itself (legends)."""
    L, R, T, B = box
    dark = (luma(rgb) < 110) & (rgb.max(2) - rgb.min(2) < 60)  # black/gray, not a dark colored series
    inner = np.zeros_like(dark)
    m = 4
    inner[int(T) + m : int(B) - m + 1, int(L) + m : int(R) - m + 1] = True
    lab, n = ndi.label(dark & inner, structure=np.ones((3, 3)))
    rects = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        if w > 0.12 * (R - L) and h > 0.08 * (B - T):
            rects.append((sl[1].start, sl[0].start, sl[1].stop, sl[0].stop))
    return rects


def detect_colors(rgb, box, ignore, tol, min_pixels):
    L, R, T, B = box
    inside = np.zeros(rgb.shape[:2], bool)
    m = 4
    inside[int(T) + m : int(B) - m + 1, int(L) + m : int(R) - m + 1] = True
    for x0, y0, x1, y1 in ignore:
        inside[int(y0) : int(y1) + 1, int(x0) : int(x1) + 1] = False
    mx, mn = rgb.max(2), rgb.min(2)
    sat = (mx - mn) / np.maximum(mx, 1)
    cand = inside & (sat > 0.4) & (mx > 80)
    pix = rgb[cand]
    colors = []
    while len(pix) >= min_pixels:
        q = (pix // 16).astype(int)
        key = q[:, 0] * 256 + q[:, 1] * 16 + q[:, 2]
        top = np.bincount(key).argmax()
        c0 = pix[key == top].mean(0)
        near = np.linalg.norm(pix - c0, axis=1) < tol
        if near.sum() < min_pixels:
            pix = pix[~near]
            continue
        c = pix[near].mean(0)
        pix = pix[~near]
        # Skip antialiasing tints, which lie on the segment between an existing color and white
        if not any(dist_to_segment(c, c0, np.full(3, 255.0)) < 25 for c0 in colors):
            colors.append(c)
    return colors


def dist_to_segment(p, a, b):
    t = np.clip(np.dot(p - a, b - a) / np.dot(b - a, b - a), 0, 1)
    return np.linalg.norm(p - (a + t * (b - a)))


def hex_to_rgb(s):
    s = s.lstrip("#")
    return np.array([int(s[i : i + 2], 16) for i in (0, 2, 4)], float)


def rgb_to_hex(c):
    return "#%02X%02X%02X" % tuple(int(round(v)) for v in c)


def parse_series(spec):
    name, rest = spec.split("=", 1)
    region = (None, None, None, None)
    if ":" in rest:
        rest, reg = rest.split(":", 1)
        vals = reg.split(",")
        if len(vals) != 4:
            sys.exit(f"Bad region in --series '{spec}': need XMIN,XMAX,YMIN,YMAX")
        region = tuple(float(v) if v.strip() else None for v in vals)
    return name.strip(), hex_to_rgb(rest.strip()), region


def make_axis(args_cal, lim, box_lo, box_hi, log, which):
    """Returns pixel -> data function for one axis."""
    if args_cal is not None:
        p1, v1, p2, v2 = args_cal
    elif lim is not None:
        (p1, p2), (v1, v2) = (box_lo, box_hi), lim
    else:
        sys.exit(f"Give --{which}lim or --{which}cal")
    if log:
        v1, v2 = np.log10(v1), np.log10(v2)
    def f(p):
        v = v1 + (np.asarray(p, float) - p1) / (p2 - p1) * (v2 - v1)
        return 10**v if log else v
    def inv(v):
        v = np.log10(v) if log else np.asarray(v, float)
        return p1 + (v - v1) / (v2 - v1) * (p2 - p1)
    return f, inv


def cross_centers(mask, size):
    """Centers of marker-like blobs: color extends `size` pixels both vertically and horizontally."""
    # Bridge thin occlusions (e.g. an axis line drawn over a marker)
    mask = ndi.binary_closing(mask, structure=np.ones((3, 3), bool))
    v = ndi.binary_erosion(mask, structure=np.ones((size, 1), bool))
    h = ndi.binary_erosion(mask, structure=np.ones((1, size), bool))
    core = ndi.binary_dilation(v, np.ones((1, 3), bool)) & ndi.binary_dilation(h, np.ones((3, 1), bool))
    lab, n = ndi.label(core, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros((0, 2))
    return np.array(ndi.center_of_mass(core, lab, range(1, n + 1)))[:, ::-1]  # (x, y) in pixels


def curve_columns(mask):
    """Mean row of the mask pixels at each pixel column (NaN where empty)."""
    cnt = mask.sum(0)
    rows = np.arange(mask.shape[0])[:, None]
    with np.errstate(invalid="ignore"):
        y = np.where(cnt > 0, (mask * rows).sum(0) / np.maximum(cnt, 1), np.nan)
    return y


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("image")
    ap.add_argument("--xlim", type=float, nargs=2)
    ap.add_argument("--ylim", type=float, nargs=2)
    ap.add_argument("--xcal", type=float, nargs=4)
    ap.add_argument("--ycal", type=float, nargs=4)
    ap.add_argument("--xlog", action="store_true")
    ap.add_argument("--ylog", action="store_true")
    ap.add_argument("--box", type=float, nargs=4, metavar=("L", "R", "T", "B"))
    ap.add_argument("--legend", type=float, nargs=4, action="append", metavar=("X0", "Y0", "X1", "Y1"))
    ap.add_argument("--series", nargs="+", default=None)
    ap.add_argument("--marker-size", type=int, default=9)
    ap.add_argument("--tol", type=float, default=60)
    ap.add_argument("--min-pixels", type=int, default=150)
    ap.add_argument("--curve", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    img = Path(args.image)
    rgb = load_rgb(img)
    H, W = rgb.shape[:2]

    # Plot box and axes calibration
    box = tuple(args.box) if args.box else detect_box(rgb)
    L, R, T, B = box
    print(f"Image {W}x{H}px; plot box: left={L:.1f} right={R:.1f} top={T:.1f} bottom={B:.1f}")
    px2x, x2px = make_axis(args.xcal, args.xlim, L, R, args.xlog, "x")
    py2y, y2py = make_axis(args.ycal, args.ylim, B, T, args.ylog, "y")

    # Regions to ignore (legends)
    ignore = list(map(tuple, args.legend)) if args.legend else detect_legend(rgb, box)
    for r in ignore:
        print("Ignoring legend/frame region (px): x %d-%d, y %d-%d" % (r[0], r[2], r[1], r[3]))

    # Series
    if args.series:
        series = [parse_series(s) for s in args.series]
    else:
        cols = detect_colors(rgb, box, ignore, args.tol, args.min_pixels)
        if not cols:
            sys.exit("No colored series found (black/gray series are not supported by auto-detection; use --series)")
        series = [(f"color{i+1}_{rgb_to_hex(c)}", c, (None,) * 4) for i, c in enumerate(cols)]
        print("Detected colors:", ", ".join(rgb_to_hex(c) for c in cols))

    # Pixels available to the series: the whole image (markers at the axes limits stick out of the box)
    # except the ignored regions
    valid = np.ones((H, W), bool)
    for x0, y0, x1, y1 in ignore:
        valid[int(y0) : int(y1) + 1, int(x0) : int(x1) + 1] = False

    xg, yg = np.meshgrid(np.arange(W), np.arange(H))
    outprefix = Path(args.out) if args.out else Path(__file__).resolve().parent / "outputs" / img.stem
    outprefix.parent.mkdir(parents=True, exist_ok=True)

    rows, curve_rows = [], []
    allpts = []
    for name, color, (xmin, xmax, ymin, ymax) in series:
        mask = (np.linalg.norm(rgb - color, axis=2) < args.tol) & valid
        # Region in data coordinates -> pixel mask
        X, Y = px2x(xg), py2y(yg)
        if xmin is not None: mask &= X >= xmin
        if xmax is not None: mask &= X <= xmax
        if ymin is not None: mask &= Y >= ymin
        if ymax is not None: mask &= Y <= ymax
        pts = cross_centers(mask, args.marker_size)
        pts = pts[np.argsort(pts[:, 0])] if len(pts) else pts
        print(f"{name} ({rgb_to_hex(color)}): {len(pts)} markers")
        for px, py in pts:
            x, y = float(px2x(px)), float(py2y(py))
            rows.append((name, x, y))
            print(f"    x = {x:.6g}, y = {y:.6g}")
        allpts.append((name, pts))
        if args.curve:
            yc = curve_columns(mask)
            for px in np.flatnonzero(~np.isnan(yc)):
                curve_rows.append((name, float(px2x(px)), float(py2y(yc[px]))))

    with open(f"{outprefix}.csv", "w") as f:
        f.write("series,x,y\n")
        for name, x, y in rows:
            f.write(f'"{name}",{x:.8g},{y:.8g}\n')
    print(f"Markers written to {outprefix}.csv")
    if args.curve:
        with open(f"{outprefix}_curve.csv", "w") as f:
            f.write("series,x,y\n")
            for name, x, y in curve_rows:
                f.write(f'"{name}",{x:.8g},{y:.8g}\n')
        print(f"Curves written to {outprefix}_curve.csv")

    # Check image: detected points (black ring + series color dot), plot box and ignored regions
    chk = Image.fromarray(rgb.astype(np.uint8))
    d = ImageDraw.Draw(chk)
    d.rectangle([L, T, R, B], outline=(0, 160, 0))
    for x0, y0, x1, y1 in ignore:
        d.rectangle([x0, y0, x1, y1], outline=(160, 0, 160))
    for (name, pts), (_, color, _) in zip(allpts, series):
        for px, py in pts:
            d.ellipse([px - 9, py - 9, px + 9, py + 9], outline=(0, 0, 0), width=2)
            d.ellipse([px - 2, py - 2, px + 2, py + 2], fill=(0, 0, 0))
    chk.save(f"{outprefix}_check.png")
    print(f"Check image written to {outprefix}_check.png")


if __name__ == "__main__":
    main()
