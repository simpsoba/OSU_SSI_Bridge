#!/usr/bin/env python3
"""
Goals
-----
Stamp lab time (model clock, same as meaSigOS) on video t = 0.

Track one bolt on the orange cylinder relative to the fixed beam, normalize,
and cross-correlate against actuator meaSigOS. Peak lag is t_lab at frame 0.

This file is self-contained (no plot/ imports) so it can run on Sherlock
from a single .py plus the video and a slim mea npz.

Usage
-----
  # 1. Slim the F04 extract (run once, anywhere with the full npz):
  python plot/SyncVideoToLabTime.py --pack-mea \\
      --mea 0821_GusBridge_rowPos1_16Core.npz \\
      --out mea_F04_slim.npz

  # 2. Write frame 0 with a pixel grid (pick bolt / beam if defaults are off):
  python plot/SyncVideoToLabTime.py --dump-frame --video 2026-08-21-09h56m.mov

  # 3. Track + correlate (Sherlock interactive):
  python plot/SyncVideoToLabTime.py \\
      --video 2026-08-21-09h56m.mov --mea mea_F04_slim.npz \\
      --bolt 227 363 --beam 309 264 --out-dir ./video_lab_sync

Files to copy to Sherlock
-------------------------
  plot/SyncVideoToLabTime.py                         this script
  2026-08-21-09h56m.mov                              ~34 MB  (Shared Drive Videos/)
  mea_F04_slim.npz                                   from --pack-mea
    source: .../mat_extract/0821_GusBridge_rowPos1_16Core.npz  (F04)

Python: numpy scipy matplotlib imageio imageio-ffmpeg
Optional: opencv-python-headless (faster template match)

Defaults below are for this clip after auto-rotate to portrait (540 x 960):
lower-left bolt on the 2x2 plate; beam = dark stain on the gray bar (right of the cylinder).

Units: lab Time in s (model); video time in s (wall = model, not prototype).
Do not multiply the phone clip by √λ. meaSigOS Time is already model scale.
Motion in px then z-score.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

# Portrait-frame defaults for 2026-08-21-09h56m.mov (imageio decode).
DEFAULT_BOLT = (227.0, 363.0)  # lower-left plate bolt (more of the head visible)
DEFAULT_BEAM = (309.0, 264.0)  # dark stain on the gray beam, right of the cylinder
PATCH = 21
SEARCH = 48
MIN_NCC = 0.35
OVERLAP_MIN_S = 30.0


# ------------------------------------------------------------
# 1. ARRAY HELPERS
# ------------------------------------------------------------


def zscore(x: np.ndarray) -> np.ndarray:
    """Zero-mean / unit-std. Constant traces stay zero.
    Args:    x
    Returns: same shape (float)
    """
    y = np.asarray(x, dtype=float)
    m = float(np.nanmean(y))
    s = float(np.nanstd(y))
    if not np.isfinite(s) or s < 1.0e-18:
        return np.zeros_like(y)
    return (y - m) / s


def ncc_1d(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson r of two equal-length finite series.
    Args:    a, b
    Returns: r in [-1, 1], or nan
    """
    mask = np.isfinite(a) & np.isfinite(b)
    if int(mask.sum()) < 8:
        return float("nan")
    aa = zscore(a[mask])
    bb = zscore(b[mask])
    return float(np.mean(aa * bb))


def gray(frame: np.ndarray) -> np.ndarray:
    """RGB/RGBA/gray → float32 luminance.
    Args:    frame  HxW or HxWxC
    Returns: HxW float32
    """
    arr = np.asarray(frame)
    if arr.dtype == np.uint16:
        arr = (arr.astype(np.float32) * (255.0 / 65535.0)).astype(np.uint8)
    elif arr.dtype != np.uint8 and arr.max() <= 1.5:
        arr = (np.clip(arr, 0.0, 1.0) * 255.0).astype(np.uint8)
    if arr.ndim == 2:
        return arr.astype(np.float32)
    rgb = arr[:, :, :3].astype(np.float32)
    return 0.299 * rgb[:, :, 0] + 0.587 * rgb[:, :, 1] + 0.114 * rgb[:, :, 2]


def maybe_portrait(frame: np.ndarray, rotate: str) -> np.ndarray:
    """Rotate so the cylinder is vertical when rotate=auto and the frame is landscape.
    Args:    frame; rotate  auto | 0 | 90 | 180 | 270  (90 = CCW)
    Returns: possibly rotated frame
    """
    k = {"0": 0, "90": 1, "180": 2, "270": 3}.get(str(rotate), None)
    if k is None:
        h, w = frame.shape[:2]
        k = 1 if w > h else 0
    if k == 0:
        return frame
    return np.rot90(frame, k)


def extract_patch(image: np.ndarray, cx: float, cy: float, half: int) -> np.ndarray | None:
    """Centered square patch; None if it would clip the frame.
    Args:    image HxW; cx cy px; half  (patch = 2*half+1)
    Returns: (2*half+1, 2*half+1) or None
    """
    x0 = int(round(cx)) - half
    y0 = int(round(cy)) - half
    x1 = x0 + 2 * half + 1
    y1 = y0 + 2 * half + 1
    if x0 < 0 or y0 < 0 or y1 > image.shape[0] or x1 > image.shape[1]:
        return None
    return image[y0:y1, x0:x1].copy()


def ncc_map(region: np.ndarray, templ: np.ndarray) -> np.ndarray:
    """Normalized cross-correlation map (valid mode).
    Args:    region, templ  2d float
    Returns: map with shape (H-h+1, W-w+1)
    """
    try:
        import cv2

        return cv2.matchTemplate(
            region.astype(np.float32),
            templ.astype(np.float32),
            cv2.TM_CCOEFF_NORMED,
        )
    except Exception:
        pass
    from scipy.signal import fftconvolve

    t = templ.astype(np.float64)
    img = region.astype(np.float64)
    t0 = t - t.mean()
    denom_t = float(np.sqrt((t0**2).sum())) + 1.0e-12
    ones = np.ones_like(t)
    loc_sum = fftconvolve(img, ones[::-1, ::-1], mode="valid")
    loc_ssq = fftconvolve(img**2, ones[::-1, ::-1], mode="valid")
    n = float(t.size)
    loc_var = np.maximum(loc_ssq - loc_sum**2 / n, 0.0)
    num = fftconvolve(img, t0[::-1, ::-1], mode="valid")
    return (num / (np.sqrt(loc_var) * denom_t + 1.0e-12)).astype(np.float32)


def subpixel_peak(ncc: np.ndarray, iy: int, ix: int) -> tuple[float, float]:
    """Parabola peak around the integer NCC maximum.
    Args:    ncc map; iy, ix  integer peak
    Returns: (dy, dx) offset from (iy, ix)
    """

    def axis(vec: np.ndarray, i: int) -> float:
        if i <= 0 or i >= len(vec) - 1:
            return 0.0
        ym, y0, yp = float(vec[i - 1]), float(vec[i]), float(vec[i + 1])
        den = 2.0 * y0 - yp - ym
        if abs(den) < 1.0e-12:
            return 0.0
        return 0.5 * (yp - ym) / den

    return axis(ncc[:, ix], iy), axis(ncc[iy, :], ix)


def match_near(
    image: np.ndarray,
    templ: np.ndarray,
    cx: float,
    cy: float,
    search: int,
) -> tuple[float, float, float] | None:
    """NCC peak near (cx, cy).
    Args:    image; templ; cx cy last px; search  half-window (px)
    Returns: (cx, cy, ncc) or None
    """
    half = templ.shape[0] // 2
    y0 = int(np.floor(cy)) - search - half
    x0 = int(np.floor(cx)) - search - half
    y1 = int(np.ceil(cy)) + search + half + 1
    x1 = int(np.ceil(cx)) + search + half + 1
    y0c, x0c = max(0, y0), max(0, x0)
    y1c = min(image.shape[0], y1)
    x1c = min(image.shape[1], x1)
    region = image[y0c:y1c, x0c:x1c]
    if region.shape[0] < templ.shape[0] or region.shape[1] < templ.shape[1]:
        return None
    ncc = ncc_map(region, templ)
    iy, ix = np.unravel_index(int(np.argmax(ncc)), ncc.shape)
    score = float(ncc[iy, ix])
    dy, dx = subpixel_peak(ncc, int(iy), int(ix))
    # valid-mode (0,0) is the top-left of the template inside region
    cy_new = y0c + iy + dy + half
    cx_new = x0c + ix + dx + half
    return float(cx_new), float(cy_new), score


# ------------------------------------------------------------
# 2. VIDEO
# ------------------------------------------------------------


def open_reader(path: Path):
    """Open an FFMPEG reader (imageio).
    Args:    path
    Returns: reader with .get_meta_data() / iterator
    """
    import imageio.v2 as iio

    return iio.get_reader(str(path), format="FFMPEG")


def video_meta(path: Path) -> dict:
    """fps, duration, size from the container (no full decode).
    Args:    path
    Returns: dict with fps, duration_s, source_size
    """
    r = open_reader(path)
    try:
        meta = r.get_meta_data()
    finally:
        r.close()
    fps = float(meta.get("fps") or 30.0)
    dur = meta.get("duration")
    duration_s = float(dur) if dur is not None else None
    return {
        "fps": fps,
        "duration_s": duration_s,
        "source_size": meta.get("source_size") or meta.get("size"),
        "plugin_size": meta.get("size"),
    }


def iter_gray_frames(path: Path, rotate: str, stride: int, max_seconds: float | None):
    """Yield (frame_index, t_s, gray HxW).
    Args:    path; rotate; stride  (>=1); max_seconds  or None
    Returns: iterator
    """
    r = open_reader(path)
    try:
        meta = r.get_meta_data()
        fps = float(meta.get("fps") or 30.0)
        for i, frame in enumerate(r):
            if max_seconds is not None and i / fps > max_seconds:
                break
            if i % stride:
                continue
            fr = maybe_portrait(np.asarray(frame), rotate)
            yield i, i / fps, gray(fr)
    finally:
        r.close()


def first_color_frame(path: Path, rotate: str) -> np.ndarray:
    """Frame 0 as uint8 RGB, after rotate.
    Args:    path; rotate
    Returns: HxWx3
    """
    r = open_reader(path)
    try:
        frame = np.asarray(next(iter(r)))
    finally:
        r.close()
    fr = maybe_portrait(frame, rotate)
    if fr.ndim == 2:
        fr = np.stack([fr, fr, fr], axis=-1)
    if fr.dtype != np.uint8:
        mx = float(fr.max()) if fr.size else 1.0
        scale = 255.0 / mx if mx > 0 else 1.0
        fr = np.clip(fr.astype(np.float32) * scale, 0, 255).astype(np.uint8)
    return fr[:, :, :3]


# ------------------------------------------------------------
# 3. meaSigOS
# ------------------------------------------------------------


def typeconv3_column(names: list[str]) -> int | None:
    """Index of typeConv3 in stateOS signalNames.
    Args:    names
    Returns: column index or None
    """
    for i, n in enumerate(names):
        if "typeconv3" in str(n).lower().replace("_", ""):
            return i
    return None


def load_mea(path: Path) -> dict:
    """Load Time + primary (and typeConv3 when present) from a slim or full extract.
    Args:    path  .npz
    Returns: dict t_lab_s, u_m, typeConv3 (or None), mat, dt_s
    """
    z = np.load(path, allow_pickle=True)
    if "meaSigOS_time" not in z.files or "meaSigOS_primary" not in z.files:
        raise SystemExit(f"SyncVideoToLabTime: {path} missing meaSigOS_time/primary")
    t = np.asarray(z["meaSigOS_time"], dtype=float)
    u = np.asarray(z["meaSigOS_primary"], dtype=float)
    conv = None
    if "stateOS_typeConv3" in z.files:
        conv = np.asarray(z["stateOS_typeConv3"], dtype=float)
    elif "stateOS_data" in z.files:
        names = [str(x) for x in z["stateOS_signalNames"].tolist()]
        j = typeconv3_column(names)
        if j is not None:
            conv = np.asarray(z["stateOS_data"][:, j], dtype=float)
    mat = str(z["mat"]) if "mat" in z.files else path.name
    dt = float(np.median(np.diff(t))) if t.size > 1 else float("nan")
    return {"t_lab_s": t, "u_m": u, "typeConv3": conv, "mat": mat, "dt_s": dt}


def pack_mea(src: Path, dst: Path) -> None:
    """Write a slim npz (Time, primary, typeConv3) for scp.
    Args:    src  full mat_extract npz; dst
    Returns: none
    """
    mea = load_mea(src)
    payload = {
        "meaSigOS_time": mea["t_lab_s"].astype(np.float32),
        "meaSigOS_primary": mea["u_m"].astype(np.float32),
        "mat": np.array(mea["mat"]),
        "source": np.array(src.name),
    }
    if mea["typeConv3"] is not None:
        payload["stateOS_typeConv3"] = mea["typeConv3"].astype(np.int8)
    dst.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dst, **payload)
    print(
        f"pack-mea {dst}  n={mea['t_lab_s'].size}  "
        f"t_lab={mea['t_lab_s'][0]:.4g}..{mea['t_lab_s'][-1]:.4g} s  "
        f"dt={mea['dt_s']:.6g} s"
    )


# ------------------------------------------------------------
# 4. TRACK
# ------------------------------------------------------------


def track_points(
    video: Path,
    bolt: tuple[float, float],
    beam: tuple[float, float],
    *,
    rotate: str,
    stride: int,
    patch: int,
    search: int,
    max_seconds: float | None,
) -> dict:
    """Track bolt and beam; u_px = bolt_x - beam_x (portrait x = likely ux).
    Args:    video; bolt/beam (x,y) px on frame 0; rotate; stride; patch (odd); search; max_seconds
    Returns: dict with t_s, u_px, bolt_xy, beam_xy, ncc_bolt, ncc_beam, fps, shape
    """
    half = patch // 2
    it = iter_gray_frames(video, rotate, stride, max_seconds)
    try:
        i0, t0, g0 = next(it)
    except StopIteration:
        raise SystemExit("SyncVideoToLabTime: no frames")
    templ_bolt = extract_patch(g0, bolt[0], bolt[1], half)
    templ_beam = extract_patch(g0, beam[0], beam[1], half)
    if templ_bolt is None or templ_beam is None:
        raise SystemExit(
            f"SyncVideoToLabTime: patch clips frame {g0.shape} "
            f"bolt={bolt} beam={beam} patch={patch}"
        )
    bx, by = float(bolt[0]), float(bolt[1])
    mx, my = float(beam[0]), float(beam[1])
    rows: list[dict] = []
    n_frames = 0
    last_i = i0
    for i, t_s, g in _chain_first((i0, t0, g0), it):
        n_frames += 1
        last_i = i
        hit_b = match_near(g, templ_bolt, bx, by, search)
        hit_m = match_near(g, templ_beam, mx, my, search)
        if hit_b is None or hit_m is None:
            rows.append(
                {
                    "i": i,
                    "t_s": t_s,
                    "bolt_x": np.nan,
                    "bolt_y": np.nan,
                    "beam_x": np.nan,
                    "beam_y": np.nan,
                    "u_px": np.nan,
                    "ncc_bolt": np.nan,
                    "ncc_beam": np.nan,
                }
            )
            continue
        bx, by, sb = hit_b
        mx, my, sm = hit_m
        if sb < MIN_NCC:
            hit_b = match_near(g, templ_bolt, bx, by, search * 2)
            if hit_b is not None:
                bx, by, sb = hit_b
        if sm < MIN_NCC:
            hit_m = match_near(g, templ_beam, mx, my, search * 2)
            if hit_m is not None:
                mx, my, sm = hit_m
        rows.append(
            {
                "i": i,
                "t_s": t_s,
                "bolt_x": bx,
                "bolt_y": by,
                "beam_x": mx,
                "beam_y": my,
                "u_px": bx - mx,
                "ncc_bolt": sb,
                "ncc_beam": sm,
            }
        )
    meta = video_meta(video)
    return {
        "rows": rows,
        "fps": meta["fps"],
        "shape": tuple(int(x) for x in g0.shape),
        "n_decoded": n_frames,
        "last_i": last_i,
        "duration_s": meta["duration_s"],
    }


def _chain_first(first, rest):
    yield first
    yield from rest


# ------------------------------------------------------------
# 5. LAG
# ------------------------------------------------------------


def mea_clock(mea: dict, clock: str) -> np.ndarray:
    """Lab Time column, or uniform DAQ-sample time.
    Args:    mea from load_mea; clock  lab | sample
    Returns: t_s same length as mea t
    """
    t = mea["t_lab_s"]
    if clock == "lab":
        return t
    dt = mea["dt_s"]
    return np.arange(t.size, dtype=float) * dt


def best_lag(
    t_vid: np.ndarray,
    u_vid: np.ndarray,
    t_mea: np.ndarray,
    u_mea: np.ndarray,
    *,
    overlap_min_s: float,
    tau_min: float | None = None,
    tau_max: float | None = None,
) -> dict:
    """Sliding-window Pearson r of video vs mea (both signs).
    Args:    t_vid, u_vid; t_mea, u_mea; overlap_min_s; tau_min/max  (lab s, optional)
    Returns: dict tau_s, r, sign, r_of_tau, tau_grid
    """
    finite = np.isfinite(t_vid) & np.isfinite(u_vid)
    t_v = t_vid[finite]
    u_v = u_vid[finite]
    if t_v.size < 16:
        raise SystemExit("SyncVideoToLabTime: too few tracked samples")
    t_v = t_v - t_v[0]
    t_span = float(t_v[-1] - t_v[0])
    t0 = float(t_mea[0])
    t1 = float(t_mea[-1])
    tau_lo = t0 - (t_span - overlap_min_s)
    tau_hi = t1 - overlap_min_s
    if tau_min is not None:
        tau_lo = max(tau_lo, float(tau_min))
    if tau_max is not None:
        tau_hi = min(tau_hi, float(tau_max))
    if tau_hi <= tau_lo:
        raise SystemExit(
            f"SyncVideoToLabTime: empty tau range [{tau_lo:.3g}, {tau_hi:.3g}] "
            f"(overlap_min={overlap_min_s}, t_span={t_span:.3g})"
        )

    def score(tau: float, sign: float) -> float:
        t_query = t_v + tau
        u_m = np.interp(t_query, t_mea, u_mea, left=np.nan, right=np.nan)
        return ncc_1d(sign * u_v, u_m)

    coarse = np.arange(tau_lo, tau_hi + 1.0e-9, 0.05)
    best = {"r": -2.0, "tau_s": 0.0, "sign": 1.0}
    r_coarse = []
    for sign in (1.0, -1.0):
        rs = np.array([score(tau, sign) for tau in coarse], dtype=float)
        r_coarse.append(rs)
        j = int(np.nanargmax(rs))
        if rs[j] > best["r"]:
            best = {"r": float(rs[j]), "tau_s": float(coarse[j]), "sign": sign}
    fine = np.arange(best["tau_s"] - 0.4, best["tau_s"] + 0.4 + 1.0e-12, 0.002)
    fine = fine[(fine >= tau_lo) & (fine <= tau_hi)]
    rf = np.array([score(tau, best["sign"]) for tau in fine], dtype=float)
    j = int(np.nanargmax(rf))
    best["tau_s"] = float(fine[j])
    best["r"] = float(rf[j])
    best["tau_grid"] = coarse
    best["r_of_tau"] = r_coarse[0] if best["sign"] > 0 else r_coarse[1]
    best["r_of_tau_pos"] = r_coarse[0]
    best["r_of_tau_neg"] = r_coarse[1]
    best["t_vid0_s"] = 0.0
    best["tau_lo"] = tau_lo
    best["tau_hi"] = tau_hi
    best["t_span_s"] = t_span
    return best


def load_track_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read t_s, u_px, ncc_bolt from a previous track.csv.
    Args:    path
    Returns: (t_s, u_px, ncc_bolt)
    """
    t: list[float] = []
    u: list[float] = []
    ncc: list[float] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            t.append(float(row["t_s"]))
            u.append(float(row["u_px"]))
            ncc.append(float(row["ncc_bolt"]))
    return np.asarray(t), np.asarray(u), np.asarray(ncc)


# ------------------------------------------------------------
# 6. PLOTS / IO
# ------------------------------------------------------------


def write_dump_frame(
    video: Path,
    out: Path,
    bolt: tuple[float, float],
    beam: tuple[float, float],
    rotate: str,
) -> None:
    """Frame 0 with a 50 px grid and bolt/beam marks (coords in this decode).
    Args:    video; out png/jpg; bolt; beam; rotate
    Returns: none
    """
    from PIL import Image, ImageDraw

    rgb = first_color_frame(video, rotate)
    im = Image.fromarray(rgb)
    draw = ImageDraw.Draw(im)
    h, w = rgb.shape[:2]
    for x in range(0, w, 50):
        draw.line([(x, 0), (x, h)], fill=(40, 40, 40), width=1)
        draw.text((x + 2, 4), str(x), fill=(255, 255, 0))
    for y in range(0, h, 50):
        draw.line([(0, y), (w, y)], fill=(40, 40, 40), width=1)
        draw.text((4, y + 2), str(y), fill=(255, 255, 0))
    bx, by = bolt
    mx, my = beam
    r = PATCH // 2
    draw.rectangle([bx - r, by - r, bx + r, by + r], outline=(0, 255, 0), width=2)
    draw.rectangle([mx - r, my - r, mx + r, my + r], outline=(0, 180, 255), width=2)
    draw.text((bx + r + 4, by - 8), "bolt", fill=(0, 255, 0))
    draw.text((mx + r + 4, my - 8), "beam", fill=(0, 180, 255))
    out.parent.mkdir(parents=True, exist_ok=True)
    im.save(out, quality=90)
    print(f"dump-frame {out}  shape={w}x{h}  bolt={bolt}  beam={beam}")


def write_track_csv(path: Path, rows: list[dict]) -> None:
    """One row per tracked frame.
    Args:    path; rows
    Returns: none
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = [
        "i",
        "t_s",
        "bolt_x",
        "bolt_y",
        "beam_x",
        "beam_y",
        "u_px",
        "ncc_bolt",
        "ncc_beam",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in keys})


def write_overlay(
    path: Path,
    t_vid: np.ndarray,
    u_vid: np.ndarray,
    t_mea: np.ndarray,
    u_mea: np.ndarray,
    lag: dict,
) -> None:
    """NCC vs tau (left) and shifted overlay (right).
    Args:    path; series; lag from best_lag
    Returns: none
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tau = lag["tau_s"]
    sign = lag["sign"]
    u_v = sign * u_vid
    t_shift = t_vid + tau
    u_m_at = np.interp(t_shift, t_mea, u_mea, left=np.nan, right=np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), layout="constrained")
    ax = axes[0]
    ax.plot(lag["tau_grid"], lag["r_of_tau_pos"], color="#455a64", lw=0.8, label="+u")
    ax.plot(lag["tau_grid"], lag["r_of_tau_neg"], color="#c45c12", lw=0.8, label="−u")
    ax.axvline(tau, color="#001F3F", lw=0.8)
    ax.set_xlabel(r"$t_{\mathrm{lab}}$ at video $t=0$ (s)")
    ax.set_ylabel("Pearson $r$")
    ax.set_title("cross-correlation")
    ax.legend(loc="lower right")

    ax = axes[1]
    ax.plot(
        t_shift,
        zscore(u_m_at),
        color="#001F3F",
        lw=0.9,
        label="meaSigOS",
        zorder=2,
    )
    ax.plot(
        t_shift,
        zscore(u_v),
        color="#c45c12",
        lw=0.8,
        alpha=0.85,
        label="bolt−beam (px)",
        zorder=3,
    )
    ax.set_xlabel(r"$t_{\mathrm{lab}}$ (s)")
    ax.set_ylabel("z-score")
    ax.set_title(rf"overlay  $\tau={tau:.3f}$ s  $r={lag['r']:.3f}$")
    ax.legend(loc="upper right")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ------------------------------------------------------------
# 7. CLI
# ------------------------------------------------------------


def parse_xy(vals: list[float], default: tuple[float, float]) -> tuple[float, float]:
    """CLI pair or default.
    Args:    vals  length 0 or 2; default
    Returns: (x, y)
    """
    if not vals:
        return default
    if len(vals) != 2:
        raise SystemExit("expected two numbers: x y")
    return float(vals[0]), float(vals[1])


def main(argv: list[str] | None = None) -> int:
    """Pack, dump-frame, or track+correlate."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pack-mea", action="store_true", help="Slim a full mat_extract npz")
    ap.add_argument("--dump-frame", action="store_true", help="Write frame 0 with a px grid")
    ap.add_argument("--video", type=Path, help="Phone .mov / .mp4")
    ap.add_argument("--mea", type=Path, help="Slim or full meaSigOS npz")
    ap.add_argument("--out", type=Path, help="pack-mea / dump-frame output path")
    ap.add_argument("--out-dir", type=Path, default=Path("video_lab_sync"))
    ap.add_argument("--bolt", type=float, nargs=2, metavar=("X", "Y"), help="Bolt px on frame 0")
    ap.add_argument("--beam", type=float, nargs=2, metavar=("X", "Y"), help="Beam px on frame 0")
    ap.add_argument("--rotate", default="auto", help="auto | 0 | 90 | 180 | 270 (90=CCW)")
    ap.add_argument("--stride", type=int, default=2, help="Use every Nth frame (default 2)")
    ap.add_argument("--patch", type=int, default=PATCH)
    ap.add_argument("--search", type=int, default=SEARCH)
    ap.add_argument("--clock", choices=("lab", "sample"), default="lab")
    ap.add_argument("--max-seconds", type=float, default=None, help="Track only the first T s (debug)")
    ap.add_argument(
        "--overlap-min",
        type=float,
        default=OVERLAP_MIN_S,
        help="Min overlap of video and mea (s). Short windows inflate Pearson r.",
    )
    ap.add_argument(
        "--tau-min",
        type=float,
        default=None,
        help="Min t_lab at video t=0 (s). Use 0 to forbid a pre-mat start.",
    )
    ap.add_argument(
        "--tau-max",
        type=float,
        default=None,
        help="Max t_lab at video t=0 (s)",
    )
    ap.add_argument(
        "--from-track",
        type=Path,
        help="Skip video decode; reuse track.csv and only recompute the lag",
    )
    args = ap.parse_args(argv)

    bolt = parse_xy(args.bolt or [], DEFAULT_BOLT)
    beam = parse_xy(args.beam or [], DEFAULT_BEAM)

    if args.pack_mea:
        if args.mea is None or args.out is None:
            raise SystemExit("SyncVideoToLabTime: --pack-mea needs --mea and --out")
        pack_mea(args.mea, args.out)
        return 0

    if args.dump_frame:
        if args.video is None:
            raise SystemExit("SyncVideoToLabTime: --dump-frame needs --video")
        out = args.out or Path("frame0_grid.jpg")
        write_dump_frame(args.video, out, bolt, beam, args.rotate)
        return 0

    if args.mea is None:
        raise SystemExit("SyncVideoToLabTime: need --mea")
    if args.from_track is None and args.video is None:
        raise SystemExit("SyncVideoToLabTime: need --video or --from-track")

    mea = load_mea(args.mea)
    t_mea = mea_clock(mea, args.clock)
    print(
        f"mea {args.mea.name}  mat={mea['mat']}  n={t_mea.size}  "
        f"t={t_mea[0]:.4g}..{t_mea[-1]:.4g} s  dt={mea['dt_s']:.6g}  clock={args.clock}"
    )
    if mea["typeConv3"] is not None:
        frac2 = float(np.mean(mea["typeConv3"] == 2))
        print(f"  typeConv3==2 sample fraction {frac2:.4g}")

    rows: list[dict] | None = None
    tr: dict | None = None
    if args.from_track is not None:
        t_vid, u_vid, ncc_b = load_track_csv(args.from_track)
        print(f"from-track {args.from_track}  n={t_vid.size}")
    else:
        print(f"track {args.video.name}  bolt={bolt}  beam={beam}  stride={args.stride}")
        tr = track_points(
            args.video,
            bolt,
            beam,
            rotate=args.rotate,
            stride=max(1, args.stride),
            patch=args.patch if args.patch % 2 else args.patch + 1,
            search=args.search,
            max_seconds=args.max_seconds,
        )
        rows = tr["rows"]
        t_vid = np.array([r["t_s"] for r in rows], dtype=float)
        u_vid = np.array([r["u_px"] for r in rows], dtype=float)
        ncc_b = np.array([r["ncc_bolt"] for r in rows], dtype=float)
        print(
            f"  frames kept {len(rows)}  shape={tr['shape']}  fps={tr['fps']}  "
            f"ncc_bolt median={np.nanmedian(ncc_b):.3f}"
        )

    lag = best_lag(
        t_vid,
        u_vid,
        t_mea,
        mea["u_m"],
        overlap_min_s=args.overlap_min,
        tau_min=args.tau_min,
        tau_max=args.tau_max,
    )
    print(
        f"  tau search [{lag['tau_lo']:.3f}, {lag['tau_hi']:.3f}] s  "
        f"video span {lag['t_span_s']:.3f} s"
    )
    print(
        f"  t_lab at video t=0: {lag['tau_s']:.4f} s  "
        f"r={lag['r']:.4f}  sign={int(lag['sign'])}"
    )
    # Checkpoint: r at +30 s lab (GM onset is ~26 s; D5–95 starts ~36 s).
    t_v = t_vid[np.isfinite(t_vid) & np.isfinite(u_vid)]
    u_v = u_vid[np.isfinite(t_vid) & np.isfinite(u_vid)]
    t_v = t_v - t_v[0]
    for probe in (30.0,):
        tq = t_v + probe
        um = np.interp(tq, t_mea, mea["u_m"], left=np.nan, right=np.nan)
        rp = ncc_1d(u_v, um)
        rn = ncc_1d(-u_v, um)
        print(f"  r at tau={probe:g} s:  +u {rp:.4f}  -u {rn:.4f}")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if rows is not None:
        write_track_csv(out_dir / "track.csv", rows)
    write_overlay(
        out_dir / "overlay.png",
        t_vid,
        u_vid,
        t_mea,
        mea["u_m"],
        lag,
    )
    result = {
        "video": str(args.video) if args.video else None,
        "from_track": str(args.from_track) if args.from_track else None,
        "mea": str(args.mea),
        "mat": mea["mat"],
        "clock": args.clock,
        "tau_min": args.tau_min,
        "tau_max": args.tau_max,
        "overlap_min_s": args.overlap_min,
        "bolt_xy": list(bolt),
        "beam_xy": list(beam),
        "stride": args.stride,
        "t_lab_at_video_t0_s": lag["tau_s"],
        "pearson_r": lag["r"],
        "sign": lag["sign"],
        "ncc_bolt_median": float(np.nanmedian(ncc_b)),
        "n_track": int(t_vid.size),
    }
    if tr is not None:
        result["fps"] = tr["fps"]
        result["frame_shape"] = list(tr["shape"])
    (out_dir / "lag.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (out_dir / "lag.txt").write_text(
        f"t_lab_at_video_t0_s  {lag['tau_s']:.6f}\n"
        f"pearson_r            {lag['r']:.6f}\n"
        f"sign                 {lag['sign']:.0f}\n",
        encoding="utf-8",
    )
    print(f"wrote {out_dir / 'lag.txt'}  {out_dir / 'overlay.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
