#!/usr/bin/env python3
"""Goals
-----
One-off preview frame: deformed EQ window (left) + running ux histories
(right). Uses the same layout as ``PlotEQ.plot_frames`` (``DO_FRAME_HIST``).

Usage
-----
  python plot/PreviewFrameHist.py [eqOutDir] [--t SEC] [--t-model SEC]

  --t         frame time, prototype scale (t_num)
  --t-model   frame time, model scale (t / sqrt(lambda)); overrides --t

Units: N, m, s.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import PlotEQ as peq
import PlotEQParallel as pep
from lab_paths import run_eq_plots_dir

PREVIEW_DPI = 120


def main() -> int:
    """Stitch one lab dump (if needed) and write a preview PNG."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "eq",
        nargs="?",
        default=str(
            Path(__file__).resolve().parents[1]
            / "OSU_SSI_BRIDGE_DATA_LOCAL"
            / "opensees_data"
            / "r-17_20260821_1649_Storm_Wave"
        ),
    )
    ap.add_argument("--t", type=float, default=None, help="Frame time, prototype scale (s)")
    ap.add_argument(
        "--t-model",
        type=float,
        default=50.0,
        help="Frame time, model scale (s); default 50",
    )
    args = ap.parse_args()
    eq = Path(args.eq).resolve()

    tmp = None
    work = eq
    if (eq / "window_nodes.txt.0").is_file():
        np_run, metas = pep.load_np(eq)
        tmp = Path(tempfile.mkdtemp(prefix="eqmp_prev_"))
        print(f"PreviewFrameHist: stitch np={np_run} -> {tmp}")
        pep.stitch(eq, tmp, np_run, metas)
        work = tmp

    try:
        meta = peq.read_meta(work)
        tags, xy = peq.read_nodes(work)
        disp_tags = peq.read_disp_nodes(work)
        lines, quads = peq.read_eles(work)
        disp_files = meta.get("dispFiles", "window_disp.out").split()
        t, ux, uy = peq.load_window_disp(work, disp_tags, disp_files)
        if peq.SUBTRACT_T0:
            ux = peq.maybe_t0(ux)
            uy = peq.maybe_t0(uy)
        idx = {tg: i for i, tg in enumerate(disp_tags)}
        js = peq.load_spring_json(meta)
        traces = peq.pick_frame_traces(disp_tags, xy, meta, idx)
        if not traces:
            print("PreviewFrameHist: no nodes for histories", file=sys.stderr)
            return 1
        if args.t is not None:
            t_target = float(args.t)
        else:
            t_target = peq.t_model_to_proto(float(args.t_model))
        k = int(__import__("numpy").argmin(__import__("numpy").abs(t - t_target)))
        mesh = peq._frame_mesh_parts(disp_tags, xy, lines, quads, js, meta)
        amp = max(float(__import__("numpy").max(__import__("numpy").abs(ux))),
                  float(__import__("numpy").max(__import__("numpy").abs(uy))), 1e-6)
        pad = peq.SCALE * amp + 0.5
        xlim = (float(mesh["X0"].min()) - pad, float(mesh["X0"].max()) + pad)
        ylim = (float(mesh["Y0"].min()) - pad, float(mesh["Y0"].max()) + pad)
        ctx = peq._create_frame_hist_figure(t, ux, mesh, traces, xlim, ylim, 1)
        peq._update_frame_hist_figure(ctx, k, 0, ux[k], uy[k])
        out = run_eq_plots_dir(eq.name) / "frame_hist_preview.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        ctx["fig"].savefig(out, dpi=PREVIEW_DPI, facecolor="white")
        plt.close(ctx["fig"])
        print(f"PreviewFrameHist: wrote {out}")
        t_now = float(t[k])
        print(
            f"  t={t_now:.3f} s proto  ({peq.t_proto_to_model(t_now):.3f} s model)  "
            f"k={k}/{len(t)-1}  traces={len(traces)}"
        )
        for lab, tg, _c in traces:
            y = xy.get(tg, (None, None))[1]
            print(f"    tag {tg:>6}  y={y}  {lab}")
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
