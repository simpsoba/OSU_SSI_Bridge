"""Offline Newmark vs MKR-alpha comparison in one figure (dissertation Ch. 6; read-only on data).

Combines plot/PlotLumpedHoldIntegratorPaper.py (displacement histories) and plot/PlotLumpedHoldHysteresisPaper.py
(hysteresis loops) without modifying them: same runs, stations, styling, NRMSE, and zoom windows.
  (a)         near-field mesh with all stations (one legend and one time step under it)
  (b)-(e)     horizontal displacement histories with two zoom columns (x.1), (x.2), as in the displacement figure
              (bold tag only; station names in the caption)
  (f)-(i)     hysteresis loops in one row (bold tag inside each axis, names in the caption, one-line NRMSE above):
              pier base zero-length section, p-y spring at the cap soffit, pile section below the cap, soil element
              near the base. Station letters on (a) are gray and regular weight, as in the two source figures.
The pile-head displacement node (c) and the p-y spring (g) are at the same point (x = 0, y = -0.99 m) and share one
mark. Writes <out_stem>.png / .pdf.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]  # OSU_SSI_Bridge
sys.path.insert(0, str(ROOT / "plot"))
import PlotLumpedHoldIntegratorPaper as disp  # noqa: E402
import PlotLumpedHoldHysteresisPaper as hyst  # noqa: E402
import PlotResponseSpectrum as prs  # noqa: E402

OUT = Path(sys.argv[1])

# Top block: the displacement-figure geometry (disp.*). Bottom row: four loop panels.
LOOP_TOP_GAP = 0.20  # history x label -> one-line NRMSE above the loops
PH = 1.0  # loop panel height (in); width PW follows the margins below
LOOP_X0 = 0.38  # left margin (y label + tick labels of the first loop)
LOOP_GAP = 0.50  # between loop panels (end y tick labels of the next one; its label sits beside its 0 tick)
LOOP_MR = 0.06  # right margin
PW = (disp.FIG_W - LOOP_X0 - 3 * LOOP_GAP - LOOP_MR) / 4
Y_LABEL_X = -0.20  # y label beside the narrow 0 tick label (axes fraction of PW); the wide end labels are above/below it
LOOP_MB = 0.42  # loop x labels
FIG_W = disp.FIG_W
FIG_H = disp.FIG_H + LOOP_TOP_GAP + PH + LOOP_MB


def mark_all(ax, xy, pts):
    """All stations on the mesh: one dot per location, letters off to the side with thin leaders."""
    pier, pile_head, pile_mid = xy[disp.PIER_TAG], xy[disp.PILE_TAGS[0]], xy[disp.PILE_TAGS[1]]
    soil = xy[disp.soil_tags(xy)[0]]
    spec = [  # (label, point, label position in data units, ha)
        ("(b)", pier, (5.6, pier[1] - 2.2), "left"),
        ("(c, g)", pile_head, (5.6, pile_head[1] - 1.6), "left"),
        ("(d)", pile_mid, (-4.6, pile_mid[1]), "right"),
        ("(e)", soil, (soil[0] - 1.8, soil[1] + 2.4), "right"),
        ("(f)", pts["hinge"], (-4.6, 2.2), "right"),
        ("(h)", pts["pile"], (5.6, -5.6), "left"),
        ("(i)", pts["quad"], (5.6, -21.0), "left"),
    ]
    for lab, (x, y), (lx, ly), ha in spec:
        ax.scatter([x], [y], s=10, c="#111111", marker="o", zorder=5, linewidths=0.3, edgecolors="white")
        ax.annotate(lab, xy=(x, y), xytext=(lx, ly), textcoords="data", ha=ha, va="center",
                    color=disp.STATION_GRAY, zorder=6, clip_on=False,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=disp.STATION_GRAY, shrinkA=1.5, shrinkB=1.5),
                    bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85))


# Space kept free of the zoom-box tags at the top left of each history: the tag only (no station name).
disp.TITLE_W_IN = 0.30


def axes_in(fig, x, y, w, h):
    return fig.add_axes((x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H))


def main():
    mkr, nm = disp.ROOT / disp.MKR_DUMP, disp.ROOT / disp.NM_DUMP
    prs.configure_font()
    xy = disp.read_xy(mkr / "window_nodes.txt")
    ele_rows = disp.read_eles(mkr / "window_eles.txt")
    layers = disp.read_layers(mkr / "window_quads.txt")
    order_m, order_n = disp.read_order(mkr / "disp_nodes.txt"), disp.read_order(nm / "disp_nodes.txt")
    soil = disp.soil_tags(xy)

    t_pm, u_pm = disp.load_pier_ux(mkr / "pier_node_5.out")
    t_pn, u_pn = disp.load_pier_ux(nm / "pier_node_5.out")
    t_pile_m, u_pile_m = disp.series_for(mkr, order_m, list(disp.PILE_TAGS))
    t_pile_n, u_pile_n = disp.series_for(nm, order_n, list(disp.PILE_TAGS))
    t_soil_m, u_soil_m = disp.series_for(mkr, order_m, soil)
    t_soil_n, u_soil_n = disp.series_for(nm, order_n, soil)
    rows = [("pier", "b", "Pier top", t_pm, u_pm, t_pn, u_pn),
            ("pile_head", "c", "Center pile", t_pile_m, u_pile_m[:, 0], t_pile_n, u_pile_n[:, 0]),
            ("pile", "d", "Center pile", t_pile_m, u_pile_m[:, 1], t_pile_n, u_pile_n[:, 1]),
            ("soil", "e", "Soil", t_soil_m, u_soil_m[:, 0], t_soil_n, u_soil_n[:, 0])]

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=150)
    top = FIG_H - disp.MT
    ax_a = axes_in(fig, disp.ML, top - disp.MESH_H, disp.MESH_W, disp.MESH_H)
    disp.draw_mesh(ax_a, xy, ele_rows, layers, [])
    pts = hyst.station_points(xy, dict(ele_rows))
    mark_all(ax_a, xy, pts)

    # ---- (b)-(e): displacement histories, as in PlotLumpedHoldIntegratorPaper.render(relative=False) ----
    hist, zoom1, zoom2 = [], [], []
    for i in range(disp.N_HIST):
        b = top - disp.HIST_H - i * (disp.HIST_H + disp.HIST_GAP)
        hist.append(axes_in(fig, disp.HIST_X, b, disp.HIST_W, disp.HIST_H))
        zoom1.append(axes_in(fig, disp.ZOOM_X, b, disp.ZOOM_W, disp.HIST_H))
        zoom2.append(axes_in(fig, disp.ZOOM2_X, b, disp.ZOOM_W, disp.HIST_H))
    group_lim = {}
    for group, _l, _r, t_m, u_m, t_n, u_n in rows:
        group_lim[group] = max(group_lim.get(group, 0.0), disp.window_peak_mm(t_m, u_m, t_n, u_n))
    group_lim = {g: disp.round_lim(pk) for g, pk in group_lim.items()}
    w1, w2 = disp.pick_windows(rows, group_lim)
    for i, (ax, az1, az2, (group, letter, rest, t_m, u_m, t_n, u_n)) in enumerate(zip(hist, zoom1, zoom2, rows)):
        last = i == len(rows) - 1
        disp.plot_pair(ax, t_m, u_m, t_n, u_n)
        disp.style_history(ax)
        ax.set_xlim(disp.T_MIN, disp.T_MAX)
        lim = group_lim[group]
        ax.set_ylim(-lim, lim)
        ax.set_yticks([-lim, 0.0, lim])
        tag1, tag2 = f"({letter}.1)", f"({letter}.2)"
        box1 = disp.draw_zoom(ax, az1, w1, t_m, u_m, t_n, u_n, tag1, last)
        box2 = disp.draw_zoom(ax, az2, w2, t_m, u_m, t_n, u_n, tag2, last)
        disp.place_box_tags(ax, [(tag1, *box1), (tag2, *box2)], t_m, u_m, t_n, u_n)
        # Bold tag only; station names are in the caption.
        ax.text(0.012, 0.97, rf"\textbf{{({letter})}}", transform=ax.transAxes, ha="left", va="top", zorder=7,
                bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none", alpha=0.85))
        e = disp.nrmse(t_m, u_m, t_n, u_n, disp.T_MIN, disp.T_MAX)
        print(f"({letter}) NRMSE {e:.3f} %")
        ax.text(0.012, 0.04, rf"NRMSE $= {e:.1f}\%$", transform=ax.transAxes, ha="left", va="bottom",
                fontsize=disp.ZOOM_TICK_FS, zorder=8,
                bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85))
        ax.set_xticks(np.arange(disp.T_MIN, disp.T_MAX + 1.0, 30.0))
    for ax in hist[:-1]:
        ax.tick_params(labelbottom=False)
    hist[-1].set_xlabel(r"Time, $t$ (s)")
    fig.text((disp.HIST_X - 0.47) / FIG_W, (top - 0.5 * disp.STACK_H) / FIG_H,
             r"Horizontal displacement, $u_x$ (mm)", rotation=90, ha="center", va="center")

    # ---- (f)-(i): hysteresis loops, as in PlotLumpedHoldHysteresisPaper.render(), one row ----
    s_m, s_n = hyst.load_station_series(mkr), hyst.load_station_series(nm)
    panels = [("f", "hinge", "Pier base section", r"$\theta$ (mrad)", r"$M$ (MN$\cdot$m)"),
              ("g", "spring", r"$p$--$y$ spring", r"$y$ (mm)", r"$p$ (kN)"),
              ("h", "pile", "Pile section", r"$\kappa$ ($10^{-3}$/m)", r"$M$ (MN$\cdot$m)"),
              ("i", "quad", "Soil element", r"$\gamma_{xy}$ (\%)", r"$\tau_{xy}$ (kPa)")]
    y0 = LOOP_MB
    for k, (letter, key, name, xlab, ylab) in enumerate(panels):
        ax = axes_in(fig, LOOP_X0 + k * (PW + LOOP_GAP), y0, PW, PH)
        xm, ym = s_m[key]
        xn, yn = s_n[key]
        hyst.plot_loop(ax, xm, ym, xn, yn)
        disp.style_history(ax)
        xl = hyst.nice_lim(max(np.abs(xm).max(), np.abs(xn).max()))
        yl = hyst.nice_lim(max(np.abs(ym).max(), np.abs(yn).max()))
        ax.set_xlim(-xl, xl)
        ax.set_ylim(-yl, yl)
        ax.set_xticks(hyst.sym_ticks(xl))
        ax.set_yticks(hyst.sym_ticks(yl))
        ax.axhline(0.0, color="0.8", lw=0.4, zorder=1)
        ax.axvline(0.0, color="0.8", lw=0.4, zorder=1)
        ax.set_xlabel(xlab, labelpad=1.5)
        ax.set_ylabel(ylab)
        ax.yaxis.set_label_coords(Y_LABEL_X * 1.0 / PW, 0.5)
        # Bold panel tag inside the axes (names go in the caption); one-line NRMSE above the axes.
        ax.text(0.03, 0.97, rf"\textbf{{({letter})}}", transform=ax.transAxes, ha="left", va="top", zorder=8,
                bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none", alpha=0.85))
        e = disp.nrmse(s_m["time"][0], xm, s_n["time"][0], xn, 0.0, disp.T_MAX)
        print(f"({letter}) {key}: NRMSE of x {e:.3f} %")
        ax.set_title(rf"NRMSE $= {e:.1f}\%$", loc="right", pad=2)

    # one legend and one time step, under the mesh
    integ = [Line2D([0], [0], color=disp.NM_COLOR, lw=disp.NM_LW, ls="-",
                    label="Newmark\n" + r"$\gamma = 0.5$" + "\n" + r"$\beta = 0.25$"),
             Line2D([0], [0], color="black", lw=disp.MKR_LW, ls=disp.MKR_LS,
                    label=r"MKR-$\alpha$" + "\n" + r"$\rho_\infty^{2} = 0.5$")]
    dt = disp.read_dt(mkr)
    assert abs(dt - disp.read_dt(nm)) < 1e-12
    leg = fig.legend(handles=integ, title=rf"$\Delta t = {1e3 * dt:.2f}$ ms", loc="upper left",
                     bbox_to_anchor=(0.0, (top - disp.MESH_H - disp.LEG_DROP) / FIG_H), bbox_transform=fig.transFigure,
                     handlelength=1.0, handletextpad=0.3, labelspacing=1.1, frameon=False, borderaxespad=0.0,
                     alignment="left")
    leg.get_title().set_fontsize(plt.rcParams["legend.fontsize"])

    fig.savefig(OUT.with_suffix(".png"), dpi=300)
    fig.savefig(OUT.with_suffix(".pdf"))
    print(f"figure {FIG_W:.2f} x {FIG_H:.2f} in; wrote {OUT.with_suffix('.png')}")


if __name__ == "__main__":
    main()
