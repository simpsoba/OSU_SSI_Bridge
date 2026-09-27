#!/usr/bin/env python3
"""
Goals
-----
Elastic response spectrum (5% damping) for the Tohoku motion used in
Parameters.tcl: PEER VT2 velocity -> acceleration -> Sd, PSv, PSa vs T.

Page: 6 x 2.2 in, one row of three panels, 9 pt. usetex hybrid (CM Sans
text + CM serif math). Spectrum CSV cached next to the figure stem.

  python plot/PlotResponseSpectrum.py
  python plot/PlotResponseSpectrum.py path/to/file.VT2 [out_stem]
  python plot/PlotResponseSpectrum.py --recompute

Default out:
  plot/out/spectrum/tohoku_FKSH19_NS1_SdSvSa.{png,pdf,csv}
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# 1. PAGE, FONT, COLORS
# ------------------------------------------------------------

FIG_W = 6.0  # in
FIG_H = 1.9  # in — one row of three panels
FONT_SIZE = 9

COLORS = ("#385F96", "#CF5921", "#9EB8DB", "#E7B800", "#800000")

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DEFAULT_VT2 = REPO / "ground-motion" / "Tohoku2011-FKSH" / "FKSH19.NS1.VT2"
OUT_DIR = HERE / "out" / "spectrum"
TEXINPUTS_DIR = HERE / "texinputs"  # type1cm.sty -> type1ec (matplotlib usetex)
TEXLIVE_BIN = Path(r"C:\texlive\2026\bin\windows")

G = 9.80665  # m/s^2
ZETA = 0.05
T_MIN = 0.02  # s
T_MAX = 10.0  # s
N_PERIOD = 200
T1 = 2.3  # s — rayleighT1 / SSI system


def _ensure_tex_path() -> None:
    """Put TeX Live on PATH and prepend plot/texinputs for the type1cm shim."""
    if TEXLIVE_BIN.is_dir():
        os.environ["PATH"] = str(TEXLIVE_BIN) + os.pathsep + os.environ.get("PATH", "")
    if TEXINPUTS_DIR.is_dir():
        prev = os.environ.get("TEXINPUTS", "")
        prefix = str(TEXINPUTS_DIR) + os.pathsep
        if not prev.startswith(prefix):
            os.environ["TEXINPUTS"] = prefix + prev


def configure_font() -> None:
    """usetex: CM Sans text, CM serif math ($S_d$, $T_1$, …)."""
    _ensure_tex_path()
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Computer Modern Sans Serif"],
        "font.size": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE,
        "text.latex.preamble": r"""
\usepackage[T1]{fontenc}
\renewcommand{\familydefault}{\sfdefault}
\usepackage{amsmath}
""",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "grid.color": "0.85",
        "axes.axisbelow": True,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.5,
        "pdf.fonttype": 42,
        "figure.autolayout": False,
    })


# ------------------------------------------------------------
# 2. GROUND MOTION
# ------------------------------------------------------------

def read_peer_vt2(path: Path) -> tuple[np.ndarray, float]:
    """
    Read a PEER NGA velocity (.VT2) file.

    Args:    path  PEER VT2 path
    Returns: vel_mps (m/s), dt (s)
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 5 or "PEER" not in lines[0].upper():
        raise ValueError(f"expected PEER header in {path}")
    # NPTS=  30000, DT=   0.0100 SEC,
    meta = lines[3]
    npts = int(meta.split("NPTS=")[1].split(",")[0].strip())
    dt = float(meta.split("DT=")[1].split("SEC")[0].strip().rstrip(","))
    vals: list[float] = []
    for line in lines[4:]:
        for tok in line.split():
            vals.append(float(tok))
    if len(vals) < npts:
        raise ValueError(f"{path}: got {len(vals)} samples, header says {npts}")
    vel_cm_s = np.asarray(vals[:npts], dtype=float)
    return vel_cm_s / 100.0, dt  # cm/s -> m/s


def accel_from_vel(vel: np.ndarray, dt: float) -> np.ndarray:
    """
    Central-difference acceleration from a velocity series.

    Args:    vel (m/s), dt (s)
    Returns: acc (m/s^2)
    """
    return np.gradient(vel, dt)


# ------------------------------------------------------------
# 3. ELASTIC SDF SPECTRUM
# ------------------------------------------------------------

def newmark_umax(acc: np.ndarray, dt: float, T: float, zeta: float) -> float:
    """
    Peak relative displacement of an elastic SDF under ground accel.

    Newmark average acceleration (gamma=1/2, beta=1/4), same integrator
    family as OpenSees sdfResponse for the elastic (Fy -> inf) case.

    Args:    acc (m/s^2), dt (s), T (s), zeta (-)
    Returns: umax (m)
    """
    if T <= 0.0:
        raise ValueError("T must be > 0")
    m = 1.0
    omega = 2.0 * np.pi / T
    k = m * omega * omega
    c = 2.0 * zeta * m * omega
    gamma = 0.5
    beta = 0.25

    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = dt * (gamma / (2.0 * beta) - 1.0)
    k_hat = k + a0 * m + a1 * c

    u = 0.0
    v = 0.0
    # relative accel at t=0: ü = -ug - (c ú + k u)/m
    a = -acc[0]
    umax = abs(u)

    for i in range(1, len(acc)):
        p_hat = -m * acc[i] + m * (a0 * u + a2 * v + a3 * a) + c * (
            a1 * u + a4 * v + a5 * a
        )
        u_new = p_hat / k_hat
        a_new = a0 * (u_new - u) - a2 * v - a3 * a
        v_new = v + dt * ((1.0 - gamma) * a + gamma * a_new)
        u, v, a = u_new, v_new, a_new
        if abs(u) > umax:
            umax = abs(u)
    return umax


def spectrum_sdfresponse(
    acc: np.ndarray, dt: float, periods: np.ndarray, zeta: float
) -> np.ndarray:
    """
    Sd via openseespy.opensees.sdfResponse (elastic: Fy huge).

    Args:    acc (m/s^2), dt, periods (s), zeta
    Returns: Sd (m) array
    """
    import openseespy.opensees as ops  # noqa: WPS433

    m = 1.0
    Fy = 1.0e16
    alpha = 0.0
    # force = -m * ugddot, one value per line
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False, encoding="utf-8"
    ) as fh:
        force_path = fh.name
        for a in acc:
            fh.write(f"{-m * a:.10e}\n")
    try:
        Sd = np.empty(len(periods), dtype=float)
        for i, T in enumerate(periods):
            omega = 2.0 * np.pi / T
            k = m * omega * omega
            umax, _u, _up, _amax, _tamax = ops.sdfResponse(
                m, zeta, k, Fy, alpha, dt, force_path, dt
            )
            Sd[i] = abs(umax)
        return Sd
    finally:
        Path(force_path).unlink(missing_ok=True)


def spectrum_newmark(
    acc: np.ndarray, dt: float, periods: np.ndarray, zeta: float
) -> np.ndarray:
    """Sd via local Newmark (fallback when OpenSeesPy is unavailable)."""
    return np.array(
        [newmark_umax(acc, dt, float(T), zeta) for T in periods], dtype=float
    )


def compute_spectrum(
    acc: np.ndarray, dt: float, periods: np.ndarray, zeta: float = ZETA
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """
    Elastic pseudo-spectral Sd, Sv, Sa.

    Returns: Sd (m), Sv (m/s), Sa (g), method label
    """
    method = "Newmark (local)"
    try:
        Sd = spectrum_sdfresponse(acc, dt, periods, zeta)
        method = "OpenSeesPy sdfResponse"
    except Exception:
        Sd = spectrum_newmark(acc, dt, periods, zeta)

    omega = 2.0 * np.pi / periods
    Sv = omega * Sd
    Sa_g = (omega * omega * Sd) / G
    return Sd, Sv, Sa_g, method


# ------------------------------------------------------------
# 4. CSV CACHE + FIGURE
# ------------------------------------------------------------

def spectrum_csv_path(out_stem: Path) -> Path:
    return out_stem.with_suffix(".csv")


def save_spectrum_csv(
    path: Path,
    periods: np.ndarray,
    Sd: np.ndarray,
    Sv: np.ndarray,
    Sa_g: np.ndarray,
    *,
    vt2: Path,
    dt: float,
    method: str,
    pga_g: float,
) -> None:
    """Write T, Sd, Sv, Sa so the figure can be redrawn without re-integrating."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        f"# vt2={vt2}\n"
        f"# dt_s={dt:g}  zeta={ZETA:g}  method={method}  pga_g={pga_g:.6f}\n"
        f"# T_min={T_MIN:g}  T_max={T_MAX:g}  n={len(periods)}\n"
        "# columns: T_s, Sd_m, Sv_mps, Sa_g\n"
    )
    data = np.column_stack([periods, Sd, Sv, Sa_g])
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(header)
        np.savetxt(f, data, delimiter=",", fmt="%.8g")


def load_spectrum_csv(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, str] | None:
    """
    Load cached spectrum.

    Returns: periods, Sd, Sv, Sa_g, pga_g, method — or None if unusable.
    """
    if not path.is_file():
        return None
    pga_g = float("nan")
    method = "csv"
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") and "pga_g=" in line:
                pga_g = float(line.split("pga_g=")[1].split()[0])
            if line.startswith("#") and "method=" in line:
                method = line.split("method=")[1].split()[0]
            if not line.startswith("#"):
                break
    arr = np.loadtxt(path, delimiter=",", comments="#", skiprows=0)
    if arr.ndim != 2 or arr.shape[1] < 4:
        return None
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], pga_g, method


def plot_spectrum(
    periods: np.ndarray,
    Sd: np.ndarray,
    Sv: np.ndarray,
    Sa_g: np.ndarray,
    out_stem: Path,
) -> None:
    """
    One row: Sd | Sv | Sa vs T (log–log). No in-figure marks.

    Args:    periods (s), Sd (m), Sv (m/s), Sa_g (g), out_stem
    Returns: none (writes png + pdf)
    """
    configure_font()
    fig, axes = plt.subplots(
        1, 3, figsize=(FIG_W, FIG_H), sharex=True, layout="constrained"
    )

    from matplotlib.ticker import LogFormatterSciNotation

    log_fmt = LogFormatterSciNotation(base=10, labelOnlyBase=True)

    series = (
        (Sd, r"$S_d$ (m)"),
        (Sv, r"$S_v$ (m/s)"),
        (Sa_g, r"$S_a$ (g)"),
    )

    for ax, (y, ylabel) in zip(axes, series):
        y = np.asarray(y, dtype=float)
        y_pos = np.maximum(y, np.nanmax(y) * 1e-6)
        ax.plot(periods, y_pos, color=COLORS[0], lw=1.5)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(r"$T$ (s)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(T_MIN, T_MAX)
        ymax = float(np.nanmax(y))
        ymin = float(np.nanmin(y[y > 0])) if np.any(y > 0) else ymax * 1e-3
        ax.set_ylim(ymin * 0.8, ymax * 1.4)
        ax.xaxis.set_major_formatter(log_fmt)
        ax.yaxis.set_major_formatter(log_fmt)
        ax.grid(True, which="major")
        ax.grid(False, which="minor")

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = out_stem.with_suffix(".pdf")
    png_path = out_stem.with_suffix(".png")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    _png_from_pdf(pdf_path, png_path, dpi=300)


def _png_from_pdf(pdf_path: Path, png_path: Path, dpi: int = 300) -> None:
    """
    Rasterize PDF -> PNG. Prefer Ghostscript (clean NewCM Type-42); else pdftoppm.

    Args:    pdf_path, png_path, dpi
    Returns: none
    """
    import shutil
    import subprocess

    gs = shutil.which("gswin64c") or shutil.which("gs")
    if gs:
        subprocess.run(
            [
                gs,
                "-dSAFER",
                "-dBATCH",
                "-dNOPAUSE",
                "-dQUIET",
                "-sDEVICE=png16m",
                f"-r{dpi}",
                f"-sOutputFile={png_path}",
                str(pdf_path),
            ],
            check=True,
        )
        return

    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm:
        stem = png_path.with_suffix("")
        subprocess.run(
            [pdftoppm, "-png", "-r", str(dpi), "-singlefile", str(pdf_path), str(stem)],
            check=True,
        )
        return

    raise RuntimeError("need Ghostscript or pdftoppm to rasterize PNG from PDF")


def pier_sdof_period() -> tuple[float, float, dict[str, float]]:
    """
    Fixed-base cantilever: tip mass = m_deck + half pier mass, k = 3EI/H^3.

    Returns: T_cracked (s), T_gross/uncracked (s), intermediates (SI).
    Cracked uses pierCrackedFactor=0.5 on transformed I (Parameters.tcl).
    """
    inch = 0.0254
    foot = 12.0 * inch
    pi = np.pi

    D_pier = 4.0 * foot
    H_pier = 3.02 * ((4.0 * 12.0) / 20.0)
    fc = 28.0e6
    Ec = 4700.0 * np.sqrt(fc / 1.0e6) * 1.0e6
    Es = 200.0e9
    dens_c, dens_s = 2400.0, 7850.0
    cover = 2.0 * inch
    db_long = (10.0 / 8.0) * inch
    db_tran = (7.0 / 8.0) * inch
    n_long = 28
    As_tot = n_long * (pi / 4.0) * db_long**2
    A_g = (pi / 4.0) * D_pier**2
    R_core = 0.5 * D_pier - cover - 0.5 * db_long - db_tran
    I_g = (pi / 64.0) * D_pier**4
    I_s = 0.5 * As_tot * R_core**2
    I_uncr = I_g + (Es / Ec - 1.0) * I_s
    cracked_factor = 0.5  # pierCrackedFactor
    I_cr = cracked_factor * I_uncr
    rhoL = dens_c * (A_g - As_tot) + dens_s * As_tot
    M_pier = rhoL * H_pier

    A_deck = 8869.0 * 6.4516e-4
    dens_deck = 22.78e3 / 9.81
    m_deck = dens_deck * A_deck * (150.0 * foot)
    M_eq = m_deck + 0.5 * M_pier
    k_cr = 3.0 * Ec * I_cr / H_pier**3
    k_uncr = 3.0 * Ec * I_uncr / H_pier**3
    T_cr = 2.0 * pi * np.sqrt(M_eq / k_cr)
    T_uncr = 2.0 * pi * np.sqrt(M_eq / k_uncr)
    return float(T_cr), float(T_uncr), {
        "H_pier": float(H_pier),
        "m_deck": float(m_deck),
        "M_pier": float(M_pier),
        "M_eq": float(M_eq),
        "I_g": float(I_g),
        "I_uncr": float(I_uncr),
        "I_cr": float(I_cr),
        "Ec": float(Ec),
        "k_cr": float(k_cr),
        "k_uncr": float(k_uncr),
        "cracked_factor": cracked_factor,
    }


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    recompute = "--recompute" in args
    if recompute:
        args.remove("--recompute")

    vt2 = Path(args[0]) if args else DEFAULT_VT2
    if not vt2.is_file():
        print(f"VT2 not found: {vt2}", file=sys.stderr)
        return 1

    if len(args) >= 2:
        out_stem = Path(args[1])
    else:
        out_stem = OUT_DIR / "tohoku_FKSH19_NS1_SdSvSa"

    csv_path = spectrum_csv_path(out_stem)
    cached = None if recompute else load_spectrum_csv(csv_path)

    if cached is not None:
        periods, Sd, Sv, Sa_g, pga_g, method = cached
        vel = None
        dt = float("nan")
        print(f"cache    {csv_path}")
    else:
        vel, dt = read_peer_vt2(vt2)
        acc = accel_from_vel(vel, dt)
        pga_g = float(np.max(np.abs(acc)) / G)
        periods = np.logspace(np.log10(T_MIN), np.log10(T_MAX), N_PERIOD)
        Sd, Sv, Sa_g, method = compute_spectrum(acc, dt, periods, ZETA)
        save_spectrum_csv(
            csv_path, periods, Sd, Sv, Sa_g, vt2=vt2, dt=dt, method=method, pga_g=pga_g
        )
        print(f"wrote    {csv_path}")

    def at_T(y: np.ndarray, T: float) -> float:
        return float(np.interp(T, periods, y))

    sa_t1_g = at_T(Sa_g, T1)
    sv_t1 = at_T(Sv, T1)
    sd_t1 = at_T(Sd, T1)

    plot_spectrum(periods, Sd, Sv, Sa_g, out_stem)

    print(f"file     {vt2}")
    if vel is not None:
        print(f"dt       {dt:g} s   npts={len(vel)}")
        print(f"PGA      {float(np.max(np.abs(accel_from_vel(vel, dt)))):.4f} m/s^2  ({pga_g:.4f} g)")
        print(f"PGV      {float(np.max(np.abs(vel))):.4f} m/s")
    else:
        print(f"PGA      {pga_g:.4f} g  (from csv)")
    print(f"method   {method}")
    print(f"Sd max   {float(np.max(Sd)):.4f} m at T={float(periods[np.argmax(Sd)]):.3f} s")
    print(f"Sa max   {float(np.max(Sa_g)):.4f} g at T={float(periods[np.argmax(Sa_g)]):.3f} s")
    print("---")
    print(f"T_1  (rayleighT1)  {T1:.3f} s")
    print(f"  Sa={sa_t1_g:.3f} g  Sv={sv_t1:.3f} m/s  Sd={sd_t1:.3f} m")
    print(f"wrote    {out_stem.with_suffix('.png')}")
    print(f"wrote    {out_stem.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
