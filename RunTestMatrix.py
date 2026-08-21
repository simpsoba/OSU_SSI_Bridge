# RunTestMatrix.py
#
# How to use:
#   python RunTestMatrix.py --row 1
#   python RunTestMatrix.py --row=-1
# writes Overrides.tcl for that Run id (from TestMatrix.csv).
#
# Run sign is the EQ constraint handler: +N = Auto, -N = Transformation.
# Same |N| is the same case. IDs start at 1. A case may omit +N (Transformation
# only; stiffer soil is -13). Use --row=-N (equals); --row -1 is parsed as a flag.
#
# How to override more knobs later:
#   1. Add a column to TestMatrix.csv (same name as the Tcl variable when possible).
#   2. Add that name to PASS_THROUGH or to a special-case block below.
#   3. Re-run this script for the row you want.
# Leave cells blank to keep the Parameters.tcl / Run*.tcl default.
#
# soilProfile / soilMesh cells may be bare numbers (1) or labeled (0 (PRODUCTION));
# both work. Wave Name lookup (Storm Wave, Big Tsunami): see WaveCatalog.csv
# (prototype vs lab scale). Not written into Overrides.tcl yet.
# Rayleigh: T/xi/offFac/stiff + region *ON columns; see analysis/RayleighDamping.tcl.
#
# This script only WRITES Overrides.tcl. You launch OpenSees yourself, e.g.:
#   mpiexec -n 8 OpenSeesMP RunParallel.tcl Overrides.tcl
#   OpenSees Run.tcl Overrides.tcl

import argparse
import csv
import datetime
import os
import re
import sys

# Folder that holds this script (= repo root)
HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(HERE, "TestMatrix.csv")
OUT_PATH = os.path.join(HERE, "Overrides.tcl")

# ---------------------------------------------------------------------------
# Columns that become plain:  set name value
# (Add a new Tcl knob name here when you add a matching CSV column.)
# ---------------------------------------------------------------------------
PASS_THROUGH = [
    "realTimeON",
    "expElementType",
    "holdPierON",
    "gmStartTime",
    "DT_FACTOR",
    "gmScaleFactor",
    "h_water",
    "eqIntegrator",
    "prePartitionSystem",
    "postPartitionSystem",
    "constraintsHandler",
    # Rayleigh (analysis/RayleighDamping.tcl). T2 filled as sqrt(cylinderSF)/20 today.
    "rayleighT1",
    "rayleighT2",
    "rayleighXi1",
    "rayleighXi2",
    "rayleighOffFac",
    "rayleighStiff",
    "rayleighSoilON",
    "rayleighBoundON",
    "rayleighSprPileON",
    "rayleighSprCapFaceON",
    "rayleighSprSoffitON",
    "rayleighPilesON",
    "rayleighCapON",
    "rayleighDeckON",
    "rayleighPierON",
    "rayleighPierHingeON",
]

# Columns we read for comments / hints / outDIR slug — not written as set lines
META_COLUMNS = ["Run", "File", "Number of Procs", "Name", "Goal"]

# Columns handled in special blocks below (not in PASS_THROUGH)
SPECIAL_COLUMNS = [
    "soilProfile",
    "soilEleType",
    "soilConstitutive",
    "soilMesh",
    "gmDir",
    "gmVelFile",
    "outDIR",
]


def clean_cell(raw):
    """Empty / missing CSV cells -> None."""
    if raw is None:
        return None
    s = str(raw).strip()
    if s == "":
        return None
    return s


def strip_label(raw):
    """Accept '1', '1 (FINE)', or '4 (SOFT)' -> the leading integer as a string."""
    s = clean_cell(raw)
    if s is None:
        return None
    m = re.match(r"^([+-]?\d+)\s*\(", s)
    if m:
        return m.group(1)
    # bare number (or anything else) — pass through as-is
    return s


def strip_paren_note(raw):
    """'UmfPack (or ProfileSPD)' -> 'UmfPack'."""
    s = clean_cell(raw)
    if s is None:
        return None
    if " (" in s:
        s = s.split(" (", 1)[0].strip()
    return s


def enum_soil_ele(raw):
    s = clean_cell(raw)
    if s is None:
        return None
    low = s.lower()
    if low in ("sspquad", "ssquad"):
        return "SSPquad"
    if low == "quad":
        return "quad"
    return s


def enum_constitutive(raw):
    s = clean_cell(raw)
    if s is None:
        return None
    return s.lower()


def tcl_quote(s):
    """Wrap a string for Tcl: set x \"...\""""
    return '"' + s.replace("\\", "/").replace('"', '\\"') + '"'


def tcl_number_or_string(s):
    """Use bare number if it looks numeric; else a quoted string."""
    if re.match(r"^[+-]?\d+(\.\d+)?([eE][+-]?\d+)?$", s):
        return s
    return tcl_quote(s)


def slug(text, max_len=24):
    """Filesystem-safe short name from Name/Goal."""
    if not text:
        return ""
    out = []
    for ch in text.strip():
        if ch.isalnum():
            out.append(ch)
        elif ch in (" ", "-", "_") and out and out[-1] != "_":
            out.append("_")
    s = "".join(out).strip("_")
    return s[:max_len]


def canon_run_id(raw):
    """'8'/'+8'/'8.0' -> '8'; '-8' -> '-8'. No run 0."""
    s = clean_cell(raw)
    if s is None:
        return None
    m = re.match(r"^([+-])?(\d+)(?:\.0+)?$", s)
    if not m:
        return None
    sign, digits = m.group(1), m.group(2)
    mag = int(digits)
    if mag == 0:
        return None
    if sign == "-":
        return "-" + str(mag)
    return str(mag)


def run_handler(canon):
    """+N -> Auto; -N -> Transformation."""
    if canon.startswith("-"):
        return "Transformation"
    return "Auto"


def run_out_tag(canon):
    """Filesystem slug: r+01, r-01, r+08, r-13."""
    if canon.startswith("-"):
        return "r-" + canon[1:].zfill(2)
    return "r+" + canon.zfill(2)


def find_row(rows, row_id):
    want = canon_run_id(row_id)
    if want is None:
        return None
    for r in rows:
        got = canon_run_id(r.get("Run"))
        if got is not None and got == want:
            return r
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Write Overrides.tcl from one TestMatrix.csv row."
    )
    parser.add_argument(
        "--row",
        type=str,
        required=True,
        help="Run id from TestMatrix.csv (1, -1, 8, -8). Use --row=-N for negatives.",
    )
    parser.add_argument(
        "--csv",
        default=CSV_PATH,
        help="Path to TestMatrix.csv (default: next to this script)",
    )
    parser.add_argument(
        "--out",
        default=OUT_PATH,
        help="Path for Overrides.tcl (default: next to this script)",
    )
    parser.add_argument(
        "--outDIR",
        default=None,
        help="Force exact outDIR (skips row/timestamp naming)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        print("ERROR: cannot find CSV:", args.csv, file=sys.stderr)
        sys.exit(1)

    with open(args.csv, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    row_id = canon_run_id(args.row)
    if row_id is None:
        print("ERROR: bad Run id", repr(args.row), "(want 1, -1, 8, -8, ...; no 0)", file=sys.stderr)
        sys.exit(1)
    row = find_row(rows, row_id)
    if row is None:
        print("ERROR: no Run ==", row_id, "in", args.csv, file=sys.stderr)
        sys.exit(1)

    # Warn about unknown columns (ignored on purpose, e.g. wave metadata)
    known = set(PASS_THROUGH + META_COLUMNS + SPECIAL_COLUMNS)
    for col in row.keys():
        if col and col not in known:
            print("note: ignoring CSV column", repr(col))

    name = clean_cell(row.get("Name")) or ""
    goal = clean_cell(row.get("Goal")) or ""
    file_tcl = clean_cell(row.get("File")) or "RunParallel.tcl"
    np_hint = clean_cell(row.get("Number of Procs")) or "?"

    # --- outDIR: always set for matrix rows ---
    if args.outDIR is not None:
        out_dir = args.outDIR
    else:
        now = datetime.datetime.now()
        stamp = now.strftime("%Y%m%d_%H%M")
        base = clean_cell(row.get("outDIR"))
        if base is None:
            base = "runs"
        tag = "%s_%s" % (run_out_tag(row_id), stamp)
        extra = slug(name) or slug(goal)
        if extra:
            tag = tag + "_" + extra
        # if CSV gave a base folder, append the tag
        out_dir = base.rstrip("/\\") + "/" + tag

    lines = []
    lines.append("# Overrides.tcl — generated by RunTestMatrix.py")
    lines.append("# row=%s  Name=%s  Goal=%s" % (row_id, name, goal))
    lines.append("# File=%s  np=%s (hint only — launch OpenSees yourself)" % (file_tcl, np_hint))
    lines.append("# Applied only if Run*.tcl has overridesON 1 and this path is argv.")
    lines.append("")

    # Pass-through knobs (skip blanks)
    for col in PASS_THROUGH:
        val = clean_cell(row.get(col))
        if val is None:
            continue
        if col in ("eqIntegrator", "prePartitionSystem", "postPartitionSystem", "constraintsHandler"):
            val = strip_paren_note(val)
        if col == "constraintsHandler":
            want = run_handler(row_id)
            if val != want:
                print("ERROR: Run %s must use constraintsHandler %s (got %s)" % (
                    row_id, want, val), file=sys.stderr)
                sys.exit(1)
        if col == "rayleighStiff":
            val = val.lower()
            if val not in ("committed", "initial"):
                print("ERROR: rayleighStiff must be committed|initial, got", repr(val), file=sys.stderr)
                sys.exit(1)
        if col == "expElementType":
            if val not in ("generic", "twoNodeLink"):
                print("ERROR: expElementType must be generic|twoNodeLink, got", repr(val), file=sys.stderr)
                sys.exit(1)
        if col == "holdPierON":
            try:
                val = "1" if float(val) > 0 else "0"
            except ValueError:
                print("ERROR: holdPierON must be a number, got", repr(val), file=sys.stderr)
                sys.exit(1)
        lines.append("set %s %s" % (col, tcl_number_or_string(val)))

    lines.append("set outDIR %s" % tcl_quote(out_dir))

    # soilProfile / soilMesh: strip "4 (SOFT)" style labels
    sp = strip_label(row.get("soilProfile"))
    if sp is not None:
        lines.append("set soilProfile %s" % sp)

    sm = strip_label(row.get("soilMesh"))
    if sm is not None:
        lines.append("set soilMesh %s" % sm)

    se = enum_soil_ele(row.get("soilEleType"))
    if se is not None:
        lines.append("set soilEleType %s" % tcl_quote(se))

    sc = enum_constitutive(row.get("soilConstitutive"))
    if sc is not None:
        lines.append("set soilConstitutive %s" % tcl_quote(sc))

    # GM paths: relative names join under gmRoot from Parameters.tcl
    gm_dir = clean_cell(row.get("gmDir"))
    gm_file = clean_cell(row.get("gmVelFile"))
    if gm_dir is not None:
        # absolute path? keep as-is; else join with gmRoot
        if os.path.isabs(gm_dir) or gm_dir.startswith("$"):
            lines.append("set gmDir %s" % tcl_quote(gm_dir))
        else:
            lines.append("set gmDir [file join $gmRoot %s]" % gm_dir)
    if gm_file is not None:
        if os.path.isabs(gm_file) or gm_file.startswith("$") or "/" in gm_file or "\\" in gm_file:
            lines.append("set gmVelFile %s" % tcl_quote(gm_file.replace("\\", "/")))
        else:
            lines.append("set gmVelFile [file join $gmDir %s]" % gm_file)

    lines.append("")

    text = "\n".join(lines) + "\n"
    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)

    print("Wrote", args.out)
    print("  outDIR =", out_dir)
    if file_tcl.lower().endswith("runparallel.tcl"):
        print("  hint: mpiexec -n %s OpenSeesMP %s %s" % (np_hint, file_tcl, os.path.basename(args.out)))
    else:
        print("  hint: OpenSees %s %s" % (file_tcl, os.path.basename(args.out)))


if __name__ == "__main__":
    main()
