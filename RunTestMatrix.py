"""
Goals
-----
Write Overrides.tcl from one TestMatrix.csv row (keyed by Test).
Leave blank cells at their Parameters.tcl or Run*.tcl defaults.

  python RunTestMatrix.py --test F07
  python RunTestMatrix.py --test S04F07
  python RunTestMatrix.py --test S04F07a
  python RunTestMatrix.py --test S04F22p

``Test`` is unique (W## / F## / Fd## / Fx## / S04F## / S04F##a / S04F22p / …).
``constraintsHandler`` comes from the CSV cell (no sign-derived rule).

To add another override:
  1. Add a TestMatrix.csv column, preferably with the Tcl variable name.
  2. Add the name to PASS_THROUGH or a special-case block below.
  3. Run this script for the desired Test.

soilProfile and soilMesh may be bare numbers or labeled values such as
``0 (BASELINE)``. Rayleigh columns map to analysis/RayleighDamping.tcl.
Lab-log metadata (DumpFolder, MatFile, Note, …) is not written as set lines.

This script writes Overrides.tcl; it does not launch OpenSees. Example launches:
  mpiexec -n 8 OpenSeesMP RunParallel.tcl Overrides.tcl
  OpenSees Run.tcl Overrides.tcl
"""

import argparse
import csv
import datetime
import os
import re
import sys

# ------------------------------------------------------------
# 1. PATHS AND CSV COLUMN GROUPS
# ------------------------------------------------------------

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
META_COLUMNS = [
    "Test",
    "File",
    "Number of Procs",
    "Name",
    "Goal",
    "DOFs",
    "DateTime",
    "DumpFolder",
    "MatFile",
    "LabTrial",
    "Note",
]

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


# ------------------------------------------------------------
# 2. CSV VALUE NORMALIZATION
# ------------------------------------------------------------


def clean_cell(raw):
    """
    Convert an empty or missing CSV cell to None.

    Args:    raw
    Returns: stripped string, or None
    """
    if raw is None:
        return None
    s = str(raw).strip()
    if s == "":
        return None
    return s


def strip_label(raw):
    """
    Remove a parenthetical label from a numeric matrix value.

    Args:    raw  e.g. "1", "1 (FINE)", or "4 (SOFT)"
    Returns: leading integer string, original value, or None
    """
    s = clean_cell(raw)
    if s is None:
        return None
    m = re.match(r"^([+-]?\d+)\s*\(", s)
    if m:
        return m.group(1)
    # bare number (or anything else) — pass through as-is
    return s


def strip_paren_note(raw):
    """
    Remove a trailing parenthetical note.

    Args:    raw  e.g. "UmfPack (or ProfileSPD)"
    Returns: base value, or None
    """
    s = clean_cell(raw)
    if s is None:
        return None
    if " (" in s:
        s = s.split(" (", 1)[0].strip()
    return s


def enum_soil_ele(raw):
    """
    Normalize accepted soil-element spellings.

    Args:    raw
    Returns: "SSPquad", "quad", original value, or None
    """
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
    """
    Normalize the soil constitutive name to lowercase.

    Args:    raw
    Returns: lowercase value, or None
    """
    s = clean_cell(raw)
    if s is None:
        return None
    return s.lower()


def tcl_quote(s):
    """
    Quote a string for a Tcl set command and normalize slashes.

    Args:    s
    Returns: Tcl double-quoted string
    """
    return '"' + s.replace("\\", "/").replace('"', '\\"') + '"'


def tcl_number_or_string(s):
    """
    Keep numeric values bare and quote all other Tcl values.

    Args:    s
    Returns: Tcl value token
    """
    if re.match(r"^[+-]?\d+(\.\d+)?([eE][+-]?\d+)?$", s):
        return s
    return tcl_quote(s)


def slug(text, max_len=24):
    """
    Build a short filesystem-safe token from Name or Goal.

    Args:    text, max_len
    Returns: slug string
    """
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


# ------------------------------------------------------------
# 3. TEST-ID CONVENTIONS
# ------------------------------------------------------------


def canon_test_id(raw):
    """
    Canonicalize a Test ID (case + digit padding).

    Args:    raw  e.g. "f7", "F07", "s04f7a", "W01", "Fd01"
    Returns: canonical Test string, or None
    """
    s = clean_cell(raw)
    if s is None:
        return None

    m = re.fullmatch(r"(W|F|Fd|Fx)(\d+)", s, flags=re.IGNORECASE)
    if m:
        prefix = {"w": "W", "f": "F", "fd": "Fd", "fx": "Fx"}[m.group(1).lower()]
        return "%s%02d" % (prefix, int(m.group(2)))

    # S04F07 / S04F07a / S04F22p / S04F22ma (month-day dry + optional letter suffix)
    m = re.fullmatch(r"(S)(\d{2})(F)(\d+)([A-Za-z]*)", s, flags=re.IGNORECASE)
    if m:
        return "S%sF%02d%s" % (m.group(2), int(m.group(4)), m.group(5).lower())

    # Unknown shape: keep stripped text with letters as typed length
    return s


def find_row_by_test(rows, test_id):
    """
    Find one TestMatrix row by canonical Test ID.

    Args:    rows, test_id
    Returns: matching row dictionary, or None
    """
    want = canon_test_id(test_id)
    if want is None:
        return None
    for r in rows:
        got = canon_test_id(r.get("Test"))
        if got is not None and got.lower() == want.lower():
            return r
    return None


def test_out_tag(test_id):
    """
    Format a Test ID for an output-directory name.

    Args:    test_id  canonical Test
    Returns: filesystem-safe tag (same as Test for current IDs)
    """
    return re.sub(r"[^\w.+-]+", "_", test_id)


# ------------------------------------------------------------
# 4. READ MATRIX AND WRITE OVERRIDES
# ------------------------------------------------------------


def main():
    """
    Parse CLI options and write Overrides.tcl for one matrix row.

    Args:    command-line arguments from argparse
    Returns: none (exits on invalid input; writes Overrides.tcl)
    """
    parser = argparse.ArgumentParser(
        description="Write Overrides.tcl from one TestMatrix.csv row (Test ID)."
    )
    parser.add_argument(
        "--test",
        type=str,
        required=True,
        help="Test id from TestMatrix.csv (F07, S04F07, S04F07a, W01, …).",
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
        help="Force exact outDIR (skips Test/timestamp naming)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        print("ERROR: cannot find CSV:", args.csv, file=sys.stderr)
        sys.exit(1)

    with open(args.csv, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    if rows and "Test" not in rows[0]:
        print(
            "ERROR: %s has no Test column (expected lab-style matrix)." % args.csv,
            file=sys.stderr,
        )
        sys.exit(1)

    test_id = canon_test_id(args.test)
    if test_id is None:
        print("ERROR: bad Test id", repr(args.test), file=sys.stderr)
        sys.exit(1)
    row = find_row_by_test(rows, test_id)
    if row is None:
        print("ERROR: no Test ==", test_id, "in", args.csv, file=sys.stderr)
        sys.exit(1)
    # Prefer the CSV spelling after canon match
    test_id = canon_test_id(row.get("Test")) or test_id

    # Warn about unknown columns (ignored on purpose, e.g. wave metadata)
    known = set(PASS_THROUGH + META_COLUMNS + SPECIAL_COLUMNS)
    for col in row.keys():
        if col and col not in known:
            print("note: ignoring CSV column", repr(col))

    name = clean_cell(row.get("Name")) or ""
    goal = clean_cell(row.get("Goal")) or ""
    file_tcl = clean_cell(row.get("File")) or "RunParallel.tcl"
    np_hint = clean_cell(row.get("Number of Procs")) or "?"
    note = clean_cell(row.get("Note")) or ""

    # --- outDIR: always set for matrix rows ---
    if args.outDIR is not None:
        out_dir = args.outDIR
    else:
        now = datetime.datetime.now()
        stamp = now.strftime("%Y%m%d_%H%M")
        base = clean_cell(row.get("outDIR"))
        if base is None:
            base = "runs"
        tag = "%s_%s" % (test_out_tag(test_id), stamp)
        extra = slug(name) or slug(goal)
        if extra:
            tag = tag + "_" + extra
        out_dir = base.rstrip("/\\") + "/" + tag

    lines = []
    lines.append("# Overrides.tcl — generated by RunTestMatrix.py")
    lines.append("# Test=%s  Name=%s  Goal=%s" % (test_id, name, goal))
    if note:
        lines.append("# Note=%s" % note)
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
            if val not in ("Auto", "Transformation"):
                print(
                    "ERROR: constraintsHandler must be Auto|Transformation, got",
                    repr(val),
                    file=sys.stderr,
                )
                sys.exit(1)
        if col == "realTimeON":
            try:
                val = "1" if float(val) > 0 else "0"
            except ValueError:
                print("ERROR: realTimeON must be a number, got", repr(val), file=sys.stderr)
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
    print("  Test =", test_id)
    print("  outDIR =", out_dir)
    if file_tcl.lower().endswith("runparallel.tcl"):
        print("  hint: mpiexec -n %s OpenSeesMP %s %s" % (np_hint, file_tcl, os.path.basename(args.out)))
    else:
        print("  hint: OpenSees %s %s" % (file_tcl, os.path.basename(args.out)))


if __name__ == "__main__":
    main()
