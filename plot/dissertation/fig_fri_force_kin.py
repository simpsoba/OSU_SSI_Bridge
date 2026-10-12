"""Friday F05 version of fig_wed_force_kin.py: interface force against measured displacement, velocity, and
acceleration over the whole record, colored by phase (earthquake, free vibration, solitary wave); read-only on data."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from fig_wed_force_kin import make  # noqa: E402

if __name__ == "__main__":
    make("F05", Path(sys.argv[1]))
