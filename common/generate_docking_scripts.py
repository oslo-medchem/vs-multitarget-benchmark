#!/usr/bin/env python3
"""
Generate Uni-Dock docking shell scripts for both naive and skill-guided protocols.

Usage:
    python generate_docking_scripts.py \
        --target-dir targets/egfr/ \
        --center-x 10.5 --center-y 20.3 --center-z 30.1

Outputs:
    target_dir/04_docking/naive/run_unidock.sh
    target_dir/04_docking/skill/run_unidock.sh
"""
import argparse
import os
import stat
from pathlib import Path

TEMPLATE = r'''#!/usr/bin/env bash
# {Protocol} protocol: Uni-Dock docking
# Box {size}A, exhaustiveness={exh}, Vina scoring
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
TARGET_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

RECEPTOR="$SCRIPT_DIR/receptor_{protocol}.pdbqt"
PDBQT_DIR="$TARGET_DIR/03_library_preparation/{protocol}/pdbqt"
RESULTS_DIR="$SCRIPT_DIR/results"
LIGAND_LIST="$SCRIPT_DIR/ligand_list_{protocol}.txt"

mkdir -p "$RESULTS_DIR"

CENTER_X={cx}
CENTER_Y={cy}
CENTER_Z={cz}
SIZE={size}
EXHAUSTIVENESS={exh}
NUM_MODES={modes}

echo "Building ligand list..."
find "$PDBQT_DIR" -name "*.pdbqt" -size +50c | sort > "$LIGAND_LIST"
N=$(wc -l < "$LIGAND_LIST")
echo "Ligands to dock: $N"

echo "Starting Uni-Dock ({protocol} protocol)..."
START_TIME=$(date +%s)

unidock \
  --receptor "$RECEPTOR" \
  --ligand_index "$LIGAND_LIST" \
  --center_x $CENTER_X \
  --center_y $CENTER_Y \
  --center_z $CENTER_Z \
  --size_x $SIZE \
  --size_y $SIZE \
  --size_z $SIZE \
  --exhaustiveness $EXHAUSTIVENESS \
  --num_modes $NUM_MODES \
  --dir "$RESULTS_DIR" \
  --scoring vina \
  2>&1 | tee "$SCRIPT_DIR/unidock_{protocol}.log"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Elapsed: ${{ELAPSED}}s" | tee -a "$SCRIPT_DIR/unidock_{protocol}.log"
echo "{Protocol} docking complete."
'''


def main():
    parser = argparse.ArgumentParser(
        description="Generate Uni-Dock docking scripts for naive and skill protocols."
    )
    parser.add_argument("--target-dir", required=True,
                        help="Path to target directory (e.g. targets/egfr/)")
    parser.add_argument("--center-x", required=True, type=float,
                        help="Box center X coordinate")
    parser.add_argument("--center-y", required=True, type=float,
                        help="Box center Y coordinate")
    parser.add_argument("--center-z", required=True, type=float,
                        help="Box center Z coordinate")
    parser.add_argument("--naive-box-size", type=int, default=20,
                        help="Box size for naive protocol (default: 20)")
    parser.add_argument("--skill-box-size", type=int, default=25,
                        help="Box size for skill protocol (default: 25)")
    parser.add_argument("--naive-exhaustiveness", type=int, default=8,
                        help="Exhaustiveness for naive protocol (default: 8)")
    parser.add_argument("--skill-exhaustiveness", type=int, default=32,
                        help="Exhaustiveness for skill protocol (default: 32)")
    parser.add_argument("--num-modes", type=int, default=9,
                        help="Number of binding modes (default: 9)")
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve()
    docking_dir = target_dir / "04_docking"

    protocols = {
        "naive": {
            "protocol": "naive",
            "Protocol": "Naive",
            "size": args.naive_box_size,
            "exh": args.naive_exhaustiveness,
        },
        "skill": {
            "protocol": "skill",
            "Protocol": "Skill",
            "size": args.skill_box_size,
            "exh": args.skill_exhaustiveness,
        },
    }

    for pname, params in protocols.items():
        out_dir = docking_dir / pname
        out_dir.mkdir(parents=True, exist_ok=True)

        script_content = TEMPLATE.format(
            protocol=params["protocol"],
            Protocol=params["Protocol"],
            cx=args.center_x,
            cy=args.center_y,
            cz=args.center_z,
            size=params["size"],
            exh=params["exh"],
            modes=args.num_modes,
        )

        script_path = out_dir / "run_unidock.sh"
        script_path.write_text(script_content)
        # Make executable
        st = os.stat(script_path)
        os.chmod(script_path, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

        print(f"Generated: {script_path}")
        print(f"  Box center: ({args.center_x}, {args.center_y}, {args.center_z})")
        print(f"  Box size: {params['size']}A, exhaustiveness: {params['exh']}, modes: {args.num_modes}")

    print("\nDone. Run with:")
    print(f"  bash {docking_dir}/naive/run_unidock.sh")
    print(f"  bash {docking_dir}/skill/run_unidock.sh")


if __name__ == "__main__":
    main()
