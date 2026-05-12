"""Fix v7 l3 pseudobulk *_rows.csv / *_cols.csv to add headers.

Reads each headerless csv, prepends "gene" or "donor", writes back.
Idempotent: skips files that already have the header.
"""
from pathlib import Path
import sys

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
ROOTS = sorted(DATA_DIR.glob("pseudobulk_v7_l3*"))
print(f"Scanning {len(ROOTS)} v7 root(s):")
for r in ROOTS:
    print(f"  {r}")

n_fixed = 0
n_skipped = 0
for ROOT in ROOTS:
 for csv in ROOT.rglob("*_rows.csv"):
    txt = csv.read_text()
    first_line = txt.split("\n", 1)[0].strip()
    if first_line == '"gene"' or first_line == "gene":
        n_skipped += 1
        continue
    # Prepend "gene" header
    new_txt = '"gene"\n' + "\n".join(f'"{x.strip()}"' for x in txt.splitlines() if x.strip()) + "\n"
    csv.write_text(new_txt)
    n_fixed += 1

 for csv in ROOT.rglob("*_cols.csv"):
    txt = csv.read_text()
    first_line = txt.split("\n", 1)[0].strip()
    if first_line == '"donor"' or first_line == "donor":
        n_skipped += 1
        continue
    new_txt = '"donor"\n' + "\n".join(f'"{x.strip()}"' for x in txt.splitlines() if x.strip()) + "\n"
    csv.write_text(new_txt)
    n_fixed += 1

print(f"Fixed {n_fixed} files, skipped {n_skipped} (already had header).")
