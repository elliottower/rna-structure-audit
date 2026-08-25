"""Rebuild the HTT CAG-repeat fragments and fold them, as panel-format records.

The deposited HTT case study recorded repeat counts, sequence lengths and MFE
values but neither the sequences nor the dot-brackets, so a re-run has to
reconstruct both. The construction is recoverable from the stored lengths: 92
nucleotides of real HTT mRNA on each side of the CAG tract, with the tract
replaced by `CAG` repeated n times. That reproduces all five deposited lengths
exactly -- 235, 247, 292, 304 and 364 for n = 17, 21, 36, 40 and 60.

Folding happens here rather than inside a GPU image, so the structures become
versioned data that is hashed with everything else, and no model image needs
ViennaRNA.

    uv run --no-project --with "ViennaRNA==2.7.0" --python 3.12 \\
        python scripts/build_htt_fragments.py
"""

import json
import re
from pathlib import Path

import RNA

REPO = Path(__file__).resolve().parents[1]
FASTA = REPO / "data" / "HTT_NM_002111.7.fasta"
OUT = REPO / "data" / "htt_fragments"

# The deposited run recorded lengths and MFEs but not its construction. Total
# flank is fixed at 184 by the lengths; the split was recovered by scanning every
# split for one reproducing all five deposited MFEs. Six do -- left 28 through 33 --
# because the extra 5' nucleotides do not pair and cannot change the energy. Of
# those, only left = 30 begins on a start codon (ATGAAGGCCTTC at position 337),
# so it is the one used. A symmetric 92/92 split matches all five lengths and
# none of the five energies.
LEFT_FLANK = 30
RIGHT_FLANK = 154
REPEAT_COUNTS = [17, 21, 36, 40, 60]
# What the deposited run recorded, as a check on the reconstruction. Lengths fix
# the total flank; energies fix the split.
DEPOSITED_LENGTHS = {17: 235, 21: 247, 36: 292, 40: 304, 60: 364}
DEPOSITED_MFE = {17: -102.3, 21: -105.5, 36: -124.0, 40: -129.0, 60: -154.0}
CLINICAL = {17: "normal", 21: "normal", 36: "pathogenic", 40: "pathogenic",
            60: "pathogenic"}


def main():
    sequence = "".join(line.strip() for line in FASTA.read_text().splitlines()
                       if not line.startswith(">"))
    tract = re.search(r"(CAG){5,}", sequence)
    if tract is None:
        raise ValueError(f"no CAG tract in {FASTA.name}")
    left = sequence[tract.start() - LEFT_FLANK:tract.start()]
    right = sequence[tract.end():tract.end() + RIGHT_FLANK]
    if len(left) != LEFT_FLANK or len(right) != RIGHT_FLANK:
        raise ValueError("the tract sits closer to an end than the flanks allow")

    OUT.mkdir(parents=True, exist_ok=True)
    for n in REPEAT_COUNTS:
        fragment = (left + "CAG" * n + right).replace("T", "U")
        if len(fragment) != DEPOSITED_LENGTHS[n]:
            raise ValueError(
                f"CAG{n} rebuilt to {len(fragment)} nt, deposited run recorded "
                f"{DEPOSITED_LENGTHS[n]}; the construction is not the one used")
        structure, mfe = RNA.fold(fragment)
        if abs(mfe - DEPOSITED_MFE[n]) > 0.05:
            raise ValueError(
                f"CAG{n} folds to {mfe:.2f}, deposited run recorded "
                f"{DEPOSITED_MFE[n]}; the construction is not the one used")
        name = f"HTT_CAG{n}_{CLINICAL[n]}"
        (OUT / f"{name}.json").write_text(json.dumps({
            "name": name,
            "sequence": fragment,
            "dot_bracket": structure,
            "n_cag": n,
            "clinical_label": CLINICAL[n],
            "mfe": round(mfe, 2),
            "source": f"NM_002111.7, {LEFT_FLANK} nt upstream and "
                      f"{RIGHT_FLANK} nt downstream of the CAG tract, folded "
                      "with ViennaRNA 2.7.0",
        }, indent=2) + "\n")
        print(f"  {name:28s} {len(fragment):4d} nt  MFE {mfe:8.2f}  "
              f"{structure.count('(')} pairs")

    print(f"\n  wrote {len(REPEAT_COUNTS)} fragments to {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
