"""Prepare Rfam family data for Phase 6 (and future local experiments).

Fetches consensus sequences and secondary structures from Rfam for
the 52 families used in Phases 1-5. Saves each as a JSON file in
data/rfam_families/.

For the 12 Phase 1 families with PDB crystal structures, uses the
curated sequences from the paper (Table 1). For the 40 Phase 2
expansion families, fetches from Rfam API.

Usage:
    uv run python scripts/prepare_rfam_families.py
"""

import json
import re
import urllib.request
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "data" / "rfam_families"

PHASE1_FAMILIES = {
    "tRNA_Phe_yeast": {
        "rfam_id": "RF00005",
        "sequence": "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
        "dot_bracket": "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....",
    },
    "tRNA_Ala_human": {
        "rfam_id": "RF00005",
        "sequence": "GGGGGAUUAGCUCAGCUGGGAGAGCGUCUGCUUAGCACGCAGAGGUCCUGGGUUCAAUCCCCAGCAUCUCCACCA",
        "dot_bracket": "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....",
    },
    "5S_rRNA_ecoli": {
        "rfam_id": "RF00001",
        "sequence": "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU",
        "dot_bracket": "(((((((...(((((...........))))).(((((...((((((((.......)))).)))).....))))).......(((((.......)))))))))))).................",
    },
    "hammerhead_ribozyme": {
        "rfam_id": "RF00163",
        "sequence": "GGAGACCGCUGCCGAAACGCGAAAGCGUCUAGCGUCCG",
        "dot_bracket": "((((...((((..........))))(((...)))..))))",
    },
    "SAM_riboswitch": {
        "rfam_id": "RF00162",
        "sequence": "GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAAC",
        "dot_bracket": "(((((((..((((..........)))).(((((.((...((((....((((....))))..........))))..)).))))).......(((((.......)))))))))))).........",
    },
    "TPP_riboswitch": {
        "rfam_id": "RF00059",
        "sequence": "GCGUUCCUAUAAUGUUGAUAUGGAUUUGAUUAACAGAAGAUCGGCGAACCUGCAUUAAAGAGAGCGGACGGAUUAAUAGUCUGAUCACGAAGUCGAGAAUUCCC",
        "dot_bracket": "..(((((((..((((.........)))).(((((.((.......)).)))))....(((((...(((((.........)))))...))))).))))))).........",
    },
    "SRP_RNA_helix8": {
        "rfam_id": "RF00017",
        "sequence": "GGGCUCAGUGGCUCGAUUUGGACUUAGAUGCUCCAACCUGUAAGCUCAGGUAGGAUGUAAACGGCUGAGCCUGGAGGCAGAAGCUGCCCC",
        "dot_bracket": "((((((.((((((((.((((((.((.....)).)))))).)))))))).))))))..(((((((((..........)))))))))........",
    },
    "mir_21_precursor": {
        "rfam_id": "RF00001",
        "sequence": "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA",
        "dot_bracket": "(((((((((((((((((((((((((((((((((..........)))))))))))))))))))))))))))))))))..",
    },
    "IRES_HCV_domainII": {
        "rfam_id": None,
        "sequence": "GCCAGCCCCCGAUUGGGGGCGACACUCCACCAUGGAUCACUCCCUGUCAGGCGUUCCGCAGUCCCCGUCCUUCCUUCUGUUAGGAC",
        "dot_bracket": "(((((.(((((((.((.((.(((((((((((((.........)))))))..........)))))))).)).))))))).))))).....  ",
    },
    "HDV_ribozyme": {
        "rfam_id": "RF00094",
        "sequence": "GGCCGGCAUGGUCCCAGCCUCCUCGCUGGCGCCGGCUGGGCAACAUUCCGAGGGGACCGUCCCCUCGGUAAUGGCGAAUGGGUCC",
        "dot_bracket": "(((((.((((((.....(((.((.((..((((..((((.......)))).)))).)).)).)))..))))))..))))).........",
    },
    "RNaseP_specificity": {
        "rfam_id": "RF00010",
        "sequence": "GAGGAAAGUCCGGGCUCCAUAGGGCAGGGUGCCAGGUAACGCCUGGGGGGGAAACCCACGACCAGUGCAACAGAGAGCAAACCGCCGAUGGCCC",
        "dot_bracket": "........((((((.((((((((..((((((........))))))...)))))))).........((((((.......)))))).))))))",
    },
    "U2_snRNA_stem": {
        "rfam_id": "RF00004",
        "sequence": "AUCCUUUGCUUUGGCUUAGAUUCAAAUUGAAAUUUGAUCAAGUGCGUAUAUCCUUAUACUAACAAUAUCGGAGAAAAGCUGAU",
        "dot_bracket": "((((((((((((..(((((.........)))))..((((.........)))).((((((......)))))).)))))))))))..",
    },
}

PHASE2_RFAM_IDS = {
    "6S_RNA": "RF00013",
    "7SK_RNA": "RF00100",
    "Bacterial_SRP": "RF00169",
    "CRISPR_leader": "RF01315",
    "Corona_5UTR": "RF03120",
    "Corona_s2m": "RF00164",
    "CrPV_IRES": "RF00458",
    "FMN_riboswitch": "RF00050",
    "Hepatitis_C_IRES_III": "RF00061",
    "Histone_3prime": "RF00032",
    "IRE_stem_loop": "RF00037",
    "SAH_riboswitch": "RF01057",
    "SECIS_element": "RF00031",
    "THF_riboswitch": "RF01831",
    "T_box_leader": "RF00230",
    "U1_snRNA": "RF00003",
    "U4_snRNA": "RF00015",
    "U5_snRNA": "RF00020",
    "U6_snRNA": "RF00026",
    "Vault_RNA": "RF00006",
    "Y_RNA": "RF00019",
    "ZMP_riboswitch": "RF01750",
    "c_di_GMP_riboswitch": "RF01051",
    "cobalamin_riboswitch": "RF00174",
    "fluoride_riboswitch": "RF01734",
    "glmS_ribozyme": "RF00234",
    "glycine_riboswitch": "RF00504",
    "group_II_intron_D5": "RF00029",
    "group_I_intron_P4P6": "RF00028",
    "hatchet_ribozyme": "RF02678",
    "lysine_riboswitch": "RF00168",
    "manganese_riboswitch": "RF01786",
    "mir_122_precursor": "RF00001",
    "mir_155_precursor": "RF00001",
    "mir_let7_precursor": "RF00001",
    "pistol_ribozyme": "RF02679",
    "preQ1_riboswitch": "RF00522",
    "purine_riboswitch": "RF00167",
    "tmRNA": "RF00023",
    "twister_ribozyme": "RF02681",
}


def fetch_rfam_consensus(rfam_id):
    """Fetch consensus sequence and structure from Rfam API."""
    url = f"https://rfam.org/family/{rfam_id}/alignment?acc={rfam_id}&format=stockholm&download=0"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode("utf-8")
    except Exception as e:
        print(f"  Failed to fetch {rfam_id}: {e}")
        return None, None

    sequence = None
    structure = None
    for line in text.split("\n"):
        if line.startswith("#=GC SS_cons"):
            structure = line.split(None, 2)[2].strip()
            structure = re.sub(r"[<{]", "(", structure)
            structure = re.sub(r"[>}]", ")", structure)
            structure = re.sub(r"[^().]", ".", structure)
        elif line.startswith("#=GC RF"):
            sequence = line.split(None, 2)[2].strip().upper()
            sequence = sequence.replace("T", "U")
            sequence = re.sub(r"[^ACGU]", "", sequence)

    return sequence, structure


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    saved = 0
    failed = []

    for name, data in PHASE1_FAMILIES.items():
        seq = data["sequence"].strip()
        db = data["dot_bracket"].strip()
        if len(seq) != len(db):
            db = db[:len(seq)]

        family = {
            "name": name,
            "sequence": seq,
            "dot_bracket": db,
            "source": "Phase 1 curated (PDB/Rfam)",
            "rfam_id": data.get("rfam_id"),
            "length": len(seq),
        }
        out_path = OUT_DIR / f"{name}.json"
        with open(out_path, "w") as f:
            json.dump(family, f, indent=2)
        saved += 1
        print(f"  {name}: {len(seq)} nt (Phase 1 curated)")

    for name, rfam_id in PHASE2_RFAM_IDS.items():
        print(f"  Fetching {name} ({rfam_id})...", end=" ")
        seq, db = fetch_rfam_consensus(rfam_id)
        if seq is None or db is None:
            print("FAILED")
            failed.append(name)
            continue

        if len(seq) != len(db):
            min_len = min(len(seq), len(db))
            seq = seq[:min_len]
            db = db[:min_len]

        family = {
            "name": name,
            "sequence": seq,
            "dot_bracket": db,
            "source": f"Rfam {rfam_id} consensus",
            "rfam_id": rfam_id,
            "length": len(seq),
        }
        out_path = OUT_DIR / f"{name}.json"
        with open(out_path, "w") as f:
            json.dump(family, f, indent=2)
        saved += 1
        print(f"{len(seq)} nt")

    print(f"\nSaved {saved} families to {OUT_DIR}")
    if failed:
        print(f"Failed: {failed}")
    print(f"Timestamp: {timestamp}")


if __name__ == "__main__":
    main()
