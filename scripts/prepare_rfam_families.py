"""Prepare Rfam family data for Phase 6 (and future local experiments).

Fetches seed alignment members and secondary structures from Rfam for
the 52 families used in Phases 1-5. Saves each as a JSON file in
data/rfam_families/.

For the 12 Phase 1 families with PDB crystal structures, uses the
curated sequences from the paper (Table 1). For the 40 Phase 2
expansion families, fetches the best seed alignment member (highest
canonical WC pair count) from the Rfam seed Stockholm alignment.

Usage:
    uv run python scripts/prepare_rfam_families.py
    uv run python scripts/prepare_rfam_families.py --dry-run
"""

import argparse
import json
import re
import ssl
import urllib.request
from collections import OrderedDict
from pathlib import Path
from datetime import datetime

from tqdm import tqdm

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


WC_PAIRS = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")}


def parse_dot_bracket(db_string):
    stack = []
    pairs = []
    for idx, char in enumerate(db_string):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if stack:
                partner = stack.pop()
                pairs.append((partner, idx))
    return sorted(pairs)


def count_wc_pairs(sequence, db_string):
    pairs = parse_dot_bracket(db_string)
    return sum(1 for i, j in pairs if (sequence[i], sequence[j]) in WC_PAIRS)


def fetch_seed_alignment(rfam_id):
    """Fetch Rfam seed Stockholm alignment and return parsed members.

    Returns list of (seqname, ungapped_sequence, projected_dot_bracket).
    """
    url = (
        f"https://rfam.org/family/{rfam_id}/alignment"
        f"?acc={rfam_id}&format=stockholm&download=0&type=seed"
    )
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        with urllib.request.urlopen(url, timeout=60, context=ctx) as resp:
            text = resp.read().decode("utf-8")
    except Exception as e:
        print(f"  Failed to fetch {rfam_id}: {e}")
        return []

    seq_chunks = OrderedDict()
    ss_chunks = []

    for line in text.split("\n"):
        line = line.rstrip()
        if not line or line.startswith("#=GF") or line.startswith("#=GR") or line == "//":
            continue
        if line.startswith("#=GC SS_cons"):
            chunk = line.split(None, 2)[2]
            ss_chunks.append(chunk)
        elif line.startswith("#"):
            continue
        else:
            parts = line.split()
            if len(parts) == 2:
                name, chunk = parts
                if name not in seq_chunks:
                    seq_chunks[name] = []
                seq_chunks[name].append(chunk)

    if not ss_chunks:
        return []

    full_ss = "".join(ss_chunks)
    full_ss = re.sub(r"[<{\[]", "(", full_ss)
    full_ss = re.sub(r"[>}\]]", ")", full_ss)
    full_ss = re.sub(r"[^().]", ".", full_ss)

    members = []
    for seqname, chunks in seq_chunks.items():
        aligned_seq = "".join(chunks)
        if len(aligned_seq) != len(full_ss):
            continue

        ungapped_seq = []
        projected_ss = []
        for col_idx, (s_char, ss_char) in enumerate(zip(aligned_seq, full_ss)):
            if s_char in (".", "-", "_"):
                continue
            nuc = s_char.upper().replace("T", "U")
            if nuc not in "ACGU":
                nuc = "N"
            ungapped_seq.append(nuc)
            projected_ss.append(ss_char)

        seq_str = "".join(ungapped_seq)
        ss_str = "".join(projected_ss)

        seq_str = seq_str.replace("N", "")
        if len(seq_str) != len(ss_str):
            clean_seq = []
            clean_ss = []
            for s, d in zip(ungapped_seq, projected_ss):
                if s != "N":
                    clean_seq.append(s)
                    clean_ss.append(d)
            seq_str = "".join(clean_seq)
            ss_str = "".join(clean_ss)

        if len(seq_str) < 30:
            continue

        stack = []
        balanced = True
        for ch in ss_str:
            if ch == "(":
                stack.append(ch)
            elif ch == ")":
                if not stack:
                    balanced = False
                    break
                stack.pop()
        if stack:
            balanced = False
        if not balanced:
            continue

        members.append((seqname, seq_str, ss_str))

    return members


def pick_best_member(members):
    """Pick the seed member with the most canonical WC pairs."""
    best = None
    best_wc = -1
    best_len = -1
    for seqname, seq, ss in members:
        wc = count_wc_pairs(seq, ss)
        if wc > best_wc or (wc == best_wc and len(seq) > best_len):
            best = (seqname, seq, ss)
            best_wc = wc
            best_len = len(seq)
    return best, best_wc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="Report pair counts without overwriting files")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    saved = 0
    failed = []
    results = []

    print("Phase 1 curated families:")
    for name, data in PHASE1_FAMILIES.items():
        seq = data["sequence"].strip()
        db = data["dot_bracket"].strip()
        if len(seq) != len(db):
            db = db[:len(seq)]
        wc = count_wc_pairs(seq, db)

        if not args.dry_run:
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

        status = "PASS" if wc >= 15 else "FAIL"
        print(f"  {name}: {len(seq)} nt, {wc} WC pairs [{status}]")
        results.append((name, len(seq), wc, status))

    print(f"\nPhase 2 expansion families (fetching seed alignments):")
    for name, rfam_id in tqdm(list(PHASE2_RFAM_IDS.items()), desc="Fetching"):
        members = fetch_seed_alignment(rfam_id)
        if not members:
            print(f"  {name} ({rfam_id}): NO VALID MEMBERS")
            failed.append(name)
            results.append((name, 0, 0, "FAILED"))
            continue

        best, wc = pick_best_member(members)
        if best is None:
            print(f"  {name} ({rfam_id}): NO VALID MEMBERS")
            failed.append(name)
            results.append((name, 0, 0, "FAILED"))
            continue

        seqname, seq, ss = best
        status = "PASS" if wc >= 15 else "FAIL"
        print(f"  {name}: {len(seq)} nt, {wc} WC pairs, "
              f"from {seqname} ({len(members)} members) [{status}]")
        results.append((name, len(seq), wc, status))

        if not args.dry_run:
            family = {
                "name": name,
                "sequence": seq,
                "dot_bracket": ss,
                "source": f"Rfam {rfam_id} seed alignment member ({seqname})",
                "rfam_id": rfam_id,
                "length": len(seq),
                "seed_member": seqname,
                "seed_members_available": len(members),
            }
            out_path = OUT_DIR / f"{name}.json"
            with open(out_path, "w") as f:
                json.dump(family, f, indent=2)
            saved += 1

    pass_count = sum(1 for _, _, wc, s in results if s == "PASS")
    fail_count = sum(1 for _, _, wc, s in results if s == "FAIL")
    fetch_fail = sum(1 for _, _, wc, s in results if s == "FAILED")

    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Pass (>= 15 WC pairs): {pass_count}")
    print(f"  Fail (< 15 WC pairs):  {fail_count}")
    print(f"  Fetch failed:          {fetch_fail}")
    print(f"  Total:                 {len(results)}")
    if not args.dry_run:
        print(f"  Saved {saved} families to {OUT_DIR}")
    else:
        print(f"  (dry run — no files written)")
    if failed:
        print(f"  Fetch failures: {failed}")
    print(f"  Timestamp: {timestamp}")


if __name__ == "__main__":
    main()
