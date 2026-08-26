"""Run inside the Modal container, once per branch, as a subprocess.

Imports must happen fresh per branch, so this cannot be a function call in the
parent process -- Python caches modules and the editable checkout changes under
it between branches.

Prints one JSON object on stdout.
"""

import argparse
import importlib
import json
import pickle
import traceback

WT = "GGCUAGCUAAGGCUAGCC"
MUT = "GGCUAGCUAUGGCUAGCC"  # single substitution at position 10

MODELS = [
    ("multimolecule/rnabert", "RnaBertModel"),
    ("multimolecule/rnafm", "RnaFmModel"),
    ("multimolecule/splicebert", "SpliceBertModel"),
]


def check(repo: str, cls_name: str) -> dict:
    import torch
    from multimolecule import RnaTokenizer
    from nnsight import NNsight

    entry: dict = {}
    multimolecule = importlib.import_module("multimolecule")
    hf = getattr(multimolecule, cls_name).from_pretrained(
        repo, attn_implementation="eager"
    )
    hf.eval()
    tok = RnaTokenizer.from_pretrained(repo)
    model = NNsight(hf)

    layers = model.encoder.layer
    cls = type(layers[0])
    entry["envoy_class"] = cls.__name__
    entry["dot_in_name"] = "." in cls.__name__
    module = importlib.import_module(cls.__module__)
    entry["resolves_back"] = getattr(module, cls.__name__, None) is cls
    entry["siblings_share_class"] = type(layers[0]) is type(layers[1])
    entry["n_distinct_layer_classes"] = len({type(l) for l in layers})

    try:
        pickle.dumps(cls)
        entry["class_picklable"] = True
    except Exception as exc:
        entry["class_picklable"] = False
        entry["pickle_error"] = f"{type(exc).__name__}: {str(exc)[:200]}"

    def hidden_states(seq):
        ids = tok(seq, return_tensors="pt")["input_ids"]
        with torch.no_grad():
            return hf(ids, output_hidden_states=True).hidden_states

    hs_wt = hidden_states(WT)
    hs_mut = hidden_states(MUT)
    shift = (hs_wt[-1] - hs_mut[-1]).norm(dim=-1)[0]
    entry["n_layers"] = len(hs_wt)
    entry["layer_shape"] = list(hs_wt[0].shape)
    entry["per_position_l2_shift"] = [round(float(x), 4) for x in shift]
    entry["argmax_position"] = int(shift.argmax())

    with model.trace(tok(WT, return_tensors="pt")["input_ids"]):
        mid = layers[len(layers) // 2].nns_output[0].save()
    entry["traced_shape"] = list(mid.shape)
    return entry


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    import nnsight

    results = {"nnsight_version": nnsight.__version__, "models": {}}
    for repo, cls_name in MODELS:
        try:
            entry = check(repo, cls_name)
            entry["status"] = "ok"
        except Exception as exc:
            entry = {
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        results["models"][repo] = entry

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(json.dumps(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
