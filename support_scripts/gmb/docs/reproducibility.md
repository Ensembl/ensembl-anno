# Reproducing a GMB run

Every run records what produced it. Given the manifest, the inputs and the software version,
the run can be repeated.

---

## What is recorded, and where

| artifact | written by | contains |
|---|---|---|
| `build/run_manifest.json` / `.tsv` | `gmb-build` | **the run**: versions, command, inputs + hashes, resolved roles/weights/policy, runtime, peak RSS, QC outcome |
| `build/resolved_config.yaml` | `gmb-build` | every setting after layering |
| `build/resolved_config_sha256` | `gmb-build` | its hash |
| `finalise/handover_manifest.json` / `.tsv` | `gmb-finalise` | **the outputs**: every handover file with size and SHA-256 |
| `preflight/preflight_report.json` | `gmb-preflight` | the input verdict and every measurement behind it |
| `build/gmb.log` | `gmb-build` | full stage-by-stage log |

## Run manifest fields

```
stage, started_at, finished_at, runtime_seconds, peak_rss_kb
hostname, platform, cpu_count, python_version, gmb_version
git.commit, git.branch, git.dirty, git.dirty_files
external_tools.{diamond, psauron, diamond_db}
command_line, preset, config_overlays
resolved_config_path, resolved_config_sha256
inputs[].{label, path, sha256, size_bytes, resolved_role, resolved_weight}
evidence_roles{}, evidence_weights{}
resolved_policy{}
outputs{}.{size_bytes, sha256}
qc.{fasta_qc_pass, failed_checks}
```

### Three fields that earn their place

**`git.dirty`** — a commit id does not identify what ran if the tree was dirty. When
`dirty: true`, the commit alone is not enough to reproduce the run; `git.dirty_files` lists
what differed.

**`inputs[].sha256`** — filenames get reused. The hash is what proves two runs saw the same
bytes. Reproducing a run means matching hashes, not paths.

**`resolved_policy`** — some policies are decided **at run time from the evidence**. The
`backbone_intron_rescue` applicability gate is the clearest case: the config says `auto`, and
whether it actually fired depends on the measured backbone quality. The manifest records the
configured mode, the resolved outcome, the reason, **and the numbers behind it**:

```json
"backbone_intron_rescue": {
  "mode": "auto",
  "enabled": true,
  "reason": "auto: backbone is 17.3% multi-exon (<= 55%) while assembled transcripts are 98.7% multi-exon (5.71x more) at 100.0% canonical splice sites ...",
  "evidence": {
    "backbone_multi_exon_fraction": 0.173,
    "assembled_multi_exon_fraction": 0.987,
    "assembled_to_backbone_ratio": 5.71,
    "assembled_canonical_splice_fraction": 1.0
  }
}
```

Recording the mode without the outcome would hide what actually happened.

---

## Reproducing a run

```bash
# 1. same software
git clone git@github.com:Ensembl/ensembl-anno.git && cd ensembl-anno
git checkout <git.commit from the manifest>
git status --short                     # must be empty if git.dirty was false
pip install -e support_scripts/gmb

# 2. same inputs — verify by hash, not by path
sha256sum /path/to/genome.fa           # must match inputs[].sha256

# 3. same command
#    manifest.command_line is the verbatim invocation

# 4. verify the configuration resolved identically
sha256sum "$OUT/build/resolved_config.yaml"   # must match resolved_config_sha256
```

If `resolved_config_sha256` matches and the input hashes match, the two runs saw identical
configuration and identical bytes.

## Determinism

GMB is deterministic given identical inputs, configuration and software version. It uses no
randomness in the selection path. `--sample-loci` takes `--seed` and is a testing aid, not part
of a production run.

What is **not** guaranteed identical across versions: generated gene and transcript IDs are
positional, so an upstream change that alters gene count shifts subsequent IDs. **Compare
structures, not IDs**, when checking two builds for equivalence.

## Comparing two builds

```bash
# structure-level equality, independent of IDs: the set of (chrom, strand, CDS-interval) keys
python - <<'PY'
import collections
def structures(path):
    mrna, cds = {}, collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] in ("mRNA", "transcript"):
            mrna[a["ID"]] = (f[0], f[6])
        elif f[2] == "CDS":
            for t in a.get("Parent", "").split(","):
                if t: cds[t].append((int(f[3]), int(f[4])))
    return {(*mrna[t], tuple(sorted(v))) for t, v in cds.items() if t in mrna}

a, b = structures("run_a/finalise/consensus.gff3"), structures("run_b/finalise/consensus.gff3")
print(f"identical {len(a & b)}  only_a {len(a - b)}  only_b {len(b - a)}")
PY
```

## Auditing a handover you were given

```bash
python - <<'PY'
import hashlib, json, os
m = json.load(open("finalise/handover_manifest.json"))
for name, meta in m["outputs"].items():
    path = os.path.join("finalise", name)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    print(f"{'OK ' if h == meta['sha256'] else 'BAD'} {name}")
PY
```
