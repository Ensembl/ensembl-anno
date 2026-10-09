"""``gmb-compare`` has moved out of GMB.

Reference-annotation comparison is evaluation tooling, not part of building
gene models, and now lives in the ``ensembl-genes`` repository as
``annotation-qc pairwise-compare``. This stub only tells a caller where it went.
"""

from __future__ import annotations

import sys

MOVED_MESSAGE = """\
gmb-compare has been removed from Gene Model Builder.

Reference comparison now lives in the ensembl-genes repository
(src/python/ensembl/genes/annotation_qc, Python >= 3.12):

    annotation-qc pairwise-compare \\
        --query     <gmb_out>/finalise/consensus.gff3 \\
        --reference <reference.gff3> \\
        --genome    <genome.fa> \\
        --evaluation-mode protein_coding \\
        --reference-transcript-biotypes protein_coding \\
        --evidence-attribution <gmb_out>/build/evidence_attribution.tsv \\
        --outdir    <comparison_dir>

It writes the same comparison_summary.{json,tsv} and comparison_details.tsv
files gmb-compare did. See support_scripts/gmb/docs/qc.md.
"""


def main(argv=None) -> int:
    sys.stderr.write(MOVED_MESSAGE)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
