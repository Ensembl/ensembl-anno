"""Gene Model Builder (GMB) — consensus gene annotation for eukaryotic genomes.

GMB is the evidence-integration and final gene-model selection stage of an
annotation pipeline. It consumes an ab initio backbone, assembled transcript
models and protein alignments, and emits a finished annotation with per-model
evidence attribution, QC and provenance.

Public interface
----------------
Command line::

    gmb-preflight   validate an evidence bundle before building
    gmb-build       integrate evidence and select gene models
    gmb-finalise    regenerate FASTA from the final GFF3, QC, canonical, handover
    gmb-compare     evaluate against a reference annotation (EVALUATION ONLY)

Python::

    from gmb import run_gene_model_builder

Everything under ``gmb.pipeline``, ``gmb.preflight`` and ``gmb.provenance`` is
internal. An orchestration layer should use ``run_gene_model_builder`` or the
CLI, never those modules directly.

Documentation: ``docs/quickstart.md``, ``docs/pipeline_integration.md``.
"""

from gmb.api import GmbResult, run_gene_model_builder  # noqa: F401

__version__ = "2.0.0"

#: Human-facing name, short name, and the pipeline stage identifier an
#: orchestration layer should use to refer to this module.
PRODUCT_NAME = "Gene Model Builder"
SHORT_NAME = "GMB"
PIPELINE_STAGE = "gene_model_builder"

__all__ = ["run_gene_model_builder", "GmbResult", "__version__",
           "PRODUCT_NAME", "SHORT_NAME", "PIPELINE_STAGE"]
