"""Run provenance for Gene Model Builder.

Records everything needed to answer "what produced this annotation, and could I
produce it again?" — software versions, the exact command, every input path with
its SHA-256, the resolved evidence roles and weights, the selection policy that
actually applied (including gates resolved at run time), and resource use.

This complements the handover manifest written by ``gmb-finalise``: that one
describes the OUTPUTS, this one describes the RUN.

Used by ``gmb-build``, which writes ``run_manifest.json`` and ``run_manifest.tsv``
into its output directory.
"""

from gmb.provenance.manifest import (  # noqa: F401
    RunManifest,
    build_run_manifest,
    write_run_manifest,
)

__all__ = ["RunManifest", "build_run_manifest", "write_run_manifest"]
