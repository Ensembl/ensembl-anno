#!/usr/bin/env python3
"""CLI entry point: prepare InterProScan review candidates.

    python -m gmb.cli.interpro_review \
        --canonical-transcripts out/canonical_selection/canonical_transcripts.tsv \
        --transcript-ranking out/canonical_selection/transcript_ranking.tsv \
        --protein-validation out/protein_validation.tsv \
        --prot-fa out/prot.fa \
        --config my_config.yaml \
        --output-dir out/interpro_review

Requires canonical_selection.interpro_resolver.enabled: true in config (the
default is false -- see gmb.pipeline.config.InterProResolverConfig). When
that section also sets run_interproscan: true, this also launches
InterProScan 6 itself via Nextflow (Mode A) using the
canonical_selection.interpro_resolver.nextflow: settings; the default
(false) is Mode B -- prepare candidates only, and consume an
already-completed InterProScan run separately via gmb.cli.interpro_resolve.
"""

from gmb.pipeline.interpro_review import main

if __name__ == "__main__":
    main()
