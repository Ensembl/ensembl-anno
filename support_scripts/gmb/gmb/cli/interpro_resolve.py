#!/usr/bin/env python3
"""CLI entry point: review ambiguous canonical choices against InterProScan output.

    python -m gmb.cli.interpro_resolve \
        --manifest out/interpro_review/interpro_review_manifest.tsv \
        --interpro-jsonl results.faa.jsonl \
        --config my_config.yaml \
        --canonical-transcripts out/canonical_selection/canonical_transcripts.tsv \
        --consensus-gff3 out/consensus.gff3 \
        --output-dir out/interpro_review

--canonical-transcripts/--consensus-gff3 are optional; when both are given,
also writes canonical_decisions.tsv and consensus.final_canonical_annotated.gff3
(see docs/interpro_resolver.md). See gmb.pipeline.config.InterProResolverConfig
(canonical_selection.interpro_resolver: in YAML) for the replacement policy's
enabled/apply_replacements/safeguard settings.
"""

from gmb.pipeline.interpro_resolver import main

if __name__ == "__main__":
    main()
