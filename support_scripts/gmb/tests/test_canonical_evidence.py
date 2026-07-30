"""Tests for gmb.pipeline.canonical_evidence -- the shared evidence-class
vocabulary and balanced-coverage helper used by both canonical_selection.py
and interpro_review.py.
"""

from gmb.pipeline.canonical_evidence import (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_LONG_READ,
    EVIDENCE_CLASS_OTHER,
    EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
    EVIDENCE_CLASS_PROTEIN_VALIDATION,
    EVIDENCE_CLASS_SHORT_READ,
    balanced_coverage,
    evidence_classes_for_transcript,
    named_source_evidence_classes,
)


class TestNamedSourceEvidenceClasses:
    def test_scallop_and_stringtie_collapse_to_one_short_read_class(self):
        classes, unknown = named_source_evidence_classes("Scallop,StringTie")
        assert classes == {EVIDENCE_CLASS_SHORT_READ}
        assert unknown == set()

    def test_long_and_short_read_stay_distinct(self):
        classes, _ = named_source_evidence_classes("Scallop,Minimap2")
        assert classes == {EVIDENCE_CLASS_SHORT_READ, EVIDENCE_CLASS_LONG_READ}

    def test_backbone_matched_against_configured_label(self):
        classes, _ = named_source_evidence_classes("Tiberius", backbone_label="Tiberius")
        assert classes == {EVIDENCE_CLASS_BACKBONE}
        # Same source string, different configured backbone -> not backbone.
        classes2, unknown2 = named_source_evidence_classes("Tiberius", backbone_label="Helixer")
        assert classes2 == {EVIDENCE_CLASS_OTHER}
        assert unknown2 == {"tiberius"}

    def test_protein_alignment_sources(self):
        classes, _ = named_source_evidence_classes("OrthoDB,GenBlast,UniProt")
        assert classes == {EVIDENCE_CLASS_PROTEIN_ALIGNMENT}

    def test_unknown_source_mapped_to_other_and_reported_never_dropped(self):
        classes, unknown = named_source_evidence_classes("BrandNewAssembler")
        assert classes == {EVIDENCE_CLASS_OTHER}
        assert unknown == {"brandnewassembler"}

    def test_empty_and_none(self):
        assert named_source_evidence_classes("") == (set(), set())
        assert named_source_evidence_classes(None) == (set(), set())

    def test_never_includes_protein_validation(self):
        # protein_validation is derived from the transcript's own validation
        # result, never from the named-source string -- even a source
        # literally named that would not collide, since sources are matched
        # case-insensitively against _SOURCE_TO_CLASS/backbone_label only.
        classes, _ = named_source_evidence_classes("Scallop,StringTie,Minimap2,OrthoDB")
        assert EVIDENCE_CLASS_PROTEIN_VALIDATION not in classes


class TestEvidenceClassesForTranscript:
    def test_adds_protein_validation_when_supported(self):
        classes, _ = evidence_classes_for_transcript(
            "Scallop,StringTie", "Helixer", has_protein_validation_support=True
        )
        assert classes == {EVIDENCE_CLASS_SHORT_READ, EVIDENCE_CLASS_PROTEIN_VALIDATION}

    def test_omits_protein_validation_when_unsupported(self):
        classes, _ = evidence_classes_for_transcript(
            "Scallop,StringTie", "Helixer", has_protein_validation_support=False
        )
        assert classes == {EVIDENCE_CLASS_SHORT_READ}

    def test_protein_validation_alone_for_source_free_transcript(self):
        classes, _ = evidence_classes_for_transcript(
            "", "Helixer", has_protein_validation_support=True
        )
        assert classes == {EVIDENCE_CLASS_PROTEIN_VALIDATION}


class TestBalancedCoverage:
    def test_uses_minimum_not_mean(self):
        # A fragment match (100% of a short query, 10% of a long target)
        # must not be flattered by averaging.
        assert balanced_coverage(100.0, 10.0) == 0.1

    def test_symmetric(self):
        assert balanced_coverage(30.0, 90.0) == balanced_coverage(90.0, 30.0)

    def test_none_when_either_missing(self):
        assert balanced_coverage(None, 90.0) is None
        assert balanced_coverage(90.0, None) is None
        assert balanced_coverage(None, None) is None

    def test_full_coverage_both_sides(self):
        assert balanced_coverage(100.0, 100.0) == 1.0
