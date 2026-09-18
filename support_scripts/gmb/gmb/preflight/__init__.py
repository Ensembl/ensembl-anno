"""Pre-build input validation for Gene Model Builder.

Runs before an expensive build and answers one question: *is this evidence
bundle fit to annotate from?*

Every check is reference-free. Preflight never reads a reference annotation and
never looks at the organism name; it inspects only the files that will actually
be handed to ``gmb-build``.

The check that earns this module its place is **splice quality**. A long-read
consensus track used for a production P. falciparum build turned out to be
16.4% canonically spliced, against 99-100% for every other track; locus-by-locus
it worsened as many gene models as it improved. Nothing in the pipeline noticed,
because every downstream QC check it touched still passed. Preflight would have
flagged it before the build.

Public API::

    from gmb.preflight import run_preflight, PreflightReport
    report = run_preflight(inputs, config)
    report.verdict          # "PASS" | "WARN" | "FAIL"
    report.to_dict()        # JSON-serialisable
"""

from gmb.preflight.checks import (  # noqa: F401
    FAIL,
    PASS,
    WARN,
    CheckResult,
    PreflightReport,
    run_preflight,
)

__all__ = ["run_preflight", "PreflightReport", "CheckResult", "PASS", "WARN", "FAIL"]
