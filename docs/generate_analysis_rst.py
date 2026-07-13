"""Generate RST files for each analysis module in the ensembl.tools.anno package."""
from pathlib import Path
import importlib

SRC = Path("src/python/ensembl/tools/anno")
OUT = Path("docs/source/analyses")

GROUPS = {
    "repeat_annotation": "Repeat Annotation",
    "simple_feature_annotation": "Simple Feature Annotation",
    "snc_rna_annotation": "Small ncRNA Annotation",
    "transcriptomic_annotation": "Transcriptomic Annotation",
    "protein_annotation": "Protein Annotation",
    "ab_initio": "Ab initio Annotation",
}

HEADER = """{title}
{underline}

{summary}

API
---

.. automodule:: {module}
    :members:
    :undoc-members:
"""

OUT.mkdir(parents=True, exist_ok=True)

for package, group in GROUPS.items():

    package_dir = SRC / package

    for module in sorted(package_dir.glob("*.py")):

        if module.name.startswith("_"):
            continue

        name = module.stem

        module_name = f"ensembl.tools.anno.{package}.{name}"

        try:
            mod = importlib.import_module(module_name)
            summary = (mod.__doc__ or "").strip().split("\n")[0]
        except Exception:#pylint: disable=broad-except
            summary = ""

        title = name.replace("_", " ").title()

        rst = HEADER.format(
            title=title,
            underline="=" * len(title),
            module=module_name,
            summary=summary,
        )

        with open(OUT / f"{name}.rst", "w", encoding="utf-8") as fh:
            fh.write(rst)