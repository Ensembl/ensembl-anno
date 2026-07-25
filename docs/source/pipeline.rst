Pipeline
========

The Ensembl Anno workflow consists of five analysis groups followed by
gene set finalisation.

Pipeline overview

.. code-block:: text
Genome
   │
   ▼
Repeat annotation
   │
   ├── Red
   ├── Dust
   ├── TRF
   └── RepeatMasker
   │
   ▼
Simple features
   │
   ├── CpG
   └── Eponine
   │
   ▼
Transcriptomic evidence
   │
   ├── STAR
   ├── StringTie
   ├── Scallop
   └── minimap2
   │
   ▼
Protein evidence
   │
   ├── GenBlast
   └── Miniprot
   │
   ▼
Gene set finalisation
   │
   ▼
Optional loading into Ensembl Core DB