Wrapper
=======

The wrapper script (ensembl_anno.py) orchestrates the complete annotation
workflow.

Basic usage

.. code-block:: bash

    python ensembl_anno.py \
        --genome_file genome.fa \
        --output_dir output

Analysis groups

* Repeat annotation
* Simple feature annotation
* Transcriptomic annotation
* Protein annotation
* Small ncRNA annotation
* Finalisation

Automatic dependency resolution
-------------------------------

The wrapper enables analyses only when the required inputs are available.

Examples

* STAR requires short-read FASTQs.
* minimap2 requires long-read FASTQs.
* cmsearch requires an Rfam accession file.
* GenBlast requires protein FASTA.