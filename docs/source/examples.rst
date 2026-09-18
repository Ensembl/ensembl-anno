Examples
========

Complete annotation

.. code-block:: bash

    python ensembl_anno.py \
        --genome_file genome.fa \
        --output_dir output \
        --run_full_annotation

Transcriptome only

.. code-block:: bash

    python ensembl_anno.py \
        --run_transcriptomic \
        --short_read_fastq_dir reads/

Protein annotation

.. code-block:: bash

    python ensembl_anno.py \
        --run_genblast \
        --protein_file proteins.fa