.. _tutorial_htseq_exon:

******************************************************
Quantifying exon-level expression with ``htseq-count``
******************************************************

.. currentmodule:: HTSeq

``htseq-count`` and ``htseq-count-barcodes`` are commonly used to quantify gene expression. They identify reads in exonic areas of each gene and, unless there is an ambiguous overlap with a second gene in the same genomic region, assign each read - or read pair - to a gene.

Once can also use these scripts to quantify how many reads overlap with **each specific exon**. This tutorial explains how to leverage that feature of ``HTSeq``.


How to use it
-------------
To run ``htseq-count`` with exon-level counts, run::

   >>> htseq-count -i gene_id -i exon_number --additional-attr gene_name --additional-atr exon_number

Because exons do not have canonical names, this call creates a table where each row (or column, depending on the output format) is identified by a string ``geneid:exon_number``, e.g. ``ENSG00000223972:1`` for the first exon of gene DDX11L1 (which has gene_id ENSG00000223972). The table will have two additional metadata columns, one with the gene name and one with the exon number, to make life a little easier for the user.

Caveats with this approach
--------------------------
Because of splicing, counting at the exon level is intrinsically more noisy than at the gene level. The molecules captures in typical assays, such as RNA-Seq, correspond to spliced isoforms or unspliced transcripts, not directly to exons. Therefore, whenever a read is spannig an exon-exon junction, an exon-level counter can be confused.

If the data comes from single-end sequencing, a single exon-exon spanning read might be considered ambiguous by ``htseq-count``. To count it 50%-50% for each exon, you can use the option ``--nonunique=fraction``. The option ``--nonunique=random`` assigns the read to either exon at random, which gives fairly similar results to using fractions while maintaining integer (not floating point) results.

If the data comes from paired-end sequencing, it's likely that at least one of the two mate reads falls into a single exon, therefore the pair should be uniquely assigned to that exon in most cases. If one read is assigned to one exon and the other to another exon, the result depends on the ``-m/--mode`` option. For ``--mode=union``, the pair is assigned to both exons and therefore a multimapper, and the same ``--nonunique`` options apply as for single-end sequencing (see above). For *intersection* modes, the script cannot found a single exon that covers the whole pair, and the read pair is assigned as ``__no_feature``.

.. note:: To be sure: this counting scheme does not enable biologically faithful reconstruction of isoform counts, however it can be used to detect up- or down-regulation of specific exons beyond standard RNA-Seq analyses in both bulk and single cell experiments. **If you have used this feature, we'd like to hear about your experience**: write us an email!
