.. _exon_example:

****************************************
Opening the black box of ``htseq-count``
****************************************

.. currentmodule:: HTSeq

``htseq-count`` and ``htseq-count-barcodes`` are popular scripts to quantify transcript abundances of a RNA-Seq (bulk or single cell) experiments. They are designed to be used from the command line and support a number of options (see :ref:`htseqcount` and :ref:`htseqcount_with_barcodes`). However, it is instructive to describe the main steps of those scripts to understand what they do and hopefully encourage power users to customize these scripts to fit their own needs.

As an example, we will look at a typical ``htseq-count`` call applied on a yeast experiment. First, you need two files: a GTF file with the genome annotations and a SAM/BAM/CRAM file with the mapped reads. For the annotations, we will use ``Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz`` from the ``example_data`` folder, while for the reads we will use ``yeast_RNASeq_excerpt.sam`` from the same folder.

Step one: loading the annotation GTF file
-----------------------------------------
The first step is to learn where the genomic features are in the genome.

We start by opening the GTF file::

   >>> import HTSeq
   >>> gtffile = HTSeq.GFF_Reader("Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz")

We then scan the GTF file to store the genomic locations (chromosome, start and end coordinates) of the features of interest::

   >>> feature_scan = HTSeq.make_feature_genomicarrayofsets(
   ...     gtffile,
   ...     id_attribute='gene_id',
   ...     feature_type='exon',
   ... )

``feature_scan`` is an efficient data structure (:class:`GenomicArrayOfSets`) that contains many intervals (:class:`GenomicInterval`) and one or more exons covering each interval. ``id_attribute`` describes which attribute of each feature we want to use to label the intervals covered by that feature.

For instance, if there is an exon at the interval ``Chromosome 1, pos 1566 to 1589`` and **no other exon covers that area**, ``feature_scan`` would connect that interval to the ``gene_id`` associated with that exon. If, however, an exon from an overlapping gene covers the interval ``Chromosome 1, pos 1580 to 1589``, there will be two intervals defined:

- ``Chromosome 1, pos 1566 to 1579`` associated with the first ``gene_id``
- ``Chromosome 1, pos 1580 to 1589`` associated with both genes

So, if a read falls onto the first interval, it will be assigned to the first gene, while if it falls onto the second interval, it will be ambiguous: we cannot know which transcript that read was originating from.

Step two: scanning the BAM file for reads or read pairs
-------------------------------------------------------
The second step is to scan the BAM file and classify each read - or, for paired end sequencing, each read pair - into one or more genomic features based on its location.

We open the SAM/BAM file (in this case, for simplicity, it's an uncompressed SAM file)::

   >>> bamfile = HTSeq.BAM_Reader("yeast_RNASeq_excerpt.sam")

.. note:: For the sake of simplicity, this SAM file is not paired-end, so we can examine one read at a time. Paired-end sequencing creates some challenges that will be mentioned at the end of this page.

We have to prepare some simple data structures where the counts will be saved::

   >>> attributes = feature_scan['attributes']
   >>> feature_attr = sorted(attributes.keys())
   >>> counts = {key: 0 for key in feature_attr}

We also prepare a few "quasi-features" that describe unmapped, ambiguous, and other not uniquely assigned reads::

   >>> counts['notaligned'] = 0
   >>> counts['no_feature'] = 0
   >>> counts['ambiguous'] = 0

(there might be a few more categories in a real-world scenario).

We then iterate over ``bamfile`` and examine each read::

   >>> for read in bamfile:
   ...     if not r.aligned:
   ...         counts['notaligned'] += 1
   ...         continue

(script continues below). These lines just increase the ``count`` dictionary whenever we find an unaligned read. Of course, the more interesting reads are the ones that do align with the genome. For those, we have to collect all genomic intervals where the read aligns, skipping over and insertions and deletion. This operation is enabled by the so-called CIGAR codes of the read::

   ...     aligned_codes = ('M', '=', 'X')
   ...     iv_read = (co.ref_iv for co in read.cigar if co.type in aligned_codes)

In practice, there are a few subtleties, such as checking for strandedness and read quality, but we'll gloss over those fine points here.

Step three: overlapping each read with the gene intervals
---------------------------------------------------------
Now for each read we know which genomic interval it covers: it's time to compare those coordinates with each GTF feature of interest to check for overlaps::

   ...     gene_ids_read = None
   ...     for iv in iv_read:
   ...         for _, gene_ids in feature_scan['features'][iv].steps():
   ...             if features_read is None:
   ...                 gene_ids_read = gene_ids.copy()
   ...             else:
   ...                 features_read.intersection(gene_ids)

This is a lot to unpack, so let's walk through it slowly. First, because of insertion/deletions and intron splicing, a single read might have two or more distinct intervals that align to the genome. Each of these intervals (in most cases there's only one) is indexed by the variable ``iv``.

For each ``iv``, we go back to the ``feature_scan`` and extract all ``gene_ids`` associated with that interval. In most cases, there's only one gene: in any case, we add them all to ``gene_ids_read``.

You may wonder what the ``intersection`` is used for. That's what ``htseq-count`` calls a "mode". In "intersection-strict" mode (this example), we exclude genes that overlap with only part of the read (i.e. not all ``.steps()``), which ensures we only count reads that are fully within that gene. If we were to choose another "mode", things would be different. For example, in the mode "union", **all** genes overlapping with the read, even by just a single base, are counted: if that's more than one gene, this read will be classified as "multimapper".

Step four: uniquely mapped VS multimapped
-----------------------------------------
The next step is: we just check if any gene overlapped with the read at all and, if so, how many::

   ...     if gene_ids_read is None or len(gene_ids_read) == 0:
   ...         counts['no_feature'] += 1
   ...         continue

Okay, this is simple: if the read didn't overlap with any gene, put it into a special bin. However, if it did overlap::

   ...     if len(gene_ids) > 1:
   ...         counts['ambiguous'] += 1
   ...         continue

This is the opposite situation: the read overlaps with more than one gene and it's impossible to know which one is the true origin of that molecule. These reads are also put into a special bin.

The most common case, in most application, is when a read is **uniquely mapped**, i.e. there was only a single gene that overlapped with its intervals. In that case, we assign the read to that gene::

   ...     gene_id = list(gene_ids)[0]
   ...     counts[gene_id] += 1

Step five: closing the files
-------------------------------
At the end of the for loop, we are done with the counting, and the only thing left for us to do is to diligently close the SAM and GTF file we opened earlier on::

   >>> bamfile.close()
   >>> gtffile.close()

Conclusion
----------
Hopefully this case study helped you understand the design and mechanics of ``htseq-count``. If you want to learn more, you should take a look at the script itself, or start coding your own custom variant straight away!

.. note::
   Note that instead of counting on a gene level (using ``gene_id``), we can count reads on an exon level, i.e. each exon will be a separate line in the output. For that purpose, you can use a different call for :func:`make_feature_genomicarrayofsets`::

     >>> feature_scan = HTSeq.make_feature_genomicarrayofsets(
     ...     gtffile,
     ...     id_attribute=['gene_id', 'exon_number'],
     ...     feature_type='exon',
     ... )

.. note:: **Paired-end** BAM files can be a little more tricky to analyze. Notably, paired-end BAM files can be unsorted (the two reads of each pair appear on subsequent lines) or sorted by position. Sorting usually separates the two reads within each pair, and the gap can consist of million of other reads. When we count genes, each read *pair* counts for 1, but since the two reads are now disjoint, a buffering mechanism is used to "delay" counting until the second ("mate") read has been found. This can slow down the script, increase memory consumption and, in extreme cases, crash the analysis. As a general rule, always try to run ``htseq-count`` on unsorted (or "sorted-by-name") BAM files.
