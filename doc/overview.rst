.. _overview:

.. currentmodule:: HTSeq

************************************************************
Overview
************************************************************

* :ref:`install`

  Download links and installation instructions can be found here

* :ref:`tour`

  The Tour shows you how to get started. It explains how to install HTSeq, and then
  demonstrates typical analysis steps with explicit examples. Read this first, and 
  then see the Reference for details.
  
* :ref:`tss`

  This chapter explains typical usage patterns for HTSeq by explaining in detail 
  three different solutions to the same programming task.

* :ref:`counting`

  This chapter explorer in detail the use case of counting the overlap of reads
  with annotation features and explains how to implement custom logic by
  writing one's own customized counting scripts


* Reference documentation

  The various classes of `HTSeq` are described here.

  * :ref:`refoverview` 

    A brief overview over all classes.

  * :ref:`sequences` 
  
    In order to represent sequences and reads (i.e., sequences with base-call quality 
    information), the classes :class:`Sequence` and :class:`SequenceWithQualities` are used.
    The classes :class:`FastaReader` and :class:`FastqReader` allow to parse FASTA and FASTQ
    files.
  
  * :ref:`genomic`
  
    The classes :class:`GenomicInterval` and :class:`GenomicPosition` represent intervals and
    positions in a genome. The class :class:`GenomicArray` is an all-purpose container
    with easy access via a genomic interval or position, and :class:`GenomicArrayOfSets`
    is a special case useful to deal with genomic features (such as genes, exons,
    etc.)
    
  * :ref:`alignments`
  
    To process the output from short read aligners in various formats (e.g., SAM),
    the classes described here are used, to represent output files and alignments,
    i.e., reads with their alignment information.

  * :ref:`features`
  
    The classes :class:`GenomicFeature` and :class:`GFF_Reader` help to deal with genomic
    annotation data.
    
  * :ref:`otherparsers`
  
    This page describes classes to parse VCF, Wiggle and BED files.

  * :ref:`misc`


* Scripts

  The following scripts can be used without any Python knowledge.
  
  * :ref:`qa`
  
    Given a FASTQ or SAM file, this script produces a PDF file with plots depicting
    the base calls and base-call qualities by position in the read. This is useful to
    assess the technical quality of a sequencing run.
  
  * :ref:`htseqcount`
  
    Given one/multiple SAM/BAM/CRAM files with alignments and a GTF file with genomic
    features, this script counts how many reads map to each feature. This script is
    especially popular for bulk and single-cell RNA-Seq analysis.

  * :ref:`htseqcount_with_barcodes`

    Similar to `htseq-count`, but for a single SAM/BAM/CRAM file containing reads with
    cell and molecular barcodes (e.g. 10X Genomics `cellranger` output). This script
    enables customization of single-cell RNA-Seq pipelines, e.g. to quantify exon-level
    expression or simply to obtain a count matrix that contains chromosome information
    additional feature metadata.

* Appendices

..

  * :ref:`history`
  
..  
  
  * :ref:`contrib`

..

  * :ref:`Table of Contents<tableofcontents>`

..

