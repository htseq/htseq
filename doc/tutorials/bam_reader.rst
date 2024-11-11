.. _tutorial_bam:

**********************************************
Tutorial: Using the SAM/BAM/CRAM parsers
**********************************************

.. currentmodule:: HTSeq

High-throughput sequencing machines commonly produce FASTQ files. For many analyses,
reads are *aligned* to a reference genome by a "mapper" or "aligner," which typically
writes the output *alignments* as SAM/BAM/CRAM files.
``HTSeq`` has a set of parsers for these file formats.

This tutorial will walk you through a few routine operations with `SAM/BAM/CRAM`_ files.
The example files used here are all in the ``example_data`` folder of the ``HTSeq``
repository on GitHub.

.. _`SAM/BAM/CRAM`: http://samtools.github.io/hts-specs/SAMv1.pdf

What's the difference?
----------------------
tl;dr:

- SAM is an uncompressed, text format
- BAM is its compressed, binary sibling
- CRAM is a newer compressed, binary format that evolved from BAM but is less widespread

``samtools`` can be used to simply convert between them from the command line.

.. note:: SAM files sometimes have a header containing the names and lengths of
   chromosomes and a few more info (e.g. the aligner name). For BAM files, the
   header is mandatory instead of optional.


Parsing a SAM file
------------------
To parse a SAM (text) file, you can use the :class:`SAM_Reader` class::

   >>> f = HTSeq.SAM_Reader('yeast_RNASeq_excerpt_withNH.sam')

You can use context management to ensure the file is closed properly. If we look
at the first 3 reads only::

   >>> with HTSeq.SAM_Reader('yeast_RNASeq_excerpt_withNH.sam') as f:
   ...     for i, read in enumerate(f):
   ...         print(read)
   ...         if i == 2:
   ...             break
   <SAM_Alignment object: Read 'HWI-EAS225:1:10:1284:1974#0/1' aligned to VIII:[394100,394136)/+>
   <SAM_Alignment object: Read 'HWI-EAS225:1:10:1284:986#0/1', not aligned>
   <SAM_Alignment object: Read 'HWI-EAS225:1:10:1284:2012#0/1' aligned to VII:[1001605,1001641)/+>

.. note:: Some SAM/BAM files derive from *paired-end* experiments. The ``for`` loop
   above iterates over reads, not read pairs. See below for a parser that specifically
   deals with read pairs.

Parsing a BAM file
--------------------
To parse a BAM (binary) file, you can use the sibling class :class:`BAM_Reader`. Using a
context::

   >>> with HTSeq.BAM_Reader('SRR001432_head.bam') as f:
   ...     for i, read in enumerate(f):
   ...         print(read)
   ...         if i == 2:
   ...             break
   <SAM_Alignment object: Read 'SRR001432.1 USI-EAS21_0008_3445:8:1:107:882 length=25', not aligned>
   <SAM_Alignment object: Read 'SRR001432.2 USI-EAS21_0008_3445:8:1:82:90 length=25', not aligned>
   <SAM_Alignment object: Read 'SRR001432.3 USI-EAS21_0008_3445:8:1:639:904 length=25', not aligned>

As for SAM files, the iteration is over each single read, not read pairs.

Paired-end reads
----------------
Illumina sequencers can sequence a DNA molecule from both ends, leading to so-called
*paired-end* reads. It is sometimes useful to examine both reads of each pair (called
mates) at the *same time*, because they are two representative of the same original
DNA molecule.

``HTSeq`` has two classes to achieve this goal, depending on whether the SAM/BAM file
is "unsorted" aka *sorted by name*, or "sorted" aka *sorted by position*. In the former
case, the alignments from the two mates are found in the BAM file on consecutive lines.

For name-sorted SAM/BAM files, you can use :func:`pair_SAM_alignments`::

   >>> with HTSeq.BAM_Reader('SRR001432_head.bam') as f:
   ...     for (read1, read2) in HTSeq.pair_SAM_alignments(f):
   ...         print(read1)
   ...         print(read2)

For position-sorted SAM/BAM files, you can use :func:`pair_SAM_alignments_with_buffer`::

   >>> with HTSeq.BAM_Reader('SRR001432_head.bam') as f:
   ...     for (read1, read2) in HTSeq.pair_SAM_alignments_with_buffer(f):
   ...         print(read1)
   ...         print(read2)

In practical terms, the main difference between these two functions is their operation
when they meet a read but the next read is not the mate: the first function assumes the
first read is an incomplete pair, while the second function waits to see if the mate will
appear later on in the file.

.. note:: The second function is used in ``htseq-count`` to provide the `-r pos` option
   for position-sorted BAM files. However, because the first read has to be stored in
   memory until the second read is found, this approach incurs a significant memory cost.
   It is recommended to call ``htseq-count`` on unsorted/name-sorted BAM files, and
   ``samtools sort`` can be used to "unsort" a BAM file.
