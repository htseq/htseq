.. _step_vs_stretch:

******************************************
Tutorial: StepVector versus StretchVector
******************************************

.. currentmodule:: HTSeq

One of the main uses of ``HTSeq`` is to represent data along a genome: coverage, genes and other features, transcription factor (TF) binding, etc. Storing long genomes in memory as a dense array (a la ``numpy``) is not only impractical on most laptop/desktop computers, but also rarely useful. For instance, if you are only interested in gene bodies, you can save a lot of memory by ignoring intergenic regions altogether. For TF binding sites, most of the genome does not harbour a compatible binding motif with any given TF. For chromatin conformation, most pairs of sites in the genome are not in contact. This observation leads to a need for "sparse" data representations.

However the same implementation of sparsity might not suit each and every application. ``HTSeq`` offers two main sparse representations for genomic data, through the classes
:class:`StepVector` and :class:`StretchVector`. This tutorial will give an overview of the concepts and potential applications for both of these interfaces.


``StepVector``: Piecewise constant values at low resolution
-----------------------------------------------------------
The first way to represent sparse genomic data, implemented by :class:`StepVector`, is useful for data that are piecewise constant. A typical example are the gene overlaps used in ``htseq-count``: perhaps the first exon starts at 1450972 and occupies 1kb, the second exon is 500bp on its own and 300bp overlapping with another gene, the third exon is 10kb downstream, and so on. :class:`StepVector` would represent this situation as four steps:

1. ``[1450972, 1451972)`` (i.e. from 1450972, 0-based and included to 1451972 excluded): exon 1
2. ``[1453972, 1454472]``: exon 2
3. ``[1454472, 1454772]``: exon 2 and overlapping gene
4. ``[1464772, ...]``: exon 3

Within each step, the data is constant and the step or interval width is much larger than 1bp. In a way, this representation is great for low-resolution information compared to 1bp, just like bedGraph files are. If you imagine this kind of data as a function of genomic coordinate, this representation is ideal for slowly fluctuating functions.

``StretchVector``: Islands of rapidly fluctuaring data
--------------------------------------------------------
The second option included in ``HTSeq`` to represent sparse genomic data is via a :class:`StretchVector`. An example would be targeted sequencing of a specific gene or a few genes: most of the genome has no data, but a lot of extremely detailed information, down to the single nucleotide level, is available for the genomic stretches covered by the sequencing.

The core idea behind :class:`StretchVector` is to store two connected lists, one containing dense numpy arrays representing the data on each stretch, and the second containing the start and end coordinate of each stretch. In a way, each stretch behaves a bit like a step in :class:`StepVector`, but the value needs not be constant within the stretch.

In its simplest form, when you need only one stretch, :class:`StretchVector` basically behaves like a numpy array with an associated start coordinate. Keeping track of the start saves you some headaches in offsetting the genomic reference coordinates all the time. If you have several stretches, this bookkeeping can become onerous and lead to bugs, so the infrastructure provided by this class is even more useful.


``ChromVector`` and ``GenomicArray``
------------------------------------
Both :class:`StepVector` and :class:`StretchVector` are relatively low-level data structures: they cover exactly *one* linear genomic stretch, have no notion of strandedness, and have their own constructor and methods.

``HTSeq`` also offers higher-level data structures that can make life a little easier. The first is :class:`ChromVector`, which implement the abstract concept of data on *one* chromosome. When you create a new :class:`ChromVector`, you can choose a storage method depending on what kind of data you have:

1. short chromosome, dense data: ``ndarray`` (a full ``numpy`` array)
2. long chromosome, dense data: ``memmap`` (i.e. a disk-stored memory map, which is slow but enables you to store data larger than your computer RAM)
3. long chromosome, piecewise constant sparsity: ``step``, which uses :class:`StepVector`
4. long chromosome, island sparsity: ``stretch``, which uses :class:`StretchVector`

The methods of this class will work independently of the storage chosen, even though their speed and memory performance will be affected by the nature of the data.

If you are juggling more than one chromosome at once, including whole genomes, your friend is :class:`GenomicArray`, which is essentially a dictionary of :class:`ChromVector` entries, one per chromosome for unstranded data and two (strands) per chromosome if the data is stranded.

Ultimately, :class:`GenomicArray` is a good starting point if you are just dipping your feet in the ``HTSeq`` API, while :class:`StepVector` and :class:`StretchVector` might become more useful as you become more familiar with the library and want to develop custom analyses.
