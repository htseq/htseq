import itertools
import pysam
import random
import HTSeq
import sys

from HTSeq.scripts.utils import (
    UnknownChrom,
    invert_strand
)

from HTSeq.scripts.count_reads_in_feature.reads_stats import ReadsStatistics

def get_read_intervals(do_invert, read_sequence):
    com = ('M', '=', 'X')
    if do_invert:
        iv_seq = (invert_strand(co.ref_iv)
                  for co in read_sequence.cigar if (co.type in com and
                                        co.size > 0))
    else:
        iv_seq = (co.ref_iv for co in read_sequence.cigar if co.type in com
                  and co.size > 0)
    return iv_seq

def check_paired_end_read(read_sequence,
                           reads_stats,
                           secondary_alignment_mode,
                           supplementary_alignment_mode,
                           minaqual,
                           multimapped_mode):

    """
    Function to check the read for paired end.

    Parameters
    ----------
    read_sequence
    reads_stats : ReadsStatistics object
        Object which stores a bunch of statistics about the read sequences.
    secondary_alignment_mode : str
        Whether to score secondary alignments (0x100 flag).
        Choices: score or ignore.
    supplementary_alignment_mode : str
        Whether to score supplementary alignments (0x800 flag).
        Choices: score or ignore.
    minaqual : int
        Value denoting the MAPQ alignment quality of reads to skip.
    multimapped_mode : str
        Whether and how to score reads that are not uniquely aligned or
        ambiguously assigned to features.
        Choices: none, all, fraction, random.

    Returns
    -------
    True if read is to be skipped. False otherwise.
    """

    if read_sequence[0] is None or not read_sequence[0].aligned:
        reads_stats.add_not_aligned_read(read_sequence)
        return True

    if secondary_alignment_mode == 'ignore':
        if (
            read_sequence[0] is not None and
            read_sequence[0].not_primary_alignment
        ):
            return True
        elif (
            read_sequence[1] is not None and
            read_sequence[1].not_primary_alignment
        ):
            return True

    if supplementary_alignment_mode == 'ignore':
        if (read_sequence[0] is not None) and read_sequence[0].supplementary:
            return True
        elif (read_sequence[1] is not None) and read_sequence[1].supplementary:
            return True

    try:
        not_unique_read = False
        for i in range(0, 2):
            if (
                read_sequence[i] is not None and
                read_sequence[i].optional_field("NH") > 1
            ):
                not_unique_read = True
        if not_unique_read:
            reads_stats.add_not_unique_read(read_sequence)
            # TODO: Add unit test case for this
            if multimapped_mode == 'none':
                return True
    except KeyError:
        pass

    low_qual_read = False
    for i in range(0, 2):
        if (
            read_sequence[i] is not None and
            read_sequence[i].aQual < minaqual
        ):
            low_qual_read = True
    if low_qual_read:
        reads_stats.add_low_quality_read(read_sequence)
        return True

    return False

def check_non_paired_end_read(read_sequence,
                               reads_stats,
                               secondary_alignment_mode,
                               supplementary_alignment_mode,
                               minaqual,
                               multimapped_mode):
    """
    Function to check the read for non paired end.

    Parameters
    ----------
    read_sequence
    reads_stats : ReadsStatistics object
        Object which stores a bunch of statistics about the read sequences.
    secondary_alignment_mode : str
        Whether to score secondary alignments (0x100 flag).
        Choices: score or ignore.
    supplementary_alignment_mode : str
        Whether to score supplementary alignments (0x800 flag).
        Choices: score or ignore.
    minaqual : int
        Value denoting the MAPQ alignment quality of reads to skip.
    multimapped_mode : str
        Whether and how to score reads that are not uniquely aligned or
        ambiguously assigned to features.
        Choices: none, all, fraction, random.

    Returns
    -------
    True if read is to be skipped. False otherwise.
    """


    if not read_sequence.aligned:
        reads_stats.add_not_aligned_read(read_sequence)
        return True

    if (
        secondary_alignment_mode == 'ignore' and read_sequence.not_primary_alignment
    ):
        return True

    if (
        supplementary_alignment_mode == 'ignore' and
        read_sequence.supplementary
    ):
        return True

    try:
        if read_sequence.optional_field("NH") > 1:
            reads_stats.add_not_unique_read(read_sequence)
            if multimapped_mode == 'none':
                return True
    except KeyError:
        pass

    if read_sequence.aQual < minaqual:
        reads_stats.add_low_quality_read(read_sequence)
        return True

    return False

def _check_read_then_get_read_intervals(read_sequence, paired_end, reads_stats,
                        secondary_alignment_mode, supplementary_alignment_mode,
                        minaqual, stranded, multimapped_mode):

    """
    Function to break down the read sequence into intervals which will
    subsequently be processed.
    Some checks will be done to ensure downstream processing can be done.

    Parameters
    ----------
    read_sequence
    paired_end : boolean
        Is this a paired-end data?
    reads_stats : ReadsStatistics object
        Object which stores a bunch of statistics about the read sequences.
    secondary_alignment_mode : str
        Whether to score secondary alignments (0x100 flag).
        Choices: score or ignore.
    supplementary_alignment_mode : str
        Whether to score supplementary alignments (0x800 flag).
        Choices: score or ignore.
    minaqual : int
        Value denoting the MAPQ alignment quality of reads to skip.
    stranded : str
        Whether the data to be aligned is from a strand-specific assay.
        Option is yes, no, reverse.
        reverse means yes with reversed strand interpretation.
    multimapped_mode : str
        Whether and how to score reads that are not uniquely aligned or
        ambiguously assigned to features.
        Choices: none, all, fraction, random.

    Returns
    -------
    None or iv_seq
        None if the read is to be skipped.
        Else, the read broken into intervals (iv_seq variable).
    """

    if not paired_end:
        skip_this_read = check_non_paired_end_read(
            read_sequence = read_sequence,
            reads_stats = reads_stats,
            secondary_alignment_mode = secondary_alignment_mode,
            supplementary_alignment_mode = supplementary_alignment_mode,
            minaqual = minaqual,
            multimapped_mode = multimapped_mode
        )
        if skip_this_read:
            return None

        iv_seq = get_read_intervals(do_invert = stranded == 'reverse',
                                read_sequence = read_sequence)

    else:
        skip_this_read = check_paired_end_read(
            read_sequence = read_sequence,
            reads_stats = reads_stats,
            secondary_alignment_mode = secondary_alignment_mode,
            supplementary_alignment_mode = supplementary_alignment_mode,
            minaqual = minaqual,
            multimapped_mode = multimapped_mode
        )
        if skip_this_read:
            return None

        # This old if statement is no longer needed as by the time it
        # gets here, r[0] won't be none and it will have been aligned.
        # See check_paired_end_read.
        # if r[0] is not None and r[0].aligned:
        iv_seq = get_read_intervals(do_invert = stranded == 'reverse',
                                    read_sequence = r[0])

        # Note: the next end can be None or not aligned it will still
        # go through!
        if read_sequence[1] is not None and read_sequence[1].aligned:
            iv_seq_next_end = get_read_intervals(
                do_invert = stranded != 'reverse',
                read_sequence = read_sequence[1])
            iv_seq = itertools.chain(iv_seq, iv_seq_next_end)

    return iv_seq


def count_reads_single_file(
        isam,
        sam_filename,
        features,
        feature_attr,
        order,
        max_buffer_size,
        stranded,
        overlap_mode,
        multimapped_mode,
        secondary_alignment_mode,
        supplementary_alignment_mode,
        feature_type,
        id_attribute,
        additional_attributes,
        verbose,
        minaqual,
        samout_format,
        samout_filename,
        ):
    """
    The function that does the counting for each input BAM/SAM file.

    Note, if you change any of my parameters, please also change the following
    in count_reads_in_feature/parameters.py:
    1. CountParameters class, update the attributes.
    2. get_args_for_count function, update what is being appended into args
    variable returned by the function.

    Refactoring todo:
    1. Split this into reading the BAM files for reads
    2. Process the aligment.

    isam : int
        input files' indexing for the purpose of parallel processing.
        This basically tell you which input file is being processed by this
        instance of function.
    sam_filename : str
        Path to the SAM/BAM file containing the mapped reads.
    features : array
        TODO check the type of this parameter.
        Supplied by HTSeq.make_feature_genomicarrayofsets
    feature_attr : array
        TODO check the type of this parameter.
        Supplied by HTSeq.make_feature_genomicarrayofsets
    order : str
        Can only be either 'pos' or 'name'. Sorting order of <alignment_file>.
    max_buffer_size : int
        The number of reads allowed to stay in memory until mates are found.
        Used when <alignment_file> is paired end sorted by position.
    stranded : str
        Whether the data to be aligned is from a strand-specific assay.
        Option is yes, no, reverse.
        reverse means yes with reversed strand interpretation.
    overlap_mode : str
        Mode to handle reads overlapping more than one feature.
        Choices: union, intersection-strict, intersection-nonempty.
    multimapped_mode : str
        Whether and how to score reads that are not uniquely aligned or
        ambiguously assigned to features.
        Choices: none, all, fraction, random.
    secondary_alignment_mode : str
        Whether to score secondary alignments (0x100 flag).
        Choices: score or ignore.
    supplementary_alignment_mode : str
        Whether to score supplementary alignments (0x800 flag).
        Choices: score or ignore.
    feature_type : str
        Feature type (3rd column in GTF file) to be used, all features of other
        type are ignored (default, suitable for Ensembl, GTF files: exon).
    id_attribute : str
        GTF attribute to be used as feature ID.
        Normally gene_id, suitable for Ensembl GTF files.
    additional_attributes : array
        Additional feature attributes.
        Commonly, gene_name is suitable for Ensembl GTF files.
    verbose : boolean
        Whether to suppress progress report.
    minaqual : int
        Value denoting the MAPQ alignment quality of reads to skip.
    samout_format : str
        Format of the output files denoted by samouts.
        Choices: SAM, BAM, sam, bam.
    samout_filename : str
        The name of SAM/BAM file to write out all SAM alignment records into.

    Returns
    -------
    Dictionary
        TODO update me when done refactoring
    """

    read_io_object = ReadsIO(
        sam_filename = sam_filename,
        supplementary_alignment_mode = supplementary_alignment_mode,
        secondary_alignment_mode = secondary_alignment_mode,
        order = order,
        samout_format = samout_format
    )

    # To store some stats about the reads
    reads_stats = ReadsStatistics(
        read_io_object = read_io_object,
        feature_attr = feature_attr
    )

    try:
        for r in read_io_object.read_seq:
            iv_seq = _check_read_then_get_read_intervals(
                read_sequence = r,
                paired_end = read_io_object.pe_mode,
                reads_stats = reads_stats,
                secondary_alignment_mode = secondary_alignment_mode,
                supplementary_alignment_mode = supplementary_alignment_mode,
                minaqual = minaqual,
                stranded = stranded,
                multimapped_mode = multimapped_mode)

            # It means this read need to be skipped
            if iv_seq is None:
                continue

            try:
                map_and_count_read(overlap_mode = overlap_mode,
                             iv_seq = iv_seq,
                             features = features,
                             read_sequence = r,
                             multimapped_mode = multimapped_mode,
                             reads_stats = reads_stats,
                             read_io_object = read_io_object)

            except UnknownChrom:
                reads_stats.add_empty_read(r)

            if not verbose:
                reads_stats.print_progress()
            reads_stats.add_num_reads_processed()

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_io_object.bam_sam_file_reader.get_line_number_string()))
        raise

    if not verbose:
        reads_stats.print_progress(force_print=True)

    reads_stats.close_samout()

    res = reads_stats.generate_output(isam=isam)
    return(res)

def _get_feature_set(overlap_mode, iv_seq, features):

    # Combine all features as long as they somewhat overlap to the read
    if overlap_mode == "union":
        fs = set()
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                raise UnknownChrom
            for iv2, fs2 in features[iv].steps():
                fs = fs.union(fs2)

    # Only combine the features which fully overlap to the read.
    # How does the strict work? It will look for features which overlap
    # with every part of the read.
    # The difference with nonempty is when one part of the read align with
    # a gene, but the other align with no gene.
    # Strict option will say no feature aligned as there is a part of the
    # read that cannot be aligned to any feature.
    # Nonempty option will say the read aligns with that gene as the other
    # part is "empty" and thus not considered.
    elif overlap_mode in ("intersection-strict",
                          "intersection-nonempty"):
        fs = None
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                raise UnknownChrom
            for iv2, fs2 in features[iv].steps():
                if ((len(fs2) > 0) or
                   (overlap_mode == "intersection-strict")):
                    if fs is None:
                        fs = fs2.copy()
                    else:
                        fs = fs.intersection(fs2)
    else:
        sys.exit("Illegal overlap mode.")

    return fs

def _classify_read_against_feature(fs, read_sequence, reads_stats, read_io_object):
    if fs is None or len(fs) == 0:
        reads_stats.add_empty_read(read_sequence)
    elif len(fs) > 1:
        reads_stats.add_ambiguous_read(
            read_sequence = read_sequence,
            assignment = "__ambiguous[" + '+'.join(sorted(fs)) + "]"
        )
    else:
        # Read only aligns with one feature
        read_io_object.write_to_samout(
            read_sequence = read_sequence,
            assignment = list(fs)[0]
            )

def _assign_count(fs, read_sequence, multimapped_mode, reads_stats):
    if multimapped_mode == 'none':
        # Read only aligns with one feature, add to count.
        if len(fs) == 1:
            reads_stats.add_to_count(feature=list(fs)[0])
    elif multimapped_mode == 'all':
        # All features are counted.
        for fsi in list(fs):
            reads_stats.add_to_count(feature=fsi)
    elif multimapped_mode == 'fraction':
        # Count is divided among the features evenly.
        value = 1.0 / len(fs)
        for fsi in list(fs):
            reads_stats.add_to_count(feature=fsi,
                                     value=value)
    elif multimapped_mode == 'random':
        # Randomly pick a feature to assign the count to.
        fsi = random.choice(list(fs))
        reads_stats.add_to_count(feature=fsi)
    else:
        sys.exit("Illegal multimap mode.")

def map_and_count_read(overlap_mode, iv_seq,
                 features, read_sequence,
                 multimapped_mode, reads_stats,
                 read_io_object):
    """
    This function process each read by assigning it to features.

    Parameters
    ----------
    overlap_mode : str
        The overlap resolution modes, i.e. what to do if a read overlap with more than 1 feature.
        Option: union, intersection-strict, intersection-nonempty.

    iv_seq : array
        Intervals of reads.

    features : array
        A list of features to align the reads to.

    read_sequence : array
        The read sequence.

    multimapped_mode : str
        How to handle read which align with multiple features.
        Option: none, all, fraction, random

    reads_stats : ReadsStatistics object
        ReadsStatistics object which stores the statistics about the reads

    read_io_object : ReadsIO object
    """

    fs = _get_feature_set(overlap_mode, iv_seq, features)
    _classify_read_against_feature(fs, read_sequence, reads_stats, read_io_object)
    if fs is not None and len(fs) > 0:
        _assign_count(fs, read_sequence, multimapped_mode, reads_stats)
