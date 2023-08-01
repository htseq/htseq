import itertools
import random
import sys

from HTSeq.scripts.utils import invert_strand, UnknownChrom
from HTSeq.scripts.count_features.reads_io_processor import ReadsIO
from HTSeq.scripts.count_features.reads_stats import ReadsStatistics


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
    quiet,
    minaqual,
    samout_format,
    samout_filename,
):
    """
    The function that does the counting for each input BAM/SAM file.
    Fixme: there are some redundant parameters here.. feature_type, id_attribute, additional_attributes

    Parameters
    ----------
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
    quiet : boolean
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
    try:
        read_io_obj = ReadsIO(
            sam_filename=sam_filename,
            samout_filename=samout_filename,
            samout_format=samout_format,
            supplementary_alignment_mode=supplementary_alignment_mode,
            secondary_alignment_mode=secondary_alignment_mode,
            order=order,
            max_buffer_size=max_buffer_size,
        )

        # If the BAM header is available, check that at least one of the
        # chromosomes is also found in the GTF/GFF file, otherwise the user
        # is probably doing something wrong (e.g. "chr1" vs "1").
        bam_chroms = read_io_obj.get_chromosome_names_header()
        if bam_chroms is not None:
            bam_chroms = set(bam_chroms)
            feature_chroms = set(features.chrom_vectors.keys())
            if not (bam_chroms & feature_chroms):
                sys.stderr.write(
                    f"The alignment file has no chromosomes in common with the GFF/GTF "
                    "file. This will result in zero feature counts. Please check if the "
                    "references match, e.g. if you are using 'chr1' or '1' as "
                    "chromosome names.\n")

    except:
        sys.stderr.write("Error occurred when reading beginning of SAM/BAM file.\n")
        raise

    try:
        read_stats = ReadsStatistics(
            feature_attr=feature_attr, read_io_object=read_io_obj
        )
    except:
        sys.stderr.write(
            "Error occurred when preparing object to store the reads' assignments\n"
        )
        raise

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ("M", "=", "X")

    try:
        for r in read_io_obj.read_seq:
            read_stats.print_progress()
            read_stats.add_num_reads_processed()

            #  get the interval/read sequence.
            if not read_io_obj.pe_mode:
                skip_read = _assess_non_pe_read(
                    read_sequence=r,
                    read_stats=read_stats,
                    secondary_alignment_mode=secondary_alignment_mode,
                    supplementary_alignment_mode=supplementary_alignment_mode,
                    multimapped_mode=multimapped_mode,
                    minaqual=minaqual,
                )

                if skip_read:
                    continue
                iv_seq = _get_iv_seq_non_pe_read(com, r, stranded)
            else:

                # NOTE: the logic here is a little arbitrary and might benefit
                # from an optional arg. If the reads are paired-end but one of
                # the two is missing, ATM we rely on the other one for info,
                # however the data is technically inconsistent and we might
                # want to let the user choose.
                skip_read = _assess_pe_read(
                    minaqual,
                    multimapped_mode,
                    r,
                    read_stats,
                    secondary_alignment_mode,
                    supplementary_alignment_mode,
                )
                if skip_read:
                    continue

                iv_seq = _get_iv_seq_pe_read(com, r, stranded)

            # this bit updates the counts obtained from aligning reads to feature sets.
            try:
                fs = _align_reads_to_feature_set(features, iv_seq, overlap_mode)

                _update_feature_set_counts(fs, multimapped_mode, r, read_stats)

            except UnknownChrom:
                read_stats.add_empty_read(read_sequence=r)

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n"
            % (read_io_obj.read_seq_file.get_line_number_string())
        )
        raise

    if not quiet:
        read_stats.print_progress(force_print=True)

    read_io_obj.close_samoutfile()

    res = read_stats.get_output(isam)
    return res


def _update_feature_set_counts(fs, multimapped_mode, read_sequence, read_stats):
    """
    Distribute the counts among the aligned feature set.

    Parameters
    ----------
    fs : array
        A list of feature set previously aligned to the read
    multimapped_mode : str
        How to handle read mapped to multiple features
    read_sequence : array
        Read sequence
    read_stats : ReadsStatistics object
        For updating bad reads

    """
    if fs is None or len(fs) == 0:
        read_stats.add_empty_read(read_sequence=read_sequence)
    elif len(fs) > 1:
        read_stats.add_ambiguous_read(
            read_sequence=read_sequence,
            assignment="__ambiguous[" + "+".join(sorted(fs)) + "]",
        )
    else:
        read_stats.add_good_read_assignment(
            read_sequence=read_sequence, assignment=list(fs)[0]
        )
    if fs is not None and len(fs) > 0:
        fs = list(fs)
        if multimapped_mode == "none":
            if len(fs) == 1:
                read_stats.add_to_count(feature=fs[0])
        elif multimapped_mode == "all":
            for fsi in fs:
                read_stats.add_to_count(feature=fsi)
        elif multimapped_mode == "fraction":
            val = 1.0 / len(fs)
            for fsi in fs:
                read_stats.add_to_count(feature=fsi, value=val)
        elif multimapped_mode == "random":
            fsi = random.choice(fs)
            read_stats.add_to_count(feature=fsi)
        else:
            sys.exit("Illegal multimap mode.")


def _align_reads_to_feature_set(features, iv_seq, overlap_mode):
    """
    Align reads to feature set.

    Parameters
    ----------
    features : array
        A set of features to align the reads to
        TODO not sure the type yet.
    iv_seq : array
        TODO not sure the type yet.
        Read (or interval?) sequence
    overlap_mode : str
        How to select the features for read that not 100% aligned to a feature.

    Returns
    -------
    fs : array
        A set of features to align the reads to
        TODO not sure the type yet.

    """
    if overlap_mode == "union":
        fs = set()
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                raise UnknownChrom
            for iv2, fs2 in features[iv].steps():
                fs = fs.union(fs2)
    elif overlap_mode in ("intersection-strict", "intersection-nonempty"):
        fs = None
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                raise UnknownChrom
            for iv2, fs2 in features[iv].steps():
                if (len(fs2) > 0) or (overlap_mode == "intersection-strict"):
                    if fs is None:
                        fs = fs2.copy()
                    else:
                        fs = fs.intersection(fs2)
    else:
        sys.exit("Illegal overlap mode.")
    return fs


def _get_iv_seq_pe_read(com, r, stranded):
    """
    Function to break down the read sequence into intervals which will
    subsequently be processed.

    Parameters
    ----------
    com : array
        CIGAR match characters (including alignment match, sequence match, and
        sequence mismatch
    r :
        todo update type
        Read sequence
    stranded : str
        Whether the data to be aligned is from a strand-specific assay.
        Option is yes, no, reverse.
        reverse means yes with reversed strand interpretation.

    Returns
    -------
    iv_seq :
        todo update type

    """
    if r[0] is not None and r[0].aligned:
        iv_seq = _get_iv_seq_pe_read_first(com, r[0], stranded)
    else:
        iv_seq = tuple()
    if r[1] is not None and r[1].aligned:
        iv_seq = _get_iv_seq_pe_read_second(com, iv_seq, r[1], stranded)
    return iv_seq


def _assess_pe_read(
    minaqual,
    multimapped_mode,
    read_sequence,
    read_stats,
    secondary_alignment_mode,
    supplementary_alignment_mode,
):
    """
    Function to check the read for paired end.

    Parameters
    ----------
    minaqual : int
        Value denoting the MAPQ alignment quality of reads to skip.
    multimapped_mode : str
        Whether and how to score reads that are not uniquely aligned or
        ambiguously assigned to features.
        Choices: none, all, fraction, random.
    read_sequence :
        todo update type
    read_stats : ReadsStatistics object
        Object which stores a bunch of statistics about the read sequences.
    secondary_alignment_mode : str
        Whether to score secondary alignments (0x100 flag).
        Choices: score or ignore.
    supplementary_alignment_mode : str
        Whether to score supplementary alignments (0x800 flag).
        Choices: score or ignore.

    Returns
    -------

    """
    # NOTE: Sometimes read1 is None or not aligned but read2 is fine, in that
    # case we should not exclude the entire pair but rather use the interval
    # of the second read
    read1_miss = (read_sequence[0] is None) or (not read_sequence[0].aligned)
    read2_miss = (read_sequence[1] is None) or (not read_sequence[1].aligned)
    if read1_miss and read2_miss:
        read_stats.add_not_aligned_read(read_sequence=read_sequence)
        return True

    if secondary_alignment_mode == "ignore":
        if (read_sequence[0] is not None) and read_sequence[0].not_primary_alignment:
            return True
        elif (read_sequence[1] is not None) and read_sequence[1].not_primary_alignment:
            return True
    if supplementary_alignment_mode == "ignore":
        if (read_sequence[0] is not None) and read_sequence[0].supplementary:
            return True
        elif (read_sequence[1] is not None) and read_sequence[1].supplementary:
            return True
    try:
        if (
            read_sequence[0] is not None and read_sequence[0].optional_field("NH") > 1
        ) or (
            read_sequence[1] is not None and read_sequence[1].optional_field("NH") > 1
        ):
            read_stats.add_not_unique_read(read_sequence=read_sequence)
            if multimapped_mode == "none":
                return True
    except KeyError:
        pass
    if (read_sequence[0] and read_sequence[0].aQual < minaqual) or (
        read_sequence[1] and read_sequence[1].aQual < minaqual
    ):
        read_stats.add_low_quality_read(read_sequence=read_sequence)
        return True
    return False


# Get GenomicInterval for each read, whether single-end or paired-end
def _get_iv_seq_non_pe_read(com, r, stranded):
    if stranded != "reverse":
        iv_seq = (co.ref_iv for co in r.cigar if co.type in com and co.size > 0)
    else:
        iv_seq = (
            invert_strand(co.ref_iv)
            for co in r.cigar
            if (co.type in com and co.size > 0)
        )
    return iv_seq


def _get_iv_seq_pe_read_first(com, read, stranded):
    if stranded != "reverse":
        iv_seq = (co.ref_iv for co in read.cigar if co.type in com and co.size > 0)
    else:
        iv_seq = (
            invert_strand(co.ref_iv)
            for co in read.cigar
            if co.type in com and co.size > 0
        )
    return iv_seq


def _get_iv_seq_pe_read_second(com, iv_seq, read, stranded):
    if stranded != "reverse":
        iv_seq = itertools.chain(
            iv_seq,
            (
                invert_strand(co.ref_iv)
                for co in read.cigar
                if co.type in com and co.size > 0
            ),
        )
    else:
        iv_seq = itertools.chain(
            iv_seq, (co.ref_iv for co in read.cigar if co.type in com and co.size > 0)
        )
    return iv_seq


def _assess_non_pe_read(
    read_sequence,
    read_stats,
    secondary_alignment_mode,
    supplementary_alignment_mode,
    multimapped_mode,
    minaqual,
):
    if not read_sequence.aligned:
        read_stats.add_not_aligned_read(read_sequence=read_sequence)
        return True
    if (secondary_alignment_mode == "ignore") and read_sequence.not_primary_alignment:
        return True
    if (supplementary_alignment_mode == "ignore") and read_sequence.supplementary:
        return True
    try:
        if read_sequence.optional_field("NH") > 1:
            read_stats.add_not_unique_read(read_sequence=read_sequence)
            if multimapped_mode == "none":
                return True
    except KeyError:
        pass
    if read_sequence.aQual < minaqual:
        read_stats.add_low_quality_read(read_sequence=read_sequence)
        return True

    return False
