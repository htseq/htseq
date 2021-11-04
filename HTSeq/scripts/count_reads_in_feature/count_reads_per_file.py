import itertools
import pysam
import random
import HTSeq
import sys

from HTSeq.scripts.utils import (
    UnknownChrom,
    invert_strand
)

def get_outputfile_and_template(
        samout_filename,
        samout_format,
        bam_sam_file_reader):
    """
    Get the template and the object for the output BAM/SAM if possible.

    Parameters
    ----------
    samout_format : str
        Format of the output files denoted by samouts.
        Choices: SAM, BAM, sam, bam.
    samout_filename : str
        The name of SAM/BAM file to write out all SAM alignment records into.
    bam_sam_file_reader : HTSeq.BAM_Reader
        Parser for SAM/BAM/CRAM files. See __init__.py for HTSeq.

    Returns
    -------
    template : pysam.AlignmentFile or None
        pysam.AlignmentFile object to be used as template.
        None if no samout_filename given (i.e. samout_filename is None) or if output
    samoutfile : pysam.AlignmentFile or None or normal file
        pysam.AlignmentFile object the output file is templated BAM or SAM file.
        None if no samout_filename given (i.e. samout_filename is None).
        Normal file if the requested output is neither SAM nor BAM file.

    """

    no_samout_filename_given = samout_filename is None
    want_templated_bam_output = samout_format in ('bam', 'BAM')
    # The latter check if we can produce a template
    want_templated_sam_output = (samout_format in ('sam', 'SAM')) and hasattr(bam_sam_file_reader, 'get_template')


    if no_samout_filename_given:
        template = None
        samoutfile = None

    # Templated BAM output
    # The 'b' qualifier in 'wb' indicate a BAM file.
    elif want_templated_bam_output:
        template = bam_sam_file_reader.get_template()
        samoutfile = pysam.AlignmentFile(
                samout_filename, 'wb',
                template=template,
                )

    # Theoretically, this and the templated BAM can be merged into 1 if, but
    # it's cleaner this way. Leave it.
    # TODO: what happen if you want untemplated SAM file? It won't reach the
    # else...
    elif want_templated_sam_output:
        template = bam_sam_file_reader.get_template()
        samoutfile = pysam.AlignmentFile(
                samout_filename, 'w',
                template=template,
                )

    else:
        template = None
        samoutfile = open(samout_filename, 'w')

    return template, samoutfile


def prepare_bam_sam_file_parser(
    sam_filename,
    supplementary_alignment_mode,
    secondary_alignment_mode,
    order
    ):
    """
    Prepare the BAM/SAM file parser.
    This will create the parser and prepare an iterator for it.
    Depending on whether we have paired-end reads or not, different iterator
    will be returned.

    Parameters
    ----------
    samout_filename : str
        The name of SAM/BAM file to write out all SAM alignment records into.
    secondary_alignment_mode : str
        Whether to score secondary alignments (0x100 flag).
        Choices: score or ignore.
    supplementary_alignment_mode : str
        Whether to score supplementary alignments (0x800 flag).
        Choices: score or ignore.
    order : str
        Can only be either 'pos' or 'name'. Sorting order of <alignment_file>.

    Returns
    -------
    bam_sam_file_reader : HTSeq.BAM_Reader
        Parser for SAM/BAM/CRAM files. See __init__.py for HTSeq.
    read_seq : itertools.chain
        Containing the very first read followed by the iterator for the bam_sam_file_reader.
        If the input file is empty, this will be an empty array.
    pe_mode : boolean
        Is this a paired-end data?
    """
    try:
        if sam_filename == "-":
            # BAM_Reader is in HTSeq __init__ file.
            bam_sam_file_reader = HTSeq.BAM_Reader(sys.stdin)
        else:
            bam_sam_file_reader = HTSeq.BAM_Reader(sam_filename)

        # See the __iter__ of the SAM_Reader object in HTSeq __init__ file.
        read_seq_iter = iter(bam_sam_file_reader)

        # Catch empty BAM files
        try:
            first_read = next(read_seq_iter)
            # Is this a paired-end data?
            pe_mode = first_read.paired_end
        # FIXME: catchall can hide subtle bugs
        except:
            # I thought it is nice to have some kind of warning..
            sys.stderr.write("Input BAM/SAM file is empty!")
            first_read = None
            pe_mode = False

        if first_read is not None:
            read_seq = itertools.chain([first_read], read_seq_iter)
        else:
            read_seq = []

        # What to do if paired-end data?
        if pe_mode:
            primary_only = supplementary_alignment_mode == 'ignore' and secondary_alignment_mode == 'ignore'

            if order == "name":
                read_seq = HTSeq.pair_SAM_alignments(
                        read_seq,
                        primary_only=primary_only)
            elif order == "pos":
                read_seq = HTSeq.pair_SAM_alignments_with_buffer(
                        read_seq,
                        max_buffer_size=max_buffer_size,
                        primary_only=primary_only)
            else:
                raise ValueError("Illegal order specified.")

    except:
        sys.stderr.write(
            "Error occured when reading beginning of SAM/BAM file.\n")
        raise

    return bam_sam_file_reader, read_seq, pe_mode

def update_progress(num_reads_processed, pe_mode):
    """
    Simple function to update the progress of reads processing.

    Parameters
    ----------
    num_reads_processed : int
        Number of alignment records/pairs processed.
    pe_mode : boolean
        Whether the data is paired-end.
    """
    if num_reads_processed > 0 and num_reads_processed % 100000 == 0:
        sys.stderr.write(
            "%d alignment record%s processed.\n" %
            (num_reads_processed, "s" if not pe_mode else " pairs"))
        sys.stderr.flush()


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
    The one that does the counting.

    You can technically get this to read the parameters object by just extending
    it to have a field that say which sam_file to process.
    But then the input and output files array will be duplicated several times,
    very wasteful...

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

    def write_to_samout(r, assignment, samoutfile, template=None):
        if samoutfile is None:
            return
        if not pe_mode:
            # pe_mode is set somewhere in the outer function.
            r = (r,)
        for read in r:
            if read is not None:
                read.optional_fields.append(('XF', assignment))
                if template is not None:
                    samoutfile.write(read.to_pysam_AlignedSegment(template))
                elif samout_format in ('SAM', 'sam'):
                    samoutfile.write(read.get_sam_line() + "\n")
                else:
                    raise ValueError(
                        'BAM/SAM output: no template and not a test SAM file',
                    )

    read_seq_file, read_seq, pe_mode = prepare_bam_sam_file_parser(
        sam_filename,
        supplementary_alignment_mode,
        secondary_alignment_mode,
        order
        )

    # This was originally inside the try and catch block within the
    # get_bam_sam_file_parser. Not sure why as the code inside there was just
    # trying to read the input SAM and BAM file rather than preparing output
    # file.
    # TODO: if both the template and the samoutfile is none, what happen then?
    template, samoutfile = get_outputfile_and_template(samout_filename,
                                            samout_format,
                                            read_seq_file)

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')
    counts = {key: 0 for key in feature_attr}

    try:
        empty = 0
        ambiguous = 0
        notaligned = 0
        lowqual = 0
        nonunique = 0
        # Just for progress update
        num_reads_processed = 0

        for r in read_seq:
            if not verbose:
                update_progress(num_reads_processed=i, pe_mode=pe_mode)
            # TODO: Move me to the bottom! Shouldn't be updating the count
            # unless the process has finished
            num_reads_processed += 1

            if not pe_mode:
                if not r.aligned:
                    notaligned += 1
                    write_to_samout(
                            r, "__not_aligned", samoutfile,
                            template)
                    continue
                if ((secondary_alignment_mode == 'ignore') and
                   r.not_primary_alignment):
                    continue
                if ((supplementary_alignment_mode == 'ignore') and
                   r.supplementary):
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        nonunique += 1
                        write_to_samout(
                                r,
                                "__alignment_not_unique",
                                samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    lowqual += 1
                    write_to_samout(
                            r, "__too_low_aQual", samoutfile,
                            template)
                    continue
                if stranded != "reverse":
                    iv_seq = (co.ref_iv for co in r.cigar if co.type in com
                              and co.size > 0)
                else:
                    iv_seq = (invert_strand(co.ref_iv)
                              for co in r.cigar if (co.type in com and
                                                    co.size > 0))
            else:
                if r[0] is not None and r[0].aligned:
                    if stranded != "reverse":
                        iv_seq = (co.ref_iv for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                    else:
                        iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                else:
                    iv_seq = tuple()
                if r[1] is not None and r[1].aligned:
                    if stranded != "reverse":
                        iv_seq = itertools.chain(
                                iv_seq,
                                (invert_strand(co.ref_iv) for co in r[1].cigar
                                if co.type in com and co.size > 0))
                    else:
                        iv_seq = itertools.chain(
                                iv_seq,
                                (co.ref_iv for co in r[1].cigar
                                 if co.type in com and co.size > 0))
                else:
                    if (r[0] is None) or not (r[0].aligned):
                        write_to_samout(
                                r, "__not_aligned", samoutfile,
                                template)
                        notaligned += 1
                        continue
                if secondary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].not_primary_alignment:
                        continue
                    elif (r[1] is not None) and r[1].not_primary_alignment:
                        continue
                if supplementary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].supplementary:
                        continue
                    elif (r[1] is not None) and r[1].supplementary:
                        continue
                try:
                    if ((r[0] is not None and r[0].optional_field("NH") > 1) or
                       (r[1] is not None and r[1].optional_field("NH") > 1)):
                        nonunique += 1
                        write_to_samout(
                                r, "__alignment_not_unique", samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if ((r[0] and r[0].aQual < minaqual) or
                   (r[1] and r[1].aQual < minaqual)):
                    lowqual += 1
                    write_to_samout(
                            r, "__too_low_aQual", samoutfile,
                            template)
                    continue

            try:
                if overlap_mode == "union":
                    fs = set()
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            fs = fs.union(fs2)
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

                if fs is None or len(fs) == 0:
                    write_to_samout(
                            r, "__no_feature", samoutfile,
                            template)
                    empty += 1
                elif len(fs) > 1:
                    write_to_samout(
                            r, "__ambiguous[" + '+'.join(sorted(fs)) + "]",
                            samoutfile,
                            template)
                    ambiguous += 1
                else:
                    write_to_samout(
                            r, list(fs)[0], samoutfile,
                            template)

                if fs is not None and len(fs) > 0:
                    if multimapped_mode == 'none':
                        if len(fs) == 1:
                            counts[list(fs)[0]] += 1
                    elif multimapped_mode == 'all':
                        for fsi in list(fs):
                            counts[fsi] += 1
                    elif multimapped_mode == 'fraction':
                        for fsi in list(fs):
                            counts[fsi] += 1.0 / len(fs)
                    elif multimapped_mode == 'random':
                        fsi = random.choice(list(fs))
                        counts[fsi] += 1
                    else:
                        sys.exit("Illegal multimap mode.")

            except UnknownChrom:
                write_to_samout(
                        r, "__no_feature", samoutfile,
                        template)
                empty += 1

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_seq_file.get_line_number_string()))
        raise

    if not verbose:
        sys.stderr.write(
            "%d %s processed.\n" %
            (i, "alignments " if not pe_mode else "alignment pairs"))
        sys.stderr.flush()

    if samoutfile is not None:
        samoutfile.close()

    return {
        'isam': isam,
        'counts': counts,
        'empty': empty,
        'ambiguous': ambiguous,
        'lowqual': lowqual,
        'notaligned': notaligned,
        'nonunique': nonunique,
    }
