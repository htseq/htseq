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

class ReadsStatistics(object):
    """
    For storing a bunch of statistics about the read sequences.
    """

    def __init__(self, samoutfile, pe_mode, template=None):
        self.samoutfile = samoutfile
        self.template = template
        self.pe_mode = pe_mode

        self.empty = 0
        self.ambiguous = 0
        self.notaligned = 0
        self.lowqual = 0
        self.nonunique = 0
        self.num_reads_processed = 0

    def add_num_reads_processed(self):
        self.num_reads_processed += 1

    def add_empty_read(self, read_sequence):
        self.empty += 1
        self.write_to_samout(read_sequence, "__no_feature")

    def add_ambiguous_read(self, read_sequence, assignment):
        self.ambiguous += 1
        self.write_to_samout(read_sequence, assignment)

    def add_low_quality_read(self, read_sequence):
        self.lowqual += 1
        self.write_to_samout(read_sequence, "__too_low_aQual")

    def add_not_unique_read(self, read_sequence):
        self.nonunique += 1
        self.write_to_samout(read_sequence, "__alignment_not_unique")

    def add_not_aligned_read(self, read_sequence):
        self.notaligned += 1
        self.write_to_samout(read_sequence, "__not_aligned")


    def write_to_samout(self, read_sequence, assignment):
        if self.samoutfile is None:
            return
        if not self.pe_mode:
            # What the heck is this?
            updated_read_sequence = (read_sequence,)
        for read in updated_read_sequence:
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
        return updated_read_sequence

    def close_samout(self):
        if self.samoutfile is not None:
            self.samoutfile.close()

    def print_progress(self, force_print=False):
        """
        Simple function to update the progress of reads processing.

        """

        if force_print:
            do_print = True
        else:
            do_print = self.num_reads_processed > 0 and self.num_reads_processed % 100000 == 0

        if do_print:
            sys.stderr.write(
                "%d alignment record%s processed.\n" %
                (self.num_reads_processed, "s" if not self.pe_mode else " pairs"))
            sys.stderr.flush()

    def generate_output(self, isam, counts):
        return {
            'isam': isam,
            'counts': counts,
            'empty': self.empty,
            'ambiguous': self.ambiguous,
            'lowqual': self.lowqual,
            'notaligned': self.notaligned,
            'nonunique': self.nonunique,
        }


# TODO: Rename me later. Sounds inappropriate now.
def process_strand(do_invert, read_sequence):
    com = ('M', '=', 'X')
    if do_invert:
        iv_seq = (invert_strand(co.ref_iv)
                  for co in read_sequence.cigar if (co.type in com and
                                        co.size > 0))
    else:
        iv_seq = (co.ref_iv for co in read_sequence.cigar if co.type in com
                  and co.size > 0)
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
    # TODO: Move me to the object
    # com = ('M', '=', 'X')
    counts = {key: 0 for key in feature_attr}

    # To store some stats about the reads
    reads_stats = ReadsStatistics(samoutfile, pe_mode, template)

    try:
        for r in read_seq:
            if not verbose:
                reads_stats.print_progress()
            # TODO: Move me to the bottom! Shouldn't be updating the count
            # unless the process has finished
            reads_stats.add_num_reads_processed()

            if not pe_mode:
                if not r.aligned:
                    reads_stats.add_not_aligned_read(r)
                    continue

                if ((secondary_alignment_mode == 'ignore') and
                   r.not_primary_alignment):
                    continue
                if ((supplementary_alignment_mode == 'ignore') and
                   r.supplementary):
                    continue

                try:
                    if r.optional_field("NH") > 1:
                        reads_stats.add_not_unique_read(r)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    reads_stats.add_low_quality_read(r)
                    continue

                iv_seq = process_strand(do_invert = stranded == 'reverse',
                                        read_sequence = r)

            else:
                if r[0] is not None and r[0].aligned:

                    iv_seq = process_strand(do_invert = stranded == 'reverse',
                                            read_sequence = r[0])
                else:
                    iv_seq = tuple()

                # Not sure what's the deal with 0 and 1 here really..
                if r[1] is not None and r[1].aligned:
                    # TODO: maybe inappropriate variable name
                    iv_seq_next_read = process_strand(
                        do_invert = stranded != 'reverse',
                        read_sequence = r[1])
                    iv_seq = itertools.chain(iv_seq, iv_seq_next_read)
                else:
                    # So this is when the previous sequence has not been aligned.
                    if (r[0] is None) or not (r[0].aligned):
                        reads_stats.add_not_aligned_read(r)
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
                        reads_stats.add_not_unique_read(r)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if ((r[0] and r[0].aQual < minaqual) or
                   (r[1] and r[1].aQual < minaqual)):
                    reads_stats.add_low_quality_read(r)
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
                    reads_stats.add_empty_read(r)
                elif len(fs) > 1:
                    reads_stats.add_ambiguous_read(
                        read_sequence = r,
                        assignment = "__ambiguous[" + '+'.join(sorted(fs)) + "]"
                    )
                else:
                    reads_stats.write_to_samout(
                        read_sequence = r,
                        assignment = list(fs)[0]
                        )

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
                reads_stats.add_empty_read(r)

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_seq_file.get_line_number_string()))
        raise

    if not verbose:
        reads_stats.print_progress(force_print=True)

    reads_stats.close_samout()

    res = reads_stats.generate_output(isam=isam,
                                counts=counts)
    return(res)
