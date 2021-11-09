import itertools
import pysam
import HTSeq
import sys


class ReadsIO(object):
    """docstring for ReadsIO."""

    def __init__(self, sam_filename, samout_filename, samout_format, supplementary_alignment_mode,
                 secondary_alignment_mode, order, max_buffer_size):

        # Set by _prepare_bam_sam_file_parser function below.
        self.pe_mode = None
        self.read_seq = None
        self.read_seq_file = None
        self.template = None
        self.samoutfile = None
        self.samout_format = samout_format

        self._prepare_samoutfile_and_read_seq(sam_filename, samout_filename, samout_format)
        self._set_read_seq(supplementary_alignment_mode, secondary_alignment_mode, order, max_buffer_size)

    def write_to_samout(self, read_sequence, assignment):
        if self.samoutfile is None:
            return
        if not self.pe_mode:
            # TODO not sure if this is good in all honesty..
            read_sequence = (read_sequence,)
        for read in read_sequence:
            if read is not None:
                read.optional_fields.append(('XF', assignment))
                if self.template is not None:
                    self.samoutfile.write(read.to_pysam_AlignedSegment(self.template))
                elif self.samout_format in ('SAM', 'sam'):
                    self.samoutfile.write(read.get_sam_line() + "\n")
                else:
                    raise ValueError(
                        'BAM/SAM output: no template and not a test SAM file',
                    )

    def close_samoutfile(self):
        if self.samoutfile is not None:
            self.samoutfile.close()

    def _prepare_samoutfile_and_read_seq(self, sam_filename, samout_filename, samout_format):
        if sam_filename == "-":
            self.read_seq_file = HTSeq.BAM_Reader(sys.stdin)
        else:
            self.read_seq_file = HTSeq.BAM_Reader(sam_filename)
        # Get template for output BAM/SAM if possible
        self._set_output_template(samout_filename, samout_format)

    def _set_read_seq(self, supplementary_alignment_mode, secondary_alignment_mode, order,
                      max_buffer_size):
        read_seq_iter = iter(self.read_seq_file)
        # Catch empty BAM files
        try:
            first_read = next(read_seq_iter)
            self.pe_mode = first_read.paired_end
        # FIXME: catchall can hide subtle bugs
        except:
            first_read = None
            self.pe_mode = False
        if first_read is not None:
            self.read_seq = itertools.chain([first_read], read_seq_iter)
        else:
            self.read_seq = []

        if self.pe_mode:
            if ((supplementary_alignment_mode == 'ignore') and
                    (secondary_alignment_mode == 'ignore')):
                primary_only = True
            else:
                primary_only = False
            if order == "name":
                self.read_seq = HTSeq.pair_SAM_alignments(
                    self.read_seq,
                    primary_only=primary_only)
            elif order == "pos":
                self.read_seq = HTSeq.pair_SAM_alignments_with_buffer(
                    self.read_seq,
                    max_buffer_size=max_buffer_size,
                    primary_only=primary_only)
            else:
                raise ValueError("Illegal order specified.")

    def _set_output_template(self, samout_filename, samout_format):
        if samout_filename is None:
            self.template = None
            self.samoutfile = None
        elif samout_format in ('bam', 'BAM'):
            self.template = self.read_seq_file.get_template()
            self.samoutfile = pysam.AlignmentFile(
                samout_filename, 'wb',
                template=self.template,
            )
        elif (samout_format in ('sam', 'SAM')) and \
                hasattr(self.read_seq_file, 'get_template'):
            self.template = self.read_seq_file.get_template()
            self.samoutfile = pysam.AlignmentFile(
                samout_filename, 'w',
                template=self.template,
            )
        else:
            self.template = None
            self.samoutfile = open(samout_filename, 'w')
