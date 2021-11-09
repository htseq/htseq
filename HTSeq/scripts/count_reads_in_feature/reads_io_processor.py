import itertools
import pysam
import HTSeq
import sys

class ReadsIO(object):
    """
    Object to store information and provide functionalities related to the
    input and output files.
    This include providing obtaining the parser for the input file, providing
    facility to write to the output SAM file, etc.

    Attributes
    ----------
    bam_sam_file_reader : HTSeq.BAM_Reader
        Parser for SAM/BAM/CRAM files. See __init__.py for HTSeq.
        Set by __prepare_bam_sam_file_parser based on the sam_filename parameter.

    read_seq : itertools.chain
        Containing the very first read followed by the iterator for the bam_sam_file_reader.
        Set by __prepare_bam_sam_file_parser.

    pe_mode : boolean
        Boolean denoting whether the input data is paired-end.
        Set by __prepare_bam_sam_file_parser.

    template : pysam.AlignmentFile or None
        pysam.AlignmentFile object to be used as template.
        None if no samout_filename given (i.e. samout_filename is None).
        Set by __get_outputfile_and_template

    samoutfile : pysam.AlignmentFile or None or normal file
        pysam.AlignmentFile object the output file is templated BAM or SAM file.
        None if no samout_filename given (i.e. samout_filename is None).
        Normal file if the requested output is neither SAM nor BAM file.
        Set by __get_outputfile_and_template
    """

    def __init__(self, sam_filename, supplementary_alignment_mode,
                 secondary_alignment_mode, order, samout_format,
                 samout_filename):

        """
        Initialise the ReadIO object.

        Parameters
        ----------
        sam_filename : str
            Path to the SAM/BAM file containing the mapped reads.

        order : str
            Can only be either 'pos' or 'name'. Sorting order of <alignment_file>.

        secondary_alignment_mode : str
            Whether to score secondary alignments (0x100 flag).
            Choices: score or ignore.

        supplementary_alignment_mode : str
            Whether to score supplementary alignments (0x800 flag).
            Choices: score or ignore.

        samout_format : str
            Format of the output files denoted by samouts.
            Choices: SAM, BAM, sam, bam.

        samout_filename : str
            The name of SAM/BAM file to write out all SAM alignment records into.
        """

        # Set by __prepare_bam_sam_file_parser function below.
        # Setting to None for now for clarity only.
        self.bam_sam_file_reader = None
        self.read_seq = None
        self.pe_mode = None

        self.__prepare_bam_sam_file_parser(sam_filename,
                                          supplementary_alignment_mode,
                                          secondary_alignment_mode,
                                          order)

        # Set by __get_outputfile_and_template function below.
        self.samoutfile = None
        self.template = None

        self.__get_outputfile_and_template(
            samout_filename,
            samout_format
        )

    def __prepare_bam_sam_file_parser(self,
                                     sam_filename,
                                     supplementary_alignment_mode,
                                     secondary_alignment_mode,
                                     order):
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
                self.bam_sam_file_reader = HTSeq.BAM_Reader(sys.stdin)
            else:
                self.bam_sam_file_reader = HTSeq.BAM_Reader(sam_filename)

            # See the __iter__ of the SAM_Reader object in HTSeq __init__ file.
            read_seq_iter = iter(self.bam_sam_file_reader)

            # Catch empty BAM files
            try:
                first_read = next(read_seq_iter)
                # Is this a paired-end data?
                self.pe_mode = first_read.paired_end
            # FIXME: catchall can hide subtle bugs
            except:
                # I thought it is nice to have some kind of warning..
                sys.stderr.write("Input BAM/SAM file is empty!")
                first_read = None
                self.pe_mode = False

            if first_read is not None:
                self.read_seq = itertools.chain([first_read], read_seq_iter)
            else:
                self.read_seq = []

            # What to do if paired-end data?
            if self.pe_mode:
                primary_only = supplementary_alignment_mode == 'ignore' and secondary_alignment_mode == 'ignore'

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

        except:
            sys.stderr.write(
                "Error occured when reading beginning of SAM/BAM file.\n")
            raise

    def write_to_samout(self, read_sequence, assignment):
        """
        Function to write the assignment of a read sequence to SAM/BAM output
        file.

        Parameters
        ----------
        read_sequence : array
            The read sequence to write out

        assignment : str
            The assignment of the read sequence.
            Options:
                empty: __no_feature,
                ambiguous: customised message!,
                not unique: __alignment_not_unique,
                low quality: __too_low_aQual,
                not aligned: __not_aligned
        """


        if self.samoutfile is None:
            return
        if not self.pe_mode:
            # I think this makes it such that paired end or not paired end can
            # be handled in the same manner.
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
        """
        Close the output SAM file
        """
        if self.samoutfile is not None:
            self.samoutfile.close()

    def __get_outputfile_and_template(self, samout_filename, samout_format):
        """
        Get the template and the object for the output BAM/SAM if possible.
        This was originally inside the try and catch block within the
        get_bam_sam_file_parser. Not sure why as the code inside there was just
        trying to read the input SAM and BAM file rather than preparing output
        file.
        TODO: if both the template and the samoutfile is none, what happen then?

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
        want_templated_sam_output = (samout_format in ('sam', 'SAM')) and hasattr(self.bam_sam_file_reader, 'get_template')


        if no_samout_filename_given:
            template = None
            samoutfile = None

        # Templated BAM output
        # The 'b' qualifier in 'wb' indicate a BAM file.
        elif want_templated_bam_output:
            template = self.bam_sam_file_reader.get_template()
            samoutfile = pysam.AlignmentFile(
                    samout_filename, 'wb',
                    template=template,
                    )

        # Theoretically, this and the templated BAM can be merged into 1 if, but
        # it's cleaner this way. Leave it.
        # TODO: what happen if you want untemplated SAM file? It won't reach the
        # else...
        elif want_templated_sam_output:
            template = self.bam_sam_file_reader.get_template()
            samoutfile = pysam.AlignmentFile(
                    samout_filename, 'w',
                    template=template,
                    )

        else:
            template = None
            samoutfile = open(samout_filename, 'w')

        return template, samoutfile
