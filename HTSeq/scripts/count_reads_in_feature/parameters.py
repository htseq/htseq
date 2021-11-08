import pysam
import HTSeq
import argparse
import sys

class CountParameters(object):
    """
    A class used to represent the arguments user pass for htseq-count command.

    Attributes
    ----------
    sam_filenames : str
        Path to the SAM/BAM files containing the mapped reads.
    gff_filename : str
        Path to the GTF file containing the features
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
        For each position i in the read, a set S(i) is defined as the set of
        all features overlapping position i. Then, consider the set S, which is
        (with i running through all position within the read or a read pair):
        1) union: the union of all the sets S(i).
        2) intersection-strict: the intersection of all the sets S(i).
        3) intersection-nonempty: the intersection of all non-empty sets S(i).
        See: https://htseq.readthedocs.io/en/master/htseqcount.html#htseqcount

    multimapped_mode : str
        Whether and how to count reads that are not uniquely aligned or
        ambiguously assigned to ONE feature.
        Equivalent to the --nonunique cmd line option.
        Choices:
        --------
        1) None: the read (or read pair) is counted as ambiguous and not counted for any features. Also, if the read (or read pair) aligns to more than one location in the reference, it is scored as alignment_not_unique.
        2) All: the read (or read pair) is counted as ambiguous and is also counted in all features to which it was assigned. Also, if the read (or read pair) aligns to more than one location in the reference, it is scored as alignment_not_unique and also separately for each location.
        3) Fraction: the read (or read pair) is counted as ambiguous and is also counted fractionally in all features to which it was assigned. For example, if the read overlaps with 3 features, it will be counted 1/3 to each of them.
        4) Random: the read (or read pair) is counted as ambiguous and is also counted uniformly at random to one of the features to which it was assigned.
        See: https://htseq.readthedocs.io/en/master/htseqcount.html#htseqcount

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
    samouts : array
        The SAM/BAM files' name to write out all SAM alignment records into.
    samout_format : str
        Format of the output files denoted by samouts.
        Choices: SAM, BAM, sam, bam.
    output_delimiter : str
        Column delimiter in the output files.
    output_filename : str
        Filename to output the counts to instead.
    output_append : boolean
        Append counts output to an existing file instead of creating a new one.
    nprocesses : int
        Number of parallel CPU processes to use.
    feature_query : str
        Restrict to features descibed in this expression.
    counts_output_sparse : boolean
        Whether to store the counts as a sparse matrix (mtx, h5ad, loom).
    """

    def __init__(self, arg_parser_obj):
        self.sam_filenames = arg_parser_obj.samfilenames
        self.gff_filename = arg_parser_obj.featuresfilename
        self.order = arg_parser_obj.order
        self.max_buffer_size = arg_parser_obj.max_buffer_size
        self.stranded = arg_parser_obj.stranded
        self.overlap_mode = arg_parser_obj.mode
        self.multimapped_mode = arg_parser_obj.nonunique
        self.secondary_alignment_mode = arg_parser_obj.secondary_alignments
        self.supplementary_alignment_mode = arg_parser_obj.supplementary_alignments
        self.feature_type = arg_parser_obj.featuretype
        self.id_attribute = arg_parser_obj.idattr
        self.additional_attributes = arg_parser_obj.additional_attr
        self.add_chromosome_info = arg_parser_obj.add_chromosome_info
        self.verbose = arg_parser_obj.quiet
        self.minaqual = arg_parser_obj.minaqual
        self.samouts = arg_parser_obj.samouts
        self.samout_format = arg_parser_obj.samout_format
        self.output_delimiter = arg_parser_obj.output_delimiter
        self.output_filename = arg_parser_obj.output_filename
        self.output_append = arg_parser_obj.output_append
        self.nprocesses = arg_parser_obj.nprocesses
        self.feature_query = arg_parser_obj.feature_query
        self.counts_output_sparse = arg_parser_obj.counts_output_sparse

        # just run these here rather than at the start of do_count_reads_in_features
        self.check_samouts_integrity()
        self.check_sam_filenames_integrity()

        # this will be filled by prepare_features
        self.features = None
        self.feature_attr = None
        self.attributes = None

    def check_sam_filenames_integrity(self):
        if (len(self.sam_filenames) != 1) or (self.sam_filenames[0] != '-'):
            for sam_filename in self.sam_filenames:
                with pysam.AlignmentFile(sam_filename, 'r') as sf:
                    pass

    def check_samouts_integrity(self):
        if self.samouts != []:
            if len(samouts) != len(sam_filenames):
                raise ValueError(
                        'Select the same number of input and output files')

            # Try to open samout files early in case any of them has issues.
            if self.samout_format.lower() == 'sam':
                for samout in self.samouts:
                    with open(samout, 'w'):
                        pass
            else:
                # We don't have a template if the input is stdin
                # The "-" check is to ensure that the input is not stdin
                # as "-" implies stdin.
                # I think this is run if the samout files are not SAM format?
                # So basically "converting" it to SAM format.
                if (len(self.sam_filenames) != 1) or (self.sam_filenames[0] != '-'):
                    for sam_filename, samout in zip(self.sam_filenames, self.samouts):
                        with pysam.AlignmentFile(sam_filename, 'r') as sf:
                            with pysam.AlignmentFile(samout, 'w', template=sf):
                                pass
        else:
            self.samouts = [None for x in self.sam_filenames]

    def prepare_features(self):
        """
        Not quite sure what this is yet. It says prepare features, but for what?
        Anyway, only calling make_feature_genomicarrayofsets...

        Parameters
        ----------
        in_param : CountParameters
            Custom object (see parameters module) which stores the
            input parameters given by user to run htseq-count.
        """
        # Prepare features
        gff = HTSeq.GFF_Reader(self.gff_filename)
        feature_scan = HTSeq.make_feature_genomicarrayofsets(
            gff,
            self.id_attribute,
            feature_type = self.feature_type,
            feature_query = self.feature_query,
            additional_attributes = self.additional_attributes,
            stranded = self.stranded != 'no',
            verbose = not self.verbose,
            add_chromosome_info = self.add_chromosome_info,
            )
        self.features = feature_scan['features']
        self.attributes = feature_scan['attributes']
        self.feature_attr = sorted(self.attributes.keys())

        if len(self.feature_attr) == 0:
            sys.stderr.write(
                "Warning: No features of type '%s' found.\n" % self.feature_type)

    def get_args_for_count(self):
        """
        This prepares all the arguments required to run the counting function.
        Each element contains the parameters to be passed to an instance of counting function.

        Returns
        -------
        args : array
            Array containing the arguments for the counting function.
        """

        args = []
        for isam_file_idx, (sam_filename, samout_filename) in enumerate(zip(self.sam_filenames, self.samouts)):
            args.append((
                isam_file_idx,
                sam_filename,
                self.features,
                self.feature_attr,
                self.order,
                self.max_buffer_size,
                self.stranded,
                self.overlap_mode,
                self.multimapped_mode,
                self.secondary_alignment_mode,
                self.supplementary_alignment_mode,
                self.feature_type,
                self.id_attribute,
                self.additional_attributes,
                self.verbose,
                self.minaqual,
                self.samout_format,
                samout_filename,
                ))
        return args

def process_cmd_line_options():
    """
    This creates an ArgumentParser object embedded with the message for --help
    argument when running htseq-count.

    If --version flag is passed, it will show the version of htseq installed
    and exit.

    Otherwise, it will dd all the command-line options available for htseq-
    count.
    If you update this with new command-line option or change any of their
    name, please update the CountParameters class above.

    Returns
    -------
    args : Object
        An object returned by ArgumentParser.parse_args() function.
        Each command-line option is assigned as the attributes of this object.
    """

    description = """This script takes one or more alignment files in SAM/BAM
    format and a feature file in GFF format and calculates for each feature the
    number of reads mapping to it. See http://htseq.readthedocs.io/en/master/
    count.html for details."""

    epilog = """Written by Simon Anders (sanders@fs.tum.de),
    European Molecular Biology Laboratory (EMBL) and Fabio Zanini
    (fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2020.
    Released under the terms of the GNU General Public License v3.
    Part of the 'HTSeq' framework, version %s.""" % HTSeq.__version__

    pa = argparse.ArgumentParser(
        usage = "%(prog)s [options] alignment_file gff_file",
        description = description,
        epilog = epilog)

    pa.add_argument(
            "--version", action="store_true",
            help='Show software version and exit')

    args, argv = pa.parse_known_args()
    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    pa.add_argument(
            "samfilenames", nargs='+', type=str,
            help="Path to the SAM/BAM files containing the mapped reads. " +
            "If '-' is selected, read from standard input")

    pa.add_argument(
            "featuresfilename", type=str,
            help="Path to the GTF file containing the features")

    pa.add_argument(
            "-f", "--format", dest="samtype",
            choices=("sam", "bam", "auto"), default="auto",
            help="Type of <alignment_file> data. DEPRECATED: " +
            "file format is detected automatically. This option is ignored.")

    pa.add_argument(
            "-r", "--order", dest="order",
            choices=("pos", "name"), default="name",
            help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
            "data must be sorted either by position or by read name, and the sorting order " +
            "must be specified. Ignored for single-end data.")

    pa.add_argument(
            "--max-reads-in-buffer", dest="max_buffer_size", type=int,
            default=30000000,
            help="When <alignment_file> is paired end sorted by position, " +
            "allow only so many reads to stay in memory until the mates are " +
            "found (raising this number will use more memory). Has no effect " +
            "for single end or paired end sorted by name")

    pa.add_argument(
            "-s", "--stranded", dest="stranded",
            choices=("yes", "no", "reverse"), default="yes",
            help="Whether the data is from a strand-specific assay. Specify 'yes', " +
            "'no', or 'reverse' (default: yes). " +
            "'reverse' means 'yes' with reversed strand interpretation")

    pa.add_argument(
            "-a", "--minaqual", type=int, dest="minaqual",
            default=10,
            help="Skip all reads with MAPQ alignment quality lower than the given " +
            "minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM " +
            "file and its usage depends on the software used to map the reads.")

    pa.add_argument(
            "-t", "--type", type=str, dest="featuretype",
            default="exon",
            help="Feature type (3rd column in GTF file) to be used, " +
            "all features of other type are ignored (default, suitable for Ensembl " +
            "GTF files: exon)")

    pa.add_argument(
            "-i", "--idattr", type=str, dest="idattr",
            default="gene_id",
            help="GTF attribute to be used as feature ID (default, " +
            "suitable for Ensembl GTF files: gene_id). All feature of the " +
            "right type (see -t option) within the same GTF attribute will " +
            "be added together. The typical way of using this option is to " +
            "count all exonic reads from each gene and add the exons " +
            "but other uses are possible as well.")

    pa.add_argument(
            "--additional-attr", type=str,
            action='append',
            default=[],
            help="Additional feature attributes (default: none, " +
            "suitable for Ensembl GTF files: gene_name). Use multiple times " +
            "for more than one additional attribute. These attributes are " +
            "only used as annotations in the output, while the determination " +
            "of how the counts are added together is done based on option -i.")

    pa.add_argument(
            "--add-chromosome-info", action='store_true',
            help="Store information about the chromosome of each feature as " +
            "an additional attribute (e.g. colunm in the TSV output file).",
            )

    pa.add_argument(
            "-m", "--mode", dest="mode",
            choices=("union", "intersection-strict", "intersection-nonempty"),
            default="union",
            help="Mode to handle reads overlapping more than one feature " +
            "(choices: union, intersection-strict, intersection-nonempty; default: union)")

    pa.add_argument(
            "--nonunique", dest="nonunique", type=str,
            choices=("none", "all", "fraction", "random"), default="none",
            help="Whether and how to score reads that are not uniquely aligned " +
            "or ambiguously assigned to features " +
            "(choices: none, all, fraction, random; default: none)")

    pa.add_argument(
            "--secondary-alignments", dest="secondary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score secondary alignments (0x100 flag)")

    pa.add_argument(
            "--supplementary-alignments", dest="supplementary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score supplementary alignments (0x800 flag)")

    pa.add_argument(
            "-o", "--samout", type=str, dest="samouts",
            action='append',
            default=[],
            help="Write out all SAM alignment records into " +
            "SAM/BAM files (one per input file needed), annotating each line " +
            "with its feature assignment (as an optional field with tag 'XF')" +
            ". See the -p option to use BAM instead of SAM.")

    pa.add_argument(
            "-p", '--samout-format', type=str, dest='samout_format',
            choices=('SAM', 'BAM', 'sam', 'bam'), default='SAM',
            help="Format to use with the --samout option."
            )

    pa.add_argument(
            "-d", '--delimiter', type=str, dest='output_delimiter',
            default='\t',
            help="Column delimiter in output (default: TAB)."
            )
    pa.add_argument(
            "-c", '--counts_output', type=str, dest='output_filename',
            default='',
            help="Filename to output the counts to instead of stdout."
            )

    pa.add_argument(
            "--counts-output-sparse", action='store_true',
            help="Store the counts as a sparse matrix (mtx, h5ad, loom)."
            )

    pa.add_argument(
            '--append-output', action='store_true', dest='output_append',
            help='Append counts output to an existing file instead of ' +
            'creating a new one. This option is useful if you have ' +
            'already creates a TSV/CSV/similar file with a header for your ' +
            'samples (with additional columns for the feature name and any ' +
            'additionl attributes) and want to fill in the rest of the file.'
            )

    pa.add_argument(
            "-n", '--nprocesses', type=int, dest='nprocesses',
            default=1,
            help="Number of parallel CPU processes to use (default: 1). " +
            "This option is useful to process several input files at once. " +
            "Each file will use only 1 CPU. It is possible, of course, to " +
            "split a very large input SAM/BAM files into smaller chunks " +
            "upstream to make use of this option."
            )

    pa.add_argument(
            '--feature-query', type=str, dest='feature_query',
            default=None,
            help='Restrict to features descibed in this expression. Currently ' +
            'supports a single kind of expression: attribute == "one attr" to ' +
            'restrict the GFF to a single gene or transcript, e.g. ' +
            '--feature-query \'gene_name == "ACTB"\' - notice the single ' +
            'quotes around the argument of this option and the double ' +
            'quotes around the gene name. Broader queries might become ' +
            'available in the future.',
            )

    pa.add_argument(
            "-q", "--quiet", action="store_true", dest="quiet",
            help="Suppress progress report")  # and warnings" )

    args = pa.parse_args()

    return args
