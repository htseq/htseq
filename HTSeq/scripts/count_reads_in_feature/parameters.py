import pysam

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
