class UserDefinedArguments(object):
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
        self.quiet = arg_parser_obj.quiet
        self.minaqual = arg_parser_obj.minaqual
        self.samouts = arg_parser_obj.samouts
        self.samout_format = arg_parser_obj.samout_format
        self.output_delimiter = arg_parser_obj.output_delimiter
        self.output_filename = arg_parser_obj.output_filename
        self.output_append = arg_parser_obj.output_append
        self.nprocesses = arg_parser_obj.nprocesses
        self.feature_query = arg_parser_obj.feature_query
        self.counts_output_sparse = arg_parser_obj.counts_output_sparse
