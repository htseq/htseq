# Anything to do with processing the options/arguments for htseq-count
# command is here.

import argparse

def create_input_argument():
    """
    This creates an ArgumentParser object embedded with the message for --help
    argument when running htseq-count.

    Returns
    -------
    pa : ArgumentParser
        ArgumentParser object embedded with message.
    """

    description = "This script takes one or more alignment files in SAM/BAM " +
    "format and a feature file in GFF format and calculates for each feature " +
    "the number of reads mapping to it. See " +
    "http://htseq.readthedocs.io/en/master/count.html for details."

    epilog = "Written by Simon Anders (sanders@fs.tum.de), " +
    "European Molecular Biology Laboratory (EMBL) and Fabio Zanini " +
    "(fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2020. " +
    "Released under the terms of the GNU General Public License v3. " +
    "Part of the 'HTSeq' framework, version %s." % HTSeq.__version__

    pa = argparse.ArgumentParser(
        usage = "%(prog)s [options] alignment_file gff_file",
        description = description,
        epilog = epilog)

    pa.add_argument(
            "--version", action="store_true",
            help='Show software version and exit')

    return pa

def add_options(pa):
    """
    Add all the arguments/options available for htseq-count.

    TODO: Refactor this so it's easier to maintain? If need be. Atm it is
    just very long.....

    If you update this with new argument or change any of the arguments'
    name, please update UserDefinedArguments class.

    Parameters
    ----------
    pa : argparse.ArgumentParser
        An ArgumentParser object. Generally the one returned by process_options.
    """

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
