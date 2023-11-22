import sys
import argparse
import operator
import itertools
import warnings
import traceback
import os.path
import multiprocessing
import numpy as np
import pysam

import HTSeq
from HTSeq.scripts.count_features.count_features_per_file import count_reads_single_file
from HTSeq.scripts.utils import (
    my_showwarning,
    _write_output,
)


def count_reads_in_features(args):
    """Count reads in features, parallelizing by file

    Args:
        args: ArgumentParser args, i.e. each argument is a property of this
          instance. Check the CLI parsing function below for a full list
          of properties, i.e. command-line options.

    This function can be conceptually split into the following steps:

    1. Load features from GTF file into memory
    2. Parse the reads from each BAM file in parallel or series
    3. Write output table

    Step 2 can be further split into two main components:

    1. Find what features overlap with each read/read pair
    2. Assign that read/pair to a feature (if unique) or a corner case
       (e.g. multimappers)
    """

    # Load feature GenomicArrayOfSets to mark overlaps
    gff = HTSeq.GFF_Reader(args.featuresfilename)
    feature_scan = HTSeq.make_feature_genomicarrayofsets(
        gff,
        args.idattr,
        feature_type=args.feature_type,
        feature_query=args.feature_query,
        additional_attributes=args.additional_attributes,
        stranded=args.stranded != "no",
        verbose=not args.quiet,
        add_chromosome_info=args.add_chromosome_info,
    )
    features = feature_scan["features"]
    attributes = feature_scan["attributes"]
    feature_attr = sorted(attributes.keys())
    if len(feature_attr) == 0:
        sys.stderr.write("Warning: No features of type '{args.feature_type}' found.\n")

    # Count reads in parallel or in series
    count_args, attributes = _prepare_args_for_counting(
        features,
        feature_attr,
        attributes,
        args.add_chromosome_info,
        args.additional_attributes,
        args.feature_query,
        args.feature_type,
        args.featuresfilename,
        args.idattr,
        args.max_buffer_size,
        args.minaqual,
        args.nonunique,
        args.order,
        args.mode,
        args.quiet,
        args.samfilenames,
        args.samout_format,
        args.samouts,
        args.secondary_alignments,
        args.stranded,
        args.supplementary_alignments,
    )
    if args.nprocesses > 1:
        with multiprocessing.Pool(args.nprocesses) as pool:
            results = pool.starmap(count_reads_single_file, count_args)
        results.sort(key=operator.itemgetter("isam"))
    else:
        results = list(itertools.starmap(count_reads_single_file, count_args))

    # Merge and write output
    _write_output(
        results,
        args.samfilenames,
        attributes,
        args.additional_attributes,
        args.output_filename,
        args.output_delimiter,
        args.output_append,
        sparse=args.counts_output_sparse,
        dtype=np.float32,
        add_tsv_header=args.with_header
    )


def _prepare_args_for_counting(
    features,
    feature_attr,
    attributes,
    add_chromosome_info,
    additional_attributes,
    feature_query,
    feature_type,
    gff_filename,
    id_attribute,
    max_buffer_size,
    minaqual,
    multimapped_mode,
    order,
    overlap_mode,
    quiet,
    sam_filenames,
    samout_format,
    samouts,
    secondary_alignment_mode,
    stranded,
    supplementary_alignment_mode,
):
    args = []
    for isam, (sam_filename, samout_filename) in enumerate(zip(sam_filenames, samouts)):
        args.append(
            (
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
            )
        )
    return args, attributes


def _check_sam_files(sam_filenames):
    if (len(sam_filenames) != 1) or (sam_filenames[0] != "-"):
        for sam_filename in sam_filenames:
            with pysam.AlignmentFile(sam_filename, "r") as sf:
                pass


def _check_samouts(sam_filenames, samout_format, samouts):
    if len(samouts) != len(sam_filenames):
        raise ValueError("Select the same number of input and output files")
    # Try to open samout files early in case any of them has issues
    if samout_format in ("SAM", "sam"):
        for samout in samouts:
            with open(samout, "w"):
                pass
    else:
        # We don't have a template if the input is stdin
        if (len(sam_filenames) != 1) or (sam_filenames[0] != "-"):
            for sam_filename, samout in zip(sam_filenames, samouts):
                with pysam.AlignmentFile(sam_filename, "r") as sf:
                    with pysam.AlignmentFile(samout, "w", template=sf):
                        pass



# Adapted from: https://github.com/python/cpython/issues/60603
class OverwriteUniqueAppendAction(argparse.Action):
    """Custom action to append unique values to a list, overwriting the default.

    When using the `append` action, the default value is not removed
    from the list. This problem is described in
    https://github.com/python/cpython/issues/60603

    This custom action aims to fix this problem by removing the default
    value when the argument is specified for the first time.

    Moreover, it only appends if the value is not already there, so the resulting
    list has unique elements.
    """

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        """Initialize the action."""
        self.called_times = 0
        self.default_value = kwargs.get("default")
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        """When the argument is specified on the commandline."""
        current_values = getattr(namespace, self.dest)

        if self.called_times == 0 and current_values == self.default_value:
            current_values = []

        # Only add if not already present (unique values)
        if values not in current_values:
            current_values.append(values)

        setattr(namespace, self.dest, current_values)
        self.called_times += 1


def _parse_sanitize_cmdline_arguments():
    pa = argparse.ArgumentParser(
        add_help=False,
    )
    pa.add_argument(
        "--version", action="store_true", help="Show software version and exit"
    )
    args, argv = pa.parse_known_args()

    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    pa = argparse.ArgumentParser(
        parents=[pa],
        description="This script takes one or more alignment files in SAM/BAM "
        + "format and a feature file in GFF format and calculates for each feature "
        + "the number of reads mapping to it. See "
        + "http://htseq.readthedocs.io/en/master/count.html for details.",
        epilog="Written by Simon Anders (sanders@fs.tum.de), "
        + "European Molecular Biology Laboratory (EMBL), Givanna Putri "
        + "(g.putri@unsw.edu.au) and Fabio Zanini "
        + "(fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2021. "
        + "Released under the terms of the GNU General Public License v3. "
        + "Please cite the following paper if you use this script: \n"
        + "    G. Putri et al. Analysing high-throughput sequencing data in "
        + "Python with HTSeq 2.0. Bioinformatics (2022). "
        + "https://doi.org/10.1093/bioinformatics/btac166.\n"
        + "Part of the 'HTSeq' framework, version %s." % HTSeq.__version__,
    )

    pa.add_argument(
        "samfilenames",
        nargs="+",
        type=str,
        help="Path to the SAM/BAM files containing the mapped reads. "
        + "If '-' is selected, read from standard input",
    )
    pa.add_argument(
        "featuresfilename",
        type=str,
        help="Path to the GTF file containing the features",
    )
    pa.add_argument(
        "-f",
        "--format",
        dest="samtype",
        choices=("sam", "bam", "auto"),
        default="auto",
        help="Type of <alignment_file> data. DEPRECATED: "
        + "file format is detected automatically. This option is ignored.",
    )
    pa.add_argument(
        "-r",
        "--order",
        dest="order",
        choices=("pos", "name"),
        default="name",
        help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing "
        + "data must be sorted either by position or by read name, and the sorting order "
        + "must be specified. Ignored for single-end data.",
    )
    pa.add_argument(
        "--max-reads-in-buffer",
        dest="max_buffer_size",
        type=int,
        default=30000000,
        help="When <alignment_file> is paired end sorted by position, "
        + "allow only so many reads to stay in memory until the mates are "
        + "found (raising this number will use more memory). Has no effect "
        + "for single end or paired end sorted by name",
    )
    pa.add_argument(
        "-s",
        "--stranded",
        dest="stranded",
        choices=("yes", "no", "reverse"),
        default="yes",
        help="Whether the data is from a strand-specific assay. Specify 'yes', "
        + "'no', or 'reverse' (default: yes). "
        + "'reverse' means 'yes' with reversed strand interpretation",
    )
    pa.add_argument(
        "-a",
        "--minaqual",
        type=int,
        dest="minaqual",
        default=10,
        help="Skip all reads with MAPQ alignment quality lower than the given "
        + "minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM "
        + "file and its usage depends on the software used to map the reads.",
    )
    pa.add_argument(
        "-t",
        "--type",
        type=str,
        dest="feature_type",
        action=OverwriteUniqueAppendAction,
        default=["exon"],
        help="Feature type (3rd column in GTF file) to be used, "
        + "all features of other type are ignored (default, suitable for"
        + "Ensembl GTF files: exon). If you can call this option multiple times, "
        + "features of all specified types will be included, e.g. to include "
        + "both genes and pseudogenes you might use -t gene -t pseudogene. "
        + "Calling this option multiple times is a rare need and might result "
        + "in excessive numbers of ambiguous counts: only use if you know what "
        + "you are doing.",
    )
    pa.add_argument(
        "-i",
        "--idattr",
        type=str,
        dest="idattr",
        action=OverwriteUniqueAppendAction,
        default=["gene_id"],
        help="GTF attribute to be used as feature ID (default, "
        + "suitable for Ensembl GTF files: gene_id). All feature of the "
        + "right type (see -t option) within the same GTF attribute will "
        + "be added together. The typical way of using this option is to "
        + "count all exonic reads from each gene and add the exons "
        + "but other uses are possible as well. You can call this option "
        + "multiple times: in that case, the combination of all attributes "
        + "separated by colons (:) will be used as a unique identifier, "
        + "e.g. for exons you might use -i gene_id -i exon_number.",
    )
    pa.add_argument(
        "--additional-attr",
        type=str,
        action='append',
        dest='additional_attributes',
        default=[],
        help="Additional feature attributes (default: none, "
        + "suitable for Ensembl GTF files: gene_name). Use multiple times "
        + "for more than one additional attribute. These attributes are "
        + "only used as annotations in the output, while the determination "
        + "of how the counts are added together is done based on option -i.",
    )
    pa.add_argument(
        "--add-chromosome-info",
        action="store_true",
        help="Store information about the chromosome of each feature as "
        + "an additional attribute (e.g. colunm in the TSV output file).",
    )
    pa.add_argument(
        "-m",
        "--mode",
        dest="mode",
        choices=("union", "intersection-strict", "intersection-nonempty"),
        default="union",
        help="Mode to handle reads overlapping more than one feature "
        + "(choices: union, intersection-strict, intersection-nonempty; default: union)",
    )
    pa.add_argument(
        "--nonunique",
        dest="nonunique",
        type=str,
        choices=("none", "all", "fraction", "random"),
        default="none",
        help="Whether and how to score reads that are not uniquely aligned "
        + "or ambiguously assigned to features "
        + "(choices: none, all, fraction, random; default: none)",
    )
    pa.add_argument(
        "--secondary-alignments",
        dest="secondary_alignments",
        type=str,
        choices=("score", "ignore"),
        default="ignore",
        help="Whether to score secondary alignments (0x100 flag)",
    )
    pa.add_argument(
        "--supplementary-alignments",
        dest="supplementary_alignments",
        type=str,
        choices=("score", "ignore"),
        default="ignore",
        help="Whether to score supplementary alignments (0x800 flag)",
    )
    pa.add_argument(
        "-o",
        "--samout",
        type=str,
        dest="samouts",
        action='append',
        default=[],
        help="Write out all SAM alignment records into "
        + "SAM/BAM files (one per input file needed), annotating each line "
        + "with its feature assignment (as an optional field with tag 'XF')"
        + ". See the -p option to use BAM instead of SAM.",
    )
    pa.add_argument(
        "-p",
        "--samout-format",
        type=str,
        dest="samout_format",
        choices=("SAM", "BAM", "sam", "bam"),
        default="SAM",
        help="Format to use with the --samout option.",
    )
    pa.add_argument(
        "-d",
        "--delimiter",
        type=str,
        dest="output_delimiter",
        default="\t",
        help="Column delimiter in output (default: TAB).",
    )
    pa.add_argument(
        "-c",
        "--counts_output",
        type=str,
        dest="output_filename",
        default="",
        help="Filename to output the counts to instead of stdout.",
    )
    pa.add_argument(
        "--counts-output-sparse",
        action="store_true",
        help="Store the counts as a sparse matrix (mtx, h5ad, loom).",
    )
    pa.add_argument(
        "--append-output",
        action="store_true",
        dest="output_append",
        help="Append counts output to an existing file instead of "
        + "creating a new one. This option is useful if you have "
        + "already creates a TSV/CSV/similar file with a header for your "
        + "samples (with additional columns for the feature name and any "
        + "additionl attributes) and want to fill in the rest of the file.",
    )
    pa.add_argument(
        "-n",
        "--nprocesses",
        type=int,
        dest="nprocesses",
        default=1,
        help="Number of parallel CPU processes to use (default: 1). "
        + "This option is useful to process several input files at once. "
        + "Each file will use only 1 CPU. It is possible, of course, to "
        + "split a very large input SAM/BAM files into smaller chunks "
        + "upstream to make use of this option.",
    )
    pa.add_argument(
        "--feature-query",
        type=str,
        dest="feature_query",
        default=None,
        help="Restrict to features descibed in this expression. Currently "
        + 'supports a single kind of expression: attribute == "one attr" to '
        + "restrict the GFF to a single gene or transcript, e.g. "
        + "--feature-query 'gene_name == \"ACTB\"' - notice the single "
        + "quotes around the argument of this option and the double "
        + "quotes around the gene name. Broader queries might become "
        + "available in the future.",
    )
    pa.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        dest="quiet",
        help="Suppress progress report",
    )  # and warnings" )

    pa.add_argument(
        "--with-header",
        action="store_true",
        dest="with_header",
        help="Whether to add a column header to the output TSV file indicating which column "
             + "corresponds to which input BAM file. Only used if output to console or tsv or csv file. "
             + "Default to False."
    )

    args = pa.parse_args()

    # Never use more CPUs than files
    args.nprocesses = min(args.nprocesses, len(args.samfilenames))

    # Check and sanitize annotated SAM/BAM outputs
    if args.samouts != []:
        _check_samouts(args.samfilenames, args.samout_format, args.samouts)
    else:
        args.samouts = [None for x in args.samfilenames]

    # Try to open samfiles to fail early in case any of them is not there
    _check_sam_files(args.samfilenames)

    return args


def main():
    '''Main loop for htseq-count'''

    args = _parse_sanitize_cmdline_arguments()

    warnings.showwarning = my_showwarning
    try:
        count_reads_in_features(args)
    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write(
            "  [Exception type: %s, raised in %s:%d]\n"
            % (
                sys.exc_info()[1].__class__.__name__,
                os.path.basename(traceback.extract_tb(sys.exc_info()[2])[-1][0]),
                traceback.extract_tb(sys.exc_info()[2])[-1][1],
            )
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
