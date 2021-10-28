import warnings
import argparse # needed?
import traceback
import os.path
import sys

import HTSeq
import HTSeq.scripts.count_reads_in_feature.param_builder as param_builder
from HTSeq.scripts.count_reads_in_feature.parameters import CountParameters
from HTSeq.scripts.count_reads_in_feature.core import do_count_reads_in_features
from HTSeq.scripts.utils import (
    UnknownChrom,
    my_showwarning,
    invert_strand,
    _write_output,
)

def run():
    """
    Call me if you want to run htseq-count!
    """

    pa = param_builder.create_input_argument()
    args, argv = pa.parse_known_args()
    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    param_builder.add_options()

    # Get the parameters
    args = pa.parse_args()

    # Show warning message?
    warnings.showwarning = my_showwarning

    parameters = CountParameters(arg_parser_obj = args)

    try:
        do_count_reads_in_features()
    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write("  [Exception type: %s, raised in %s:%d]\n" %
                         (sys.exc_info()[1].__class__.__name__,
                          os.path.basename(traceback.extract_tb(
                              sys.exc_info()[2])[-1][0]),
                          traceback.extract_tb(sys.exc_info()[2])[-1][1]))
        sys.exit(1)

if __name__ == "__main__":
    run()
