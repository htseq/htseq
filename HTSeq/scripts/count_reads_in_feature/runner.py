import warnings
import argparse # needed?

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

if __name__ == "__main__":
    run()
