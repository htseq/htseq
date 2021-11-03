import warnings
import argparse # needed?
import traceback
import os.path
import sys

import HTSeq
import HTSeq.scripts.count_reads_in_feature.param_builder as param_builder
from HTSeq.scripts.count_reads_in_feature.parameters import CountParameters
from HTSeq.scripts.count_reads_in_feature.core import do_count_reads_in_features
from HTSeq.scripts.utils import my_showwarning

def main():
    """
    Call me if you want to run htseq-count!

    Where are the parameters?
    Within count_reads_in_feature there is param_builder which handles all the user defined parameters. If you need to add or modify a given parameter, please go there.
    See the comment for each step below if you want to understand better how it works.
    """

    # This will create an ArgumentParser object pre-filled with the default description of HTSeq.
    # param_builder is located in count_reads_in_feature/param_builder.py
    pa = param_builder.create_input_argument()
    args, argv = pa.parse_known_args()
    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    # This will add all the parameters/arguments you are exposing to the users.
    # So when you run htseq-count command, all the arguments/flags available to the users are defined in the method below.
    # Hence, if you need to add/modify/remove a flag/argument, modify this.
    # NOTE! If you change any argument, make sure you modify CountParameters object as well. See below.
    param_builder.add_options(pa)

    # Get the arguments/flags/parameters value given by user.
    args = pa.parse_args()

    # Show warning message?
    warnings.showwarning = my_showwarning

    # What do I do? I take all the user defined arguments/parameters/flags and turn them into a self-contained object used by the actual count function.
    # CountParameters is defined in count_reads_in_feature/parameters.py
    parameters = CountParameters(arg_parser_obj = args)

    try:
        # I'm the one who does the actual count function.
        # I'm within count_reads_in_feature/core.py script.
        do_count_reads_in_features(parameters)
    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write("  [Exception type: %s, raised in %s:%d]\n" %
                         (sys.exc_info()[1].__class__.__name__,
                          os.path.basename(traceback.extract_tb(
                              sys.exc_info()[2])[-1][0]),
                          traceback.extract_tb(sys.exc_info()[2])[-1][1]))
        sys.exit(1)

if __name__ == "__main__":
    main()
