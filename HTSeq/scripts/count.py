import warnings
import argparse # needed?
import traceback
import os.path
import sys

import HTSeq
import HTSeq.scripts.count_reads_in_feature.param_builder as param_builder
from HTSeq.scripts.count_reads_in_feature.parameters import (
    CountParameters,
    process_cmd_line_options
)
from HTSeq.scripts.count_reads_in_feature.core import do_count_reads_in_features
from HTSeq.scripts.utils import my_showwarning

def main():
    """
    Main function which runs htseq-count. This is the main function that is
    called by the htseq-count command line.
    """

    # A lot is happening in this function, but it essentially:
    # 1. Create an ArgumentParser object pre-filled with the default description of HTSeq.
    # 2. Display htseq version installed if --version flag is used.
    # 3. Process all the command-line options supplied to htseq-count command.
    # 4. Return these options as the attributes of an object (args below.)
    # See count_reads_in_feature/parameters.py file for more info.
    # Note, if you want to remove/modify/add a command-line option, please
    # modify this function.
    args = process_cmd_line_options()

    # Show warning message?
    warnings.showwarning = my_showwarning

    # What do I do? I take all the command-line options and turn them into a
    # self-contained object used by the actual count function.
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
