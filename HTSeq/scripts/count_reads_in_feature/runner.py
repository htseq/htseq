import warnings
import argparse # needed?

import HTSeq
import HTSeq.scripts.count_reads_in_feature.process_options as process_options


# Call me if you want to run htseq-count!
def main():

    pa = process_options.get_help_message()
    args, argv = pa.parse_known_args()
    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    process_options.add_options()

    args = pa.parse_args()

    warnings.showwarning = my_showwarning
