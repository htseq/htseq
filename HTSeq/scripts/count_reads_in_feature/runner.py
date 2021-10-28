import warnings
import argparse # needed?

import HTSeq
import HTSeq.scripts.count_reads_in_feature.arguments_processor as ap


# Call me if you want to run htseq-count!
def main():

    pa = ap.get_help_message()
    args, argv = pa.parse_known_args()
    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    ap.add_options()

    args = pa.parse_args()

    warnings.showwarning = my_showwarning
