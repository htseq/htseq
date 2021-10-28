# Anything to do with processing the options/arguments for htseq-count
# command is here.

import argparse

def get_help_message():
    # This will show the short opening message

    # TODO: Maybe move this to a file?
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
