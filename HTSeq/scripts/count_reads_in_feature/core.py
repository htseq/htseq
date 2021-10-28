# Needed? Don't think so!
import HTSeq
from HTSeq.scripts.count_reads_in_feature.parameters import CountParameters

def do_count_reads_in_features(in_param):
    """
    I am the new version of count_reads_in_features.

    Count reads in features, parallelizing by file

    Parameters
    ----------
    in_param : CountParameters
        Custom object (see parameters module) which stores the
        input parameters given by user to run htseq-count.
    """

    nprocesses = min(in_param.nprocesses, len(in_param.sam_filenames))

    in_param.prepare_features()
    
