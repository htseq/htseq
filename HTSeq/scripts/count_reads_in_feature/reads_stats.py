import HTSeq
import sys

from HTSeq.scripts.count_reads_in_feature.reads_io_processor import ReadsIO

class ReadsStatistics(object):
    """
    For storing a bunch of statistics about the reads.

    read_io_obj : ReadsIO object

    """

    def __init__(self, feature_attr, read_io_object):
        self.read_io_obj = read_io_object

        self.empty = 0
        self.ambiguous = 0
        self.notaligned = 0
        self.lowqual = 0
        self.nonunique = 0
        self.num_reads_processed = 0
        self.counts = {key: 0 for key in feature_attr}

    def add_num_reads_processed(self):
        self.num_reads_processed += 1

    def add_empty_read(self, read_sequence):
        self.empty += 1
        self.read_io_obj.write_to_samout(read_sequence, "__no_feature")

    def add_ambiguous_read(self, read_sequence, assignment):
        self.ambiguous += 1
        self.read_io_obj.write_to_samout(read_sequence, assignment)

    def add_low_quality_read(self, read_sequence):
        self.lowqual += 1
        self.read_io_obj.write_to_samout(read_sequence, "__too_low_aQual")

    def add_not_unique_read(self, read_sequence):
        self.nonunique += 1
        self.read_io_obj.write_to_samout(read_sequence,
                                         "__alignment_not_unique")

    def add_not_aligned_read(self, read_sequence):
        self.notaligned += 1
        self.read_io_obj.write_to_samout(read_sequence, "__not_aligned")

    def add_to_count(self, feature, value=1):
        self.counts[feature] += value

    def print_progress(self, force_print=False):
        """
        Simple function to update the progress of reads processing.

        Parameters
        ----------
        force_print : boolean
            Whether to print progress regardless of number of reads processed
            or not.

        """

        if force_print:
            do_print = True
        else:
            do_print = self.num_reads_processed > 0 and self.num_reads_processed % 100000 == 0

        if do_print:
            sys.stderr.write(
                "%d alignment record%s processed.\n" %
                (self.num_reads_processed, "s" if not self.read_io_obj.pe_mode else " pairs"))
            sys.stderr.flush()

    def generate_output(self, isam):
        return {
            'isam': isam,
            'counts': self.counts,
            'empty': self.empty,
            'ambiguous': self.ambiguous,
            'lowqual': self.lowqual,
            'notaligned': self.notaligned,
            'nonunique': self.nonunique,
        }
