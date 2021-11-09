import itertools
import random
import sys

import pysam

import HTSeq
from HTSeq.scripts.utils import invert_strand, UnknownChrom
from HTSeq.scripts.count_features.reads_io_processor import ReadsIO
from HTSeq.scripts.count_features.reads_stats import ReadsStatistics


def count_reads_single_file(
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
):
    try:
        read_io_obj = ReadsIO(sam_filename=sam_filename,
                              samout_filename=samout_filename,
                              samout_format=samout_format,
                              supplementary_alignment_mode=supplementary_alignment_mode,
                              secondary_alignment_mode=secondary_alignment_mode,
                              order=order,
                              max_buffer_size=max_buffer_size)
    except:
        sys.stderr.write(
            "Error occurred when reading beginning of SAM/BAM file.\n")
        raise

    try:
        read_stats = ReadsStatistics(feature_attr=feature_attr,
                                     read_io_object=read_io_obj)
    except:
        sys.stderr.write(
            "Error occurred when preparing object to store the reads' assignments\n")
        raise

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')

    try:

        for r in read_io_obj.read_seq:

            read_stats.print_progress()
            read_stats.add_num_reads_processed()

            # todo can move this into a function, but not necessary.
            if not read_io_obj.pe_mode:
                skip_read = _assess_non_pe_read(read_sequence=r,
                                                read_stats=read_stats,
                                                secondary_alignment_mode=secondary_alignment_mode,
                                                supplementary_alignment_mode=supplementary_alignment_mode,
                                                multimapped_mode=multimapped_mode,
                                                minaqual=minaqual)

                if skip_read:
                    continue
                iv_seq = _get_iv_seq_non_pe_read(com, r, stranded)
            else:

                # todo these assessor used to be at the bottom, after creating the iv_seq and checking whether
                #  the first element of the paired end is aligned. Kind of nuts really as it wastes time?
                #  need more testing though
                skip_read = _assess_pe_read(minaqual, multimapped_mode, r, read_stats,
                                            secondary_alignment_mode, supplementary_alignment_mode)
                if skip_read:
                    continue

                iv_seq = _get_iv_seq_pe_read(com, r, stranded)

            try:
                fs = _align_reads_to_feature_set(features, iv_seq, overlap_mode)

                if fs is None or len(fs) == 0:
                    read_stats.add_empty_read(read_sequence=r)
                elif len(fs) > 1:
                    read_stats.add_ambiguous_read(read_sequence=r,
                                                  assignment="__ambiguous[" + '+'.join(sorted(fs)) + "]")
                else:
                    read_stats.add_good_read_assignment(read_sequence=r, assignment=list(fs)[0])

                if fs is not None and len(fs) > 0:
                    fs = list(fs)
                    if multimapped_mode == 'none':
                        if len(fs) == 1:
                            read_stats.add_to_count(feature=fs[0])
                    elif multimapped_mode == 'all':
                        for fsi in fs:
                            read_stats.add_to_count(feature=fsi)
                    elif multimapped_mode == 'fraction':
                        val = 1.0 / len(fs)
                        for fsi in fs:
                            read_stats.add_to_count(feature=fsi, value=val)
                    elif multimapped_mode == 'random':
                        fsi = random.choice(fs)
                        read_stats.add_to_count(feature=fsi)
                    else:
                        sys.exit("Illegal multimap mode.")

            except UnknownChrom:
                read_stats.add_empty_read(read_sequence=r)

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_io_obj.read_seq_file.get_line_number_string()))
        raise

    if not quiet:
        read_stats.print_progress(force_print=True)

    read_io_obj.close_samoutfile()

    res = read_stats.get_output(isam)
    return res


def _align_reads_to_feature_set(features, iv_seq, overlap_mode):
    if overlap_mode == "union":
        fs = set()
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                raise UnknownChrom
            for iv2, fs2 in features[iv].steps():
                fs = fs.union(fs2)
    elif overlap_mode in ("intersection-strict",
                          "intersection-nonempty"):
        fs = None
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                raise UnknownChrom
            for iv2, fs2 in features[iv].steps():
                if ((len(fs2) > 0) or
                        (overlap_mode == "intersection-strict")):
                    if fs is None:
                        fs = fs2.copy()
                    else:
                        fs = fs.intersection(fs2)
    else:
        sys.exit("Illegal overlap mode.")
    return fs


def _get_iv_seq_pe_read(com, r, stranded):
    if r[0] is not None and r[0].aligned:
        iv_seq = _get_iv_seq_pe_read_first(com, r, stranded)
    else:
        iv_seq = tuple()
    if r[1] is not None and r[1].aligned:
        iv_seq = _get_iv_seq_pe_read_second(com, iv_seq, r, stranded)
    return iv_seq


def _assess_pe_read(minaqual, multimapped_mode, r, read_stats, secondary_alignment_mode,
                    supplementary_alignment_mode):
    if (r[0] is None) or not (r[0].aligned):
        read_stats.add_not_aligned_read(read_sequence=r)
        return True

    if secondary_alignment_mode == 'ignore':
        if (r[0] is not None) and r[0].not_primary_alignment:
            return True
        elif (r[1] is not None) and r[1].not_primary_alignment:
            return True
    if supplementary_alignment_mode == 'ignore':
        if (r[0] is not None) and r[0].supplementary:
            return True
        elif (r[1] is not None) and r[1].supplementary:
            return True
    try:
        if ((r[0] is not None and r[0].optional_field("NH") > 1) or
                (r[1] is not None and r[1].optional_field("NH") > 1)):
            read_stats.add_not_unique_read(read_sequence=r)
            if multimapped_mode == 'none':
                return True
    except KeyError:
        pass
    if ((r[0] and r[0].aQual < minaqual) or
            (r[1] and r[1].aQual < minaqual)):
        read_stats.add_low_quality_read(read_sequence=r)
        return True
    return False


def _get_iv_seq_pe_read_second(com, iv_seq, r, stranded):
    if stranded != "reverse":
        iv_seq = itertools.chain(
            iv_seq,
            (invert_strand(co.ref_iv) for co in r[1].cigar
             if co.type in com and co.size > 0))
    else:
        iv_seq = itertools.chain(
            iv_seq,
            (co.ref_iv for co in r[1].cigar
             if co.type in com and co.size > 0))
    return iv_seq


def _get_iv_seq_pe_read_first(com, r, stranded):
    if stranded != "reverse":
        iv_seq = (co.ref_iv for co in r[0].cigar
                  if co.type in com and co.size > 0)
    else:
        iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar
                  if co.type in com and co.size > 0)
    return iv_seq


def _get_iv_seq_non_pe_read(com, r, stranded):
    # TODO: rename me. I'm not sure what this is doing..
    if stranded != "reverse":
        iv_seq = (co.ref_iv for co in r.cigar if co.type in com
                  and co.size > 0)
    else:
        iv_seq = (invert_strand(co.ref_iv)
                  for co in r.cigar if (co.type in com and
                                        co.size > 0))
    return iv_seq


def _assess_non_pe_read(read_sequence, read_stats, secondary_alignment_mode, supplementary_alignment_mode,
                        multimapped_mode, minaqual):
    if not read_sequence.aligned:
        read_stats.add_not_aligned_read(read_sequence=read_sequence)
        return True
    if ((secondary_alignment_mode == 'ignore') and
            read_sequence.not_primary_alignment):
        return True
    if ((supplementary_alignment_mode == 'ignore') and
            read_sequence.supplementary):
        return True
    try:
        if read_sequence.optional_field("NH") > 1:
            read_stats.add_not_unique_read(read_sequence=read_sequence)
            if multimapped_mode == 'none':
                return True
    except KeyError:
        pass
    if read_sequence.aQual < minaqual:
        read_stats.add_low_quality_read(read_sequence=read_sequence)
        return True

    return False
