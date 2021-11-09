import itertools
import random
import sys

import pysam

import HTSeq
from HTSeq.scripts.utils import invert_strand, UnknownChrom
from HTSeq.scripts.count_features.reads_io_processor import ReadsIO


def write_to_samout(r, assignment, samoutfile, pe_mode, samout_format, template=None):
    if samoutfile is None:
        return
    if not pe_mode:
        r = (r,)
    for read in r:
        if read is not None:
            read.optional_fields.append(('XF', assignment))
            if template is not None:
                samoutfile.write(read.to_pysam_AlignedSegment(template))
            elif samout_format in ('SAM', 'sam'):
                samoutfile.write(read.get_sam_line() + "\n")
            else:
                raise ValueError(
                    'BAM/SAM output: no template and not a test SAM file',
                )


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
        # pe_mode, read_seq, read_seq_file, samoutfile, template = _prepare_samoutfile_and_read_seq(sam_filename,
        #                                                                                           samout_filename,
        #                                                                                           samout_format)
        read_io_obj = ReadsIO(sam_filename, samout_filename, samout_format, supplementary_alignment_mode,
                              secondary_alignment_mode, order, max_buffer_size)
    except:
        sys.stderr.write(
            "Error occured when reading beginning of SAM/BAM file.\n")
        raise

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')
    counts = {key: 0 for key in feature_attr}

    try:

        empty = 0
        ambiguous = 0
        notaligned = 0
        lowqual = 0
        nonunique = 0
        i = 0
        for r in read_io_obj.read_seq:
            if i > 0 and i % 100000 == 0 and not quiet:
                sys.stderr.write(
                    "%d alignment record%s processed.\n" %
                    (i, "s" if not read_io_obj.pe_mode else " pairs"))
                sys.stderr.flush()

            i += 1
            if not read_io_obj.pe_mode:
                if not r.aligned:
                    notaligned += 1
                    write_to_samout(r=r, assignment="__not_aligned", samoutfile=read_io_obj.samoutfile,
                                    pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                                    template=read_io_obj.template)
                    continue
                if ((secondary_alignment_mode == 'ignore') and
                        r.not_primary_alignment):
                    continue
                if ((supplementary_alignment_mode == 'ignore') and
                        r.supplementary):
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        nonunique += 1
                        write_to_samout(
                            r=r, assignment="__alignment_not_unique", samoutfile=read_io_obj.samoutfile,
                            pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                            template=read_io_obj.template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    lowqual += 1
                    write_to_samout(
                        r=r, assignment="__too_low_aQual", samoutfile=read_io_obj.samoutfile,
                        pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                        template=read_io_obj.template)
                    continue
                if stranded != "reverse":
                    iv_seq = (co.ref_iv for co in r.cigar if co.type in com
                              and co.size > 0)
                else:
                    iv_seq = (invert_strand(co.ref_iv)
                              for co in r.cigar if (co.type in com and
                                                    co.size > 0))
            else:
                if r[0] is not None and r[0].aligned:
                    if stranded != "reverse":
                        iv_seq = (co.ref_iv for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                    else:
                        iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                else:
                    iv_seq = tuple()
                if r[1] is not None and r[1].aligned:
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
                else:
                    if (r[0] is None) or not (r[0].aligned):
                        write_to_samout(r=r, assignment="__not_aligned", samoutfile=read_io_obj.samoutfile,
                                        pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                                        template=read_io_obj.template)
                        notaligned += 1
                        continue
                if secondary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].not_primary_alignment:
                        continue
                    elif (r[1] is not None) and r[1].not_primary_alignment:
                        continue
                if supplementary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].supplementary:
                        continue
                    elif (r[1] is not None) and r[1].supplementary:
                        continue
                try:
                    if ((r[0] is not None and r[0].optional_field("NH") > 1) or
                            (r[1] is not None and r[1].optional_field("NH") > 1)):
                        nonunique += 1
                        write_to_samout(r=r, assignment="__alignment_not_unique", samoutfile=read_io_obj.samoutfile,
                                        pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                                        template=read_io_obj.template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if ((r[0] and r[0].aQual < minaqual) or
                        (r[1] and r[1].aQual < minaqual)):
                    lowqual += 1
                    write_to_samout(
                        r=r, assignment="__too_low_aQual", samoutfile=read_io_obj.samoutfile,
                        pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                        template=read_io_obj.template)
                    continue

            try:
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

                if fs is None or len(fs) == 0:
                    write_to_samout(r=r, assignment="__no_feature", samoutfile=read_io_obj.samoutfile,
                                    pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                                    template=read_io_obj.template)
                    empty += 1
                elif len(fs) > 1:
                    write_to_samout(
                        r=r, assignment="__ambiguous[" + '+'.join(sorted(fs)) + "]", samoutfile=read_io_obj.samoutfile,
                        pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                        template=read_io_obj.template)
                    ambiguous += 1
                else:
                    write_to_samout(
                        r=r, assignment=list(fs)[0], samoutfile=read_io_obj.samoutfile,
                        pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                        template=read_io_obj.template)

                if fs is not None and len(fs) > 0:
                    if multimapped_mode == 'none':
                        if len(fs) == 1:
                            counts[list(fs)[0]] += 1
                    elif multimapped_mode == 'all':
                        for fsi in list(fs):
                            counts[fsi] += 1
                    elif multimapped_mode == 'fraction':
                        for fsi in list(fs):
                            counts[fsi] += 1.0 / len(fs)
                    elif multimapped_mode == 'random':
                        fsi = random.choice(list(fs))
                        counts[fsi] += 1
                    else:
                        sys.exit("Illegal multimap mode.")

            except UnknownChrom:
                write_to_samout(r=r, assignment="__no_feature", samoutfile=read_io_obj.samoutfile,
                                pe_mode=read_io_obj.pe_mode, samout_format=samout_format,
                                template=read_io_obj.template)
                empty += 1

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_io_obj.read_seq_file.get_line_number_string()))
        raise

    if not quiet:
        sys.stderr.write(
            "%d %s processed.\n" %
            (i, "alignments " if not read_io_obj.pe_mode else "alignment pairs"))
        sys.stderr.flush()

    if read_io_obj.samoutfile is not None:
        read_io_obj.samoutfile.close()

    return {
        'isam': isam,
        'counts': counts,
        'empty': empty,
        'ambiguous': ambiguous,
        'lowqual': lowqual,
        'notaligned': notaligned,
        'nonunique': nonunique,
    }


def _prepare_samoutfile_and_read_seq(sam_filename, samout_filename, samout_format):
    if sam_filename == "-":
        read_seq_file = HTSeq.BAM_Reader(sys.stdin)
    else:
        read_seq_file = HTSeq.BAM_Reader(sam_filename)
    # Get template for output BAM/SAM if possible
    if samout_filename is None:
        template = None
        samoutfile = None
    elif samout_format in ('bam', 'BAM'):
        template = read_seq_file.get_template()
        samoutfile = pysam.AlignmentFile(
            samout_filename, 'wb',
            template=template,
        )
    elif (samout_format in ('sam', 'SAM')) and \
            hasattr(read_seq_file, 'get_template'):
        template = read_seq_file.get_template()
        samoutfile = pysam.AlignmentFile(
            samout_filename, 'w',
            template=template,
        )
    else:
        template = None
        samoutfile = open(samout_filename, 'w')
    read_seq_iter = iter(read_seq_file)
    # Catch empty BAM files
    try:
        first_read = next(read_seq_iter)
        pe_mode = first_read.paired_end
    # FIXME: catchall can hide subtle bugs
    except:
        first_read = None
        pe_mode = False
    if first_read is not None:
        read_seq = itertools.chain([first_read], read_seq_iter)
    else:
        read_seq = []
    return pe_mode, read_seq, read_seq_file, samoutfile, template
