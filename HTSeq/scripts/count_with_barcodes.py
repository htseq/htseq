import sys
import argparse
from collections import Counter, defaultdict
import operator
import itertools
import warnings
import traceback
import os.path
import multiprocessing
import numpy as np
import pysam

import HTSeq
from HTSeq.scripts.utils import (
    UnknownChrom,
    my_showwarning,
    invert_strand,
    _write_output,
)


def correct_barcodes(counts, hamming=1):
    '''Correct barcodes, usually UMIs.

    Notice: This function does not use recursive correction. Recursion sounds
    great in theory, but due to experimental corner cases it can lead to
    overcorrection and loss of signal.
    '''
    if hamming == 0:
        return

    # Count reads from all feature per barcode, and prepare to oder
    n_reads = Counter({key: sum(val.values()) for key, val in counts.items()})

    # Order by counts, from most to least
    order = [umi for (umi, _) in n_reads.most_common()]

    # Get close Hamming distances, and aggregate them into higher count ones
    umi_vectors = np.array([list(x) for x in order])
    idx_left = list(range(len(order)))
    while idx_left:
        i = idx_left.pop(0)
        vi = ''.join(umi_vectors[i])
        # Distance from all remaining barcodes
        dis = (umi_vectors[i] != umi_vectors[idx_left]).sum(axis=1)
        # Get indices of barcodes within reach
        idx = (dis <= hamming).nonzero()[0]
        js = [idx_left[idxi] for idxi in idx]
        for j in js:
            # Merge barcode counts into higher count one
            vj = ''.join(umi_vectors[j])
            counts[vi].update(counts.pop(vj))
        # Shorten list of remaining UMIs
        idx_left = [j for j in idx_left if j not in js]


def count_reads_with_barcodes(
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
        cb_tag,
        ub_tag,
        correct_ub_distance,
        ):

    def write_to_samout(r, assignment, samoutfile, template=None):
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

    def identify_barcodes(r):
        '''Identify barcode from the read or pair (both must have the same)'''
        if not pe_mode:
            r = (r,)

        # If either cell or UMI barcode doesn't exist, just raise exception
        has_cb_tag = False
        has_ub_tag = False
        for read in r:
            if read is not None:
                # If tags have not been found, then try to find it
                if not has_cb_tag:
                    has_cb_tag = read.has_optional_field(cb_tag)
                if not has_ub_tag:
                    has_ub_tag = read.has_optional_field(ub_tag)
        if not has_cb_tag or not has_ub_tag:
            raise Exception("Missing cell or UMI barcode")

        # cell, UMI
        barcodes = [None, None]
        nbar = 0
        for read in r:
            if read is not None:
                for tag, val in read.optional_fields:
                    if tag == cb_tag:
                        barcodes[0] = val
                        nbar += 1
                        if nbar == 2:
                            return barcodes
                    elif tag == ub_tag:
                        barcodes[1] = val
                        nbar += 1
                        if nbar == 2:
                            return barcodes
        return barcodes

    try:
        if sam_filename == "-":
            read_seq_file = HTSeq.BAM_Reader(sys.stdin)
        else:
            read_seq_file = HTSeq.BAM_Reader(sam_filename)

        # Get template for output BAM
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
    except:
        sys.stderr.write(
            "Error occured when reading beginning of SAM/BAM file.\n")
        raise

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')

    try:
        if pe_mode:
            if ((supplementary_alignment_mode == 'ignore') and
               (secondary_alignment_mode == 'ignore')):
                primary_only = True
            else:
                primary_only = False
            if order == "name":
                read_seq = HTSeq.pair_SAM_alignments(
                        read_seq,
                        primary_only=primary_only)
            elif order == "pos":
                read_seq = HTSeq.pair_SAM_alignments_with_buffer(
                        read_seq,
                        max_buffer_size=max_buffer_size,
                        primary_only=primary_only)
            else:
                raise ValueError("Illegal order specified.")

        # The nesting is cell barcode, UMI, feature
        counts = defaultdict(lambda: defaultdict(Counter))
        i = 0
        for r in read_seq:
            if i > 0 and i % 100000 == 0 and not quiet:
                sys.stderr.write(
                    "%d alignment record%s processed.\n" %
                    (i, "s" if not pe_mode else " pairs"))
                sys.stderr.flush()

            i += 1

            try:
                cb, ub = identify_barcodes(r)
            except:
                # Happens when cb or ub is not found
                write_to_samout(
                    r, "__too_low_aQual", samoutfile,
                    template)
                continue

            if not pe_mode:
                if not r.aligned:
                    counts[cb][ub]['__not_aligned'] += 1
                    write_to_samout(
                            r, "__not_aligned", samoutfile,
                            template)
                    continue
                if ((secondary_alignment_mode == 'ignore') and
                   r.not_primary_alignment):
                    continue
                if ((supplementary_alignment_mode == 'ignore') and
                   r.supplementary):
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        counts[cb][ub]['__alignment_not_unique'] += 1
                        write_to_samout(
                                r,
                                "__alignment_not_unique",
                                samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    counts[cb][ub]['__too_low_aQual'] += 1
                    write_to_samout(
                            r, "__too_low_aQual", samoutfile,
                            template)
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
                        write_to_samout(
                                r, "__not_aligned", samoutfile,
                                template)
                        counts[cb][ub]['__not_aligned'] += 1
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
                        write_to_samout(
                                r, "__alignment_not_unique", samoutfile,
                                template)
                        counts[cb][ub]['__alignment_not_unique'] += 1
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if ((r[0] and r[0].aQual < minaqual) or
                   (r[1] and r[1].aQual < minaqual)):
                    write_to_samout(
                            r, "__too_low_aQual", samoutfile,
                            template)
                    counts[cb][ub]['__too_low_aQual'] += 1
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
                    write_to_samout(
                            r, "__no_feature", samoutfile,
                            template)
                    counts[cb][ub]['__no_feature'] += 1
                elif len(fs) > 1:
                    write_to_samout(
                            r, "__ambiguous[" + '+'.join(fs) + "]",
                            samoutfile,
                            template)
                    counts[cb][ub]['__ambiguous'] += 1
                else:
                    write_to_samout(
                            r, list(fs)[0], samoutfile,
                            template)

                if fs is not None and len(fs) > 0:
                    if multimapped_mode == 'none':
                        if len(fs) == 1:
                            counts[cb][ub][list(fs)[0]] += 1
                    elif multimapped_mode == 'all':
                        for fsi in list(fs):
                            counts[cb][ub][fsi] += 1
                    else:
                        sys.exit("Illegal multimap mode.")


            except UnknownChrom:
                write_to_samout(
                        r, "__no_feature", samoutfile,
                        template)
                counts[cb][ub]['__no_feature'] += 1

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_seq_file.get_line_number_string()))
        raise

    if not quiet:
        sys.stderr.write(
            "%d %s processed.\n" %
            (i, "alignments " if not pe_mode else "alignment pairs"))
        sys.stderr.flush()

    if samoutfile is not None:
        samoutfile.close()

    # A UMI could be mapped to more than one feature. We count the feature
    # with the highest number of reads. In case of ties, we discard the whole
    # UMI to be on the safe side (it should not happen anyway).
    cbs = sorted(counts.keys())
    counts_noumi = {}
    for cb in cbs:
        counts_cell = Counter()

        # Correct barcodes within a certain Hamming distance
        correct_barcodes(counts[cb], hamming=correct_ub_distance)

        for ub, udic in counts.pop(cb).items():
            # In case of a tie, do not increment either feature
            top = udic.most_common(2)
            if (len(top) == 2) and (top[0][1] == top[1][1]):
                continue
            counts_cell[top[0][0]] += 1
        counts_noumi[cb] = counts_cell

    return {
        'cell_barcodes': cbs,
        'counts': counts_noumi,
        }


def count_reads_in_features(
        sam_filename,
        gff_filename,
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
        add_chromosome_info,
        quiet,
        minaqual,
        samout,
        samout_format,
        output_delimiter,
        output_filename,
        counts_output_sparse,
        cb_tag,
        ub_tag,
        correct_ub_distance,
        ):
    '''Count reads in features, parallelizing by file'''

    if samout is not None:
        # Try to open samout file early in case any of them has issues
        if samout_format in ('SAM', 'sam'):
            with open(samout, 'w'):
                pass
        else:
            # We don't have a template if the input is stdin
            if sam_filename != '-':
                with pysam.AlignmentFile(sam_filename, 'r') as sf:
                    with pysam.AlignmentFile(samout, 'w', template=sf):
                        pass

    # Try to open samfiles to fail early in case any of them is not there
    if sam_filename != '-':
        with pysam.AlignmentFile(sam_filename, 'r') as sf:
            pass

    # Prepare features
    gff = HTSeq.GFF_Reader(gff_filename)
    feature_scan = HTSeq.make_feature_genomicarrayofsets(
        gff,
        id_attribute,
        feature_type=feature_type,
        additional_attributes=additional_attributes,
        stranded=stranded != 'no',
        verbose=not quiet,
        add_chromosome_info=add_chromosome_info,
        )
    features = feature_scan['features']
    attributes = feature_scan['attributes']
    feature_attr = sorted(attributes.keys())

    if len(feature_attr) == 0:
        sys.stderr.write(
            "Warning: No features of type '%s' found.\n" % feature_type)

    # Count reads
    results = count_reads_with_barcodes(
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
        samout,
        cb_tag,
        ub_tag,
        correct_ub_distance,
        )

    # Write output
    _write_output(
        results,
        results['cell_barcodes'],
        attributes,
        additional_attributes,
        output_filename,
        output_delimiter,
        False,
        sparse=counts_output_sparse,
        dtype=np.float32,
    )


def main():

    pa = argparse.ArgumentParser(
        add_help=False,
    )

    pa.add_argument(
            "--version", action="store_true",
            help='Show software version and exit')
    args, argv = pa.parse_known_args()

    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    pa = argparse.ArgumentParser(
        parents=[pa],
        description="This script takes one alignment file in SAM/BAM " +
        "format and a feature file in GFF format and calculates for each feature " +
        "the number of reads mapping to it, accounting for barcodes. See " +
        "http://htseq.readthedocs.io/en/master/count.html for details.",
        epilog="Written by Simon Anders (sanders@fs.tum.de), " +
        "European Molecular Biology Laboratory (EMBL) and Fabio Zanini " +
        "(fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2020. " +
        "Released under the terms of the GNU General Public License v3. " +
        "Please cite the following paper if you use this script: \n" +
        "    G. Putri et al. Analysing high-throughput sequencing data in " +
        "Python with HTSeq 2.0. Bioinformatics (2022). " +
        "https://doi.org/10.1093/bioinformatics/btac166.\n" +
        "Part of the 'HTSeq' framework, version %s." % HTSeq.__version__,
    )

    pa.add_argument(
            "samfilename", type=str,
            help="Path to the SAM/BAM file containing the barcoded, mapped " +
            "reads. If '-' is selected, read from standard input")

    pa.add_argument(
            "featuresfilename", type=str,
            help="Path to the GTF file containing the features")

    pa.add_argument(
            "-f", "--format", dest="samtype",
            choices=("sam", "bam", "auto"), default="auto",
            help="Type of <alignment_file> data. DEPRECATED: " +
            "file format is detected automatically. This option is ignored.")

    pa.add_argument(
            "-r", "--order", dest="order",
            choices=("pos", "name"), default="name",
            help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
            "data must be sorted either by position or by read name, and the sorting order " +
            "must be specified. Ignored for single-end data.")

    pa.add_argument(
            "--max-reads-in-buffer", dest="max_buffer_size", type=int,
            default=30000000,
            help="When <alignment_file> is paired end sorted by position, " +
            "allow only so many reads to stay in memory until the mates are " +
            "found (raising this number will use more memory). Has no effect " +
            "for single end or paired end sorted by name")

    pa.add_argument(
            "-s", "--stranded", dest="stranded",
            choices=("yes", "no", "reverse"), default="yes",
            help="Whether the data is from a strand-specific assay. Specify 'yes', " +
            "'no', or 'reverse' (default: yes). " +
            "'reverse' means 'yes' with reversed strand interpretation")

    pa.add_argument(
            "-a", "--minaqual", type=int, dest="minaqual",
            default=10,
            help="Skip all reads with MAPQ alignment quality lower than the given " +
            "minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM " +
            "file and its usage depends on the software used to map the reads.")

    pa.add_argument(
            "-t", "--type", type=str, dest="featuretype",
            action="append", default=["exon"],
            help="Feature type (3rd column in GTF file) to be used, all "
            + "features of other type are ignored (default, suitable for"
            + "Ensembl GTF files: exon). You can call this option multiple "
            + "times. Features of all specified types will be included. "
            + "E.g. to include both genes and pseudogenes you might use "
            + "-t gene -t pseudogene")

    pa.add_argument(
            "-i", "--idattr", type=str, dest="idattr",
            default="gene_id",
            help="GTF attribute to be used as feature ID (default, " +
            "suitable for Ensembl GTF files: gene_id)")

    pa.add_argument(
            "--additional-attr", type=str,
            action='append',
            default=[],
            help="Additional feature attributes (default: none, " +
            "suitable for Ensembl GTF files: gene_name). Use multiple times " +
            "for each different attribute")

    pa.add_argument(
            "--add-chromosome-info", action='store_true',
            help="Store information about the chromosome of each feature as " +
            "an additional attribute (e.g. colunm in the TSV output file).",
            )

    pa.add_argument(
            "-m", "--mode", dest="mode",
            choices=("union", "intersection-strict", "intersection-nonempty"),
            default="union",
            help="Mode to handle reads overlapping more than one feature " +
            "(choices: union, intersection-strict, intersection-nonempty; default: union)")

    pa.add_argument(
            "--nonunique", dest="nonunique", type=str,
            choices=("none", "all"), default="none",
            help="Whether to score reads that are not uniquely aligned " +
            "or ambiguously assigned to features")

    pa.add_argument(
            "--secondary-alignments", dest="secondary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score secondary alignments (0x100 flag)")

    pa.add_argument(
            "--supplementary-alignments", dest="supplementary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score supplementary alignments (0x800 flag)")

    pa.add_argument(
            "-o", "--samout", type=str, dest="samout",
            default=None,
            help="Write out all SAM alignment records into a" +
            "SAM/BAM file, annotating each line " +
            "with its feature assignment (as an optional field with tag 'XF')" +
            ". See the -p option to use BAM instead of SAM.")

    pa.add_argument(
            "-p", '--samout-format', type=str, dest='samout_format',
            choices=('SAM', 'BAM', 'sam', 'bam'), default='SAM',
            help="Format to use with the --samout option."
            )

    pa.add_argument(
            "-d", '--delimiter', type=str, dest='output_delimiter',
            default='\t',
            help="Column delimiter in output (default: TAB)."
            )
    pa.add_argument(
            "-c", '--counts_output', type=str, dest='output_filename',
            default='',
            help="TSV/CSV filename to output the counts to instead of stdout."
            )

    pa.add_argument(
            "--counts_output_sparse", action='store_true',
            help="Store the counts as a sparse matrix (mtx, h5ad, loom)."
            )

    pa.add_argument(
            '--cell-barcode', type=str, dest='cb_tag',
            default='CB',
            help='BAM tag used for the cell barcode (default compatible ' +
            'with 10X Genomics Chromium is CB).',
            )

    pa.add_argument(
            '--UMI', type=str, dest='ub_tag',
            default='UB',
            help='BAM tag used for the unique molecular identifier, also ' +
            'known as molecular barcode (default compatible ' +
            'with 10X Genomics Chromium is UB).',
            )

    pa.add_argument(
            '--correct-UMI-distance',
            type=int,
            choices=[0, 1, 2],
            dest='correct_ub_distance',
            default=0,
            help='Correct for sequencing errors in the UMI tag, based on ' +
            'Hamming distance. For each UMI, if another UMI with more reads ' +
            'within 1 or 2 mutations is found, merge this UMI\'s reads into ' +
            'the more popular one. The default is to not correct UMIs.',
    )

    pa.add_argument(
            "-q", "--quiet", action="store_true", dest="quiet",
            help="Suppress progress report")  # and warnings" )

    args = pa.parse_args()

    warnings.showwarning = my_showwarning
    try:
        count_reads_in_features(
                args.samfilename,
                args.featuresfilename,
                args.order,
                args.max_buffer_size,
                args.stranded,
                args.mode,
                args.nonunique,
                args.secondary_alignments,
                args.supplementary_alignments,
                args.featuretype,
                args.idattr,
                args.additional_attr,
                args.add_chromosome_info,
                args.quiet,
                args.minaqual,
                args.samout,
                args.samout_format,
                args.output_delimiter,
                args.output_filename,
                args.counts_output_sparse,
                args.cb_tag,
                args.ub_tag,
                args.correct_ub_distance,
                )
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
