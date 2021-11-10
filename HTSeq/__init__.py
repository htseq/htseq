"""HTSeq is a package to process high-throughput sequencing data.

See htseq.readthedocs.io/en/master/index.html for documentation.
"""

import itertools
import warnings
import os
import shlex
import sys

import HTSeq
from HTSeq._HTSeq import *
from HTSeq.utils import FileOrSequence
from HTSeq.features import *
from HTSeq.StretchVector import StretchVector
from HTSeq._version import __version__


#########################
# GenomicArray
#########################

def read_chrom_lens(filename, delimiter="\t"):
    return dict(
        ((chrom, int(len))
         for chrom, len in csv.reader(open(filename), delimiter=delimiter)))


#########################
# Sequence readers
#########################

_re_fasta_header_line = re.compile(r'>\s*(\S+)\s*(.*)')


class FastaReader(FileOrSequence):
    """A Fasta_Reader is associated with a FASTA file or an open connection
    to a file-like object with content in FASTA format.
    It can generate an iterator over the sequences.
    """

    def __init__(self, file_, raw_iterator=False):
        FileOrSequence.__init__(self, file_)
        self.raw_iterator = raw_iterator

    def __iter__(self):
        seq = None
        name = None
        descr = None
        for line in FileOrSequence.__iter__(self):
            if line.startswith(">"):
                if seq:
                    if self.raw_iterator:
                        s = (seq, name, descr)
                    else:
                        s = Sequence(seq.encode(), name)
                        s.descr = descr
                    yield s
                mo = _re_fasta_header_line.match(line)
                name = mo.group(1)
                descr = mo.group(2)
                seq = ""
            else:
                assert seq is not None, "FASTA file does not start with '>'."
                seq += line[:-1]
        if seq is not None:
            if self.raw_iterator:
                s = (seq, name, descr)
            else:
                s = Sequence(seq.encode(), name)
                s.descr = descr
            yield s

    def get_sequence_lengths(self):
        seqname = None
        length = 0
        seqlengths = {}
        for line in FileOrSequence.__iter__(self):
            if line.startswith(">"):
                if seqname is not None:
                    seqlengths[seqname] = length
                mo = _re_fasta_header_line.match(line)
                seqname = mo.group(1)
                length = 0
            else:
                assert seqname is not None, "FASTA file does not start with '>'."
                length += len(line.rstrip())
        if seqname is not None:
            seqlengths[seqname] = length
        return seqlengths

    @staticmethod
    def _import_pysam():
        global pysam
        try:
            import pysam
        except ImportError:
            sys.stderr.write(
                "Please install the 'pysam' package to be able to use the Fasta indexing functionality.")
            raise

    def build_index(self, force=False):
        self._import_pysam()
        if not isinstance(self.fos, str):
            raise TypeError(
                "This function only works with FastaReader objects " +
                "connected to a fasta file via file name")
        index_filename = self.fos + ".fai"
        if os.access(index_filename, os.R_OK):
            if (not force) and os.stat(self.filename_or_sequence).st_mtime <= \
                 os.stat(index_filename).st_mtime:
                # index is up to date
                return
        pysam.faidx(self.fos)
        if not os.access(index_filename, os.R_OK):
            raise SystemError(
                "Building of Fasta index failed due to unknown error.")

    def __getitem__(self, iv):
        if not isinstance(iv, GenomicInterval):
            raise TypeError("GenomicInterval expected as key.")
        if not isinstance(self.fos, str):
            raise TypeError(
                "This function only works with FastaReader objects " +
                "connected to a fasta file via file name")
        self._import_pysam()
        fasta = pysam.faidx(
                self.fos,
                "%s:%d-%d" % (iv.chrom, iv.start, iv.end - 1))
        ans = list(FastaReader(fasta))
        assert len(ans) == 1
        ans[0].name = str(iv)
        if iv.strand != "-":
            return ans[0]
        else:
            return ans[0].get_reverse_complement()


class FastqReader(FileOrSequence):
    """A Fastq object is associated with a FASTQ self.file. When an iterator
    is requested from the object, the FASTQ file is read.

    qual_scale is one of "phred", "solexa", "solexa-old".
    """

    def __init__(self, file_, qual_scale="phred", raw_iterator=False):
        FileOrSequence.__init__(self, file_)
        self.qual_scale = qual_scale
        if qual_scale not in ("phred", "solexa", "solexa-old"):
            raise ValueError("Illegal quality scale.")
        self.raw_iterator = raw_iterator

    def __iter__(self):
        fin = FileOrSequence.__iter__(self)
        il = 0
        id1 = None
        id2 = None
        seq = None
        qual = None
        for line in fin:
            if il == 0:
                id1 = line
                il += 1
                continue
            elif il == 1:
                seq = line
                il += 1
                continue
            elif il == 2:
                id2 = line
                il += 1
                continue

            qual = line
            il = 0

            if qual == "":
                if id1 != "":
                    warnings.warn(
                        "Number of lines in FASTQ file is not "
                        "a multiple of 4. Discarding the last, "
                        "incomplete record")
                break

            if not qual.endswith("\n"):
                qual += "\n"
            if not id1.startswith("@"):
                raise ValueError(
                    "Primary ID line in FASTQ file does "
                    "not start with '@'. Either this is not FASTQ data or the "
                    "parser got out of sync.")
            if not id2.startswith("+"):
                raise ValueError(
                    "Secondary ID line in FASTQ file does"
                    "not start with '+'. Maybe got out of sync.")
            if len(id2) > 2 and id1[1:] != id2[1:]:
                raise ValueError(
                    "Primary and secondary ID line in FASTQ"
                    "disagree.")

            if self.raw_iterator:
                s = (seq[:-1], id1[1:-1], qual[:-1], self.qual_scale)
            else:
                s = SequenceWithQualities(
                        seq[:-1].encode(), id1[1:-1],
                        qual[:-1].encode(),
                        self.qual_scale)
            yield s


class BowtieReader(FileOrSequence):
    """A BowtieFile object is associated with a Bowtie output file that
    contains short read alignments. It can generate an iterator of Alignment
    objects."""

    def __iter__(self):
        for line in FileOrSequence.__iter__(self):
            try:
                algnt = BowtieAlignment(line)
            except ValueError:
                if line.startswith("Reported "):
                    continue
                warnings.warn(
                    "BowtieReader: Ignoring the following line, which could "
                    "not be parsed:\n%s\n" % line,
                    RuntimeWarning)
            yield algnt


def bundle_multiple_alignments(sequence_of_alignments):
    """Some alignment programs, e.g., Bowtie, can output multiple alignments,
    i.e., the same read is reported consecutively with different alignments.
    This function takes an iterator over alignments and bundles consecutive
    alignments regarding the same read to a list of Alignment objects and
    returns an iterator over these.
    """
    alignment_iter = iter(sequence_of_alignments)
    algnt = next(alignment_iter)
    ma = [algnt]
    for algnt in alignment_iter:
        if algnt.read.name != ma[0].read.name:
            yield ma
            ma = [algnt]
        else:
            ma.append(algnt)
    yield ma


class SolexaExportAlignment(Alignment):
    """Iterating over SolexaExportReader objects will yield SoelxaExportRecord
    objects. These have four fields:
       read          - a SequenceWithQualities object
       aligned       - a boolean, indicating whether the object was aligned
       iv            - a GenomicInterval giving the alignment (or None, if not aligned)
       passed_filter - a boolean, indicating whether the object passed the filter
       nomatch_code  - a code indicating why no match was found (or None, if the
         read was aligned)

    As long as 'aligned' is True, a SolexaExportRecord can be treated as an
    Alignment object.
    """

    def __init__(self):
        # Data is filled in by SolexaExportRecord
        pass

    def __repr__(self):
        if self.aligned:
            return "< %s object: Read '%s', aligned to %s >" % (
              self.__class__.__name__, self.read.name, self.iv)
        else:
            return "< %s object: Non-aligned read '%s' >" % (
              self.__class__.__name__, self.read.name)


class SolexaExportReader(FileOrSequence):
    """Parser for *_export.txt files from the SolexaPipeline software.

    Iterating over a SolexaExportReader yields SolexaExportRecord objects.
    """

    def __init__(self, filename_or_sequence, solexa_old=False):
        FileOrSequence.__init__(self, filename_or_sequence)
        if solexa_old:
            self.qualscale = "solexa-old"
        else:
            self.qualscale = "solexa"

    @classmethod
    def parse_line_bare(dummy, line):
        if line[-1] == "\n":
            line = line[:-1]
        res = {}
        (res['machine'],
         res['run_number'],
         res['lane'],
         res['tile'],
         res['x_coord'],
         res['y_coord'],
         res['index_string'],
         res['read_nbr'],
         res['read_seq'],
         res['qual_str'],
         res['chrom'],
         res['contig'],
         res['pos'],
         res['strand'],
         res['match_descr'],
         res['single_read_algnt_score'],
         res['paired_read_algnt_score'],
         res['partner_chrom'],
         res['partner_contig'],
         res['partner_offset'],
         res['partner_strand'],
         res['passed_filtering']) = line.split("\t")
        return res

    def __iter__(self):
        for line in FileOrSequence.__iter__(self):
            record = SolexaExportAlignment()
            fields = SolexaExportReader.parse_line_bare(line)
            if fields['read_nbr'] != "1":
                warnings.warn(
                    "Paired-end read encountered. PE is so far supported only "
                    "for SAM files, not yet for SolexaExport. All PE-related "
                    "fields are ignored.")
            record.read = SequenceWithQualities(
                fields['read_seq'],
                "%s:%s:%s:%s:%s#0" % (fields['machine'],
                                      fields['lane'],
                                      fields['tile'],
                                      fields['x_coord'],
                                      fields['y_coord']),
                fields['qual_str'], self.qualscale)
            if fields['passed_filtering'] == 'Y':
                record.passed_filter = True
            elif fields['passed_filtering'] == 'N':
                record.passed_filter = False
            else:
                raise ValueError(
                    "Illegal 'passed filter' value in Solexa export data: '%s'." % fields['passed_filtering'])
            record.index_string = fields['index_string']
            if fields['pos'] == '':
                record.iv = None
                record.nomatch_code = fields['chrom']
            else:
                if fields['strand'] == 'F':
                    strand = '+'
                elif fields['strand'] == 'R':
                    strand = '-'
                else:
                    raise ValueError(
                        "Illegal strand value in Solexa export data.")
                start = int(fields['pos'])
                chrom = fields['chrom']
                if fields['chrom'] == "":
                    chrom = fields['contig']
                record.iv = GenomicInterval(
                    chrom, start,
                    start + len(fields['read_seq']), strand)
            yield record


class GenomicArrayOfSets(GenomicArray):
    """A GenomicArrayOfSets is a specialization of GenomicArray that allows to store
    sets of objects. On construction, the step vectors are initialized with empty sets.
    By using the 'add_value' method, objects can be added to intervals. If an object
    is already present in the set(s) at this interval, an the new object is added to
    the present set, and the set is split if necessary.
    """

    def __init__(self, chroms, stranded=True, storage='step', memmap_dir=""):
        GenomicArray.__init__(self, chroms, stranded, 'O', storage, memmap_dir)

    def add_chrom(self, chrom, length=sys.maxsize, start_index=0):
        GenomicArray.add_chrom(self, chrom, length, start_index)
        for cv in list(self.chrom_vectors[chrom].values()):
            cv[:] = set()
            cv.is_vector_of_sets = True


###########################
# paired-end handling
###########################

def pair_SAM_alignments(
        alignments,
        bundle=False,
        primary_only=False):
    '''Iterate over SAM aligments, name-sorted paired-end

    Args:
        alignments (iterator of SAM/BAM alignments): the alignments to wrap
        bundle (bool): if True, bundle all alignments from one read pair into a
            single yield. If False (default), each pair of alignments is
            yielded separately.
        primary_only (bool): for each read, consider only the primary line
            (SAM flag 0x900 = 0). The SAM specification requires one and only
            one of those for each read.

    Yields:
        2-tuples with each pair of alignments or, if bundle==True, each bundled
        list of alignments.
    '''

    mate_missing_count = [0]

    def process_list(almnt_list):
        '''Transform a list of alignment with the same read name into pairs

        Args:
            almnt_list (list): alignments to process

        Yields:
            each pair of alignments.

        This function is needed because each line of a BAM file is not a read
        but an alignment. For uniquely mapped and unmapped reads, those two are
        the same. For multimapped reads, however, there can be more than one
        alignment for each read. Also, it is normal for a mapper to uniquely
        map one read and multimap its mate.

        This function goes down the list of alignments for a given read name
        and tries to find the first mate. So if read 1 is uniquely mapped but
        read 2 is mapped 4 times, only (read 1, read 2 - first occurrence) will
        yield; the other 3 alignments of read 2 are ignored.
        '''

        while len(almnt_list) > 0:
            a1 = almnt_list.pop(0)
            # Find its mate
            for a2 in almnt_list:
                if a1.pe_which == a2.pe_which:
                    continue
                if a1.aligned != a2.mate_aligned or a1.mate_aligned != a2.aligned:
                    continue
                if not (a1.aligned and a2.aligned):
                    break
                if a1.iv.chrom == a2.mate_start.chrom and a1.iv.start == a2.mate_start.pos and \
                   a2.iv.chrom == a1.mate_start.chrom and a2.iv.start == a1.mate_start.pos:
                    break
            else:
                if a1.mate_aligned:
                    mate_missing_count[0] += 1
                    if mate_missing_count[0] == 1:
                        warnings.warn(
                            "Read " + a1.read.name + " claims to have an aligned mate " +
                            "which could not be found in an adjacent line.")
                a2 = None
            if a2 is not None:
                almnt_list.remove(a2)
            if a1.pe_which == "first":
                yield (a1, a2)
            else:
                assert a1.pe_which == "second"
                yield (a2, a1)

    almnt_list = []
    current_name = None
    for almnt in alignments:
        if not almnt.paired_end:
            raise ValueError(
                "'pair_alignments' needs a sequence of paired-end alignments")
        if almnt.pe_which == "unknown":
            raise ValueError(
                "Paired-end read found with 'unknown' 'pe_which' status.")

        # FIXME: almnt.not_primary_alignment currently means secondary
        if primary_only and (almnt.not_primary_alignment or almnt.supplementary):
            continue

        if almnt.read.name == current_name:
            almnt_list.append(almnt)
        else:
            if bundle:
                yield list(process_list(almnt_list))
            else:
                for p in process_list(almnt_list):
                    yield p
            current_name = almnt.read.name
            almnt_list = [almnt]
    if bundle:
        yield list(process_list(almnt_list))
    else:
        for p in process_list(almnt_list):
            yield p
    if mate_missing_count[0] > 1:
        warnings.warn("%d reads with missing mate encountered." %
                      mate_missing_count[0])


def pair_SAM_alignments_with_buffer(
        alignments,
        max_buffer_size=30000000,
        primary_only=False):
    '''Iterate over SAM aligments with buffer, position-sorted paired-end

    Args:
        alignments (iterator of SAM/BAM alignments): the alignments to wrap
        max_buffer_size (int): maxmal numer of alignments to keep in memory.
        primary_only (bool): for each read, consider only the primary line
            (SAM flag 0x900 = 0). The SAM specification requires one and only
            one of those for each read.

    Yields:
        2-tuples with each pair of alignments.
    '''

    almnt_buffer = {}
    ambiguous_pairing_counter = 0
    for almnt in alignments:
        if not almnt.paired_end:
            raise ValueError(
                "Sequence of paired-end alignments expected, but got single-end alignment.")
        if almnt.pe_which == "unknown":
            raise ValueError(
                "Cannot process paired-end alignment found with 'unknown' 'pe_which' status.")
        # FIXME: almnt.not_primary_alignment currently means secondary
        if primary_only and (almnt.not_primary_alignment or almnt.supplementary):
            continue

        matekey = (
            almnt.read.name,
            "second" if almnt.pe_which == "first" else "first",
            almnt.mate_start.chrom if almnt.mate_aligned else None,
            almnt.mate_start.pos if almnt.mate_aligned else None,
            almnt.iv.chrom if almnt.aligned else None,
            almnt.iv.start if almnt.aligned else None,
            -almnt.inferred_insert_size if almnt.aligned and almnt.mate_aligned else None)

        if matekey in almnt_buffer:
            if len(almnt_buffer[matekey]) == 1:
                mate = almnt_buffer[matekey][0]
                del almnt_buffer[matekey]
            else:
                mate = almnt_buffer[matekey].pop(0)
                if ambiguous_pairing_counter == 0:
                    ambiguous_pairing_first_occurance = matekey
                ambiguous_pairing_counter += 1
            if almnt.pe_which == "first":
                yield (almnt, mate)
            else:
                yield (mate, almnt)
        else:
            almntkey = (
                almnt.read.name, almnt.pe_which,
                almnt.iv.chrom if almnt.aligned else None,
                almnt.iv.start if almnt.aligned else None,
                almnt.mate_start.chrom if almnt.mate_aligned else None,
                almnt.mate_start.pos if almnt.mate_aligned else None,
                almnt.inferred_insert_size if almnt.aligned and almnt.mate_aligned else None)
            if almntkey not in almnt_buffer:
                almnt_buffer[almntkey] = [almnt]
            else:
                almnt_buffer[almntkey].append(almnt)
            if len(almnt_buffer) > max_buffer_size:
                raise ValueError(
                    "Maximum alignment buffer size exceeded while pairing SAM alignments.")

    if len(almnt_buffer) > 0:
        warnings.warn(
            "Mate records missing for %d records; first such record: %s." %
            (len(almnt_buffer), str(list(almnt_buffer.values())[0][0])))
        for almnt_list in list(almnt_buffer.values()):
            for almnt in almnt_list:
                if almnt.pe_which == "first":
                    yield (almnt, None)
                else:
                    yield (None, almnt)

    if ambiguous_pairing_counter > 0:
        warnings.warn(
            "Mate pairing was ambiguous for %d records; mate key for first such record: %s." %
            (ambiguous_pairing_counter, str(ambiguous_pairing_first_occurance)))


###########################
# variant calls
###########################


_re_vcf_meta_comment = re.compile("^##([a-zA-Z]+)\=(.*)$")

_re_vcf_meta_descr = re.compile(
    'ID=[^,]+,?|Number=[^,]+,?|Type=[^,]+,?|Description="[^"]+",?')

_re_vcf_meta_types = re.compile("[INFO|FILTER|FORMAT]")

_vcf_typemap = {
    "Integer": int,
    "Float": float,
    "String": str,
    "Flag": bool
}


class VariantCall:
    '''Class representing a variant call, close to VCF format'''

    def __init__(
            self,
            chrom=None,
            pos=None,
            identifier=None,
            ref=None,
            alt=None,
            qual=None,
            filtr=None,
            info=None):
        '''Class representing a variant call.

        Arguments:
           chrom (str): Chromosome
           pos (int): Position on the chromosome
           identifier (str): ID of the variant
           ref (str): Reference allele
           alt (str): Alternate allele
           qual (str): Quality of the variant
           filtr (str): Filter flag indicating if the variant passed QC.
           info (str): Additional info on the variant
        '''
        self.chrom = chrom
        self.pos = pos
        self.id = identifier
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filtr
        self.info = info
        self._original_line = None

    @classmethod
    def fromdict(cls, dictionary):
        '''Create a VariantCall instance from a dict of properties'''
        ret = cls()
        ret.chrom = dictionary["chrom"]
        ret.pos = dictionary["pos"]
        ret.id = dictionary["id"]
        ret.ref = dictionary["ref"]
        ret.alt = dictionary["alt"]
        ret.qual = dictionary["qual"]
        ret.filter = dictionary["filter"]
        ret.info = dictionary["info"]
        ret._original_line = None

    @classmethod
    def fromline(cls, line, nsamples=0, sampleids=[]):
        '''Create a VariantCall instance from a VCF line'''
        ret = cls()
        if nsamples == 0:
            ret.format = None
            ret.chrom, ret.pos, ret.id, ret.ref, ret.alt, ret.qual, ret.filter, ret.info = line.rstrip("\n").split("\t", 7)
        else:
            lsplit = line.rstrip("\n").split("\t")
            ret.chrom, ret.pos, ret.id, ret.ref, ret.alt, ret.qual, ret.filter, ret.info = lsplit[:8]
            ret.format = lsplit[8].split(":")
            ret.samples = {}
            spos = 9
            for sid in sampleids:
                ret.samples[sid] = dict((name, value) for (
                    name, value) in zip(ret.format, lsplit[spos].split(":")))
                spos += 1
        ret.pos = GenomicPosition(ret.chrom, int(ret.pos))
        ret.alt = ret.alt.split(",")
        ret._original_line = line
        return ret

    def infoline(self):
        if self.info.__class__ == dict:
            return ";".join(map((lambda key: str(key) + "=" + str(self.info[key])), self.info))
        else:
            return self.info

    def get_original_line(self):
        warnings.warn(
            "Original line is empty, probably this object was created from scratch and not from a line in a .vcf file!")
        return self._original_line

    def sampleline(self):
        if self.format == None:
            sys.stderr.write("No samples in this variant call!\n")
            return ""
        keys = self.format
        ret = [":".join(keys)]
        for sid in self.samples:
            tmp = []
            for k in keys:
                if k in self.samples[sid]:
                    tmp.append(self.samples[sid][k])
            ret.append(":".join(tmp))
        return "\t".join(ret)

    def to_line(self):
        '''Convert into a VCF line'''
        if self.format == None:
            return "\t".join(map(str, [self.pos.chrom, self.pos.pos, self.id, self.ref, ",".join(self.alt), self.qual, self.filter, self.infoline()])) + "\n"
        else:
            return "\t".join(map(str, [self.pos.chrom, self.pos.pos, self.id, self.ref, ",".join(self.alt), self.qual, self.filter, self.infoline(), self.sampleline()])) + "\n"

    def __descr__(self):
        return "<VariantCall at %s, ref '%s', alt %s >" % (str(self.pos).rstrip("/."), self.ref, str(self.alt).strip("[]"))

    def __str__(self):
        return "%s:'%s'->%s" % (str(self.pos).rstrip("/."), self.ref, str(self.alt).strip("[]"))

    def unpack_info(self, infodict):
        tmp = {}
        for token in self.info.strip(";").split(";"):
            if re.compile("=").search(token):
                token = token.split("=")
                if token[0] in infodict:
                    tmp[token[0]] = list(
                        map(infodict[token[0]], token[1].split(",")))
                else:
                    tmp[token[0]] = token[1].split(",")
                if len(tmp[token[0]]) == 1:
                    tmp[token[0]] = tmp[token[0]][0]
            else:  # Flag attribute found
                tmp[token] = True
        diff = set(infodict.keys()).difference(set(tmp.keys()))
        for key in diff:
            if infodict[key] == bool:
                tmp[key] = False
        self.info = tmp


class VCF_Reader(FileOrSequence):
    '''Reader for VCF files.

    This class parses text VCF files from scratch, independently of pysam.
    '''

    def __init__(self, filename_or_sequence):
        FileOrSequence.__init__(self, filename_or_sequence)
        self.metadata = {}
        self.info = {}
        self.filters = {}
        self.formats = {}
        self.nsamples = 0
        self.sampleids = []

    def make_info_dict(self):
        self.infodict = {}
        for key in self.info.keys():
            self.infodict[key] = _vcf_typemap[self.info[key]["Type"]]

    def parse_meta(self, header_filename=None):
        if header_filename is None:
            the_iter = FileOrSequence.__iter__(self)
        else:
            the_iter = open(header_filename, "r")

        for line in the_iter:
            if line.startswith('#'):
                if line.startswith("##"):
                    mo = _re_vcf_meta_comment.match(line)
                    if mo:
                        value = mo.group(2)
                        if mo.group(1) == "INFO":
                            value = dict(e.rstrip(",").split("=", 1)
                                         for e in _re_vcf_meta_descr.findall(value))
                            key = value["ID"]
                            del value["ID"]
                            self.info[key] = value
                        elif mo.group(1) == "FILTER":
                            value = dict(e.rstrip(",").split("=", 1)
                                         for e in _re_vcf_meta_descr.findall(value))
                            key = value["ID"]
                            del value["ID"]
                            self.filters[key] = value
                        elif mo.group(1) == "FORMAT":
                            value = dict(e.rstrip(",").split("=", 1)
                                         for e in _re_vcf_meta_descr.findall(value))
                            key = value["ID"]
                            del value["ID"]
                            self.formats[key] = value
                        else:
                            self.metadata[mo.group(1)] = mo.group(2)
                else:
                    self.sampleids = line.rstrip("\t\n").split("\t")[9:]
                    self.nsamples = len(self.sampleids)
                continue
            else:
                break

    def meta_info(self, header_filename=None):
        ret = []
        if header_filename is None:
            the_iter = FileOrSequence.__iter__(self)
        else:
            the_iter = open(header_filename, "r")

        for line in the_iter:
            if line.startswith('#'):
                ret.append(line)
            else:
                break
        return ret

    def __iter__(self):
        for line in FileOrSequence.__iter__(self):
            if line == "\n" or line.startswith('#'):
                continue
            vc = VariantCall.fromline(line, self.nsamples, self.sampleids)
            yield vc


class WiggleReader(FileOrSequence):

    def __init__(self, filename_or_sequence, verbose=True):
        FileOrSequence.__init__(self, filename_or_sequence)
        self.attributes = {}
        self.stepType = 'none'
        self.verbose = verbose

    def __iter__(self):
        span = 1
        pos = None
        step = None
        chrom = None
        for line in FileOrSequence.__iter__(self):
            if line.startswith('track'):
                fields = shlex.split(line)[1:]
                self.attributes = dict([(p[0], p[1].strip('"'))
                                        for p in [x.split("=") for x in fields]])
            elif line.startswith('fixedStep'):  # do fixed step stuff
                self.stepType = 'fixed'
                fields = shlex.split(line)[1:]
                declarations = dict([(p[0], p[1].strip('"'))
                                     for p in [x.split("=") for x in fields]])
                pos = int(declarations['start'])
                step = int(declarations['step'])
                chrom = declarations['chrom']
                if 'span' in declarations:
                    span = int(declarations['span'])
                else:
                    span = 1
            elif line.startswith('variableStep'):  # do variable step stuff
                self.stepType = 'variable'
                fields = shlex.split(line)[1:]
                declarations = dict([(p[0], p[1].strip('"'))
                                     for p in [x.split("=") for x in fields]])
                chrom = declarations['chrom']
                if 'span' in declarations:
                    span = int(declarations['span'])
                else:
                    span = 1
            elif line.startswith('browser') or line.startswith('#'):  # Comment or ignored
                if self.verbose:
                    print("Ignored line:", line)
                continue
            else:
                if self.stepType == 'fixed':
                    yield (GenomicInterval(chrom, pos, pos + span, '.'), float(line.strip()))
                    pos += step
                elif self.stepType == 'variable':
                    tmp = line.strip().split(" ")
                    pos = int(tmp[0])
                    yield (GenomicInterval(chrom, pos, pos + span, '.'), float(tmp[1]))


class BAM_Reader:
    '''Parser for SAM/BAM/CRAM files.

    This is a thin wrapper on top of pysam.AlignmentFile. It detects
    automatically whether the input file is text (SAM) or binary (BAM/CRAM) via
    the HTSlib library.
    '''

    def __init__(
            self,
            filename,
            check_sq=True):
        '''Parser for SAM/BAM/CRAM files, a thin layer over pysam.AlignmentFile.

        Arguments:
           filename (str, Path): The path to the input file to read
           check_sq (bool): check if SQ entries are present in header
        '''

        global pysam
        self.filename = filename
        self.sf = None
        self.record_no = -1
        self.check_sq = check_sq
        try:
            import pysam
        except ImportError:
            sys.stderr.write(
                "Please install pysam to use the BAM_Reader class")
            raise
        self._open_file()

    def _open_file(self):
        self.sf = pysam.AlignmentFile(
                self.filename,
                check_sq=self.check_sq,
                )

    def close(self):
        """Close the BAM file for clean up"""
        self.sf.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        self.record_no = 0
        for pa in self.sf:
            yield SAM_Alignment.from_pysam_AlignedSegment(pa, self.sf)
            self.record_no += 1

    def fetch(self, reference=None, start=None, end=None, region=None):
        self.record_no = 0
        try:
            for pa in self.sf.fetch(reference, start, end, region):
                yield SAM_Alignment.from_pysam_AlignedRead(pa, self.sf)
                self.record_no += 1
        except ValueError as e:
            if e.message == "fetch called on bamfile without index":
                print("Error: ", e.message)
                print(
                    "Your bam index file is missing or wrongly named, convention is that file 'x.bam' has index file 'x.bam.bai'!")
            else:
                raise
        except:
            raise

    def get_line_number_string(self):
        if self.record_no == -1:
            return "unopened file %s" % (self.filename)
        else:
            return "record #%d in file %s" % (self.record_no, self.filename)

    def __getitem__(self, iv):
        if not isinstance(iv, GenomicInterval):
            raise TypeError(
                "Use a HTSeq.GenomicInterval to access regions within .bam-file!")
        if self.sf is None:
            self._open_file()

        if (hasattr(self.sf, '_hasIndex') and (not self.sf._hasIndex())) or (not self.sf.has_index()):
            raise ValueError(
                "The .bam-file has no index, random-access is disabled!")
        for pa in self.sf.fetch(iv.chrom, iv.start + 1, iv.end):
            yield SAM_Alignment.from_pysam_AlignedRead(pa, self.sf)

    def get_header_dict(self):
        return self.sf.header

    def get_template(self):
        return self.sf


# NOTE: this will be deprecated
SAM_Reader = BAM_Reader


class BAM_Writer:
    '''Writer for SAM/BAM/CRAM files, a thin layer over pysam.AlignmentFile'''
    def __init__(
            self,
            filename,
            template=None,
            referencenames=None,
            referencelengths=None,
            text=None,
            header=None):
        try:
            import pysam
        except ImportError:
            sys.stderr.write(
                "Please Install pysam to use the BAM_Writer Class")
            raise

        self.filename = filename
        self.template = template
        self.referencenames = referencenames
        self.referencelengths = referencelengths
        self.text = text
        self.header = header
        self.sf = pysam.AlignmentFile(
                self.filename,
                mode="wb",
                template=self.template,
                referencenames=self.referencenames,
                referencelengths=self.referencelengths,
                text=self.text,
                header=self.header)

    @classmethod
    def from_BAM_Reader(cls, fn, br):
        return BAM_Writer(filename=fn, header=br.get_header_dict())

    def write(self, alnmt):
        self.sf.write(alnmt.to_pysam_AlignedSegment(self.sf))

    def close(self):
        self.sf.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()


class BED_Reader(FileOrSequence):
    '''Reader for BED files.

    This class simply parses the BED as a text file and converts the various
    columns into HTSeq objects. For each row it extracts:

    - a GenomicInterval with chromosome equal to the first column, start and
      end to the second and third columns, and strandedness equal to the sixth
      column;

    - a GenomicFeature with name equal to the fourth column (or "unnamed"),
      type set to "BED line", score equal to the fifth column and interval set
      as the GenomicInterval above.

    - If the "thick" start and end are provided in the BED file as seventh and
      eight columns, they are stored as a GenomicInterval in the "thick"
      attribute of the GenomicFeature above (i.e. feature.thick).

    - If the itemRgb color line is provided in the BED file as ninth column, it
      is stored as "itemRgb" attribute in the GenomicFeature above (i.e.
      feature.itemRgb).

    - The blockCount, blockStart, blockSizes columns (tenth to twelth) are
      currently ignored, this might change in the future.

    Rows starting with "track" are skipped.
    '''

    def __init__(self, filename_or_sequence):
        FileOrSequence.__init__(self, filename_or_sequence)

    def __iter__(self):
        for line in FileOrSequence.__iter__(self):
            if line.startswith("track"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError("BED file line contains less than 3 fields")
            if len(fields) > 12:
                raise ValueError("BED file line contains more than 12 fields")
            iv = GenomicInterval(
                fields[0],
                int(fields[1]),
                int(fields[2]),
                fields[5] if len(fields) > 5 else ".")
            f = GenomicFeature(
                fields[3] if len(fields) > 3 else "unnamed",
                "BED line",
                iv)
            f.score = float(fields[4]) if len(fields) > 4 else None
            f.thick = GenomicInterval(
                iv.chrom,
                int(fields[6]),
                int(fields[7]),
                iv.strand) if len(fields) > 7 else None
            f.itemRgb = [int(a) for a in fields[8].split(",")
                         ] if len(fields) > 8 else None
            yield(f)


class BigWig_Reader:
    '''A simple reader for BigWig files (using pyBigWig)'''

    def __init__(self, filename):
        '''Parser for BigWig files, a thin layer over pyBigWig.

        Arguments:
           filename (str, Path): The path to the input file to read
        '''
        global pyBigWig

        try:
            import pyBigWig
        except ImportError:
            sys.stderr.write(
                "Please Install pyBigWig to use the BigWig_Reader Class")
            raise

        self.filename = filename
        self.sf = pyBigWig.open(filename)

    def close(self):
        self.sf.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def chroms(self):
        '''Return the list of chromosomes and their lengths, as a dictionary.

        Example:

        bw.chroms() -> {'chr1': 4568999, 'chr2': 87422, ...}
        '''
        return self.sf.chroms()

    def intervals(self, chrom, strand='.', raw=False):
        '''Lazy iterator over genomic intervals

        Args:
            chrom (str): The chromosome/scaffold to find intervals for.
            strand ('.', '+', or '-'): Strandedness of the yielded
              GenomicInterval. If raw=True, this argument is ignored.
            raw (bool): If True, return the raw triplet from pyBigWig. If False,
              return the result wrapped in a GenomicInterval with the
              appropriate strandedness.
        '''
        for (chrom, start, end) in self.sf.intervals(chrom):
            if raw:
                yield (chrom, start, end)
            else:
                yield GenomicInterval(chrom, start, end, strand=strand)


# TODO: make a BigWig_Writer class with buffered write operations, i.e. move it
# from the .pyx file. One would probably want to lazy out the header by element
# as well.
