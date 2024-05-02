import pytest
import sys
import os
import glob
import sysconfig
from pathlib import Path

build_dir = "build/lib.%s-%s" % (sysconfig.get_platform(), sys.version[0:3])

sys.path.insert(0, os.path.join(os.getcwd(), build_dir))
import HTSeq


def test_fasta_parser(data_folder):
    print("Test Fasta parser")
    for seq in HTSeq.FastaReader(data_folder+'fastaExLong.fa'):
        pass
    print("Test passed")
    print("Test Fasta parser (raw iterator)")
    for seq in HTSeq.FastaReader(data_folder+'fastaExLong.fa',
                                 raw_iterator=True):
        pass
    print("Test passed")

    print('Test Fasta parser (with statement)')
    with HTSeq.FastaReader(data_folder+'fastaExLong.fa') as f:
        for seq in f:
            pass
    print('Test passed')

    print('Test Fasta parser (with statement and file handle)')
    with open(data_folder+'fastaExLong.fa') as fraw:
        f = HTSeq.FastaReader(fraw)
        for seq in f:
            pass
    print('Test passed')


def test_fastq_parser(data_folder):
    print("Test Fastq parser")
    for seq in HTSeq.FastqReader(data_folder+'fastqEx.fastq'):
        pass
    print("Test passed")
    print("Test Fastq parser on gzip input")
    for seq in HTSeq.FastqReader(data_folder+'fastqExgzip.fastq.gz'):
        pass
    print("Test passed")
    print("Test Fastq parser on gzip input (raw iterator)")
    for seq in HTSeq.FastqReader(data_folder+'fastqExgzip.fastq.gz',
                                 raw_iterator=True):
        pass
    print("Test passed")

    print('Test Fastq parser (with statement)')
    with HTSeq.FastqReader(data_folder+'fastqExgzip.fastq.gz') as f:
        for seq in f:
            pass
    print('Test passed')

    print('Test Fastq parser (with statement and file handle)')
    import gzip
    with gzip.open(data_folder+'fastqExgzip.fastq.gz', 'rt') as fraw:
        f = HTSeq.FastqReader(fraw)
        for seq in f:
            pass
    print('Test passed')


def test_path_like_reading(data_folder):
    data_file = Path(data_folder) / 'yeast_RNASeq_excerpt.sam'
    print('Test BAM reader with pathlib.Path')
    bamfile = HTSeq.BAM_Reader(data_file)
    for read in bamfile:
        pass
    print('Test passed')

    print('Test BAM reader (with statement) with pathlib.Path')
    with HTSeq.BAM_Reader(data_file) as f:
        for read in f:
            pass
    print('Test passed')


def test_open_handle_reading(data_folder):
    data_file = Path(data_folder) / 'yeast_RNASeq_excerpt.sam'
    print('Test BAM reader with open handle')
    with open(data_file) as f:
        reader = HTSeq.BAM_Reader(f)
        for read in reader:
            pass
    print('Test passed')


def test_unreadable_object_reading():
    unreadable_object = object()

    print('Test BAM reader with unreadable object')
    with pytest.raises(TypeError):
        reader = HTSeq.BAM_Reader(unreadable_object)
        for read in reader:
            pass
    print('Test passed')

    print('Test BAM reader with unreadable object, context manager')
    with pytest.raises(TypeError):
        with HTSeq.BAM_Reader(unreadable_object) as reader:
            for read in reader:
                pass
    print('Test passed')


def test_bam_reader(data_folder):
    print('Test BAM reader')
    bamfile = HTSeq.BAM_Reader(data_folder+"yeast_RNASeq_excerpt.sam")
    for read in bamfile:
        pass
    print('Test passed')

    print('Test BAM reader (with statement)')
    with HTSeq.BAM_Reader(data_folder+"yeast_RNASeq_excerpt.sam") as f:
        for read in f:
            pass
    print('Test passed')


def test_bam_inconsistent_mate(data_folder):
    print('Test inconsistent BAM file')
    bamfile = HTSeq.BAM_Reader(data_folder+"inconsistent_mate.bam")
    for read in bamfile:
        pass
    print("Test passed")


def test_bam_optional_field(data_folder):
    print('Test optional fields in BAM alignment')
    bamfile = HTSeq.BAM_Reader(data_folder+"inconsistent_mate.bam")
    for read in bamfile:
        read.has_optional_field('NO')
        break
    print("Test passed")


def test_GFF3(data_folder):
    print('Test GFF3 reader')
    with HTSeq.GFF_Reader(data_folder+'GCF_000001405.39_GRCh38.p13_genomic_subsample.gff.gz') as reader:
        for line in reader:
            pass
    print("Test passed")


def test_pickle():
    import pickle

    print('Test pickling and inpickling')
    pickles = [
            {'name': 'HTSeq.Sequence',
             'object': HTSeq.Sequence(b'ACTG', 'sequence'),
             'assert_properties': ('seq', 'name', 'descr')},
            ]

    for pic in pickles:
        print('Pickling '+pic['name'])
        pickled = pickle.dumps(pic['object'])
        print('Done')

        print('Unpickling '+pic['name'])
        unpick = pickle.loads(pickled)
        print('Done')

        if 'assert_properties' in pic:
            print('Checking serialized/deserialized')
            for prop in pic['assert_properties']:
                assert getattr(pic['object'], prop) == getattr(unpick, prop)
            print('Done')
    print("Test passed")


def test_bamfile_nosq(data_folder):
    print('Test parsing BAM file with no SQ field (e.g. PacBio)')
    bamfile = HTSeq.BAM_Reader(
            data_folder+"short_test_ccs.bam",
            check_sq=False)
    for read in bamfile:
        pass
    print("Test passed")
