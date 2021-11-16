import os
import subprocess as sp
import unittest
import numpy as np
import pytest
import conftest

try:
    import scipy
except ImportError:
    scipy = None

try:
    import anndata
except ImportError:
    anndata = None

try:
    import loompy
except ImportError:
    loompy = None


data_folder = conftest.get_data_folder()


def load_result_file(filename):
    sfx = filename.split('.')[-1]
    if sfx in ('csv', 'tsv'):
        with open(filename, 'r') as f:
            result = f.read()
    elif sfx == 'mtx':
        from scipy.io import mmread
        result = mmread(filename)
    elif sfx == 'h5ad':
        result = anndata.read_h5ad(filename)
    elif sfx == 'loom':
        result = loompy.connect(filename)
    else:
        raise ValueError(f'File extension not supported: {sfx}')

    return {'result': result, 'fmt': sfx}


def close_file(filename, resultd):
    fmt = resultd['fmt']
    if fmt == 'loom':
        resultd['result'].close()


class HTSeqCountBase(unittest.TestCase):
    def _customAssertEqual(self, outputd, expectedd):
        output_fmt = outputd['fmt']
        expected_fmt = expectedd['fmt']
        self.assertEqual(output_fmt, expected_fmt)

        fmt = output_fmt
        output = outputd['result']
        expected = expectedd['result']

        if fmt in ('tsv', 'csv'):
            self.assertEqual(output, expected)
        elif fmt == 'mtx':
            self.assertIsNone(
                np.testing.assert_array_equal(
                    output, expected,
                ))
        elif fmt == 'loom':
            self.assertIsNone(
                np.testing.assert_array_equal(
                    output[:, :], expected[:, :],
                ))
            #TODO: metadata and filenames
        elif fmt == 'h5ad':
            self.assertIsNone(
                np.testing.assert_array_equal(
                    output.X, output.X,
                ))
            #TODO: metadata and filenames
        else:
            raise ValueError(f'Format not supported: {fmt}')

    def _run(self, t):
        expected_fn = t.get('expected_fn', None)
        call = t['call']

        # Replace with injected variable
        call = [x.replace(f'{data_folder}/', data_folder) for x in call]
        if expected_fn is not None:
            expected_fn = expected_fn.replace(f'{data_folder}/', data_folder)

        ## local testing
        #if call[0] == 'htseq-count':
        #    call = ['python', 'HTSeq/scripts/count.py'] + call[1:]
        #else:
        #    call = ['python', 'HTSeq/scripts/count_with_barcodes.py'] + call[1:]

        print(' '.join(call))
        output = sp.check_output(
                ' '.join(call),
                shell=True,
        ).decode()

        if '-c' in call:
            output_fn = call[call.index('-c') + 1]
            output = load_result_file(output_fn)
        else:
            output = {'result': output, 'fmt': 'tsv'}
            output_fn = None

        if expected_fn is None:
            if '--version' in call:
                print('version:', output['result'])
            return

        if not os.path.isfile(expected_fn):
            print('Missing output file, creating one in current folder')
            out_fn = os.path.basename(expected_fn)
            if output_fn is None:
                with open(out_fn, 'wt') as f:
                    f.write(output['result'])
            else:
                import shutil
                shutil.copy(output_fn, out_fn)
            pytest.fail(
                'Expected filename not found, output filename copied in {out_fn}',
            )

        expected = load_result_file(expected_fn)

        try:
            self._customAssertEqual(output, expected)
        finally:
            if output_fn is not None:
                close_file(output_fn, output)
                close_file(expected_fn, expected)
                # FIXME
                if True:#output['fmt'] not in ['h5ad', 'loom']:
                    os.remove(output_fn)

class HTSeqCount(HTSeqCountBase):
    cmd = 'htseq-count'

    def test_version(self):
        self._run({
            'call': [
                self.cmd,
                '--version'],
            })

    def test_simple(self):
        self._run({
            'call': [
                self.cmd,
                f'{data_folder}/bamfile_no_qualities.sam',
                f'{data_folder}/bamfile_no_qualities.gtf',
            ],
            'expected_fn': f'{data_folder}/bamfile_no_qualities.tsv',
            })

    def test_output_tsv(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                f'{data_folder}/bamfile_no_qualities.sam',
                f'{data_folder}/bamfile_no_qualities.gtf',
                ],
            'expected_fn': f'{data_folder}/bamfile_no_qualities.tsv',
            })

    @unittest.skipIf(scipy is None, "test case depends on scipy")
    def test_output_mtx(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.mtx',
                f'{data_folder}/bamfile_no_qualities.sam',
                f'{data_folder}/bamfile_no_qualities.gtf',
                ],
            'expected_fn': f'{data_folder}/bamfile_no_qualities.mtx',
            })

    @unittest.skipIf(anndata is None, "test case depends on anndata")
    def test_output_h5ad(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.h5ad',
                f'{data_folder}/bamfile_no_qualities.sam',
                f'{data_folder}/bamfile_no_qualities.gtf',
                ],
            'expected_fn': f'{data_folder}/bamfile_no_qualities.h5ad',
            })

    @unittest.skipIf(loompy is None, "test case depends on loompy")
    def test_output_loom(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.loom',
                f'{data_folder}/bamfile_no_qualities.sam',
                f'{data_folder}/bamfile_no_qualities.gtf',
                ],
            'expected_fn': f'{data_folder}/bamfile_no_qualities.loom',
            })

    # Testing multiple cores on travis makes a mess
    #{'call': [
    #    'htseq-count',
    #    '-n', '2',
    #    f'{data_folder}/bamfile_no_qualities.sam',
    #    f'{data_folder}/bamfile_no_qualities.gtf',
    #    ],
    # 'expected_fn': f'{data_folder}/bamfile_no_qualities.tsv'},

    def test_no_qualities(self):
        self._run({
            'call': [
                self.cmd,
                f'{data_folder}/bamfile_no_qualities.bam',
                f'{data_folder}/bamfile_no_qualities.gtf',
            ],
            'expected_fn': f'{data_folder}/bamfile_no_qualities.tsv',
            })

    def test_some_missing_sequences(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                f'{data_folder}/yeast_RNASeq_excerpt_some_empty_seqs.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_some_empty_seqs.tsv',
            })

    def test_intersection_nonempty(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts.tsv',
            })

    def test_feature_query(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                '--feature-query', '\'gene_id == "YPR036W-A"\'',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_YPR036W-A.tsv',
            })

    def test_additional_attributes(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                '--additional-attr', 'gene_name',
                '--additional-attr', 'exon_number',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_additional_attributes.tsv',
            })

    def test_multiple_and_additional_attributes(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '-i', 'gene_id',
                '-i', 'exon_number',
                '--additional-attr', 'gene_name',
                '--additional-attr', 'exon_number',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_exon_level_and_additional_attributes.tsv',
            })

    def test_additional_attributes_chromosome_info(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                '--additional-attr', 'gene_name',
                '--additional-attr', 'exon_number',
                '--add-chromosome-info',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_additional_attributes_chromosome_info.tsv',
            })

    def test_nonunique_fraction(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'fraction',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_nonunique_fraction.tsv',
            })

    def test_withNH(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'all',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_nonunique.tsv',
            })

    def test_twocolumns(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '-i', 'gene_id',
                '--additional-attr', 'gene_name',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_twocolumns.tsv',
            })

    def test_ignore_secondary(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'ignore',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withNH.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withNH_counts_ignore_secondary.tsv',
            })


class HTSeqCountBarcodes(HTSeqCountBase):
    cmd = 'htseq-count-barcodes'

    def test_version(self):
        self._run({
            'call': [
                self.cmd,
                '--version'],
            })
    def test_simple(self):
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.tsv',
            })

    def test_output_tsv(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.tsv',
            })

    def test_output_tsv_chromosome_info(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                '-m', 'intersection-nonempty',
                '--add-chromosome-info',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes_chromosome_info.tsv',
            })

    def test_output_tsv_correct_UMI(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                '--correct-UMI-distance', '1',
                f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes_correctUMI_1.tsv',
            })

    @unittest.skipIf(anndata is None, "test case depends on anndata")
    def test_output_h5ad(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.h5ad',
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.h5ad',
            })

    @unittest.skipIf(loompy is None, "test case depends on loompy")
    def test_output_loom(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.loom',
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_withbarcodes.loom',
            })


if __name__ == '__main__':

    suite = HTSeqCount()
    suite.test_version()
    suite.test_simple()
    suite.test_output_tsv()
    suite.test_no_qualities()
    suite.test_intersection_nonempty()
    suite.test_feature_query()
    suite.test_additional_attributes()
    suite.test_nonunique_fraction()
    suite.test_withNH()
    suite.test_twocolumns()
    suite.test_ignore_secondary()

    suite = HTSeqCountBarcodes()
    suite.test_version()
    suite.test_barcodes()
