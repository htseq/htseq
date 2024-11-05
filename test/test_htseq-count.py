import os
import subprocess as sp
import unittest
import numpy as np
import pysam
import pytest
import conftest
import pandas as pd
import math

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

    def _run(self, t, remove_res_files=True):
        expected_fn = t.get('expected_fn', None)
        expected_stderr = t.get('expected_stderr', None)
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
        p = sp.run(
            ' '.join(call),
            shell=True,
            check=True,
            capture_output=True,
        )
        if expected_stderr:
            actual_stderr_str = p.stderr.decode()
            self.assertEqual(actual_stderr_str, expected_stderr)
        output = p.stdout.decode()

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
                if remove_res_files:#output['fmt'] not in ['h5ad', 'loom']:
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

    def test_output_tsv_header(self):
        # Header must be there
        self._run({
            'call': [
                self.cmd,
                '--with-header', ' ',
                '-c', f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
                f'{data_folder}/10x_pbmc1k/3_cells/cell1.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell2.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell3.bam',
                f'{data_folder}/10x_pbmc1k/HomoSapiens.GRCh38-2020-A_subsampled.gtf',
            ],
            'expected_fn': f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
        }, remove_res_files=False)

        try:
            dat = pd.read_csv(f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv', delimiter='\t')
            self.assertEqual(len(dat.columns), 4)
            self.assertTrue(f'{data_folder}10x_pbmc1k/3_cells/cell1.bam' in dat.columns)
            self.assertTrue(f'{data_folder}10x_pbmc1k/3_cells/cell2.bam' in dat.columns)
            self.assertTrue(f'{data_folder}10x_pbmc1k/3_cells/cell3.bam' in dat.columns)
        finally:
            os.remove(f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv')

    def test_output_tsv_no_header(self):
        # Backwards compatibility of no header
        self._run({
            'call': [
                self.cmd,
                '-c', f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
                f'{data_folder}/10x_pbmc1k/3_cells/cell1.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell2.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell2.bam',
                f'{data_folder}/10x_pbmc1k/HomoSapiens.GRCh38-2020-A_subsampled.gtf',
            ],
            'expected_fn': f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
        }, remove_res_files=False)
        try:
            dat = pd.read_csv(f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv', delimiter='\t', header=None)
            self.assertEqual(dat.shape[0], 7)
            self.assertEqual(dat.loc[0,0], 'ENSG00000188976')
        finally:
            os.remove(f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv')

    def test_output_tsv_header_and_append(self):
        # Make sure append still works properly
        self._run({
            'call': [
                self.cmd,
                '--with-header', ' ',
                '-c', f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
                f'{data_folder}/10x_pbmc1k/3_cells/cell1.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell2.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell3.bam',
                f'{data_folder}/10x_pbmc1k/HomoSapiens.GRCh38-2020-A_subsampled.gtf',
            ],
            'expected_fn': f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
        }, remove_res_files=False)

        self._run({
            'call': [
                self.cmd,
                '--append-output', ' ',
                '-c', f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
                f'{data_folder}/10x_pbmc1k/3_cells/cell1.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell2.bam',
                f'{data_folder}/10x_pbmc1k/3_cells/cell3.bam',
                f'{data_folder}/10x_pbmc1k/HomoSapiens.GRCh38-2020-A_subsampled.gtf',
            ],
            'expected_fn': f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv',
        }, remove_res_files=False)
        try:
            dat = pd.read_csv(f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv', delimiter='\t')
            self.assertEqual(dat.shape[0], 14)
            self.assertTrue(f'{data_folder}10x_pbmc1k/3_cells/cell1.bam' in dat.columns)
            self.assertTrue(f'{data_folder}10x_pbmc1k/3_cells/cell2.bam' in dat.columns)
            self.assertTrue(f'{data_folder}10x_pbmc1k/3_cells/cell3.bam' in dat.columns)
        finally:
            os.remove(f'{data_folder}/10x_pbmc1k/test_tsv_header.tsv')


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

    @unittest.skipIf(anndata is None, "test case depends on anndata")
    def test_additional_attributes_h5ad(self):
        # Get gene name as additional attr but output as h5ad file.
        self._run({
            'call': [
                self.cmd,
                '-m', 'intersection-nonempty',
                '--nonunique', 'none',
                '--secondary-alignments', 'score',
                '--supplementary-alignments', 'score',
                '--additional-attr', 'gene_name',
                '--additional-attr', 'exon_number',
                '-c', f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes_counts.h5ad',
                f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes.sam',
                f'{data_folder}/10x_pbmc1k/HomoSapiens.GRCh38-2020-A_subsampled.gtf',
                ],
            'expected_fn': f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes_counts.h5ad',
            }, remove_res_files=False)

        # Check gene name and exon number is there
        try:
            dat = anndata.read_h5ad(f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes_counts.h5ad')
            self.assertEqual(dat.var.loc['ENSG00000188976', 'gene_name'], 'NOC2L')
            self.assertEqual(dat.var.loc['ENSG00000251562', 'gene_name'], 'MALAT1')
            self.assertEqual(dat.var.loc['ENSG00000188976', 'exon_number'], '2')
            self.assertEqual(dat.var.loc['ENSG00000251562', 'exon_number'], '2')

            htseq_specific_features = ['__no_feature', '__ambiguous', '__too_low_aQual',
                                       '__not_aligned', '__alignment_not_unique']
            for feat in htseq_specific_features:
                self.assertTrue(math.isnan(dat.var.loc[feat, 'gene_name']))
                self.assertTrue(math.isnan(dat.var.loc[feat, 'exon_number']))

        finally:
            # clean up
            os.remove(data_folder + "/10x_pbmc1k/subsampled_with_missing_barcodes_counts.h5ad")

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

    def test_no_contig_overlap_warning(self):
        """ Ensures that no contig warning is present in stderr
            @see https://github.com/htseq/htseq/issues/63 """

        # We print full file path in the error message, so have to use template & replace
        bam_file = os.path.join(data_folder, 'SRR001432_head_sorted.bam')
        stderr_template = open(f'{data_folder}/no_contig_overlap_warning_stderr.txt').read()
        expected_stderr = stderr_template.format(FILENAME=bam_file)
        self._run({
            'call': [
                self.cmd,
                bam_file,
                f'{data_folder}/bamfile_no_qualities.gtf',
                ],
                'expected_stderr': expected_stderr,
            })

    def test_contig_overlap_no_warning(self):
        """ Ensures that warning is NOT present in stderr when there are contig overlaps
            @see https://github.com/htseq/htseq/issues/63 """

        bam_file = os.path.join(data_folder, 'bamfile_no_qualities.sam')
        expected_stderr = open(f'{data_folder}/contig_overlap_no_warning_stderr.txt').read()
        self._run({
            'call': [
                self.cmd,
                bam_file,
                f'{data_folder}/bamfile_no_qualities.gtf',
                ],
                'expected_stderr': expected_stderr,
            })

    def test_default_feature_type(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                f'{data_folder}/yeast_RNASeq_excerpt.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56_multiple_types.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_default_type.tsv',
            })

    def test_non_default_feature_type(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                '--type', 'pseudogene',
                f'{data_folder}/yeast_RNASeq_excerpt.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56_multiple_types.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_non_default_type.tsv',
            })

    def test_multiple_feature_types(self):
        self._run({
            'call': [
                self.cmd,
                '-c', 'test_output.tsv',
                '--type', 'exon',
                '--type', 'pseudogene',
                f'{data_folder}/yeast_RNASeq_excerpt.sam',
                f'{data_folder}/Saccharomyces_cerevisiae.SGD1.01.56_multiple_types.gtf.gz',
                ],
            'expected_fn': f'{data_folder}/yeast_RNASeq_excerpt_multiple_types.tsv',
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

    def test_missing_barcodes(self):
        """
        When reads are missing cell or UMI barcodes, the counts should exclude these reads
        """
        self._run({
            'call': [
                self.cmd,
                '--mode', 'union',
                '--secondary-alignments', 'ignore',
                '--supplementary-alignments', 'ignore',
                '--stranded', 'yes',
                '--counts_output_sparse',
                '-c', f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes_counts.csv',
                '--samout', f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes_counts.sam',
                f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes.sam',
                f'{data_folder}/10x_pbmc1k/HomoSapiens.GRCh38-2020-A_subsampled.gtf',
                ],
            'expected_fn': f'{data_folder}/10x_pbmc1k/subsampled_with_missing_barcodes_counts.csv',
            }, remove_res_files = False)

        # Check the result. The count csv should not register any reads with low quality
        # The SAM file should assign low quality to reads without cell or UMI barcode

        samfile = None
        try:
            df = pd.read_csv(data_folder + '/10x_pbmc1k/subsampled_with_missing_barcodes_counts.csv',
                             header=None,
                             delimiter='\t')
            df.columns = ['gene_ids', 'count']
            self.assertEqual(int(df.loc[df['gene_ids'] == 'ENSG00000188976']['count'].iloc[0]), 1)
            self.assertEqual(int(df.loc[df['gene_ids'] == 'ENSG00000251562']['count'].iloc[0]), 1)

            other_genes = ['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']

            for g in other_genes:
                self.assertEqual(int(df.loc[df['gene_ids'] == g]['count'].iloc[0]), 0)

            samfile = pysam.AlignmentFile(data_folder + "/10x_pbmc1k/subsampled_with_missing_barcodes_counts.sam")
            read_assignments = {
                "A00228:279:HFWFVDMXX:2:1385:2085:18975": "__too_low_aQual",
                "A00228:279:HFWFVDMXX:1:2158:25464:1313": "ENSG00000251562",
                "A00228:279:HFWFVDMXX:2:1168:11153:29168": "ENSG00000188976",
                "A00228:279:HFWFVDMXX:1:1425:2781:12665": "__too_low_aQual",
                "A00228:279:HFWFVDMXX:2:2322:8621:9157": "__too_low_aQual",
                "A00228:279:HFWFVDMXX:1:2119:24270:19351": "__too_low_aQual",
                "A00228:279:HFWFVDMXX:2:2249:1045:2707": "__too_low_aQual",
                "A00228:279:HFWFVDMXX:2:2401:14172:19210": "__too_low_aQual"
            }

            for read in samfile.fetch(until_eof=True):
                qname_str = read.qname
                tag_val = read.get_tag("XF")
                self.assertEqual(tag_val, read_assignments[qname_str])

        finally:
            # clean up
            if samfile is not None:
                samfile.close()
            os.remove(data_folder + '/10x_pbmc1k/subsampled_with_missing_barcodes_counts.csv')
            os.remove(data_folder + "/10x_pbmc1k/subsampled_with_missing_barcodes_counts.sam")


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
    suite.test_missing_barcodes()
