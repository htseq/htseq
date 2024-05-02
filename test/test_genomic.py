import pytest
import sys
import os
import glob
import sysconfig
from pathlib import Path
import unittest
import pytest
import conftest

build_dir = "build/lib.%s-%s" % (sysconfig.get_platform(), sys.version[0:3])

sys.path.insert(0, os.path.join(os.getcwd(), build_dir))
import HTSeq

try:
    import pyBigWig
except ImportError:
    pyBigWig = None


data_folder = conftest.get_data_folder()


class TestGenomicInterval(unittest.TestCase):
    def test_init(self):
        iv = HTSeq.GenomicInterval(
            'chr1', 100, 200, ".",
        )


class TestChromVector(unittest.TestCase):
    def test_init_step(self):
        # Autoallocation
        cv = HTSeq.ChromVector.create(
            HTSeq.GenomicInterval('chr1', 100, 200, "."),
            typecode='d',
            storage='step',
        )

    def test_init_stretch(self):
        # Autoallocation
        cv = HTSeq.ChromVector.create(
            HTSeq.GenomicInterval('chr1', 100, 200, "."),
            typecode='d',
            storage='stretch',
        )

    def test_steps_stretch(self):
        cv = HTSeq.ChromVector.create(
            HTSeq.GenomicInterval('chr1', 100, 200, "."),
            typecode='d',
            storage='stretch',
        )
        cv[120:130] = 45
        cv[129:140] = 89
        cv[180: 182] = 1
        steps = list(cv.steps())
        self.assertEqual(steps[0][0], HTSeq.GenomicInterval('chr1', 120, 129))
        self.assertEqual(steps[1][0], HTSeq.GenomicInterval('chr1', 129, 140))
        self.assertEqual(steps[2][0], HTSeq.GenomicInterval('chr1', 180, 182))
        self.assertEqual(steps[0][1], 45)
        self.assertEqual(steps[1][1], 89)
        self.assertEqual(steps[2][1], 1)


class TestGenomicArray(unittest.TestCase):
    def test_init(self):
        # Autoallocation
        ga = HTSeq.GenomicArray("auto")

        # Infinite length chromosomes
        ga = HTSeq.GenomicArray(['1', '2'])

        # Fixed chromosomes
        ga = HTSeq.GenomicArray({
            '1': 5898,
            '2': 4876,
        })

        # Store: ndarray
        ga = HTSeq.GenomicArray({
            '1': 5898,
            '2': 4876,
            },
            storage='ndarray',
        )

        # Store: memmap
        ga = HTSeq.GenomicArray({
            '1': 5898,
            '2': 4876,
            },
            storage='memmap',
            memmap_dir='.',
        )

    def test_steps(self):
        for storage in ['step', 'ndarray']:
            ga = HTSeq.GenomicArray({
                '1': 5898,
                '2': 4876,
                },
                storage=storage,
            )
            steps = ga.steps()
            steps_exp = [
                (HTSeq.GenomicInterval('1', 0, 5898, strand='+'), 0),
                (HTSeq.GenomicInterval('1', 0, 5898, strand='-'), 0),
                (HTSeq.GenomicInterval('2', 0, 4876, strand='+'), 0),
                (HTSeq.GenomicInterval('2', 0, 4876, strand='-'), 0),
            ]
            for step, step_exp in zip(steps, steps_exp):
                self.assertEqual(step, step_exp)

    def test_access_out_of_range(self):
        """ Ensure chromosomes are made with infinite size, which can be accessed  """

        def _get_arrays():
            return {
                "Provided chrom": HTSeq.GenomicArray(["1"], typecode='O'),
                "Auto chrom": HTSeq.GenomicArray("auto", typecode='O'),
            }

        # Test we can access unknown regions without error
        unknown_iv = HTSeq.GenomicInterval('1', 200, 300, "+")
        for name, genomic_array in _get_arrays().items():
            step = list(genomic_array[unknown_iv].steps())[0]
            unknown_value = step[1]
            self.assertIsNone(unknown_value, msg="Access unknown in " + name)

        # Test accessing unknown regions works the same even if we call setter first
        known_iv = HTSeq.GenomicInterval('1', 0, 100, strand='+')
        for name, genomic_array in _get_arrays().items():
            # Call setter first before getter
            genomic_array[known_iv] = "test"
            step = list(genomic_array[unknown_iv].steps())[0]
            unknown_value = step[1]
            self.assertIsNone(unknown_value, msg="Access unknown after calling setter in " + name)

    def test_bedgraph(self):
        def compare_bedgraph_line(line1, line2):
            fields1 = line1.split()
            fields2 = line2.split()
            # Chromosome
            self.assertEqual(fields1[0], fields2[0])
            # Start-end
            self.assertEqual(int(fields1[1]), int(fields2[1]))
            self.assertEqual(int(fields1[2]), int(fields2[2]))
            # Value
            self.assertEqual(float(fields1[3]), float(fields2[3]))

        ga = HTSeq.GenomicArray.from_bedgraph_file(
            data_folder+'example_bedgraph.bedgraph',
            strand='.',
        )

        steps = []
        for iv, value in ga.steps():
            steps.append((iv.chrom, iv.start, iv.end, value))

        steps_exp = [
            ('chr19', 49302000, 49302300, -1.0),
            ('chr19', 49302300, 49302600, -0.75),
            ('chr19', 49302600, 49302900, -0.50),
            ('chr19', 49302900, 49303200, -0.25),
            ('chr19', 49303200, 49303500, 0.0),
            ('chr19', 49303500, 49303800, 0.25),
            ('chr19', 49303800, 49304100, 0.50),
            ('chr19', 49304100, 49304400, 0.75),
            ('chr19', 49304400, 49304700, 1.00),
        ]
        self.assertEqual(steps, steps_exp)

        ga.write_bedgraph_file(
            'test_output.bedgraph',
            track_options='name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20',
            separator=' ',
        )
        with open(data_folder+'example_bedgraph.bedgraph') as f1, \
             open('test_output.bedgraph') as f2:
            header_found = False
            for line1, line2 in zip(f1, f2):
                if not header_found:
                    self.assertEqual(line1, line2)
                else:
                    compare_bedgraph_line(line1, line2)
                if 'track type' in line1:
                    header_found = True

    @unittest.skipIf(pyBigWig is None, "test case depends on pyBigWig")
    def test_bigwig(self):
        ga = HTSeq.GenomicArray.from_bigwig_file(
            data_folder+'example_bigwig.bw',
        )
        ga.write_bigwig_file(
            'test_output.bw',
        )

        import pyBigWig
        with pyBigWig.open(data_folder+'example_bigwig.bw') as bw1, \
                pyBigWig.open('test_output.bw') as bw2:
            self.assertEqual(bw1.chroms(), bw2.chroms())
            for chrom in bw1.chroms():
                self.assertEqual(
                    bw1.intervals(chrom),
                    bw2.intervals(chrom),
                )

    def test_access_out_of_range(self):
        """ Ensure chromosomes are made with infinite size, which can be accessed  """

        def _get_arrays():
            return {
                "Provided chrom": HTSeq.GenomicArray(["1"], typecode='O'),
                "Auto chrom": HTSeq.GenomicArray("auto", typecode='O'),
            }

        # Test we can access unknown regions without error
        unknown_iv = HTSeq.GenomicInterval('1', 200, 300, "+")
        for name, genomic_array in _get_arrays().items():
            step = list(genomic_array[unknown_iv].steps())[0]
            unknown_value = step[1]
            self.assertIsNone(unknown_value, msg="Access unknown in " + name)

        # Test accessing unknown regions works the same even if we call setter first
        known_iv = HTSeq.GenomicInterval('1', 0, 100, strand='+')
        for name, genomic_array in _get_arrays().items():
            # Call setter first before getter
            genomic_array[known_iv] = "test"
            step = list(genomic_array[unknown_iv].steps())[0]
            unknown_value = step[1]
            self.assertIsNone(unknown_value, msg="Access unknown after calling setter in " + name)


if __name__ == '__main__':

    suite = TestGenomicArray()
    suite.test_init()
