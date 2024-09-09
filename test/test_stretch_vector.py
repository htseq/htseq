import pytest
import numpy as np
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


data_folder = conftest.get_data_folder()


class TestStretchVector(unittest.TestCase):
    def test_init(self):
        sv = HTSeq.StretchVector(typecode='d')
        self.assertEqual(sv.ivs, [])
        self.assertEqual(sv.stretches, [])

    def test_setitem_number(self):
        sv = HTSeq.StretchVector(typecode='d')

        # Set initial stretch
        sv[560] = 4.5
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            np.ones(1, np.float32) * 4.5,
        )

        # Overwrite
        sv[560] = 4
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            np.ones(1, np.float32) * 4,
        )

    def test_setitem_slice(self):
        sv = HTSeq.StretchVector(typecode='d')

        # Set initial stretch
        sv[100: 300] = 4.5
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            np.ones(200, np.float32) * 4.5,
        )

        # Set overlapping stretch
        sv[50:250] = 3
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        np.testing.assert_almost_equal(
            sv.stretches[0][:200],
            np.ones(200, np.float32) * 3,
        )
        np.testing.assert_almost_equal(
            sv.stretches[0][200:],
            np.ones(50, np.float32) * 4.5,
        )

        # Set new stretch
        sv[400: 450] = np.arange(50)
        self.assertEqual(len(sv.ivs), 2)
        self.assertEqual(len(sv.stretches), 2)
        np.testing.assert_almost_equal(
            sv.stretches[1],
            np.arange(50).astype(np.float32),
        )

        # Set overlapping stretch
        sv[430: 450] = np.arange(20)
        self.assertEqual(len(sv.ivs), 2)
        self.assertEqual(len(sv.stretches), 2)
        np.testing.assert_almost_equal(
            sv.stretches[1][-20:],
            np.arange(20).astype(np.float32),
        )

    def test_getitem_number(self):
        sv = HTSeq.StretchVector(typecode='d')

        # Set initial stretch
        sv[560] = 4.5
        self.assertEqual(sv[560], 4.5)
        self.assertEqual(sv[580], None)

        sv[400: 450] = np.arange(50)
        res = sv[350: 430]
        self.assertEqual(len(res.ivs), 1)
        self.assertEqual(len(res.stretches), 1)
        np.testing.assert_almost_equal(
            res.stretches[0],
            np.arange(30).astype(np.float32),
        )

    def test_todense(self):
        sv = HTSeq.StretchVector(typecode='d')
        sv[450: 455] = 6.7
        sv[460: 465] = 1.7
        res = sv.todense()
        np.testing.assert_almost_equal(
            res,
            np.array([6.7] * 5 + [np.nan] * 5 + [1.7] * 5).astype(np.float32),
        )

    def test_from_dense(self):
        array = np.empty(20, np.float32)
        array[:] = np.nan

        # All nans
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 0)
        self.assertEqual(len(sv.stretches), 0)

        # All good
        array[:] = 78
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)

        # Single flips
        array[:5] = np.nan
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        self.assertEqual(sv.ivs[0].start, 305)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            78 * np.ones(15).astype(np.float32),
        )

        array[:] = 78
        array[15:] = np.nan
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        self.assertEqual(sv.ivs[0].start, 300)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            78 * np.ones(15).astype(np.float32),
        )

        # Double flip at edges
        array[:] = 78
        array[:4] = np.nan
        array[15:] = np.nan
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 1)
        self.assertEqual(len(sv.stretches), 1)
        self.assertEqual(sv.ivs[0].start, 304)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            78 * np.ones(11).astype(np.float32),
        )

        # Double flip in the middle
        array[:] = 78
        array[2:4] = np.nan
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 2)
        self.assertEqual(len(sv.stretches), 2)
        self.assertEqual(sv.ivs[0].start, 300)
        self.assertEqual(sv.ivs[0].end, 302)
        self.assertEqual(sv.ivs[1].start, 304)
        self.assertEqual(sv.ivs[1].end, 320)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            78 * np.ones(2).astype(np.float32),
        )
        np.testing.assert_almost_equal(
            sv.stretches[1],
            78 * np.ones(16).astype(np.float32),
        )

        # Multiple flips
        array[:] = 78
        array[2:4] = np.nan
        array[8:11] = np.nan
        array[17:19] = np.nan
        sv = HTSeq.StretchVector.from_dense(array, offset=300)
        self.assertEqual(len(sv.ivs), 4)
        self.assertEqual(len(sv.stretches), 4)
        self.assertEqual(sv.ivs[0].start, 300)
        self.assertEqual(sv.ivs[0].end, 302)
        self.assertEqual(sv.ivs[1].start, 304)
        self.assertEqual(sv.ivs[1].end, 308)
        self.assertEqual(sv.ivs[2].start, 311)
        self.assertEqual(sv.ivs[2].end, 317)
        self.assertEqual(sv.ivs[3].start, 319)
        self.assertEqual(sv.ivs[3].end, 320)
        np.testing.assert_almost_equal(
            sv.stretches[0],
            78 * np.ones(2).astype(np.float32),
        )
        np.testing.assert_almost_equal(
            sv.stretches[1],
            78 * np.ones(4).astype(np.float32),
        )
        np.testing.assert_almost_equal(
            sv.stretches[2],
            78 * np.ones(6).astype(np.float32),
        )
        np.testing.assert_almost_equal(
            sv.stretches[3],
            78 * np.ones(1).astype(np.float32),
        )


if __name__ == '__main__':

    suite = TestStretchVector()
    suite.test_init()
