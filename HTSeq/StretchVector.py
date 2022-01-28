from collections import namedtuple
import numpy as np
import warnings


Interval = namedtuple('Interval', ['start', 'end'])


class StretchVector:
    """Sparse representation for 'island' of dense data on a long line.

    This class is the most basic object representing 'windows of data' in a sea
    of unknowns of a linear structure, e.g. a chromosome. The basic design is to
    store a list of (start, end) intervals, each of them associated with a
    "stretch", i.e. a dense numpy array.

    In its simplest form, a StretchVector can be thought of as an array with
    an arbitratily large offset. This can be useful e.g. to plot coverage or
    signal enrichments while keeping real genomic coordinates.

    The class does more heavy lifting, however, when several separate arrays
    are used, e.g. to represent "peaks" of signals (ATAC-Seq, ChIP-Seq, etc.).
    In that situation, the bookkeeping involved in manipulating all offsets
    correctly can become burdensome quickly, so StretchVector takes care of it.
    In addition, extension of data to flanks is easy, so if you want to expand
    an existing data window by 1kb on each side, that can be done directly:

    >>> sv = HTSeq.StretchVector(typecode="d")
    >>> sv[6789: 8900] = 56.8
    >>> sv[6000: 7000] = 20 # <-- left extension
    >>> sv[8000: 9000] = 10 # <-- right extension

    You can use StretchVector as a storage option for higher level objects such
    as ChromVector and GenomicArray. Those classes support strandedness, unlike
    StretchVector itself.
    """
    _typecode_dict = {
        "d": np.float32,
        "i": np.int32,
        "l": np.int64,
        "O": object,
    }

    def __init__(self, typecode):
        """Create an empty StretchVector of a given type

        Args:
            typecode ("d", "i", "l", or "O"): The dtype of the stored data. Can
              be "d" (double, i.e. np.float32), "i" (np.int32), "l" (np.int64),
              or "O" (generic object). Note that np.float64 data will be recast
              to np.float32, which might lead to a loss of machine precision.

        Returns:
            A StretchVector instance of the chosen type.
        """
        self.typecode = typecode
        self.ivs = []
        self.stretches = []

    def _in_stretch(self, index):
        if len(self.stretches) == 0:
            return -1

        if isinstance(index, int):
            if (index < self.ivs[0].start) or (index >= self.ivs[-1].end):
                return -1
            for i, iv in enumerate(self.ivs):
                if index < iv.start:
                    return -1
                if index < iv.end:
                    return i

    def _get_interval(self, start, end):
        if len(self.stretches) == 0:
            return self

        ivs = []
        stretches = []
        for i, iv in enumerate(self.ivs):
            # They end before the start, skip
            if iv.end <= start:
                continue

            # They start after the end, skip all remaining
            if iv.start >= end:
                break

            # This interval overlap with start-end
            if (iv.start <= start) and (iv.end <= end):
                new_iv = Interval(start, iv.end)
                new_stretch = self.stretches[i][start - iv.start:]
                ivs.append(new_iv)
                stretches.append(new_stretch)
                continue

            if (iv.start <= start) and (iv.end > end):
                new_iv = Interval(start, end)
                new_stretch = self.stretches[i][start - iv.start:-(iv.end - end)]
                ivs.append(new_iv)
                stretches.append(new_stretch)
                break

            if (iv.start > start) and (iv.end <= end):
                new_iv = Interval(iv.start, iv.end)
                new_stretch = self.stretches[i]
                ivs.append(new_iv)
                stretches.append(new_stretch)
                continue

            if (iv.start > start) and (iv.end > end):
                new_iv = Interval(iv.start, end)
                new_stretch = self.stretches[i][:-(iv.end - end)]
                ivs.append(new_iv)
                stretches.append(new_stretch)
                break

        new_cls = self.__class__(self.typecode)
        new_cls.ivs = ivs
        new_cls.stretches = stretches
        return new_cls

    def _set_interval(self, start, end, values):
        if len(self.stretches) == 0:
            self.ivs.append(Interval(start, end))
            self.stretches.append(
                np.zeros(end - start, self._typecode_dict[self.typecode])
            )
            self.stretches[-1][:] = values
            return len(self.ivs) - 1

        # For each end, there are two possibilities, inside or outside an
        # existing stretch
        idx_start = self._in_stretch(start)
        idx_end = self._in_stretch(end - 1)

        # Neither is in, make a new stretch and delete existing stretches
        if idx_start == idx_end == -1:
            new_iv = Interval(start, end)
            new_stretch = np.zeros(end - start, self._typecode_dict[self.typecode])
            new_stretch[:] = values
            new_ivs = []
            new_stretches = []
            new_added = False
            for i, iv in enumerate(self.ivs):
                # Stretches before
                if start >= iv.end:
                    new_ivs.append(iv)
                    new_stretches.append(self.stretches[i])
                    continue

                # Add new stretch
                if not new_added:
                    new_ivs.append(new_iv)
                    new_stretches.append(new_stretch)

                # Skip overlapping stretches
                if (start <= iv.start) and (end >= iv.end):
                    continue

                # Stretches after
                new_ivs.append(iv)
                new_stretches.append(self.stretches[i])

            # Add new stretch if still missing
            if not new_added:
                new_ivs.append(new_iv)
                new_stretches.append(new_stretch)

        # Start is in a stretch, end is not
        elif (idx_start != -1) and (idx_end == -1):
            new_iv = Interval(self.ivs[idx_start].start, end)
            new_stretch = np.zeros(
                end - self.ivs[idx_start].start,
                self._typecode_dict[self.typecode],
            )
            l1 = start - self.ivs[idx_start].start
            new_stretch[:l1] = self.stretches[idx_start][:l1]
            new_stretch[l1:] = values

            new_ivs = self.ivs[:idx_start] + [new_iv]
            new_stretches = self.stretches[:idx_start] = [new_stretch]

            for i, iv in enumerate(self.ivs[idx_start:], idx_start):
                # Skip the first one
                if i == idx_start:
                    continue

                # Skip overlapping stretches
                if iv.end < end:
                    continue

                # Stretches after
                new_ivs.append(iv)
                new_stretches.append(self.stretches[i])

        # Start is not in a stretch, end is
        elif (idx_start == -1) and (idx_end != -1):
            new_iv = Interval(start, self.ivs[idx_end].end)
            new_stretch = np.zeros(
                self.ivs[idx_end].end - start,
                self._typecode_dict[self.typecode],
            )
            l2 = self.ivs[idx_end].end - end
            if l2 == 0:
                l2 = -len(self.stretches[idx_end])
            new_stretch[:-l2] = values
            if -l2 != len(self.stretches[idx_end]):
                new_stretch[-l2:] = self.stretches[idx_end][-l2:]

            new_ivs = []
            new_stretches = []
            for i, iv in enumerate(self.ivs):
                # Stretches before
                if start >= iv.end:
                    new_ivs.append(iv)
                    new_stretches.append(self.stretches[i])
                    continue
                break

            # New stretch
            new_ivs.append(new_iv)
            new_stretches.append(new_stretch)

            # If there are stretches left, add them
            if idx_end != len(self.ivs) - 1:
                new_ivs.extend(self.ivs[idx_end+1:])
                new_stretches.extend(self.stretches[idx_end+1:])

        # Both start and end are in the same stretch
        elif idx_start == idx_end:
            l1 = start - self.ivs[idx_start].start
            l2 = self.ivs[idx_end].end - end
            if l2 == 0:
                l2 = -len(self.stretches[idx_end])
            new_ivs = self.ivs
            new_stretches = self.stretches
            new_stretches[idx_start][l1:-l2] = values

        # They are in different stretches
        else:
            new_iv = Interval(
                    self.ivs[idx_start].start,
                    self.ivs[idx_end].end,
                )
            new_stretch = np.zeros(
                self.ivs[idx_end].end - self.ivs[idx_start].start,
                self._typecode_dict[self.typecode],
            )
            l1 = start - self.ivs[idx_start].start
            l2 = self.ivs[idx_end].end - end
            if l2 == 0:
                l2 = -len(self.stretches[idx_end])
            new_stretch[:l1] = self.stretches[idx_start][:l1]
            new_stretch[l1:-l2] = values
            if -l2 != len(self.stretches[idx_end]):
                new_stretch[-l2:] = self.stretches[idx_end][-l2:]

            new_ivs = self.ivs[:idx_start] + [iv]
            if idx_end != len(self.ivs) - 1:
                new_ivs.extend(self.ivs[idx_end+1:])
                new_stretches.extend(self.stretches[idx_end+1:])

        self.ivs = new_ivs
        self.stretches = new_stretches

    def _add_stretch(self, start, end, i_start=0):
        for i, iv in enumerate(self.ivs, i_start):
            if start < iv.start:
                self.ivs.insert(
                    i,
                    Interval(start, end),
                )
                self.stretches.insert(
                    i,
                    np.zeros(end - start, self._typecode_dict[self.typecode])
                )
                return i

        self.ivs.append(
            Interval(start, end),
        )
        self.stretches.append(
            np.zeros(end - start, self._typecode_dict[self.typecode])
        )
        return len(self.ivs) - 1

    def __getitem__(self, index):
        """Get a view of a portion of the StretchVector


        Args:
            index (int, slice, or GenomicInterval): Coordinate or interval to
              extract. For slices, the stretches from this intervals are
              *views* of the original array, so changes in them will be
              reflected in the parent StretchVector as well.
        Returns:
            A number of index is an int, containing the value at that site.
            A StretchVector with adapted stretches if index is anything else.

        """
        from HTSeq import GenomicInterval

        if isinstance(index, int):
            idx_iv = self._in_stretch(index)
            if idx_iv == -1:
                return None
            return self.stretches[idx_iv][index - self.ivs[idx_iv].start]

        elif isinstance(index, slice):
            if index.step is not None and index.step != 1:
                raise ValueError(
                        "Striding slices (i.e., step != 1) are not supported")
            if index.start is None:
                index.start = 0
            if index.stop is None:
                if len(self.ivs) == 0:
                    raise IndexError('No stretches, cannot find end')
                index.stop = self.ivs[-1].end

            return self._get_interval(index.start, index.stop)

        elif isinstance(index, GenomicInterval):
            return self.__getitem__(slice(index.start, index.end))

    def __setitem__(self, index, values):
        """Set/reset values within or outside the stretch

        Args:
            index (int, slice, or GenomicInterval): Coordinate or interval to
              set. These can be within current stretches, between them, outside,
              or partially overlapping.
            values (numpy.ndarray or convertible sequence): Values to set or
              reset at those locations. This will be recast as a numpy array
              of the appropriate dtype within this function.

        Returns: None
        """
        from HTSeq import GenomicInterval

        # Leave dtype out for now for speed, it will be taken care of later on
        values = np.asarray(values)

        if isinstance(index, int):
            idx_iv = self._in_stretch(index)
            if idx_iv == -1:
                idx_iv = self._add_stretch(index, index + 1)
            self.stretches[idx_iv][index - self.ivs[idx_iv].start] = values
            return

        elif isinstance(index, slice):
            if index.step is not None and index.step != 1:
                raise ValueError(
                        "Striding slices (i.e., step != 1) are not supported")
            if index.start is None:
                index.start = 0
            if index.stop is None:
                if len(self.ivs) == 0:
                    raise IndexError('No stretches, cannot find end')
                index.stop = self.ivs[-1].end

            self._set_interval(index.start, index.stop, values)

        elif isinstance(index, GenomicInterval):
            return self.__setitem__(slice(index.start, index.end), values)


    def todense(self):
        """Dense numpy array of the whole stretch, using NaNs for missing data"""
        if len(self.ivs) == 0:
            return np.empty(0, self._typecode_dict[self.typecode])

        if len(self.ivs) == 1:
            return self.stretches[0].copy()

        # At least two stretches, have to stitch them
        start = self.ivs[0].start
        res = np.empty(
            self.ivs[-1].end - start,
            self.stretches[0].dtype,
        )
        res[:] = np.nan

        for i, iv in enumerate(self.ivs):
            res[iv.start - start: iv.end - start] = self.stretches[i]

        return res

    @classmethod
    def from_dense(cls, array, offset=0):
        """Create from dense array with NaNs

        Args:
            array (numpy.array): dense array containing NaNs at positions to
              be skipped.
            offset (int): Start of the initial interval.
        """
        if array.dtype == np.float64:
            warnings.warn("np.float64 array converted to np.float32")
            array = array.astype(np.float32)

        for typecode, dtype in cls._typecode_dict.items():
            if dtype == array.dtype:
                sv = cls(typecode=typecode)
                break
        else:
            raise TypeError('Typecode not found for dtype: '+str(array.dtype))

        if len(array) == 0:
            return sv

        flips = np.diff(np.isnan(array)).nonzero()[0]

        # No flips: either all good or all skip
        if flips.sum() == 0:
            if np.isnan(array[0]):
                return sv
            sv.ivs.append(Interval(offset, offset + len(array)))
            sv.stretches.append(array.copy())
            return sv

        # If we start with nan, just increase the offset
        if np.isnan(array[0]):
            add_offset = flips[0] + 1
            offset += add_offset
            array = array[add_offset:]

            # Single flip means keep the whole rest
            if len(flips) > 1:
                flips = flips[1:]
                flips -= add_offset
            else:
                sv.ivs.append(Interval(offset, offset + len(array)))
                sv.stretches.append(array.copy())
                return sv

        # Now we have at least one flip left, and we start with a number/object
        # If we have an odd number of flips, we can forget the last block
        if len(flips) % 2:
            end = flips[-1] + 1
            array = array[:end]
            flips = flips[:-1]

        # No flip left, all good
        if len(flips) == 0:
            sv.ivs.append(Interval(offset, offset + len(array)))
            sv.stretches.append(array.copy())
            return sv

        # Now we start and end with a number/object, and there are at least
        # two flips. Initial stretch
        sv.ivs.append(Interval(offset, offset + flips[0] + 1))
        sv.stretches.append(
                array[:flips[0] + 1]
        )
        # Intermediate stretches
        for i in range((len(flips) // 2) - 1):
            new_iv = Interval(offset + flips[i * 2 + 1] + 1, offset + flips[i * 2 + 2] + 1)
            new_stretch = array[flips[i * 2 + 1] + 1: flips[i * 2 + 2] + 1]
            sv.ivs.append(new_iv)
            sv.stretches.append(new_stretch)
        # Final stretch
        new_iv = Interval(offset + flips[-1] + 1, offset + len(array))
        new_stretch = array[flips[-1] + 1:]
        sv.ivs.append(new_iv)
        sv.stretches.append(new_stretch)

        return sv

    def __iter__(self):
        """Iterate over intervals and stretches ("islands")."""
        return zip(self.ivs, self.stretches)

    def copy(self):
        """Make a copy the StretchVector and of all its stretches"""
        sv_new = StretchVector(typecode=self.typecode)
        for iv, stretch in self:
            sv_new.ivs.append(Interval(iv.start, iv.end))
            sv_new.stretches.append(stretch.copy())
        return sv_new

    def shift(self, offset):
        """Shift all stretch intervals by a constant number

        Args:
            offset (int): Shift the start and end coordinates of each stretch
              by this amount.
        Returns:
            None. The function acts in place.
        """
        for i, iv in self.ivs:
            self.ivs[i] = Interval(iv.start + offset, iv.end + offset)
