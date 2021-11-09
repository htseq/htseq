# Internal HTSeq functions, not part of the API
import HTSeq
import numpy


def GenomicInterval_range(gi, step):
    for pos in range(gi.start, gi.end, step):
        yield HTSeq.GenomicPosition(gi.chrom, pos, gi.strand)


def GenomicInterval_xranged(gi, step):
    if gi.strand == "-":
        step *= -1
    for pos in range(gi.start_d, gi.end_d, step):
        yield HTSeq.GenomicPosition(gi.chrom, pos, gi.strand)


def ChromVector_steps(cv):
    '''Steps over a ChromVector


    NOTE: ChromVectors use an offset which is also iv.start compared to their
          storage objects that start at 0.
    '''
    # "Steps" of an ndarray (or memmap?)-storaged ChromVector
    if isinstance(cv.array, numpy.ndarray):
        start = cv.iv.start
        prev_val = None
        for i in range(cv.iv.start, cv.iv.end):
            val = cv.array[i - cv.offset]
            if prev_val is None or val != prev_val:
                if prev_val is not None:
                    yield (HTSeq.GenomicInterval(cv.iv.chrom, start, i, cv.iv.strand), prev_val)
                prev_val = val
                start = i
        yield (HTSeq.GenomicInterval(
            cv.iv.chrom, start, cv.iv.end, cv.iv.strand), prev_val,
            )

    # Steps of a StepVector-storaged ChromVector
    elif isinstance(cv.array, HTSeq.StepVector.StepVector):
        for start, stop, value in cv.array[
                cv.iv.start - cv.offset: cv.iv.end - cv.offset].get_steps():
            yield (HTSeq.GenomicInterval(
                cv.iv.chrom, start + cv.offset, stop + cv.offset, cv.iv.strand), value,
                )

    # Steps in a StretchVector behave similar to a full numpy array, but uses
    # np.nan for None and treats each stretch as independent, of course
    # NOTE: one could optimize this by using np.diff and flips, as we have done
    # in the StretchVector methods. For now, we leave it like this because
    # the whole point of StretchVector is to be used for stretches, not steps.
    elif isinstance(cv.array, HTSeq.StretchVector):
        for iv, stretch in cv.array:
            start = cv.offset + iv.start
            prev_val = None
            for i, val in enumerate(stretch, cv.offset + iv.start):
                # Subsequent NaNs, ignore
                if (prev_val is not None) and numpy.isnan(prev_val) and numpy.isnan(val):
                    continue
                if prev_val is None or val != prev_val:
                    # Delay yield of the first item until you meet the first
                    # unequal item, i.e. until you know the end of the step
                    if prev_val is not None:
                        yield (HTSeq.GenomicInterval(
                                    cv.iv.chrom,
                                    start,
                                    i,
                                    cv.iv.strand),
                                prev_val)
                    prev_val = val
                    start = i
            yield (HTSeq.GenomicInterval(
                cv.iv.chrom, start, cv.offset + iv.end, cv.iv.strand), prev_val,
                )
    else:
        raise SystemError("Unknown array type.")


def GenomicArray_steps(ga):
    """Steps of a GenomicArray are just the chained steps of each ChromVector"""
    for chrom, chromstrand_dict in ga.chrom_vectors.items():
        for strand, chrom_vector in chromstrand_dict.items():
            for iv, val in chrom_vector.steps():
                yield iv, val
