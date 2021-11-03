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
    else:
        raise SystemError("Unknown array type.")


def GenomicArray_steps(ga):
    for chrom, chromstrand_dict in ga.chrom_vectors.items():
        for strand, chrom_vector in chromstrand_dict.items():
            for iv, val in chrom_vector.steps():
                yield iv, val
