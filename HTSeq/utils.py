import itertools
import warnings
import os
import shlex
import sys

import HTSeq
from HTSeq._HTSeq import *
from HTSeq._version import __version__


class FileOrSequence:
    """ The construcutor takes one argument, which may either be a string,
    which is interpreted as a file name (possibly with path), or a
    connection, by which we mean a text file opened for reading, or
    any other object that can provide an iterator over strings
    (lines of the file).

    The advantage of passing a file name instead of an already opened file
    is that if an iterator is requested several times, the file will be
    re-opened each time. If the file is already open, its lines can be read
    only once, and then, the iterator stays exhausted.

    Furthermore, if a file name is passed that end in ".gz" or ".gzip"
    (case insensitive), it is transparently gunzipped.
    """

    def __init__(self, filename_or_sequence):
        self.fos = filename_or_sequence
        self.should_close = False
        self.line_no = None
        self.lines = None

        try:
            os.fspath(self.fos)
            self.fos_is_path = True
        except TypeError:
            # assumed to be a file handle
            self.fos_is_path = False

    def __enter__(self):
        if self.fos_is_path:
            fos = os.fspath(self.fos)
            self.should_close = True
            if fos.lower().endswith((".gz", ".gzip")):
                lines = gzip.open(self.fos, 'rt')
            else:
                lines = open(self.fos)
        else:
            lines = self.fos

        self.lines = lines
        return self

    def __exit__(self, type, value, traceback):
        if self.should_close:
            self.lines.close()
        self.lines = None

    def __iter__(self):
        self.line_no = 1
        if self.lines is None:
            call_exit = True
            self.__enter__()
        else:
            call_exit = False
        lines = self.lines

        try:
            for line in lines:
                yield line
                self.line_no += 1
        finally:
            if call_exit:
                self.__exit__(None, None, None)
        self.line_no = None

    def __repr__(self):
        if isinstance(self.fos, str):
            return "<%s object, connected to file name '%s'>" % (
                self.__class__.__name__, self.fos)
        else:
            return "<%s object, connected to %s >" % (
                self.__class__.__name__, repr(self.fos))

    def get_line_number_string(self):
        if self.line_no is None:
            if isinstance(self.fos, str):
                return "file %s closed" % self.fos
            else:
                return "file closed"
        if isinstance(self.fos, str):
            return "line %d of file %s" % (self.line_no, self.fos)
        else:
            return "line %d" % self.line_no
