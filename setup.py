#!/usr/bin/env python
import os
import sys
import numpy
from setuptools import (
    setup,
    Extension,
)
from setuptools.command.build_py import build_py


# Check OS-specific quirks
def get_library_dirs_cpp():
    """OSX 10.14 and later messed up C/C++ library locations"""
    if sys.platform == "darwin":
        return ["/usr/X11R6/lib"]
    else:
        return []


def get_extra_args_cpp():
    """OSX 101.14 and later refuses to use libstdc++"""
    if sys.platform == "darwin":
        return ["-stdlib=libc++"]
    else:
        return []


class BuildWithPreprocess(build_py):
    """Cython and SWIG preprocessing followed by normal build"""

    def preprocess_cython_swig(self):
        from shutil import copy
        import subprocess as sp

        def check(x):
            return sp.run(x, shell=True, check=True)

        def announce(x):
            return self.announce(x, level=2)

        # CYTHON
        announce("cythonizing")
        cython = os.getenv("CYTHON", "cython")
        try:
            check(cython + " --version")
        except sp.SubprocessError:
            if os.path.isfile("src/_HTSeq.c"):
                announce("Cython not found, but transpiled file found")
            else:
                raise
        else:
            check(cython + " -3 src/HTSeq/_HTSeq.pyx -o src/_HTSeq.c")

        # SWIG
        announce("SWIGging")
        swig = os.getenv("SWIG", "swig")
        pyswigged = "src/StepVector.py"
        try:
            check(swig + " -Wall -c++ -python -py3 src/StepVector.i")
            announce("Files transpiled")
        except sp.SubprocessError:
            if os.path.isfile("src/StepVector_wrap.cxx") and os.path.isfile(
                "src/StepVector.py"
            ):
                announce("SWIG not found, but transpiled files found")
            else:
                announce(
                    "swig not found and traspiled files not found.\n"
                    + "Install SWIG via your package manager (linux) or "
                    + 'via "brew install swig" (OSX - via homebrew)'
                )
                raise
        announce("moving swigged .py module")
        copy(pyswigged, "HTSeq/StepVector.py")

        announce("done")

    def run(self):
        self.preprocess_cython_swig()
        super().run()


setup(
    ext_modules=[
        Extension(
            "HTSeq._HTSeq",
            ["src/_HTSeq.c"],
            include_dirs=[numpy.get_include()],
            extra_compile_args=["-w"],
        ),
        Extension(
            "HTSeq._StepVector",
            ["src/StepVector_wrap.cxx"],
            library_dirs=get_library_dirs_cpp(),
            extra_compile_args=["-w"] + get_extra_args_cpp(),
            extra_link_args=get_extra_args_cpp(),
        ),
    ],
    cmdclass={
        "build_py": BuildWithPreprocess,
    },
)
