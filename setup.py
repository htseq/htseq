#!/usr/bin/env python
import os
import sys
from setuptools import setup, Command, Extension
from setuptools.command.build_py import build_py
import numpy


this_directory = os.path.abspath(os.path.dirname(__file__))


def update_version():
    with open(os.path.join(this_directory, "HTSeq", "_version.py"), "rt") as fversion:
        version = fversion.read().strip().split("=")[1].strip().strip('"')

    return version


version = update_version()


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


class Preprocess_command(Command):
    """Cython and SWIG preprocessing"""

    description = "preprocess Cython and SWIG files for HTSeq"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        self.swig_and_cython()

    def swig_and_cython(self):
        import os
        from shutil import copy
        from subprocess import check_call
        from subprocess import SubprocessError

        def c(x):
            return check_call(x, shell=True)

        def p(x):
            return self.announce(x, level=2)

        # CYTHON
        p("cythonizing")
        cython = os.getenv("CYTHON", "cython")
        try:
            c(cython + " --version")
        except SubprocessError:
            if os.path.isfile("src/_HTSeq.c"):
                p("Cython not found, but transpiled file found")
            else:
                raise
        else:
            c(cython + " -3 src/HTSeq/_HTSeq.pyx -o src/_HTSeq.c")

        # SWIG
        p("SWIGging")
        swig = os.getenv("SWIG", "swig")
        pyswigged = "src/StepVector.py"
        try:
            c(swig + " -Wall -c++ -python -py3 src/StepVector.i")
            p("Files transpiled")
        except SubprocessError:
            if os.path.isfile("src/StepVector_wrap.cxx") and os.path.isfile(
                "src/StepVector.py"
            ):
                p("SWIG not found, but transpiled files found")
            else:
                p(
                    "swig not found and traspiled files not found.\n"
                    + "Install SWIG via your package manager (linux) or "
                    + 'via "brew install swig" (OSX - via homebrew)'
                )
                raise
        p("moving swigged .py module")
        copy(pyswigged, "HTSeq/StepVector.py")

        p("done")


class Build_with_preprocess(build_py):
    def run(self):
        self.run_command("preprocess")
        build_py.run(self)


def lazy_numpy_include_dir():
    """Lazily obtain NumPy include directory."""
    try:
        import numpy

        return os.path.join(os.path.dirname(numpy.__file__), "core", "include")
    except ImportError:
        sys.stderr.write(
            "Failed to import 'numpy'. It is required for building HTSeq.\n"
        )
        sys.exit(1)


setup(
    version=version,
    ext_modules=[
        Extension(
            "HTSeq._HTSeq",
            ["src/_HTSeq.c"],
            # include_dirs=[lazy_numpy_include_dir()],#+get_include_dirs(),
            include_dirs=[numpy.get_include()],
            extra_compile_args=["-w"],
        ),
        Extension(
            "HTSeq._StepVector",
            ["src/StepVector_wrap.cxx"],
            # include_dirs=get_include_dirs(cpp=True),
            library_dirs=get_library_dirs_cpp(),
            extra_compile_args=["-w"] + get_extra_args_cpp(),
            extra_link_args=get_extra_args_cpp(),
        ),
    ],
    cmdclass={
        "preprocess": Preprocess_command,
        "build_py": Build_with_preprocess,
    },
)
