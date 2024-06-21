#!/usr/bin/env python
import os
import sys
from setuptools import setup, Command, Extension
from setuptools.command.build_py import build_py
import numpy


this_directory = os.path.abspath(os.path.dirname(__file__))


def update_version():
    import subprocess as sp
    # Try getting the latest release tag
    version = None
    try:
        output = sp.check_output(
                'git describe --tags --abbrev=0',
                shell=True,
                )
        if not isinstance(output, str):
            output = output.decode().strip('\n')
        if output.startswith('release_'):
            version = output.split('_')[1]
        print('VERSION updated: '+version)
    except:
        pass

    # Fallback is reading from VERSION file
    if version is not None:
        with open(os.path.join(this_directory, 'VERSION'), 'wt') as fversion:
            fversion.write(version+'\n')
    else:
        with open(os.path.join(this_directory, 'VERSION')) as fversion:
            version = fversion.readline().rstrip()

    # Update version from VERSION file into module
    with open(os.path.join(this_directory, 'HTSeq', '_version.py'), 'wt') as fversion:
        fversion.write('__version__ = "'+version+'"')

    return version


if ((sys.version_info[0] == 2) or
   (sys.version_info[0] == 3 and sys.version_info[1] < 7)):
    sys.stderr.write("Error in setup script for HTSeq:\n")
    sys.stderr.write("HTSeq requires Python 3.7+.")
    sys.exit(1)


# Get version
version = update_version()

# Get README file content
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

# Check OS-specific quirks
# NOTE: setuptools < 18.0 has issues with Cython as a dependency
# NOTE: old setuptools < 18.0 has issues with extras
kwargs = dict(
    setup_requires=[
          'Cython',
          'numpy',
          'pysam',
    ],
    install_requires=[
        'numpy',
        'pysam',
    ],
    extras_require={
        'htseq-qa': ['matplotlib>=1.4'],
        'test': [
            'scipy>=1.5.0',
            'pytest>=6.2.5',
            'pandas>=1.1.0',
            'matplotlib>=1.4',
        ],
    },
  )

def get_library_dirs_cpp():
    '''OSX 10.14 and later messed up C/C++ library locations'''
    if sys.platform == 'darwin':
        return ['/usr/X11R6/lib']
    else:
        return []


def get_extra_args_cpp():
    '''OSX 101.14 and later refuses to use libstdc++'''
    if sys.platform == 'darwin':
        return ['-stdlib=libc++']
    else:
        return []


class Preprocess_command(Command):
    '''Cython and SWIG preprocessing'''
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

        def c(x): return check_call(x, shell=True)
        def p(x): return self.announce(x, level=2)

        # CYTHON
        p('cythonizing')
        cython = os.getenv('CYTHON', 'cython')
        try:
            c(cython+' --version')
        except SubprocessError:
            if os.path.isfile('src/_HTSeq.c'):
                p('Cython not found, but transpiled file found')
            else:
                raise
        else:
            c(cython+' -3 src/HTSeq/_HTSeq.pyx -o src/_HTSeq.c')

        # SWIG
        p('SWIGging')
        swig = os.getenv('SWIG', 'swig')
        pyswigged = 'src/StepVector.py'
        try:
            c(swig+' -Wall -c++ -python -py3 src/StepVector.i')
            p('Files transpiled')
        except SubprocessError:
            if (os.path.isfile('src/StepVector_wrap.cxx') and
                    os.path.isfile('src/StepVector.py')):
                p('SWIG not found, but transpiled files found')
            else:
                p('swig not found and traspiled files not found.\n' +
                  'Install SWIG via your package manager (linux) or ' +
                  'via "brew install swig" (OSX - via homebrew)')
                raise
        p('moving swigged .py module')
        copy(pyswigged, 'HTSeq/StepVector.py')

        p('done')


class Build_with_preprocess(build_py):
    def run(self):
        self.run_command('preprocess')
        build_py.run(self)


def lazy_numpy_include_dir():
    """Lazily obtain NumPy include directory."""
    try:
        import numpy
        return os.path.join(os.path.dirname(numpy.__file__), 'core', 'include')
    except ImportError:
        sys.stderr.write("Failed to import 'numpy'. It is required for building HTSeq.\n")
        sys.exit(1)


setup(name='HTSeq',
      version=version,
      author='Simon Anders, Fabio Zanini',
      author_email='fabio.zanini@unsw.edu.au',
      maintainer='Fabio Zanini',
      maintainer_email='fabio.zanini@unsw.edu.au',
      url='https://github.com/htseq',
      description="A framework to process and analyze data from " +
                  "high-throughput sequencing (HTS) assays",
      long_description=long_description,
      long_description_content_type='text/markdown',
      license='GPL3',
      classifiers=[
         'Development Status :: 5 - Production/Stable',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'License :: OSI Approved :: GNU General Public License (GPL)',
         'Operating System :: POSIX',
         'Programming Language :: Python'
      ],
      ext_modules=[
         Extension(
             'HTSeq._HTSeq',
             ['src/_HTSeq.c'],
             #include_dirs=[lazy_numpy_include_dir()],#+get_include_dirs(),
             include_dirs=[numpy.get_include()],
             extra_compile_args=['-w']),
         Extension(
             'HTSeq._StepVector',
             ['src/StepVector_wrap.cxx'],
             #include_dirs=get_include_dirs(cpp=True),
             library_dirs=get_library_dirs_cpp(),
             extra_compile_args=['-w'] + get_extra_args_cpp(),
             extra_link_args=get_extra_args_cpp(),
             ),
      ],
      py_modules=[
         'HTSeq._HTSeq_internal',
         'HTSeq.StepVector',
         'HTSeq.StretchVector',
         'HTSeq._version',
         'HTSeq.scripts.qa',
         'HTSeq.scripts.count',
         'HTSeq.scripts.count_with_barcodes',
         'HTSeq.scripts.utils',
         'HTSeq.utils',
         'HTSeq.features',
         'HTSeq.scripts.count_features.count_features_per_file',
         'HTSeq.scripts.count_features.reads_io_processor',
         'HTSeq.scripts.count_features.reads_stats',
      ],
      scripts=[
         'scripts/htseq-qa',
         'scripts/htseq-count',
         'scripts/htseq-count-barcodes',
      ],
      cmdclass={
          'preprocess': Preprocess_command,
          'build_py': Build_with_preprocess,
          },
      **kwargs
      )
