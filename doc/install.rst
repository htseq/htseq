.. _install:

****************************
Installation
****************************

HTSeq is available from the `Python Package Index (PyPI)`_:

To use HTSeq, you need Python_ 3.7 or above with:

- NumPy_, a commonly used Python package for numerical calculations
- Pysam_, a Python interface to samtools_.
- To make plots you will need matplotlib_, a plotting library. 

At the moment, HTSeq supports Linux and OSX but not Windows operating systems,
because one of the key dependencies, Pysam_, lacks automatic support and none
of the HTSeq authors have access to such a machine. However, it *might* work
with some work, if you need support for this open an issue on our Github_ page.

HTSeq follows install conventions of many Python packages. In the best case, it
should install from PyPI like this::

 pip install HTSeq

If this does not work, please open an issue on Github_ and also try the instructions
below.

.. note::

  ``pip`` will occasionally try to build HTSeq from source. If that's the case, you
  need two additional dependencies: Cython_ and SWIG_. See the Github_ page for more info.

.. _`Python Package Index (PyPI)`: http://pypi.python.org/pypi/HTSeq
.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _Pysam: https://github.com/pysam-developers/pysam
.. _samtools: http://www.htslib.org/
.. _matplotlib: http://matplotlib.org/
.. _Cython: https://cython.org/
.. _SWIG: http://swig.org/
.. _Github: https://github.com/htseq/htseq


Installation on Linux
=====================

You can choose to install HTSeq via your distribution packages or via `pip`. The former
is generally recommended but might be updated less often than the `pip` version.

Distribution package manager
----------------------------

- Ubuntu::

   sudo apt-get install build-essential python3.6-dev python-numpy python-matplotlib python-pysam python-htseq

- Arch (e.g. using ``aura``, you can grab the AUR packages manually otherwise)::

    sudo pacman -S python python-numpy python-matplotlib
    sudo aura -A python-pysam python-htseq

PIP
---
PIP should take care of the requirements for you::

  pip install HTSeq


Installing from GIT
-------------------
If you want to install a development version, just clone the git repository, switch to the branch/commit
you wish, and use ``setuptools``::

   python setup.py build
   python setup.py install

Typical setuptools options are available (e.g. ``--prefix``, ``--user``).

To test the installation, change to another director than the build directory, start Python
(by typing ``python`` or ``python2.7``) and then try whether typing ``import HTSeq`` causes an error meesage.

Remember that Python can only import from your ``PYTHONPATH`` list of folders and from direct subfolders of the current working directory: if you have a folder called ``HTSeq`` in the current directory, your import will likely fail: move somewhere else.

Installation on MacOS X
=======================
PIP
---
Try pip first::

  pip install HTSeq

Installing from source
----------------------
Mac users should install NumPy as explained here_ in the NumPy/SciPy documentation. Note that you need
to install Xcode to be able to compile NumPy. Due to the
mess that Apple recently made out of Xcode, the whole process may be a slight bit more cumbersome than necessary, especially if you work with OSX Lion, so read the instructions carefully.

.. _here: http://www.scipy.org/Installing_SciPy/Mac_OS_X

If you want to produce plots or use htseq-qa, you will also need matplotlib. (For htseq-count, it
is not required.) There seems to be a binary package (a "Python egg") available on the matplotlib
SourceForge page.

To install HTSeq itself, download the *source* package from the `HTSeq PyPI page`_, unpack the tarball,
go into the directory with the unpacked files and type there:

.. _`HTSeq PyPI page`: http://pypi.python.org/pypi/HTSeq

::

   python setup.py build

to compile HTSeq. If you get an error regarding the availability of a C compiler, you may need to
set environment variables to point Python to it. The NumPy/SciPy installation instructions above cover this topic well and
apply here, too, so simply do the same as you did to install NumPy.

Once building has been successful, use::

   python setup.py --user

to install HTSeq for the current users. To make HTSeq available to all users, use instead::

   python setup.py build
   sudo python setup.py install

To test the installation, change to another director than the build directory, start Python
(by typing ``python``) and then try whether typing ``import HTSeq`` causes an error meesage.


MS Windows
==========
If you have not yet installed Python, do so first. You can find an automatic installer 
for Windows on the `Python download page`_. Make sure to use Python 3.5 or above. 

.. _`Python download page`: http://www.python.org/getit/

Then install the newest version of NumPy. Look on `NumPy's PyPI page`_ for the automatic installer.

.. _`NumPy's PyPI page`: https://pypi.python.org/pypi/numpy

If you want to produce plots or use htseq-qa, you will also need matplotlib. (For htseq-count, it
is not required.) Follow the installation instructions on their web page.

Installation from source
------------------------
**Installation on Windows is not currently supported.** The notes below are left as a historical
record in case a generous soul wants to help us bring back HTSeq to this OS.

Old notes
---------
To install HTSeq itself, simply download the Windows installer from the `HTSeq download page`_
and run it.

.. _`HTSeq download page`: http://pypi.python.org/pypi/HTSeq

To test your installation, start Python and then try whether typing ``import HTSeq`` 
causes an error meesage.

If you get the error message "ImportError: DLL load failed", you are most likely
missing the file MSVCR110.DLL on your system, which you can get by downloading and
installing the file "vcredist_x86.exe" from `this page`_.

.. _`this page`: http://www.microsoft.com/en-us/download/details.aspx?id=30679

