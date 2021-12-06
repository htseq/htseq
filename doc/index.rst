.. _tableofcontents:

************************************************************
HTSeq: High-throughput sequence analysis in Python
************************************************************

:Author: Fabio Zanini, Simon Anders and contributors
:Date: |today|
:Version: |version|

HTSeq is a Python package for analysis of high-throughput sequencing data.

* For a high-level description of the package, see the :ref:`overview`.

* For downloads and installation instructions, see :ref:`install`.

* For a thorough example, see :ref:`tour`.

* For tutorials about specific analyses, see :ref:`tutorials`.

* For documentation on `htseq-count`, see :ref:`htseqcount`.
  
* Reference API documentation is available on the other pages. 

Citation
========

**If you use HTSeq in your research, please cite this paper:**

  | G Putri, S Anders, PT Pyl, JE Pimanda, F Zanini
  | **Analysing high-throughput sequencing data with HTSeq 2.0**
  | `arXiv`_ and submitted (2021)

.. _`arXiv`: https://arxiv.org/abs/2112.00939

.. note:: bioRxiv rejected this preprint saying it's not proper research. Thankfully,
   the arXiv is a little more supportive of open source and open science.

The original package is described in:

  | Simon Anders, Paul Theodor Pyl, Wolfgang Huber
  | *HTSeq --- A Python framework to work with high-throughput sequencing data*
  | Bioinformatics (2014), in print, online at `doi:10.1093/bioinformatics/btu638`_

.. _`doi:10.1093/bioinformatics/btu638`: https://doi.org/10.1093/bioinformatics/btu638


Indices and tables
==================
   
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Authors
=======

HTSeq is currently developed by:

* `Givanna Putri`_ at UNSW Sydney (g *dot* putri *at* unsw *dot* edu *dot* au)
* `Simon Anders`_ (anders *at* embl *dot* de) at `EMBL Heidelberg`_ (`Genome Biology Unit`_).
* `Fabio Zanini`_ at `UNSW Sydney`_ (fabio *dot* zanini *at* unsw *dot* edu *dot* au)

License
=======

HTSeq is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The full text of the GNU General Public License, version 3, can be found
here: http://www.gnu.org/licenses/gpl-3.0-standalone.html


.. _`Givanna Putri`: https://fabilab.org/pages/people/html
.. _`Fabio Zanini`: https://fabilab.org
.. _`Simon Anders`: https://www.embl.de/research/units/genome_biology/huber/members/index.php?s_personId=6001
.. _`UNSW Sydney`: https://www.unsw.edu.au/
.. _`EMBL Heidelberg`: https://www.embl.de/
.. _`Genome Biology Unit`: https://www.embl.de/research/units/genome_biology/index.html


Sitemap
=======

.. toctree::
   :glob:
   :maxdepth: 2
   :titlesonly:
   
   Home <self>
   Overview <overview>
   install
   tour
   tutorials
   counting
   refoverview
   sequences
   genomic
   alignments
   features
   otherparsers
   misc
   htseqcount      
   htseqcount_with_barcodes
   qa
   history
   contrib


