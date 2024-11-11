.. _features:

********
Features
********


.. currentmodule:: HTSeq

.. doctest::
   :hide:

   >>> import HTSeq

The easiest way to work with annotation is to use :class:`GenomicArray` with ``typecode=='O'``
or :class:`GenomicArrayOfSets`. If you have your annotation in a flat file, with each
line describing a feature and giving its coordinates, you can read in the file line for line,
parse it (see the standard Python module ``csv``), use the information on chromosome, start,
end and strand to create a :class:`GenomicInterval` object and then store the data from the line
in the genomic array at the place indicated by the genomic interval.

For example, if you have data in a tab-separated file as follows:

.. doctest::

   >>> for line in open("feature_list.txt"):  #doctest:+NORMALIZE_WHITESPACE
   ...     print(line)
   chr2  100	300	+	"gene A"
   chr2	200	400	-	"gene B"
   chr3	150	270	+	"gene C"

Then, you could load this information as follows::

   >>> import csv
   >>> genes = HTSeq.GenomicArray(["chr1", "chr2", "chr3"], typecode='O')
   >>> for (chrom, start, end, strand, name) in \
   ...        csv.reader(open("feature_list.txt"), delimiter="\t"):
   ...     iv = HTSeq.GenomicInterval(chrom, int(start), int(end), strand)
   ...     genes[iv] = name

Now, to see whether there is a feature at a given :class:`GenomicPosition`, you just query the
genomic array::

   >>> print(genes[HTSeq.GenomicPosition("chr3", 100, "+")])
   None
   >>> print(genes[HTSeq.GenomicPosition("chr3", 200, "+")])
   gene C

See :class:`GenomicArray` and :class:`GenomicArrayOfSets` for more sophisticated use.


``GFF_Reader`` and ``GenomicFeature``
=====================================

One of the most common format for annotation data is GFF_ (which includes GTF_ as
a sub-type). Hence, a parse for GFF files is included in HTSeq.

.. _GFF: http://www.sanger.ac.uk/resources/software/gff/spec.html
.. _GTF: http://mblab.wustl.edu/GTF22.html

As usual, there is a parser class, called **GFF_Reader**, that can generate an
iterator of objects describing the features. These objects are of type :class`GenomicFeature`
and each describes one line of a GFF file. See Section :ref:`tour` for an example.

.. class:: GFF_Reader(filename_or_sequence, end_included=True)

   As a subclass of :class:`FileOrSequence`, GFF_Reader can be initialized either
   with a file name or with an open file or another sequence of lines.

   When requesting an iterator, it generates objects of type :class:`GenomicFeature`.

   The GFF specification is unclear on whether the end coordinate marks the last
   base-pair of the feature (closed intervals, ``end_included=True``) or the one
   after (half-open intervals, ``end_included=False``). The default, True, is
   correct for Ensembl GTF files. If in doubt, look for a CDS or stop_codon
   feature in you GFF file. Its length should be divisible by 3. If "end-start"
   is divisible by 3, you need ``end_included=False``. If "end-start+1" is
   divisible by 3, you need ``end_included=True``.

   GFF_Reader will convert the coordinates from GFF standard (1-based, end
   maybe included) to HTSeq standard (0-base, end not included) by subtracting
   1 from the start position. This is also Python's indexing standard. If
   ``end_included=False``, the end was one-after already in the GFF, so HTSeq
   will also subtract 1 from the end position.

      .. attribute:: GFF_Reader.metadata

         GFF_Reader skips all lines starting with a single '#' as this marks
         a comment. However, lines starting with '##' contain meta data (at least
         according to the Sanger Institute's version of the GFF standard.) Such meta
         data has the format ``##key value``. When a metadata line is encountered,
         it is added to the ``metadata`` dictionary.


.. class:: GenomicFeature(name, type_, interval)

   A GenomicFeature object always contains the following attributes:

      .. attribute:: GenomicFeature.name

         A name of ID for the feature. As the GFF format does not have a dedicated
         field for this, the value of the first attribute in the *attributes* column is
         assumed to be the name of ID.

      .. attribute:: GenomicFeature.type

         The type of the feature, i.e., a string like ``"exon"`` or ``"gene"``. For GFF
         files, the 3rd column (*feature*) is taken as the type.

      .. attribute:: GenomicFeature.interval

         The interval that the feature covers on the genome. For GFF files, this information is taken
         from the first (*seqname*), the forth (*start*), the fifth (*end*), and the seventh (*strand*)
         column.

   When created by a :class:`GFF_Reader` object, the following attributes are also present, with the information
   from the remaining GFF columns:

      .. attribute:: GenomicFeature.source

         The 2nd column, denoted *source* in the specification, and intended to specify the
         data source.

      .. attribute:: GenomicFeature.frame

         The 8th column (*frame*), giving the reading frame in case of a coding feature. Its value
         is an integer (0, 1, or 2), or the string ``'.'`` in case that a frame is not specified or would not make sense.

      .. attribute:: GenomicFeature.score

         The 6th column (*score*), giving some numerical score for the feature. Its value
         is a float, or ``'.'`` in case that a score is not specified or would not make sense

      .. attribute:: GenomicFeature.attr

         The last (9th) column of a GFF file contains *attributes*, i.e. a list of name/value pairs.
         These are transformed into a dict, such that, e.g., ``gf.attr['gene_id']`` gives the value
         of the attribute ``gene_id`` in the feature described by ``GenomicFeature`` object ``gf``.
         The parser for the attribute field is reasonably flexible to deal with format variations
         (it was never clearly established whether name and value should be separated by a colon or an
         equals sign, and whether quotes need to be used) and also does a URL style decoding, as is often
         required.

      .. attribute:: GenomicFeature.attr_tuples

		Same as ``attr`` above, but a list of tuples, so as to be able to handle multiple entries with the
		same key.

   In order to write a GFF file from a sequence of features, this method is provided:

      .. method:: GenomicFeature.get_gff_line(with_equal_sign=False)

         Returns a line to describe the feature in the GFF format. This works even if the optional
         attributes given above are missing. Call this for each of your ``GenomicFeature`` objects
         and write the lines into a file to get a GFF file.

.. function:: parse_GFF_attribute_string(attrStr, extra_return_first_value=False)

   This is the function that :class:`GFF_Reader` uses to parse the attribute column. (See :attr:`GenomicFeature.attr`.)
   It returns a dict, or, if requested, a pair of the dict and the first value.

.. function::
   make_feature_genomicarrayofsets(feature_sequence, id_attribute, feature_type=None, feature_query=None, additional_attributes=None, stranded=False, verbose=False, add_chromosome_info=False)

   Organize a sequence of Feature objects into a GenomicArrayOfSets.

   :param feature_sequence: A sequence of features, e.g. as obtained from GFF_reader('myfile.gtf')
   :type feature_sequence: iterable of Feature

   :param id_attribute:
     An attribute to use to identify the feature in the output data structures (e.g.
     'gene_id'). If this is a list, the combination of all those attributes, separated
     by colons (:), will be used as an identifier. For instance,
     ['gene_id', 'exon_number'] uniquely identifies specific exons.
   :type id_attribute: string or sequence of strings

   :param feature_type:
     If None, collect all features. If a list of strings, restrict to only
     those types of features e.g. 'exon' or 'gene' and 'pseudogene'.
   :type feature_type: list or None

   :param feature_query:
     If None, all features of the selected types will be collected. If a string, it has
     to be in the format: <feature_attribute> == <attr_value>, e.g. 'gene_id == "Fn1"'
     (note the double quotes inside). Then only that feature will be collected. Using
     this argument is more efficient than collecting all features and then pruning it
     down to a single one.
   :type feature_query: string or None

   :param additional_attributes:
     A list of additional attributes to be collected into a separate dict for the same
     features, for instance ['gene_name']
   :type additional_attributes: list or None

   :param bool stranded: Whether to keep strandedness information
   :param bool verbose:  Whether to output progress and error messages
   :param bool add_chromosome_info:
      Whether to add chromosome information for each feature. If this option is True,
      the function appends at the end of the "additional_attributes" list a
      "Chromosome" attribute.

   :return: A dict with two keys, 'features' with the GenomicArrayOfSets populated
     with the features, and 'attributes' which is itself a dict with
     the id_attribute as keys and the additional attributes as values.

   Example: Let's say you load the C. elegans GTF file from Ensembl and make a
   feature dict:

   >>> gff = HTSeq.GFF_Reader("Caenorhabditis_elegans.WS200.55.gtf.gz")  #doctest: +SKIP
   >>> worm_features = HTSeq.make_feature_genomicarrayofsets(gff, 'gene_id')  #doctest: +SKIP

   (This command may take a few minutes to deal with the 430,000 features
   in the GTF file. Note that you may need a lot of RAM if you have millions
   of features.)

   This function is related but distinct from ``HTSeq.make_feature_dict``. This
   function is used in htseq-count and its barcoded twin to count gene
   expression because the output GenomicArrayOfSets is very efficient. You
   can use it in performance-critical scans of GFF files.
