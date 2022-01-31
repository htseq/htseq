'''GFF format utilities'''
import itertools
import warnings
import os
import shlex
import sys

import HTSeq
from HTSeq._HTSeq import *
from HTSeq.utils import FileOrSequence


# GFF regular expressions for cache
_re_attr_main = re.compile("\s*([^\s\=]+)[\s=]+(.*)")
_re_attr_empty = re.compile("^\s*$")
_re_gff_meta_comment = re.compile("##\s*(\S+)\s+(\S*)")



class GenomicFeature(object):
    """A genomic feature, i.e., an interval on a genome with metadata.

    At minimum, the following information should be provided by slots:

      name: a string identifying the feature (e.g., a gene symbol)
      type: a string giving the feature type (e.g., "gene", "exon")
      iv: a GenomicInterval object specifying the feature locus
    """

    def __init__(self, name, type_, interval):
        self.name = name
        self.type = sys.intern(type_)
        self.iv = interval

    def __repr__(self):
        return "<%s: %s '%s' at %s: %d -> %d (strand '%s')>" % \
            (self.__class__.__name__, self.type, self.name,
             self.iv.chrom, self.iv.start_d, self.iv.end_d, self.iv.strand)

    def __eq__(self, other):
        if not isinstance(other, GenomicFeature):
            return False
        return self.name == other.name and self.type == other.type and \
            self.iv == other.iv

    def __neq__(self, other):
        if not isinstance(other, GenomicFeature):
            return True
        return not self.__eq__(other)

    def __hash__(self):
        return (self.name, self.type, self.iv).__hash__()

    def get_gff_line(self, with_equal_sign=False):
        try:
            source = self.source
        except AttributeError:
            source = "."
        try:
            score = self.score
        except AttributeError:
            score = "."
        try:
            frame = self.frame
        except AttributeError:
            frame = "."
        try:
            attr = self.attr
        except AttributeError:
            attr = {'ID': self.name}
        if with_equal_sign:
            sep = "="
        else:
            sep = " "
        attr_str = '; '.join(
            ['%s%s\"%s\"' % (ak, sep, attr[ak]) for ak in attr])
        return "\t".join(str(a) for a in (self.iv.chrom, source,
                         self.type, self.iv.start + 1, self.iv.end, score,
                         self.iv.strand, frame, attr_str)) + "\n"







class GFF_Reader(FileOrSequence):
    """Parse a GFF file

    Pass the constructor either a file name or an iterator of lines of a
    GFF files. If a file name is specified, it may refer to a gzip compressed
    file.

    Iterating over the object then yields GenomicFeature objects.

    Args:
        filename_or_sequence: input file or iterator of lines
        end_included: whether the end coordinate of intervals is included in
          the interval itself. This is common in GTF but not the Python
          standard, hence this argument.
        gff_version: Which version of the GFF format to use (2 or 3). The None
          default has the following meaning. If the input is a filename, use
          version 2 if it ends with gtf or gtf.gz (case insensitive), else use
          version 3. If an iterator, use version 2 by default. Notice that GFF3
          does not use quotes in gene names et similia, while GTF does.
    """

    def __init__(self, filename_or_sequence, end_included=True, gff_version=None):
        super().__init__(filename_or_sequence)
        self.end_included = end_included
        self.metadata = {}
        if gff_version is None:
            self._guess_gff_version()

    def _guess_gff_version(self):
        if not self.fos_is_path:
            gff_version = 2
        else:
            fos = os.fspath(self.fos)
            if fos.lower().endswith((".gtf.gz", ".gtf.gzip", ".gtf")):
                gff_version = 2
            else:
                gff_version = 3
        self.gff_version = gff_version

    def __iter__(self):
        for line in super().__iter__():
            if isinstance(line, bytes):
                line = line.decode()
            if line == "\n":
                continue
            if line.startswith('#'):
                if line.startswith("##"):
                    mo = _re_gff_meta_comment.match(line)
                    if mo:
                        self.metadata[mo.group(1)] = mo.group(2)
                continue
            (seqname, source, feature, start, end, score,
             strand, frame, attributeStr) = line.split("\t", 8)
            (attr, name) = self.parse_GFF_attribute_string(
                    attributeStr,
                    True,
                    self.gff_version,
                    )
            iv = GenomicInterval(
                    seqname,
                    int(start) - 1, int(end) - 1 + int(self.end_included),
                    strand)
            f = GenomicFeature(name, feature, iv)
            if score != ".":
                score = float(score)
            if frame != ".":
                frame = int(frame)
            f.source = source
            f.score = score
            f.frame = frame
            f.attr = attr
            yield f

    @staticmethod
    def parse_GFF_attribute_string(
            attrStr,
            extra_return_first_value=False,
            gff_version=2,
        ):
        """Parses a GFF attribute string and returns it as a dictionary.
    
        If 'extra_return_first_value' is set, a pair is returned: the dictionary
        and the value of the first attribute. This might be useful if this is the
        ID.

        Args:
            attrStr: the GFF attribute string to parse
            extra_return_first_value: whether to return the pair explained above
            gff_version: which GFF format rules to use (2 or 3)
        """
        if attrStr.endswith("\n"):
            attrStr = attrStr[:-1]
        d = {}
        first_val = "_unnamed_"

        if gff_version == 2:
            iterator = quotesafe_split(attrStr.encode())
        else:
            # GFF3 does not care about quotes
            iterator = attrStr.encode().split(b';')

        for (i, attr) in enumerate(iterator):
            attr = attr.decode()
            if _re_attr_empty.match(attr):
                continue

            if (gff_version == 2) and (attr.count('"') not in (0, 2)):
                raise ValueError(
                    "The attribute string seems to contain mismatched quotes.")
            mo = _re_attr_main.match(attr)
            if not mo:
                raise ValueError("Failure parsing GFF attribute line")
            val = mo.group(2)

            # GFF3 does not split quotes
            if (gff_version == 2) and val.startswith('"') and val.endswith('"'):
                val = val[1:-1]
            d[sys.intern(mo.group(1))] = sys.intern(val)
            if extra_return_first_value and i == 0:
                first_val = val
        if extra_return_first_value:
            return (d, first_val)
        else:
            return d


def _parse_feature_query(feature_query):
    if '"' not in feature_query:
        raise ValueError('Invalid feature query')
    if '==' not in feature_query:
        raise ValueError('Invalid feature query')

    idx_quote1 = feature_query.find('"')
    idx_quote2 = feature_query.rfind('"')
    attr_name = feature_query[idx_quote1+1: idx_quote2]

    idx_equal = feature_query[:idx_quote1].find('==')
    attr_cat = feature_query[:idx_equal].strip()

    return {
        'attr_cat': attr_cat,
        'attr_name': attr_name,
        }


def make_feature_dict(
        feature_sequence,
        feature_type=None,
        feature_query=None,
        ):
    """Organize a sequence of Feature objects into a nested dictionary.

    Args:
        feature_sequence (iterable of Feature): A sequence of features, e.g. as
            obtained from GFF_reader('myfile.gtf')
        feature_type (string or None): If None, collect all features. If a
            string, restrict to only one type of features, e.g. 'exon'.
        feature_query (string or None): If None, all features of the selected
            types will be collected. If a string, it has to be in the format:

        <feature_attribute> == <attr_value>

        e.g.

        'gene_id == "Fn1"'

        (note the double quotes inside).

        Then only that feature will be collected. Using this argument is more
        efficient than collecting all features and then pruning it down to a
        single one.

    Returns:
        dict with all the feature types as keys. Each value is again a dict,
        now of feature names. The values of this dict is a list of features.

    Example: Let's say you load the C. elegans GTF file from Ensembl and make a
    feature dict:

    >>> gff = HTSeq.GFF_Reader("Caenorhabditis_elegans.WS200.55.gtf.gz")
    >>> worm_features_dict = HTSeq.make_feature_dict(gff)

    (This command may take a few minutes to deal with the 430,000 features
    in the GTF file. Note that you may need a lot of RAM if you have millions
    of features.)

    Then, you can simply access, say, exon 0 of gene "F08E10.4" as follows:
    >>> worm_features_dict['exon']['F08E10.4'][0]
    <GenomicFeature: exon 'F08E10.4' at V: 17479353 -> 17479001 (strand '-')>
    """

    if feature_query is not None:
        feature_qdic = _parse_feature_query(feature_query)

    features = {}
    for f in feature_sequence:
        if feature_type in (None, f.type):
            if f.type not in features:
                features[f.type] = {}
            res_ftype = features[f.type]

            if feature_query is not None:
                # Skip the features that don't even have the right attr
                if feature_qdic['attr_cat'] not in f.attr:
                    continue
                # Skip the ones with an attribute with a different name
                # from the query (e.g. other genes)
                if f.attr[feature_qdic['attr_cat']] != feature_qdic['attr_name']:
                    continue

            if f.name not in res_ftype:
                res_ftype[f.name] = [f]
            else:
                res_ftype[f.name].append(f)
    return features


def make_feature_genomicarrayofsets(
        feature_sequence,
        id_attribute,
        feature_type=None,
        feature_query=None,
        additional_attributes=None,
        stranded=False,
        verbose=False,
        add_chromosome_info=False,
        ):
    """Organize a sequence of Feature objects into a GenomicArrayOfSets.

    Args:
        feature_sequence (iterable of Feature): A sequence of features, e.g. as
            obtained from GFF_reader('myfile.gtf')
        id_attribute (string or sequence of strings): An attribute to use to
            identify the feature in the output data structures (e.g.
            'gene_id'). If this is a list, the combination of all those
            attributes, separated by colons (:), will be used as an identifier.
            For instance, ['gene_id', 'exon_number'] uniquely identifies
            specific exons.
        feature_type (string or None): If None, collect all features. If a
            string, restrict to only one type of features, e.g. 'exon'.
        feature_query (string or None): If None, all features of the selected
            types will be collected. If a string, it has to be in the format:

        <feature_attribute> == <attr_value>

        e.g.

        'gene_id == "Fn1"'

        (note the double quotes inside).

        Then only that feature will be collected. Using this argument is more
        efficient than collecting all features and then pruning it down to a
        single one.

        additional_attributes (list or None): A list of additional attributes
            to be collected into a separate dict for the same features, for
            instance ['gene_name']
        stranded (bool): Whether to keep strandedness information
        verbose (bool): Whether to output progress and error messages
        add_chromosome_info (bool): Whether to add chromosome information for
            each feature. If this option is True, the fuction appends at the
            end of the "additional_attributes" list a "Chromosome" attribute.

    Returns:
        dict with two keys, 'features' with the GenomicArrayOfSets populated
        with the features, and 'attributes' which is itself a dict with
        the id_attribute as keys and the additional attributes as values.

    Example: Let's say you load the C. elegans GTF file from Ensembl and make a
    feature dict:

    >>> gff = HTSeq.GFF_Reader("Caenorhabditis_elegans.WS200.55.gtf.gz")
    >>> worm_features = HTSeq.make_feature_genomicarrayofsets(gff)

    (This command may take a few minutes to deal with the 430,000 features
    in the GTF file. Note that you may need a lot of RAM if you have millions
    of features.)

    This function is related but distinct from HTSeq.make_feature_dict. This
    function is used in htseq-count and its barcoded twin to count gene
    expression because the output GenomicArrayofSets is very efficient. You
    can use it in performance-critical scans of GFF files.
    """

    def get_id_attr(f, id_attribute):
        '''Get feature id with a single or multiple attributes'''
        if isinstance(id_attribute, str):
            try:
                feature_id = f.attr[id_attribute]
            except KeyError:
                raise ValueError(
                        "Feature %s does not contain a '%s' attribute" %
                        (f.name, id_attribute))
        else:
            feature_id = []
            for id_attr in id_attribute:
                try:
                    feature_id.append(f.attr[id_attr])
                except KeyError:
                    raise ValueError(
                            "Feature %s does not contain a '%s' attribute" %
                            (f.name, id_attr))
            feature_id = ':'.join(feature_id)
        return feature_id

    if additional_attributes is None:
        additional_attributes = []

    if feature_query is not None:
        feature_qdic = _parse_feature_query(feature_query)

    features = HTSeq.GenomicArrayOfSets("auto", stranded)
    attributes = {}
    i = 0
    try:
        for f in feature_sequence:
            if feature_type in (None, f.type):
                feature_id = get_id_attr(f, id_attribute)

                if stranded and f.iv.strand == ".":
                    raise ValueError(
                            "Feature %s at %s does not have strand information but you are "
                            "using stranded mode. Try with unstrnded mode." %
                            (f.name, f.iv))

                if feature_query is not None:
                    # Skip the features that don't even have the right attr
                    if feature_qdic['attr_cat'] not in f.attr:
                        continue
                    # Skip the ones with an attribute with a different name
                    # from the query (e.g. other genes)
                    if f.attr[feature_qdic['attr_cat']] != feature_qdic['attr_name']:
                        continue

                features[f.iv] += feature_id
                attributes[feature_id] = [
                        f.attr[attr] if attr in f.attr else ''
                        for attr in additional_attributes]
                if add_chromosome_info:
                    attributes[feature_id] += [f.iv.chrom]

            i += 1
            if i % 100000 == 0 and verbose:
                if hasattr(feature_sequence, 'get_line_number_string'):
                    msg = "{:d} GFF lines processed.".format(i)
                else:
                    msg = "{:d} features processed.".format(i)
                sys.stderr.write(msg+'\n')
                sys.stderr.flush()
    except(KeyError, ValueError):
        if verbose:
            if hasattr(feature_sequence, 'get_line_number_string'):
                msg = "Error processing GFF file ({:}):".format(
                    feature_sequence.get_line_number_string())
            else:
                msg = "Error processing feature sequence ({:}):".format(
                    str(i+1))
            sys.stderr.write(msg+'\n')
        raise

    if verbose:
        if hasattr(feature_sequence, 'get_line_number_string'):
            msg = "{:d} GFF lines processed.".format(i)
        else:
            msg = "{:d} features processed.".format(i)
        sys.stderr.write(msg+"\n")
        sys.stderr.flush()

    if add_chromosome_info:
        additional_attributes.append('Chromosome')

    return {
        'features': features,
        'attributes': attributes,
        }
