import sys
import numpy as np


class UnknownChrom(Exception):
    pass


def my_showwarning(message, category, filename, lineno=None, file=None,
                   line=None):
    sys.stderr.write("Warning: %s\n" % message)


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


def _merge_counts(
        results,
        attributes,
        additional_attributes,
        sparse=False,
        dtype=np.float32,
        ):
    barcodes = 'cell_barcodes' in results

    if barcodes:
        cbs = results['cell_barcodes']
        counts = results['counts']

    feature_attr = sorted(attributes.keys())
    other_features = [
        ('__no_feature', 'empty'),
        ('__ambiguous', 'ambiguous'),
        ('__too_low_aQual', 'lowqual'),
        ('__not_aligned', 'notaligned'),
        ('__alignment_not_unique', 'nonunique'),
        ]

    fea_names = [fea for fea in feature_attr] + [fea[0] for fea in other_features]
    L = len(fea_names)
    if barcodes:
        n = len(cbs)
    else:
        n = len(results)
    if not sparse:
        table = np.zeros(
            (n, L),
            dtype=dtype,
        )
    else:
        from scipy.sparse import lil_matrix
        table = lil_matrix((n, L), dtype=dtype)

    if not barcodes:
        fea_ids = [fea for fea in feature_attr] + [fea[1] for fea in other_features]
        for j, r in enumerate(results):
            for i, fn in enumerate(fea_ids):
                if i < len(feature_attr):
                    countji = r['counts'][fn]
                else:
                    countji = r[fn]
                if countji > 0:
                    table[j, i] = countji
    else:
        for j, cb in enumerate(cbs):
            for i, fn in enumerate(fea_names):
                countji = counts[cb][fn]
                if countji > 0:
                    table[j, i] = countji

    if sparse:
        table = table.tocsr()

    feature_metadata = {
        'id': fea_names,
    }
    for iadd, attr in enumerate(additional_attributes):
        feature_metadata[attr] = [attributes[fn][iadd] for fn in feature_attr]

    return {
        'feature_metadata': feature_metadata,
        'table': table,
    }


def _count_results_to_tsv(
        results,
        samples_name,
        attributes,
        additional_attributes,
        output_filename,
        output_delimiter,
        output_append=False,
        add_tsv_header=False
        ):

    barcodes = 'cell_barcodes' in results

    pad = ['' for attr in additional_attributes]

    if barcodes:
        cbs = results['cell_barcodes']
        counts = results['counts']

        # Print or write header
        fields = [''] + pad + cbs
        line = output_delimiter.join(fields)
        if output_filename == '':
            print(line)
        else:
            with open(output_filename, 'w') as f:
                f.write(line)
                f.write('\n')

    elif add_tsv_header:
        # Write the header.
        # Only get here if we don't have cell barcodes, i.e. this is not called by htseq-count-barcode,
        # and user wants the tsv header
        file_header = output_delimiter.join([''] + pad + samples_name)

        if output_filename == '':
            print(file_header)
        else:
            # If append to existing file, then open as a
            file_open_opt = 'a' if output_append else 'w'

            with open(output_filename, file_open_opt) as f:
                f.write(file_header)
                f.write('\n')

    # Each feature is a row with feature id, additional attrs, and counts
    feature_attr = sorted(attributes.keys())
    for ifn, fn in enumerate(feature_attr):
        if not barcodes:
            fields = [fn] + attributes[fn] + [str(r['counts'][fn]) for r in results]
        else:
            fields = [fn] + attributes[fn] + [str(counts[cb][fn]) for cb in cbs]

        line = output_delimiter.join(fields)
        if output_filename == '':
            print(line)
        else:
            omode = 'a' if output_append or (ifn > 0) or barcodes or add_tsv_header else 'w'
            with open(output_filename, omode) as f:
                f.write(line)
                f.write('\n')

    # Add other features (unmapped, etc.)
    other_features = [
        ('__no_feature', 'empty'),
        ('__ambiguous', 'ambiguous'),
        ('__too_low_aQual', 'lowqual'),
        ('__not_aligned', 'notaligned'),
        ('__alignment_not_unique', 'nonunique'),
        ]
    for title, fn in other_features:
        if not barcodes:
            fields = [title] + pad + [str(r[fn]) for r in results]
        else:
            fields = [title] + pad + [str(counts[cb][title]) for cb in cbs]
        line = output_delimiter.join(fields)
        if output_filename == '':
            print(line)
        else:
            with open(output_filename, 'a') as f:
                f.write(line)
                f.write('\n')


def _count_table_to_mtx(
        filename,
        table,
        feature_metadata,
        samples,
        ):
    if not str(filename).endswith('.mtx'):
        raise ValueError('Matrix Marker filename should end with ".mtx"')

    try:
        from scipy.io import mmwrite
    except ImportError:
        raise ImportError('Install scipy for mtx support')

    filename_pfx = str(filename)[:-4]
    filename_feature_meta = filename_pfx+'_features.tsv'
    filename_samples = filename_pfx+'_samples.tsv'

    # Write main matrix (features as columns)
    mmwrite(
        filename,
        table,
    )

    # Write input filenames
    with open(filename_samples, 'wt') as fout:
        for fn in samples:
            fout.write(fn+'\n')

    # Write feature metadata (ids and additional attributes)
    with open(filename_feature_meta, 'wt') as fout:
        nkeys = len(feature_metadata)
        for ik, key in enumerate(feature_metadata):
            if ik != nkeys - 1:
                fout.write(key+'\t')
            else:
                fout.write(key+'\n')
        nfeatures = len(feature_metadata[key])
        for i in range(nfeatures):
            for ik, key in enumerate(feature_metadata):
                if ik != nkeys - 1:
                    fout.write(feature_metadata[key][i]+'\t')
                else:
                    fout.write(feature_metadata[key][i]+'\n')


def _count_table_to_h5ad(
        filename,
        table,
        feature_metadata,
        samples,
        ):
    try:
        import anndata
    except ImportError:
        raise ImportError('Install the anndata package for h5ad support')

    # If they have anndata, they have scipy and pandas too
    import pandas as pd

    # We don't have additional attribute (e.g. gene name) for htseq specific features like __no_feature.
    # Hence the trick is to convert the array to series so the value for htseq specific features like __no_feature
    # column is set NaN.
    # See: https://stackoverflow.com/questions/19736080/creating-dataframe-from-a-dictionary-where-entries-have-different-lengths
    feature_metadata = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in feature_metadata.items()]))
    feature_metadata.set_index(feature_metadata.columns[0], inplace=True)

    adata = anndata.AnnData(
        X=table,
        obs=pd.DataFrame([], index=samples),
        var=feature_metadata,
    )
    adata.write_h5ad(filename)


def _count_table_to_loom(
        filename,
        table,
        feature_metadata,
        samples,
        ):

    try:
        import loompy
    except ImportError:
        raise ImportError('Install the loompy package for loom support')

    # Loom uses features as rows...
    layers = {'': table.T}
    row_attrs = feature_metadata
    col_attrs = {'_index': samples}
    loompy.create(
        filename,
        layers=layers,
        row_attrs=row_attrs,
        col_attrs=col_attrs,
    )


def _write_output(
    results,
    samples,
    attributes,
    additional_attributes,
    output_filename,
    output_delimiter,
    output_append,
    sparse=False,
    dtype=np.float32,
    add_tsv_header=False
    ):

    """
    Export the gene counts as tsv/csv, mtx, loom, h5ad files.
    Note, need to update the parameter documentations.

    Parameters
    ----------
    results : list
        List of dictionaries with each element representing the counts for an input BAM file.
        Note, the list is in order of the samples parameter. So the first element in the list corresponds to
        the first file in samples parameter.
    samples : list
        List of input BAM files.

    """

    # Write output to stdout or TSV/CSV
    if output_filename == '':
        _count_results_to_tsv(
            results,
            samples,
            attributes,
            additional_attributes,
            output_filename,
            output_delimiter,
            output_append=False,
            add_tsv_header=add_tsv_header
        )
        return

    # Get file extension/format
    output_sfx = output_filename.split('.')[-1].lower()

    if output_sfx in ('csv', 'tsv'):
        _count_results_to_tsv(
            results,
            samples,
            attributes,
            additional_attributes,
            output_filename,
            output_delimiter,
            output_append,
            add_tsv_header=add_tsv_header
        )
        return

    # Make unified object of counts and feature metadata
    output_dict = _merge_counts(
        results,
        attributes,
        additional_attributes,
        sparse=sparse,
        dtype=dtype,
    )

    if output_sfx == 'mtx':
        _count_table_to_mtx(
            output_filename,
            output_dict['table'],
            output_dict['feature_metadata'],
            samples,
        )
        return

    if output_sfx == 'loom':
        _count_table_to_loom(
            output_filename,
            output_dict['table'],
            output_dict['feature_metadata'],
            samples,
        )
        return

    if output_sfx == 'h5ad':
        _count_table_to_h5ad(
            output_filename,
            output_dict['table'],
            output_dict['feature_metadata'],
            samples,
        )
        return

    raise ValueError(
        f'Format not recognized for output count file: {output_sfx}')
