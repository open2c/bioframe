from __future__ import division, print_function
from functools import partial
from six.moves.urllib.parse import urljoin, urlencode
import urllib
import posixpath as pp
import os.path as op
import pandas as pd
import requests
import glob

from .formats import (
    read_table,
    read_chromsizes,
    read_gapfile,
    read_ucsc_mrnafile,
    extract_centromeres,
)


def check_connectivity(reference='http://www.google.com'):
    try:
        urllib.request.urlopen(reference, timeout=1)
        return True
    except urllib.request.URLError:
        return False


def fetch_chromsizes(db, **kwargs):
    """
    Download chromosome sizes from UCSC as a ``pandas.Series``, indexed by
    chromosome label.
    
    Parameters
    ----------
    db : str 
        The name of a UCSC genome assembly.
    name_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
        Each corresponding set of records will be sorted in natural order.
        Default is (r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$').
    all_names : bool, optional
        Whether to return all contigs listed in the file. Default is
        ``False``.

    """
    return read_chromsizes(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/chromInfo.txt.gz'.format(db),
        **kwargs)


def fetch_ucsc_mrna(db, **kwargs):
    return read_ucsc_mrnafile(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/all_mrna.txt.gz'.format(db),
        **kwargs)


def fetch_gaps(db, **kwargs):
    return read_gapfile(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/gap.txt.gz'.format(db),
        **kwargs)


def fetch_cytoband(db, ideo=True, **kwargs):
    if ideo:
        return read_table(
            'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/cytoBandIdeo.txt.gz'.format(db),
            schema='cytoband',
            **kwargs)
    else:
        return read_table(
            'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/cytoBand.txt.gz'.format(db),
            schema='cytoband',
            **kwargs)


def fetch_centxt(db,**kwargs):
    return read_table(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/centromeres.txt.gz'.format(db),
        schema='centxt',
        **kwargs)


def fetch_centromeres(db, merge=True, verbose=False):

    # the priority goes as
    # - Local
    # - centromeres.txt
    # - cytoBandIdeo
    # - cytoBand
    # - gap.txt

    # if db in CENTROMERES:
    #     return CENTROMERES[db]

    if not check_connectivity('http://www.google.com'):
        raise ConnectionError('No internet connection!')

    if not check_connectivity('http://hgdownload.cse.ucsc.edu'):
        raise ConnectionError('No connection to the genome database at hgdownload.cse.ucsc.edu!')

    fetchers = [
        ('centxt', fetch_centxt),
        ('cytoband', fetch_cytoband),
        ('cytoband', partial(fetch_cytoband, ideo=True)),
        ('gap', fetch_gaps)
    ]

    for schema, fetcher in fetchers:
        try:
            df = fetcher(db)
            break
        except urllib.error.HTTPError:
            pass
    else:
        raise ValueError('No source for centromere data found.')

    return extract_centromeres(df, schema=schema, merge=merge)


class EncodeClient:
    BASE_URL = 'http://www.encodeproject.org/'

    def __init__(self, cachedir, assembly, metadata=None):
        self.cachedir = op.join(cachedir, assembly)
        if not op.isdir(self.cachedir):
            raise OSError("Directory doesn't exist: '{}'".format(cachedir))
        if metadata is None:
            metadata_path = op.join(cachedir, 'metadata.tsv')
            self._meta = pd.read_table(metadata_path, low_memory=False)
            self._meta = self._meta[self._meta['Assembly'] == assembly].copy()
            self._meta = self._meta.set_index('File accession')
        else:
            self._meta = metadata

    def _batch_download(self, args):
        params = urlencode(args)
        url = pp.join('batch_download', params)
        url = urljoin(self.BASE_URL, url)
        r = requests.get(url)
        r.raise_for_status()
        return r

    def _metadata(self, args):
        params = urlencode(args)
        url = pp.join('metadata', params, 'metadata.tsv')
        url = urljoin(self.BASE_URL, url)
        r = requests.get(url)
        r.raise_for_status()
        return r

    @property
    def meta(self):
        return self._meta.copy()

    def info(self, accession, width=850, height=450):
        from IPython.display import HTML
        url = urljoin(self.BASE_URL, pp.join('experiments', accession))
        return HTML(
            '<iframe width="{}px" height="{}px" src={}></iframe>'.format(
            width, height, url))

    def fetch(self, accession):
        url = self.meta.loc[accession, 'File download URL']
        #sig = self.meta.loc[accession, 'md5sum']
        filename = op.split(url)[1]
        path = op.join(self.cachedir, filename)
        if op.exists(path):
            print('File "{}" available'.format(filename))
        else:
            print('Downloading "{}"'.format(filename))
            r = requests.get(url)
            r.raise_for_status()
            with open(path, 'wb') as f:
                f.write(r.content)
        return path

    def fetch_all(self, accessions):
        return list(map(self.fetch, accessions))


class FDNClient:
    BASE_URL = 'https://data.4dnucleome.org/'

    def __init__(self, cachedir, assembly, metadata=None):
        self.cachedir = op.join(cachedir, assembly)
        if not op.isdir(self.cachedir):
            raise OSError("Directory doesn't exist: '{}'".format(cachedir))
        if metadata is None:
            metadata_paths = sorted(glob.glob(op.join(cachedir, 'metadata*.tsv')))
            metadata_path = metadata_paths[-1]
            self._meta = pd.read_table(metadata_path, low_memory=False, comment='#')
            #self._meta = self._meta[self._meta['Assembly'] == assembly].copy()
            self._meta = self._meta.set_index('File Accession')
        else:
            self._meta = metadata

    @property
    def meta(self):
        return self._meta.copy()

    def info(self, accession, width=850, height=450):
        from IPython.display import HTML
        url = urljoin(self.BASE_URL, pp.join('experiments', accession))
        return HTML(
            '<iframe width="{}px" height="{}px" src={}></iframe>'.format(
            width, height, url))

    def fetch(self, accession):
        url = self.meta.loc[accession, 'File Download URL']
        #sig = self.meta.loc[accession, 'md5sum']
        filename = op.split(url)[1]
        path = op.join(self.cachedir, filename)
        if op.exists(path):
            print('File "{}" available'.format(filename))
        else:
            print('Downloading "{}"'.format(filename))
            r = requests.get(url)
            r.raise_for_status()
            with open(path, 'wb') as f:
                f.write(r.content)
        return path

    def fetch_all(self, accessions):
        return list(map(self.fetch, accessions))
