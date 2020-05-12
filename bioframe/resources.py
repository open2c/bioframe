from __future__ import division, print_function
from functools import partial
from six.moves.urllib.parse import urljoin, urlencode
import urllib
import posixpath as pp
import os.path as op
import pandas as pd
import requests
import base64
import glob

import pkg_resources

from .schemas import SCHEMAS
from .formats import (
    read_table,
    read_chromsizes,
    read_gapfile,
    read_ucsc_mrnafile,
    extract_centromeres,
)


LOCAL_CHROMSIZES = {
    path[:-12]: lambda : read_chromsizes(
        pkg_resources.resource_filename(__name__, 'data/'+path))
    for path in pkg_resources.resource_listdir(__name__, 'data/')
    if path.endswith('.chrom.sizes')
}


LOCAL_CENTROMERES = {
    path[:-12]: lambda : read_table(
        pkg_resources.resource_filename(__name__, 'data/'+path),
        SCHEMAS['bed3'])
    for path in pkg_resources.resource_listdir(__name__, 'data/')
    if path.endswith('.centromeres')
}


def _check_connectivity(reference='http://www.google.com'):
    try:
        urllib.request.urlopen(reference, timeout=1)
        return True
    except urllib.request.URLError:
        return False


def fetch_chromsizes(db, provider=None, **kwargs):
    if provider == 'local' or db in LOCAL_CHROMSIZES:
        pass

    if provider == 'ucsc' or provider is None:
        return UCSCClient(db).fetch_chromsizes(**kwargs)
    else:
        raise ValueError("Unknown provider '{}'".format(provider))


def fetch_centromeres(db, provider=None, merge=True, verbose=False):

    # the priority goes as
    # - Local
    # - centromeres.txt
    # - cytoBandIdeo
    # - cytoBand
    # - gap.txt

    # if db in CENTROMERES:
    #     return CENTROMERES[db]

    if not _check_connectivity('http://www.google.com'):
        raise ConnectionError('No internet connection!')

    if not _check_connectivity('http://hgdownload.cse.ucsc.edu'):
        raise ConnectionError('No connection to the genome database at hgdownload.cse.ucsc.edu!')

    client = UCSCClient(db)
    fetchers = [
        ('centromeres', client.fetch_centromeres),
        ('cytoband', client.fetch_cytoband),
        ('cytoband', partial(client.fetch_cytoband, ideo=True)),
        ('gap', client.fetch_gaps)
    ]

    for schema, fetcher in fetchers:
        try:
            df = fetcher()
            break
        except urllib.error.HTTPError:
            pass
    else:
        raise ValueError('No source for centromere data found.')

    return extract_centromeres(df, schema=schema, merge=merge)


class UCSCClient:
    BASE_URL = 'http://hgdownload.cse.ucsc.edu/'

    def __init__(self, db):
        self._db = db
        self._db_url = urljoin(
            self.BASE_URL, 'goldenPath/{}/database/'.format(db)
        )

    def fetch_chromsizes(self, **kwargs):
        url = urljoin(self._db_url, 'chromInfo.txt.gz')
        return read_chromsizes(url, **kwargs)

    def fetch_centromeres(self, **kwargs):
        url = urljoin(self._db_url, 'centromeres.txt.gz')
        return read_table(url, schema='centromeres')

    def fetch_gaps(self, **kwargs):
        url = urljoin(self._db_url, 'gap.txt.gz')
        return read_gapfile(url, **kwargs)

    def fetch_cytoband(self, ideo=False, **kwargs):
        if ideo:
            url = urljoin(self._db_url, 'cytoBandIdeo.txt.gz')
        else:
            url = urljoin(self._db_url, 'cytoBand.txt.gz')
        return read_table(url, schema='cytoband')

    def fetch_mrna(self, **kwargs):
        url = urljoin(self._db_url, 'all_mrna.txt.gz')
        return read_ucsc_mrnafile(url, **kwargs)


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
            pass
            # print('File "{}" available'.format(filename))
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

    def __init__(self, cachedir, assembly, metadata=None, key_id=None, key_secret=None):
        self.cachedir = op.join(cachedir, assembly)
        if not op.isdir(self.cachedir):
            raise OSError("Directory doesn't exist: '{}'".format(cachedir))
        if metadata is None:
            metadata_paths = sorted(glob.glob(op.join(cachedir, 'metadata*.tsv')))
            metadata_path = metadata_paths[-1]
            self._meta = pd.read_table(metadata_path, low_memory=False, comment='#')
            if assembly == 'GRCh38':
                self._meta = self._meta[self._meta['Organism'] == 'human'].copy()
            self._meta = self._meta.set_index('File Accession')
        else:
            self._meta = metadata
        if key_id is not None:
            credential = (key_id + ':' + key_secret).encode('utf-8')
            self._token = base64.b64encode(credential)
        else:
            self._token = None

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
            pass
            # print('File "{}" available'.format(filename))
        else:
            print('Downloading "{}"'.format(filename))
            if self._token:
                headers = {
                    'Authorization': b'Basic ' + self._token
                }
            else:
                headers = None
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            with open(path, 'wb') as f:
                f.write(r.content)
        return path

    def fetch_all(self, accessions):
        return list(map(self.fetch, accessions))
