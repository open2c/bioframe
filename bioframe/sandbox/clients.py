import base64
import glob
import os
import os.path as op
import posixpath as pp
from urllib.parse import urlencode, urljoin

import pandas as pd
import requests


class EncodeClient:
    BASE_URL = "http://www.encodeproject.org/"

    # 2020-05-15 compatible with ENCODE Metadata at:
    METADATA_URL = "https://www.encodeproject.org/metadata/type=Experiment&status=released/metadata.tsv"

    KNOWN_ASSEMBLIES = (
        "GRCh38",
        "GRCh38-minimal",
        "ce10",
        "ce11",
        "dm3",
        "dm6",
        "hg19",
        "mm10",
        "mm10-minimal",
        "mm9",
    )

    def __init__(self, cachedir, assembly, metadata=None):
        if assembly not in self.KNOWN_ASSEMBLIES:
            raise ValueError("assembly must be in:", self.KNOWN_ASSEMBLIES)

        self.cachedir = op.join(cachedir, assembly)
        if not op.isdir(self.cachedir):
            os.makedirs(self.cachedir, exist_ok=True)

        if metadata is None:
            metadata_path = op.join(cachedir, "metadata.tsv")

            if not op.exists(metadata_path):
                print(
                    "getting metadata from ENCODE, please wait while "
                    "(~240Mb) file downloads"
                )
                with requests.get(self.METADATA_URL, stream=True) as r:
                    r.raise_for_status()
                    with open(metadata_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)

            self._meta = pd.read_table(metadata_path, low_memory=False)
            table_assemblies = sorted(
                self._meta["File assembly"].dropna().unique().tolist()
            )

            if not set(table_assemblies).issubset(set(self.KNOWN_ASSEMBLIES)):
                raise ValueError(
                    "Table assemblies do not match known assemblies, "
                    "check ENCODE metadata version"
                )
            self._meta = self._meta[self._meta["File assembly"] == assembly].copy()
            self._meta = self._meta.set_index("File accession")

        else:
            self._meta = metadata

    def _batch_download(self, args):
        params = urlencode(args)
        url = pp.join("batch_download", params)
        url = urljoin(self.BASE_URL, url)
        r = requests.get(url)
        r.raise_for_status()
        return r

    def _metadata(self, args):
        params = urlencode(args)
        url = pp.join("metadata", params, "metadata.tsv")
        url = urljoin(self.BASE_URL, url)
        r = requests.get(url)
        r.raise_for_status()
        return r

    @property
    def meta(self):
        return self._meta.copy()

    def info(self, accession, width=850, height=450):
        from IPython.display import HTML

        url = urljoin(self.BASE_URL, pp.join("experiments", accession))
        return HTML(
            f'<iframe width="{width}px" height="{height}px" src={url}></iframe>'
        )

    def fetch(self, accession):
        url = self.meta.loc[accession, "File download URL"]
        # sig = self.meta.loc[accession, 'md5sum']
        filename = op.split(url)[1]
        path = op.join(self.cachedir, filename)
        if op.exists(path):
            pass
            # print('File "{}" available'.format(filename))
        else:
            print(f'Downloading "{filename}"')
            r = requests.get(url)
            r.raise_for_status()
            with open(path, "wb") as f:
                f.write(r.content)
        return path

    def fetch_all(self, accessions):
        return list(map(self.fetch, accessions))


class FDNClient:
    BASE_URL = "https://data.4dnucleome.org/"

    def __init__(self, cachedir, assembly, metadata=None, key_id=None, key_secret=None):
        self.cachedir = op.join(cachedir, assembly)
        if not op.isdir(self.cachedir):
            raise OSError(f"Directory doesn't exist: '{cachedir}'")
        if metadata is None:
            metadata_paths = sorted(glob.glob(op.join(cachedir, "metadata*.tsv")))
            metadata_path = metadata_paths[-1]
            self._meta = pd.read_table(metadata_path, low_memory=False, comment="#")
            if assembly == "GRCh38":
                self._meta = self._meta[self._meta["Organism"] == "human"].copy()
            self._meta = self._meta.set_index("File Accession")
        else:
            self._meta = metadata
        if key_id is not None:
            credential = (key_id + ":" + key_secret).encode("utf-8")
            self._token = base64.b64encode(credential)
        else:
            self._token = None

    @property
    def meta(self):
        return self._meta.copy()

    def info(self, accession, width=850, height=450):
        from IPython.display import HTML

        url = urljoin(self.BASE_URL, pp.join("experiments", accession))
        return HTML(
            f'<iframe width="{width}px" height="{height}px" src={url}></iframe>'
        )

    def fetch(self, accession):
        url = self.meta.loc[accession, "File Download URL"]
        # sig = self.meta.loc[accession, 'md5sum']
        filename = op.split(url)[1]
        path = op.join(self.cachedir, filename)
        if op.exists(path):
            pass
            # print('File "{}" available'.format(filename))
        else:
            print(f'Downloading "{filename}"')
            if self._token:
                headers = {"Authorization": b"Basic " + self._token}
            else:
                headers = None
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            with open(path, "wb") as f:
                f.write(r.content)
        return path

    def fetch_all(self, accessions):
        return list(map(self.fetch, accessions))
