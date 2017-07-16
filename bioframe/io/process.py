from __future__ import division, print_function
from io import StringIO
import subprocess
import tempfile
import os

import pandas as pd


def tsv(df, **kwargs):
    """
    Write ``pandas.DataFrame`` to a temporary tab-delimited file.
    Works in a ``with`` block (file is deleted at context teardown).

    >>> with tsv(df1) as f1, tsv(df2) as f2:
    ...    # something that requires tsv file input (use f or f.name)

    """
    fh = tempfile.NamedTemporaryFile(mode='w+t')
    df.to_csv(fh, sep=str('\t'), index=False, header=False, **kwargs)
    fh.flush()  # DON'T FORGET TO FLUSH!!!
    fh.seek(0)
    return fh


def run(cmd, input=None, raises=True, print_cmd=False, max_msg_len=1000):
    if print_cmd:
        print(subprocess.list2cmdline(cmd))

    if input is not None:
        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate(input.encode('utf-8'))
    else:
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()


    if raises and p.returncode != 0:
        if len(out) > max_msg_len:
            out = out[:max_msg_len] + b'... [truncated]'
        raise OSError("process failed: %d\n%s\n%s" % (p.returncode,  out.decode('utf-8'), err.decode('utf-8')))

    return out.decode('utf-8')


def cmd_exists(cmd):
    return any(os.access(os.path.join(path, cmd), os.X_OK)
               for path in os.environ["PATH"].split(os.pathsep))

def to_dataframe(text, columns=None):
    # To convert decoded stdout into a dataframe
    return pd.read_csv(StringIO(text), sep='\t', header=None, names=columns)
