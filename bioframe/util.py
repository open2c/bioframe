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
    fh = tempfile.NamedTemporaryFile(mode="w+t")
    df.to_csv(fh, sep=str("\t"), index=False, header=False, na_rep="nan", **kwargs)
    fh.flush()  # DON'T FORGET TO FLUSH!!!
    fh.seek(0)
    return fh


def run(cmd, input=None, raises=True, print_cmd=False, max_msg_len=1000):
    if print_cmd:
        print(subprocess.list2cmdline(cmd))

    if input is not None:
        p = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = p.communicate(input.encode("utf-8"))
    else:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

    if raises and p.returncode != 0:
        if len(out) > max_msg_len:
            out = out[:max_msg_len] + b"... [truncated]"
        raise OSError(
            "process failed: %d\n%s\n%s"
            % (p.returncode, out.decode("utf-8"), err.decode("utf-8"))
        )

    return out.decode("utf-8")


def cmd_exists(cmd):
    return any(
        os.access(os.path.join(path, cmd), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )


def to_dataframe(text, columns=None):
    # To convert decoded stdout into a dataframe
    return pd.read_csv(StringIO(text), sep="\t", header=None, names=columns)


if cmd_exists("bedtools"):

    def _register(name):
        # Wrapper calls bedtools
        def wrapper(**kwargs):
            columns = kwargs.pop("_schema", None)

            run_kws = {}
            pandas_inputs = {}
            for arg in list(kwargs.keys()):
                if arg.startswith("_"):
                    run_kws[arg[1:]] = kwargs.pop(arg)
                elif isinstance(kwargs[arg], pd.DataFrame):
                    tmp_file = tsv(kwargs[arg])
                    pandas_inputs[arg] = tmp_file
                    kwargs[arg] = tmp_file.name

            cmd = ["bedtools", name]
            for k, v in kwargs.items():
                if isinstance(v, bool):
                    if not v:
                        continue
                    cmd.append("-{}".format(k))
                else:
                    cmd.append("-{}".format(k))
                    cmd.append(str(v))

            try:
                out = run(cmd, **run_kws)
            finally:
                for tmp_file in pandas_inputs.values():
                    tmp_file.close()

            if not len(out):
                return pd.DataFrame(columns=columns)
            return to_dataframe(out, columns=columns)

        # Call once to generate docstring from usage text
        p = subprocess.Popen(
            ["bedtools", name, "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        _, err = p.communicate()
        wrapper.__doc__ = err.decode("utf-8")
        wrapper.__name__ = str(name)

        return staticmethod(wrapper)

else:

    def _register(name):
        def wrapper(**kwargs):
            raise RuntimeError(
                "`bedtools` command not found. To use the bedtools wrapper, "
                "install the bedtools command suite and reload bioframe"
            )

        return staticmethod(wrapper)


class bedtools(object):
    intersect = _register("intersect")
    window = _register("window")
    closest = _register("closest")
    coverage = _register("coverage")
    map = _register("map")
    genomecov = _register("genomecov")
    merge = _register("merge")
    cluster = _register("cluster")
    complement = _register("complement")
    subtract = _register("subtract")
    slop = _register("slop")
    flank = _register("flank")
    sort = _register("sort")
    random = _register("random")
    shuffle = _register("shuffle")
    annotate = _register("annotate")
    jaccard = _register("jaccard")
