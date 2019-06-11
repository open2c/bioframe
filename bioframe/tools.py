from __future__ import absolute_import, division, print_function
from ._process import cmd_exists, run, to_dataframe, tsv


if cmd_exists("bedtools"):

    def _register(name):
        # Wrapper calls bedtools
        def wrapper(**kwargs):
            columns = kwargs.pop('_schema', None)
            run_kws = {kw[1:]: kwargs.pop(kw) for kw in list(kwargs.keys())
                           if kw.startswith('_')}

            cmd = ['bedtools', name]
            for k, v in kwargs.items():
                if isinstance(v, bool):
                    if not v: continue
                    cmd.append('-{}'.format(k))
                else:
                    cmd.append('-{}'.format(k))
                    cmd.append(str(v))
            out = run(cmd, **run_kws)
            if not len(out):
                return pd.DataFrame(columns=columns)
            return to_dataframe(out, columns=columns)

        # Call once to generate docstring from usage text
        p = subprocess.Popen(
                ['bedtools', name, '-h'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        _, err = p.communicate()
        wrapper.__doc__ = err.decode('utf-8')
        wrapper.__name__ = str(name)

        return staticmethod(wrapper)


    class bedtools(object):
        intersect = _register('intersect')
        window = _register('window')
        closest = _register('closest')
        coverage = _register('coverage')
        map = _register('map')
        genomecov = _register('genomecov')
        merge = _register('merge')
        cluster = _register('cluster')
        complement = _register('complement')
        subtract = _register('subtract')
        slop = _register('slop')
        flank = _register('flank')
        sort = _register('sort')
        random = _register('random')
        shuffle = _register('shuffle')
        annotate = _register('annotate')
        jaccard = _register('jaccard')


def intersect(bed1, bed2, overlap=True, outer_join=False, v=False, sort=False,
              suffixes=('_x', '_y')):

    # hacky, but we don't want to use suffixes when using -v mode
    if v:
        suffixes = ('',)

    bed1_extra = bed1[bed1.columns.difference(['chrom', 'start', 'end'])]
    bed2_extra = bed2[bed2.columns.difference(['chrom', 'start', 'end'])]

    left = bed1[['chrom', 'start', 'end']].copy()
    left['index'] = left.index

    right = bed2[['chrom', 'start', 'end']].copy()
    right['index'] = right.index

    bt_kwargs = {
        'v': v,
        'nonamecheck': False,
    }

    if outer_join:
        if overlap:
            bt_kwargs['wao'] = True
            bt_kwargs['loj'] = False
        else:
            bt_kwargs['wao'] = False
            bt_kwargs['loj'] = True
    else:
        if overlap:
            bt_kwargs['wo'] = True

    with tsv(left) as a, tsv(right) as b:
        out = bedtools.intersect(a=a.name, b=b.name, **bt_kwargs)

    bed1_extra_out = bed1_extra.iloc[out[3]].reset_index(drop=True)

    if v:
        out_final = pd.concat([out, bed1_extra_out], axis=1)
    else:
        if outer_join:
            out[4] = out[4].where(out[4] != '.')
            out[7] = out[7].where(out[7] != '.', -1).astype(int)

            bed2_extra_out = pd.DataFrame.from_items([
                (name, pd.Series(data=None, index=out.index, dtype=series.dtype))
                for name, series in bed2_extra.iteritems()])
            mask = (out[7] != -1)
            bed2_extra_out.loc[mask, :] = bed2_extra.iloc[out[7][mask]].values
        else:
            bed2_extra_out = bed2_extra.iloc[out[7]].reset_index(drop=True)
        out_final = pd.concat([out, bed1_extra_out, bed2_extra_out], axis=1)


    outcols = [c + suffixes[0] for c in ['chrom', 'start', 'end', 'index']]
    if not v:
        outcols += [c + suffixes[1] for c in ['chrom', 'start', 'end', 'index']]

    if overlap and not v:
        outcols += ['overlap']

    outcols += [c + suffixes[0] for c in bed1_extra_out.columns]
    if not v:
        outcols += [c + suffixes[1] for c in bed2_extra_out.columns]

    out_final.columns = outcols
    return out_final
