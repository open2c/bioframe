from __future__ import division, print_function, unicode_literals

import os
import subprocess

from .io import cmd_exists, run

if cmd_exists("bedtools"):

    def _register(name):
        # Wrapper calls bedtools
        def wrapper(**kwargs):
            columns = kwargs.pop('_schema', None)
            run_kws = {kw[1:]: kwargs.pop(kw) for kw in list(kwargs.keys())
                           if kw.startswith('_')}

            cmd = [os.path.join(settings['bedtools_path'], 'bedtools'), name]
            for k, v in kwargs.items():
                if isinstance(v, bool):
                    if not v: continue
                    cmd.append('-{}'.format(k))
                else:
                    cmd.append('-{}'.format(k))
                    cmd.append(str(v))
            out = run(cmd, **run_kws)
            return to_dataframe(out, columns=columns)

        # Call once to generate docstring from usage text
        p = subprocess.Popen(
                [os.path.join(settings['bedtools_path'],'bedtools'), name, '-h'],
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

