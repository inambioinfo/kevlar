#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import re
import shutil
import tempfile
import kevlar


@pytest.mark.long
def test_trio2(capsys):
    from sys import stdout, stderr
    tempdir = tempfile.mkdtemp()
    findoutfiles = ['{}/out{}'.format(tempdir, i) for i in range(4)]
    findouts = [open(fn, 'w') for fn in findoutfiles]
    for i in range(4):
        args = type('', (), {})()
        args.controls = glob.glob('tests/data/trio2/ctrl[1,2].fq.gz')
        args.ctrl_max = 1
        args.case_min = 8
        args.ksize = 31
        args.memory = 2e5
        args.max_fpr = 0.2
        args.out = findouts[i]
        args.flush = False
        args.collapse = False
        args.batch = [4, i+1]
        args.upint = 10000
        args.logfile = stderr
        args.cases = ['tests/data/trio2/case1.fq.gz']

        kevlar.find.main(args)
        findouts[i].close()

    args = type('', (), {})()
    args.memory = 5e3
    args.ksize = 31
    args.out = stdout
    args.max_fpr = 0.02
    args.ignore = False
    args.debug = False
    args.minabund = 8
    args.collapse = True
    args.find_output = findoutfiles
    args.logfile = stderr
    kevlar.collect.main(args)

    out, err = capsys.readouterr()
    assert '1 collapsed linear paths' in err
    lastline = out.strip().split('\n')[-1]
    contig, contigrc = lastline.split(',')[0:2]
    assert 'AGCCTCTG' in contig or 'AGCCTCTG' in contigrc

    shutil.rmtree(tempdir)
