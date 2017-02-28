#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import kevlar


@pytest.fixture
def bogusargs(capsys):
    args = type('', (), {})()
    args.seqid = None
    args.genomemask = None
    args.refr = 'tests/data/bogus-genome/refr.fa'
    args.reads = 'tests/data/bogus-genome/reads.bam'
    return args


def test_basic(bogusargs, capsys):
    from sys import stdout, stderr
    bogusargs.out = stdout
    bogusargs.logfile = stderr
    kevlar.dump.main(bogusargs)
    out, err = capsys.readouterr()
    outputlines = out.strip().split('\n')
    assert len(outputlines) == 5 * 4  # 5 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read4' in outputlines[4]
    assert 'read6' in outputlines[8]
    assert 'read7' in outputlines[12]
    assert 'read8' in outputlines[16]


def test_seqid_filter(bogusargs, capsys):
    from sys import stdout, stderr
    bogusargs.out = stdout
    bogusargs.logfile = stderr
    bogusargs.seqid = 'bogus-genome-chr1'
    kevlar.dump.main(bogusargs)
    out, err = capsys.readouterr()
    outputlines = out.strip().split('\n')
    assert len(outputlines) == 1 * 4  # 1 record, 4 lines per record
    assert 'read2' in outputlines[0]


def test_genomemask_filter(bogusargs, capsys):
    from sys import stdout, stderr
    bogusargs.out = stdout
    bogusargs.logfile = stderr
    bogusargs.genomemask = 'tests/data/bogus-genome/mask-chr1.fa'
    bogusargs.maskmemory = 5e7
    bogusargs.mask_k = 13
    kevlar.dump.main(bogusargs)
    out, err = capsys.readouterr()
    outputlines = out.strip().split('\n')
    assert len(outputlines) == 3 * 4  # 3 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read7' in outputlines[4]
    assert 'read8' in outputlines[8]


def test_seqid_genomemask_filters(bogusargs, capsys):
    from sys import stdout, stderr
    bogusargs.out = stdout
    bogusargs.logfile = stderr
    bogusargs.seqid = 'bogus-genome-chr1'
    bogusargs.genomemask = 'tests/data/bogus-genome/mask-chr1.fa'
    bogusargs.maskmemory = 5e7
    bogusargs.mask_k = 13
    kevlar.dump.main(bogusargs)
    out, err = capsys.readouterr()
    outputlines = out.strip().split('\n')
    assert len(outputlines) == 3 * 4  # 3 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read7' in outputlines[4]
    assert 'read8' in outputlines[8]


def test_indels(bogusargs, capsys):
    from sys import stdout, stderr
    bogusargs.out = stdout
    bogusargs.logfile = stderr
    bogusargs.genomemask = 'tests/data/bogus-genome/mask-chr2.fa'
    bogusargs.maskmemory = 5e7
    bogusargs.mask_k = 13
    bogusargs.reads = 'tests/data/bogus-genome/reads-indels.bam'
    kevlar.dump.main(bogusargs)
    out, err = capsys.readouterr()
    outputlines = out.strip().split('\n')
    assert len(outputlines) == 2 * 4  # 2 records, 4 lines per record
    assert 'read2' in outputlines[0]
    assert 'read3' in outputlines[4]


def test_suffix(bogusargs, capsys):
    from sys import stdout, stderr
    bogusargs.out = stdout
    bogusargs.logfile = stderr
    bogusargs.reads = 'tests/data/nopair.sam'
    kevlar.dump.main(bogusargs)
    out, err = capsys.readouterr()
    outputlines = out.strip().split('\n')
    assert len(outputlines) == 4
    assert outputlines[0].endswith('/1') or outputlines[0].endswith('/2')
