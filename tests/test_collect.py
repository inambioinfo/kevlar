#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import khmer
import kevlar


def test_load_all_inputs():
    cg = khmer.Countgraph(13, 2e6 / 4, 4)
    vs = kevlar.VariantSet()
    minabund = 8
    maxfpr = 0.2

    batches = [3, 4, 5, 7, 8]
    filepattern = 'tests/data/trio1/novel_1_1,2_batch{:d}.txt'
    infilenames = [filepattern.format(batch) for batch in batches]

    kevlar.collect.load_all_inputs(infilenames, cg, vs, minabund, maxfpr)
    assert vs.nkmers == 6
    assert vs.nkmerinst == 49
    assert vs.nreads == 10


def test_load_all_inputs_alpha():
    cg = khmer.Countgraph(19, 4e4 / 4, 4)
    vs = kevlar.VariantSet()
    minabund = 8
    maxfpr = 0.2
    infilenames = ['tests/data/collect.alpha.txt']

    kevlar.collect.load_all_inputs(infilenames, cg, vs, minabund, maxfpr)
    assert 'CAGGCCAGGGATCGCCGTG' not in vs.kmers and \
           kevlar.revcom('CAGGCCAGGGATCGCCGTG') not in vs.kmers

    kmers = ['TAGGGGCGTGACTTAATAA', 'AGGGGCGTGACTTAATAAG',
             'GGGGCGTGACTTAATAAGG', 'GGGCGTGACTTAATAAGGT']
    for kmer in kmers:
        assert kmer in vs.kmers or kevlar.revcom(kmer) in vs.kmers


def test_load_all_inputs_beta():
    cg = khmer.Countgraph(19, 4e4 / 4, 4)
    vs = kevlar.VariantSet()
    minabund = 8
    maxfpr = 0.2
    infilenames = glob.glob('tests/data/collect.beta.?.txt')

    kevlar.collect.load_all_inputs(infilenames, cg, vs, minabund, maxfpr)
    assert vs.nreads == 8

    readseq = 'TTAACTCTAGATTAGGGGCGTGACTTAATAAGGTGTGGGCCTAAGCGTCT'
    for kmer in cg.get_kmers(readseq):
        assert cg.get(kmer) == 8
