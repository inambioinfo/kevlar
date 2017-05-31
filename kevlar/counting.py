#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/dib-lab/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re
import sys

import khmer
import kevlar


class KevlarSampleIOError(ValueError):
    pass


def load_sample_seqfile(seqfiles, ksize, memory, maxfpr=0.2,
                        masks=None, maskmaxabund=1, numbands=None, band=None,
                        outfile=None, logfile=sys.stderr):
    """
    asdf
    """
    message = 'loading sample from ' + ','.join(seqfiles)
    print('[kevlar::counting]    ', message, file=logfile)

    sketch = khmer.Counttable(ksize, memory / 4, 4)
    n, nkmers = 0, 0
    for n, read in enumerate(kevlar.multi_file_iter_khmer(seqfiles), 1):
        for subseq in kevlar.clean_subseqs(read.sequence, ksize):
            for kmer in sketch.get_kmers(subseq):
                if numbands:
                    khash = sketch.hash(kmer)
                    if khash & (numbands - 1) != band - 1:
                        continue
                if masks:
                    for mask in masks:
                        if mask.get(kmer) > maskmaxabund:
                            break
                    else:
                        sketch.add(kmer)
                        nkmers += 1
                else:
                    sketch.add(kmer)
                    nkmers += 1

    message = 'done loading reads'
    if numbands:
        message += ' (band {:d}/{:d})'.format(band, numbands)
    fpr = kevlar.calc_fpr(sketch)
    message += '; {:d} reads processed'.format(n)
    message += ', {:d} k-mers stored'.format(nkmers)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    if fpr > maxfpr:
        message += ' (FPR too high, bailing out!!!)'
        raise SystemExit(message)
    else:
        if outfile:
            if not outfile.endswith(('.ct', '.counttable')):
                outfile += '.counttable'
            sketch.save(outfile)
            message += '; saved to "{:s}"'.format(outfile)
        print('[kevlar::counting]    ', message, file=logfile)

    return sketch


def load_samples_with_dilution(samplelists, ksize, memory, memfraction=0.1,
                               maxfpr=0.2, maxabund=1, masks=None,
                               numbands=None, band=None, skipsave=False,
                               logfile=sys.stderr):
    """
    asdf
    """
    numsamples = len(samplelists)
    message = 'computing k-mer abundances for {:d} samples'.format(numsamples)
    print('[kevlar::counting]    ', message, file=logfile)

    sketches = list()
    for samplelist in samplelists:
        if len(samplelist) < 2:
            message = 'must specify an output file and at least one input file'
            raise KevlarSampleIOError(message)
        outfile = samplelist[0]
        seqfiles = samplelist[1:]
        if masks:
            mymasks = masks
            sketchmem = memory * memfraction
        elif len(sketches) == 0:
            mymasks = None
            sketchmem = memory
        else:
            mymasks = sketches
            sketchmem = memory * memfraction
        sketch = load_sample_seqfile(
            seqfiles, ksize, sketchmem, maxfpr=maxfpr, masks=mymasks,
            maskmaxabund=maxabund, numbands=numbands, band=band,
            outfile=outfile, logfile=logfile
        )
        sketches.append(sketch)
    return sketches


def load_samples_sketchfiles(sketchfiles, maxfpr=0.2, logfile=sys.stderr):
    """
    asdf
    """
    sketches = list()
    for sketchfile in sketchfiles:
        message = 'loading sketchfile ' + sketchfile,
        print('[kevlar::counting]    ', message, end='', file=logfile)
        sketch = kevlar.sketch_autoload(sketchfile)
        fpr = kevlar.calc_fpr(sketch)
        message = 'done! estimated false positive rate is {:1.3f}'.format(fpr)
        if fpr > maxfpr:
            message += ' (FPR too high, bailing out!!!)'
            raise SystemExit(message)
        else:
            print(message, file=logfile)
        sketches.append(sketch)
    return sketches
