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
from khmer import khmer_args
import kevlar


def load_input(filelist, ksize, memory, maxfpr, logfile=sys.stderr):
    variants = kevlar.VariantSet()
    countgraph = khmer.Countgraph(ksize, memory / 4, 4)
    nr, nk = 0, 0
    for filename in filelist:
        with kevlar.open(filename, 'r') as infile:
            for record in kevlar.parse_augmented_fastx(infile):
                k = countgraph.consume(record.sequence)
                nk += k
                nr += 1
                for ikmer in record.ikmers:
                    variants.add_kmer(ikmer.sequence, record.name)

    fpr = kevlar.sketch.estimate_fpr(countgraph)
    message = '    {:d} reads'.format(nr)
    message += ' and {:d} k-mers consumed'.format(nk)
    message += '; found {:d} instances'.format(variants.nkmerinst)
    message += ' of {:d} unique novel k-mers'.format(variants.nkmers)
    message += '; estimated false positive rate is {:1.3f}'.format(fpr)
    print(message, file=logfile)
    if fpr > maxfpr:
        print('[kevlar::collect] FPR too high, bailing out', fpr, file=logfile)
        sys.exit(1)  # FIXME exception
    return variants, countgraph


def assemble_contigs(countgraph, variants, kmers_to_ignore=None,
                     collapse=True, debug=False, logfile=sys.stderr):
    asm = khmer.JunctionCountAssembler(countgraph)
    print('[kevlar::collect] Assembling contigs', file=logfile)
    for kmer in variants.kmers:
        if kmers_to_ignore and kmer in kmers_to_ignore:
            continue
        if debug:
            print('[kevlar::collect]     DEBUG kmer:', kmer, file=logfile)
        contigs = asm.assemble(kmer)
        for i, contig in enumerate(contigs):
            if contig == '':
                print('    WARNING: no contig found for k-mer', kmer,
                      file=logfile)
                continue
            variants.add_contig(contig, kmer)
    print('    {:d} contigs'.format(variants.ncontigs), file=logfile)

    if collapse:
        print('[kevlar::collect] Collapsing contigs', file=logfile)
        variants.collapse()
        print('    {:d} collapsed contigs'.format(variants.ncontigs),
              file=logfile)


def main(args):
    variants, countgraph = load_input(
        args.novel_output, args.ksize, args.memory, args.max_fpr, args.logfile
    )
    assemble_contigs(countgraph, variants, args.ignore, args.collapse,
                     args.debug, args.logfile)
    variants.write_fasta(outstream=args.out)
