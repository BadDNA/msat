#!/usr/bin/env python
# encoding: utf-8
"""
msatCommand.py

Created by Brant Faircloth on 2009-03-29.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

#import pdb
import sys
import time
import motif
import optparse
import seqsearch
from Bio import SeqIO

def interface():
    usage = "usage: %prog [options] source_file.fa"
    p = optparse.OptionParser(usage)
    p.add_option('--min-length', '-m', dest='min_length', default = 6, help='The minimum length required for repeats to be considered', metavar='MIN_LENGTH')
    p.add_option('--output-file', '-o', dest='output_file', default= None, help='The name of the output file for results [default=None]', metavar='FILE')
    p.add_option('--scan-type', '-s', default='all', help="The type of scan to run.  Valid options are: \n(1)'all' to search for everything\n '1' | '2' | '3 | '4' | '5' | '6' to search for only mononuc, dinuc, trinuc, etc.\n'2+' | '3+' | '4+' | '5+' | '6+' to search for >= dinuc, >= trinuc, etc. [default=all]")
    p.add_option('--perfect-repeats', '-p', dest='perfect', action='store_true', default=False, help="Locate only perfect repeats.  The default of false only means that when you have a repeat like 'ACACACACACACNNAC', where ambiguous bases are equal to the exact length of the repeat motif (in this case 2), ambiguous bases will be counted as part of the repeat.  This is true IFF the ambiguous bases are followed by >=1 repeat motif, also as in the example [default=False]")
    (options,arg) = p.parse_args()
    if not arg:
        p.print_help()
        sys.exit(2)
    return options, arg

def createMotifInstances(motif, min_length, perfect):
    return seqsearch.MicrosatelliteMotif(motif, min_length, perfect)

def genMotifCollection(options):
    possible_motifs = (motif.mononucleotide, motif.dinucleotide, motif.trinucleotide, motif.tetranucleotide, motif.pentanucleotide, motif.hexanucleotide)
    motif_collection = ()
    if options.scan_type == 'all':
        for m in possible_motifs:
            motif_collection += (createMotifInstances(m, options.min_length, options.perfect),)
    elif '+' in options.scan_type:
        # subtracting 1 so that we get >= options.scan_type
        scan = int(options.scan_type[0]) - 1
        for m in possible_motifs[scan:]:
            motif_collection += (createMotifInstances(m, options.min_length, options.perfect),)
    else:
        # no iteration here because tuple != nested
        scan = int(options.scan_type[0]) - 1
        motif_collection += (createMotifInstances(possible_motifs[scan], options.min_length, options.perfect),)
    return motif_collection

def stdOut(name, matches):
    if matches:
        for msat in matches:
            for match in matches[msat]:
                start = match[0]
                end = match[1]
                length = (end - start) / len(msat)
                print '%s:\t(%s)^%s\trepeat between bases\t%s\tand\t%s' % (name, msat, length, start, end)
    else:
        print '%s:\t No microsatellites found (or repeat region < min_length)' % name
                
def main():
    start_time = time.time()
    options,arg = interface()
    motif_collection = genMotifCollection(options)
    handle = open(arg[0],'rU')
    for record in SeqIO.parse(handle,'fasta'):
        search = seqsearch.RegionSearch(record.seq)
        for m in motif_collection:
            search.microsatellite(m)
        if not options.output_file:
            stdOut(record.id, search.matches)
    end_time = time.time()
    if not options.output_file:
        print '\nTime for execution: ', end_time - start_time, 'seconds'
        

if __name__ == '__main__':
    main()