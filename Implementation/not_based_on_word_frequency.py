#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import zlib
import lzma

def parse():
    parser = argparse.ArgumentParser(description = 'Methods that do not require resolving the sequence with fixed word (oligomer) length segments', add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', help='input file', required=True, type=argparse.FileType('r'))
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help="show this help message and exit")
    try:
        return parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))

def get_sequences_and_ids(input_file):
    (seq, sequences_list, ids_list)=("",[],[])
    for line in open(input_file):
        if line.startswith('>'):
            if seq!="": sequences_list.append(seq)
            line = line.replace('>','').split()
            ids_list.append(line[0])
            seq = ""
        else:
            seq+= line.rstrip()
    sequences_list.append(seq)
    return (np.vstack(sequences_list),np.vstack(ids_list))

def normalized_compression_distance_zlib(sequences_list,seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
            (x,y) = (float(len(zlib.compress(sequences_list[i],9))), float(len(zlib.compress(sequences_list[j],9))))
            matrix[i][j] =  (float(len(zlib.compress(np.core.defchararray.add(sequences_list[i],sequences_list[j]),9))) - min(x,y))/max(x,y)
    return matrix

def normalized_compression_distance_lzma(sequences_list,seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
            (x,y) = (float(len(lzma.compress(sequences_list[i]))), float(len(lzma.compress(sequences_list[j]))))
            matrix[i][j] =  (float(len(lzma.compress(np.core.defchararray.add(sequences_list[i],sequences_list[j])))) - min(x,y))/max(x,y)
    return matrix
    
def main():
    arguments = parse()
    sequences_list, ids_list = get_sequences_and_ids(arguments.input.name)
    seqs_number = len(ids_list)
    print "\nNormalized Compression Distance - zlib\n", normalized_compression_distance_zlib(sequences_list,seqs_number)
    #LZMA - wyczytałem, że lepiej kompresuje dane, dla wielkich plików strasznie długo trzeba czekać
    print "\nNormalized Compression Distance - lzma\n", normalized_compression_distance_lzma(sequences_list,seqs_number)

if __name__ == '__main__':
    main()
