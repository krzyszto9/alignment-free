#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import zlib
import os.path

methods = {'ncd_zlib':'normalized compression distance using zlib','ncd_lzma': 'normalized compression distance using lzma','all': 'use all above methods'}

def parse():
    parser = argparse.ArgumentParser(description = 'Methods that do not require resolving the sequence with fixed word (oligomer) length segments', add_help=False,formatter_class=argparse.RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', help='input file', required=True, type=argparse.FileType('r'))
    required.add_argument('--methods', '-m', help='choose one or more from: ' +"'"+"', '".join(methods.keys()) +"'.\n" + '\n'.join('{}: {}'.format(k, v) for k, v in methods.items()) +'\nIf lzma module is not installed, normalized compression distance using lzma method is unavailable.', required=True, nargs='+')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help='show this help message and exit')
    optional.add_argument('--output', '-o', action='store_true', help='write output to proper file(s)')
    try:
        pars = parser.parse_args()
        pars.methods = map(str.lower,pars.methods)
        if 'all' in pars.methods: 
            del methods['all']
            pars.methods = methods.keys()
        elif len([i for i in pars.methods if i in methods.keys()]) != len(pars.methods):
            parser.error( 'unrecignized method name')
        return pars
    except IOError, msg:
        parser.error(str(msg))

def get_seqs_and_ids(input_file):
    (seq, sequences_list, ids_list)=("",[],[])
    for line in open(input_file):
        if line.startswith('>'):
            if seq!='': sequences_list.append(seq)
            line = line.replace('>','').split()
            ids_list.append(line[0])
            seq = ''
        else:
            seq+= line.rstrip()
    sequences_list.append(seq)
    return sequences_list, ids_list

def ncd_zlib(seqs_list):
    seqs_number = len(seqs_list)
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
            (x,y) = (float(len(zlib.compress(seqs_list[i],9))), float(len(zlib.compress(seqs_list[j],9))))
            matrix[i][j] = (float(len(zlib.compress(seqs_list[i]+ seqs_list[j],9))) - min(x,y))/max(x,y)
    return matrix

def ncd_lzma(seqs_list):
    seqs_number = len(seqs_list)
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
            (x,y) = (float(len(lzma.compress(seqs_list[i]))), float(len(lzma.compress(seqs_list[j]))))
            matrix[i][j] =  (float(len(lzma.compress(seqs_list[i]+seqs_list[j]))) - min(x,y))/max(x,y)
    return matrix

def main():
    arguments = parse()
    sequences_list, ids_list = get_seqs_and_ids(arguments.input.name)
    arguments.methods = list(set(arguments.methods))
    for name in arguments.methods:
        result = eval(name +'(sequences_list)')
        print '\n' + methods[name] + '\n', result 
        if arguments.output:
            np.savetxt(name+'.mat', result,fmt='%.10f', header= " ".join(ids_list))
            if os.path.exists(name+'.mat'): print 'File ' + name +'.mat has been saved successfully\n'

if __name__ == '__main__':
    try:
        import lzma
    except ImportError:
        del methods['ncd_lzma']
    main()
