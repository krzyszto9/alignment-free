#!/usr/bin/python
import argparse
import numpy as np
import zlib
import os.path
from collections import OrderedDict

methods = 	(('ncd_zlib','normalized compression distance using zlib'),
		('ncd_lzma','normalized compression distance using lzma'),
		('all','use all above methods'))
methods = OrderedDict(methods)

def parse():
    parser = argparse.ArgumentParser(description = 'Methods that do not require resolving the sequence with fixed word (oligomer) length segments.\nAvailable methods: \n' + '\n'.join('{:<20}: {} '.format(k, v) for k, v in methods.items())+'\nIf lzma module is not installed, normalized compression distance using lzma method is unavailable.', add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--file', '-f', help='read sequence data from FILE', required=True, type=argparse.FileType('r'))
    required.add_argument('--method', '-m', required=True, nargs='+', choices = methods.keys(), help='choose one or more from: ' +"'"+"', '".join(methods.keys()) +"'", metavar='METHOD')    
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help='show this help message and exit')
    optional.add_argument('--output', '-o', action='store_true', help='write output to proper file(s)')
    optional.add_argument('--quiet', '-q', action='store_true', help='do not show any results or messages, use this only if option -o/--output is set')
    try:
        pars = parser.parse_args()
        if 'all' in pars.method: 
            del methods['all']
            pars.method = methods.keys()
        if pars.quiet and not pars.output:  parser.error( 'option -q/--quiet can only be used if option -o/--output is set')
        return pars
    except IOError, msg:
        parser.error(str(msg))

def get_seqs_and_ids(input_file):
    (seq, sequences_list, ids_list)=('',[],[])
    for line in open(input_file):
        if line.startswith('>'):
            sequences_list.append(seq)
            line = line.replace('>','').split()
            ids_list.append(line[0])
            seq = ''
        else:
            seq+= line.rstrip()
    sequences_list.append(seq)
    del sequences_list[0]
    return sequences_list, ids_list

def get_sequences(input_file):
    return get_seqs_and_ids(input_file)[0]

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
    sequences_list, ids_list = get_seqs_and_ids(arguments.file.name)
    for name in arguments.method:
        result = eval(name +'(sequences_list)')
        if not arguments.quiet: print '\n' + methods[name] + '\n', result 
        if arguments.output:
            np.savetxt(name+'.mat', result,fmt='%.10f', header= " ".join(ids_list))
            if os.path.exists(name+'.mat') and not arguments.quiet: print 'File ' + name +'.mat has been saved successfully'

if __name__ == '__main__':
    try:
        import lzma
    except ImportError:
        del methods['ncd_lzma']
    main()
