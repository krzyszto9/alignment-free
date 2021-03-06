#!/usr/bin/python
import argparse
import numpy as np
import zlib
import os.path
import re
import itertools
import math
from collections import OrderedDict

methods = 	(('ncd_zlib','normalized compression distance using zlib'),
		('ncd_lzma','normalized compression distance using lzma'),
                ('usm','universal sequence maps'),
		('all','use all above methods'))
methods = OrderedDict(methods)

def parse():
    parser = argparse.ArgumentParser(description = 'Methods that do not require resolving the sequence with fixed word (oligomer) length segments.\nAvailable methods: \n' + '\n'.join('{:<20}: {} '.format(k, v) for k, v in methods.items())+'\nIf lzma module is not installed, normalized compression distance using lzma method is unavailable.', add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--file', '-f', help='read sequence data from FILE', required=True, type=argparse.FileType('r'))
    required.add_argument('--method', '-m', required=True, nargs='+', choices = methods.keys(), help='choose one or more from: ' +"'"+"', '".join(methods.keys()) +"'", metavar='METHOD')    
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help='show this help message and exit')
    optional.add_argument('--output', '-o', nargs = '?', const = '\n', help='write output')
    optional.add_argument('--pairwise', '-pw', action='store_true', help='write output - pairwise')
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

def calc_usm_d(seqi_f,seqi_b,seqj_f,seqj_b):
    matrix = np.zeros([len(seqi_f),len(seqj_f)])
    limit = pow(10,-308)
    for i in range(0,len(seqi_f)):
        for j in range(0,len(seqj_f)):
            dist_f = -np.log2(max(limit,max(abs(seqi_f[i,:]-seqj_f[j,:]))))
            dist_b = -np.log2(max(limit,max(abs(seqi_b[i,:]-seqj_b[j,:]))))
            matrix[i,j] = dist_f+dist_b
    return matrix.mean()

def usm(seqs_list):
    d = {}
    unique = list(set(''.join(seqs_list)))
    unique.sort()
    number_of_bits=int(math.ceil(np.log2(len(unique)))) 
    for number, char in enumerate(unique): d[char] = number
    mat=[]
    for i,seq in enumerate(seqs_list):
        list_b = [np.binary_repr(numb, width=number_of_bits) for numb in [d[char] for char in seqs_list[i]]]
        matrix_usmc = np.zeros([len(seq)+2,number_of_bits*2])
        #Forward Coordinates
        matrix_usmc[0,:number_of_bits] = matrix_usmc[len(seq)+1,number_of_bits:]= np.random.rand(number_of_bits)
        for j in range(1,len(seq)+1):
            binary = np.array([int(x) for x in str(list_b[j-1])])
            matrix_usmc[j,:number_of_bits] = matrix_usmc[j-1,:number_of_bits]+(0.5*(1-matrix_usmc[j-1,:number_of_bits]))*binary -(0.5*matrix_usmc[j-1,:number_of_bits])*(1-binary)
        #Backward Coordinates
        for j in reversed(range(1,len(seq)+1)):
            binary = np.array([int(x) for x in str(list_b[j-1])])
            matrix_usmc[j,number_of_bits:] = matrix_usmc[j+1,number_of_bits:]+(0.5*(1-matrix_usmc[j+1,number_of_bits:]))*binary - (0.5*matrix_usmc[j+1,number_of_bits:])*(1-binary)
        matrix_usmc = matrix_usmc[1:-1,:]
        mat.append(matrix_usmc)
    mat= np.array(mat)
    matrix = np.zeros([len(seqs_list), len(seqs_list)])
    for i, j in itertools.combinations(range(0,len(seqs_list)),2):
          matrix[i][j] = matrix[j][i] = calc_usm_d(mat[i][:,:number_of_bits],mat[i][:,number_of_bits:],mat[j][:,:number_of_bits],mat[j][:,number_of_bits:])
    return matrix
          

def write_to_file(filename,result,ids_list, quiet='',pairwise=''):
    if pairwise:
        result= np.asmatrix('\n'.join('{} {} {}'.format(i+1,j+1,result[i,j]) for i, j in itertools.product(range(0,len(ids_list)),repeat=2)))
        result = result.reshape(result.shape[1]/3,3)
        filename = '{}_pw.mat'.format(filename)
        np.savetxt(filename,result,fmt='%.f %.f\t%.10f', header= " ".join(ids_list)) 
    else: 
        filename = '{}.mat'.format(filename)
        np.savetxt(filename,result,fmt='%.10f', header= " ".join(ids_list)) 
    if os.path.exists(filename) and not quiet: print 'File',filename,'has been saved successfully' 

def main():
    arguments = parse()
    sequences_list, ids_list = get_seqs_and_ids(arguments.file.name)
    for name in arguments.method:
        result = eval(name +'(sequences_list)')
        if not arguments.quiet: print '\n' + methods[name] + '\n', result 
        if arguments.output:
            name_f = re.sub("[.][a-zA-Z]*","",arguments.file.name)
            if arguments.output != '\n': name_f = arguments.output
            filename  = '{}_{}'.format(name_f,name)
            write_to_file(filename,result,ids_list,arguments.quiet,arguments.pairwise)

if __name__ == '__main__':
    try:
        import lzma
    except ImportError:
        del methods['ncd_lzma']
    main()
