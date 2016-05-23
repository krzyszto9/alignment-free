#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import itertools
import os.path
import numpy as np
import inspect
from collections import OrderedDict

#np.set_printoptions(threshold=np.nan) #wyswietlanie calych macierzy
np.seterr(divide='ignore',invalid='ignore')
np.set_printoptions(suppress=True)

methods = 	(('euclidean',['eucliudean distance','occurrences_list']),
		('squaredeuclidean',['squared eucliudean distance','occurrences_list']),
		('seuclidean',['standardized eucliudean distance','occurrences_list']),
		('r_seuclidean',['standardized eucliudean distance with different resolutions','occurrences_list']),
		('weuclidean',['weighted eucliudean distance','occurrences_list']),
		('r_weuclidean',['weighted eucliudean distance with different resolutions','occurrences_list']),
		('lcc',['pearson product-moment correlation coefficient - normalized and scaled','frequencies_list',]),
		('kl',['Kullback–Leibler divergence','frequencies_list']),
		('cosine',['cosine distance','occurrences_list']),
		('evol',['evolutionary distance','occurrences_list']),
		('all',['use all above methods','']))

methods = OrderedDict(methods)


def parse():
    args = argparse.ArgumentParser(description = 'Methods based on word (oligomer) frequency.\nAvailable methods: \n' + '\n'.join('{:<20}: {} {} {}'.format(k, v[0],"|",v[1]) for k, v in methods.items()), add_help=False,formatter_class=argparse.RawDescriptionHelpFormatter)
    required = args.add_argument_group('Required arguments')
    required.add_argument('--file', '-f', help='read sequence data from FILE', required=True, type=argparse.FileType('r'))
    required.add_argument('--length', '-l', help='word (oligomer) length', required=True, type=int)
    required.add_argument('--method', '-m', required=True, nargs='+', choices = methods.keys(), help='choose one or more from: ' +"'"+"', '".join(methods.keys()) +"'", metavar='METHOD')

    required_by_methods = args.add_argument_group('Arguments required by some methods')
    required_by_methods.add_argument('--maximum_length', '-ml', help="maximum length of word (oligomer) length, must be greater than -l/--length. required by 'r_seuclidean' and 'r_weuclidean'", type=int)
    required_by_methods.add_argument('--type', '-t', choices = ['p','prot','n','nucl'], help="type of sequences. 'p' or 'prot' if protein sequence, 'n' or 'nucl' if nucleotide. required by 'weuclidean' and 'r_weuclidean",metavar="TYPE")

    optional = args.add_argument_group('Optional arguments')
    optional.add_argument('--calculations', '-c', choices = ['o','occurrences','f','frequencies'], help="type of list on which calculations are based. 'o' or 'occurrences' if occurrences list,'f' or 'frequencies' if frequencies list",metavar='CALCULATIONS')
    optional.add_argument('--help', '-h', action='help', help="show this help message and exit")
    optional.add_argument('--output', '-o', action='store_true', help='write output')
    optional.add_argument('--quiet', '-q', action='store_true', help='do not show any results or messages, use this only if option -o/--output is set')
    try:
        pars = args.parse_args()
        if pars.quiet and not pars.output: args.error('option -q/--quiet can only be used if option -o/--output is set')
        if pars.calculations in ('f','frequencies'):  
            for name, value in methods.items(): value[1] = 'frequencies_list'
        elif pars.calculations in ('o','occurrences'):
            for name, value in methods.items(): value[1] = 'occurrences_list'
        if 'all' in pars.method: 
            del methods['all']
            pars.method = methods
        else:
            d=OrderedDict()
            for name in pars.method: d[name]=methods[name]
            pars.method = d
        if (not pars.type and ('r_weuclidean' in pars.method or 'weuclidean' in pars.method)): args.error( "methods 'weuclidean' and 'r_weuclidean' require option -t/--type to be set")
        if (not pars.maximum_length and ('r_weuclidean' in pars.method or 'r_seuclidean' in pars.method)): args.error( "methods 'r_seuclidean' and 'r_weuclidean' require option -ml/--maximum_length to be set")
        elif (pars.maximum_length and pars.length>=pars.maximum_length): args.error("maximum length of word (oligomer), must be greater than -l/--length")
        return pars
    except IOError, msg:
        args.error(str(msg))

def get_sequences(input_file):
    (seq, sequences_list)=('',[])
    for line in open(input_file):
        if line.startswith('>'):
            sequences_list.append(seq)
            seq = ''
        else:
            seq+= line.rstrip()
    sequences_list.append(seq)
    del sequences_list[0]
    return sequences_list

def get_motifs(length,sequences_list):
    (d,index) = ({}, 0)
    for seq in sequences_list:
        for i in range(0, len(seq)-length+1):
            word = seq[i:i+length]
            if word not in d:
                d[word] = index
                index += 1 
    return d

def get_headers(input_file):
    list_of_ids=[]
    for line in open(input_file):
        if line.startswith('>'):
            line = line.replace('>','').split()
            list_of_ids.append(line[0])
    return list_of_ids

def calculate_occurrences(length, input_file):
    sequences_list = get_sequences(input_file)
    d = get_motifs(length,sequences_list)   
    rows_num = len(d)
    cols_num = len(sequences_list)
    data = np.zeros(shape=(rows_num,cols_num))
    for col_idx, seq in enumerate(sequences_list):
        for i in range(0, len(seq)-length+1):
            word = seq[i:i+length]
            row_idx = d[word]
            data[row_idx, col_idx] += 1
    return data

def calculate_frequencies(occurrences_list,seqs_number):
    frequencies_list =[]
    for i in range(0,seqs_number):
        frequencies_list.append(occurrences_list[:,i]/np.sum(occurrences_list[:,i]))
    return (np.vstack(frequencies_list)).T

def calculate_weights(type_, seqs_number,input_file,length):
    if type_ in ('nucl','n'):
        sequences_list = get_sequences(input_file)
        d = get_motifs(length,sequences_list)
        return np.vstack([0.25 ** length] * len(d))
    elif type_ in ('p','prot'):
        sequences_list = get_sequences(input_file)
        d = get_motifs(length,sequences_list)
        prot_a = {'A':0.082 ,'R':0.055,'N':0.04,'D':0.054,'C':0.013,'Q':0.039,'E':0.067,'G':0.07,'H':0.022,'I':0.059,'L':0.096,'K':0.058,'M':0.024,'F':0.038,
'P':0.047,'S':0.065,'T':0.053,'W':0.01,'Y':0.029,'V':0.068,'Z':0.053,'B':0.047,'X':0.05}
        weight_list = np.zeros([len(d),1])
        for key, value in d.items():
            weight = 1
            for amino_acid in key:
                weight *= prot_a[amino_acid]
            weight_list[value] = weight
        return weight_list

def squaredeuclidean(list_,seqs_number):
    list_ = np.array(list_)
    matrix = np.zeros([seqs_number, seqs_number])
    for i, j in itertools.permutations(range(0,seqs_number),2):
         matrix[i][j]= matrix [j][i] = np.sum((list_[:,i] - list_[:,j])**2)
    return matrix

def euclidean(list_,seqs_number):
    return np.sqrt(squaredeuclidean(list_,seqs_number))

def seuclidean(list_, seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    std = np.nanvar(list_.T,axis=0,ddof=1)
    for i, j in itertools.permutations(range(0,seqs_number),2):
         matrix[i][j]= matrix[j][i] = np.sqrt(np.sum(np.where(std!=0,((list_[:,i] - list_[:,j])**2/std),0)))
    return matrix

def weuclidean(list_,w_list,seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for i, j in itertools.permutations(range(0,seqs_number),2):
         matrix[i][j]= matrix [j][i] = np.sum(((list_[:,i] - list_[:,j])**2) * w_list[:,0])
    return matrix

def r_seuclidean(length,maximum_length,input_file,seqs_number,string):
    matrix = np.zeros([seqs_number, seqs_number])
    for l in range(length,maximum_length+1):
        list_ = calculate_occurrences(l,input_file)
        if string.startswith('f'): list_ = calculate_frequencies(list_, seqs_number)
        matrix+=seuclidean(list_,seqs_number)
    return matrix

def r_weuclidean(length,maximum_length,input_file,seqs_number,type_,string):
    matrix = np.zeros([seqs_number, seqs_number])
    for l in range(length,maximum_length+1):
        list_ = calculate_occurrences(l,input_file)
        if string.startswith('f'): list_ = calculate_frequencies(list_, seqs_number)
        weights_list = calculate_weights(type_, seqs_number,input_file,l)
        matrix+=weuclidean(list_,weights_list,seqs_number)
    return matrix

def lcc(list_, seqs_number):
    matrix = np.ones([seqs_number, seqs_number])
    for i, j in itertools.permutations(range(0,seqs_number),2):
        m = np.vstack((list_[:,i], list_[:,j]))
        matrix[i][j] = matrix[j][i] = (np.corrcoef(m[:, np.apply_along_axis(np.count_nonzero, 0, m) >= 1])[0,1]).round(10) + 0.0
    return (matrix + 1)/2

def kl(list_, seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for i, j in itertools.permutations(range(0,seqs_number),2):
        log2 = np.log2((list_[:,i]+1)/(list_[:,j]+1))
        matrix[i][j] = (np.sum(log2 * (list_[:,i]+1))).round(10) + 0.0
    return matrix

def angle_cos(i,j,list_):
    return np.sum(list_[:,i] * list_[:,j])/(np.sqrt(np.sum(list_[:,i]**2)) * np.sqrt(np.sum(list_[:,j]**2)))

def cosine(list_, seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for i, j in itertools.permutations(range(0,seqs_number),2):
         matrix[i][j] = matrix[j][i] = (1 - angle_cos(i,j,list_)).round(10) + 0.0
    return matrix

def evol(list_, seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for i, j in itertools.permutations(range(0,seqs_number),2):
          matrix[i][j] = matrix[j][i] =(-np.log((1+angle_cos(i,j,list_))/2)).round(10) + 0.0
    return matrix

def main():
    arguments = parse()
    occurrences_list = calculate_occurrences(arguments.length,arguments.file.name)
    seqs_number = occurrences_list.shape[1]
    if 'frequencies_list' in { v2 for v1, v2 in arguments.method.values()}: 
        frequencies_list = calculate_frequencies(occurrences_list, seqs_number)
    if 'weuclidean' in arguments.method: 
        weights_list = calculate_weights(arguments.type, seqs_number,arguments.file.name, arguments.length)
    if arguments.output: ids_list = get_headers(arguments.file.name)
    for name, value in arguments.method.items():
        nargs = len(inspect.getargspec(eval(name))[0])
        if nargs == 2: result = eval(name +'('+value[1]+','+str(seqs_number)+')')
        elif nargs == 3: result = eval(name +'('+value[1]+',weights_list,'+str(seqs_number)+')')
        elif nargs == 5: result = eval(name +'('+str(arguments.length)+','+str(arguments.maximum_length)+',"'+str(arguments.file.name)+'",'+str(seqs_number)+',"'+value[1]+'")')
        else : result = eval(name +'('+str(arguments.length)+','+str(arguments.maximum_length)+',"'+str(arguments.file.name)+'",'+str(seqs_number)+',"'+str(arguments.type)+'","'+value[1]+'")')
        if not arguments.quiet: print '\n' + value[0] + '\n', result 
        if arguments.output:
            np.savetxt(name+'_'+value[1]+'.mat', result,fmt='%.10f', header= " ".join(ids_list))
            if os.path.exists(name+'_'+value[1]+'.mat') and not arguments.quiet: print 'File ' + name +'_'+value[1]+'.mat has been saved successfully'

if __name__ == '__main__':
    main()
