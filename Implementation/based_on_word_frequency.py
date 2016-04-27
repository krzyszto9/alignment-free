#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import subprocess
import numpy as np
import math
import os.path

np.set_printoptions(threshold=np.nan) #potem do usuniecia - wyświetlanie całych macierzy
np.seterr(divide='ignore',invalid='ignore')
np.set_printoptions(suppress=True)

def parse():
    parser = argparse.ArgumentParser(description = 'Methods based on word (oligomer) frequency', add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', help='input file', required=True, type=argparse.FileType('r'))
    required.add_argument('--lenght', '-l', help='word (oligomer) lenght', required=True, type=int)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help="show this help message and exit")
    optional.add_argument('--maximum_lenght', '-ml', help="maximum lenght of word (oligomer) length, must be greater than -l/--lenght", type=int)
    try:
        pars = parser.parse_args()
        if(pars.maximum_lenght and pars.lenght>=pars.maximum_lenght):
            parser.error( "maximum lenght of word (oligomer) length, must be greater than -l/--lenght")
        return pars
    except IOError, msg:
        parser.error(str(msg))

def teiresias_run(lenght,input_file):
    cmd = 'teiresias -l{0} -w{0} -k1 -i{1} -p -s'.format(lenght,input_file)
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None,shell=True).wait()
    cmd = "grep -c '^>' {0}".format(input_file)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
    if os.path.exists('output.txt'):
        return int(p.stdout.readline())
    else:
        print "error: There were some errors while running script"
        exit(1)

def teiresias_patterns(seqs_number):
    occurences_list = []
    for line in open('output.txt'):
        if line[0].isdigit():
            line = line.rstrip('\n').split()
            l=[0]*seqs_number
            for i in range(3, len(line),2):
                l[int(line[i])]+=1.0
            occurences_list.append(l)
    return np.vstack(occurences_list)

def calculate_frequencies(seqs_number, occurences_list):
    frequences_list =[]
    for i in range(0,seqs_number):
        frequences_list.append(occurences_list[:,i]/np.sum(occurences_list[:,i]))
    return np.vstack(frequences_list)

def identify_headers(input_file):
    list_of_ids=[]
    for line in open(input_file):
        if line.startswith('>'):
            line = line.replace('>','').split()
            list_of_ids.append(line[0])
    return np.vstack(list_of_ids)

def calculate_weights(seqs_number):
    (seq,weight_list) = ("",[])
    # co z B i Z?
    d = {'A':0.082 ,'R':0.055,'N':0.04,'D':0.054,'C':0.013,'Q':0.039,'E':0.067,'G':0.07,'H':0.022,'I':0.059,'L':0.096,'K':0.058,'M':0.024,'F':0.038,'P':0.047,'S':0.065,'T':0.053,'W':0.01,'Y':0.029,'V':0.068}
    for line in open('output.txt'):
        if line[0].isdigit():
            line = line.rstrip('\n').split()
            weight = 1
            for a in line[2]:
                weight *= d[a]
            weight_list.append(weight)
    return np.vstack(weight_list)

def euclidean_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
             matrix[i][j]= matrix [j][i] = np.sum((occurences_list[:,i] - occurences_list[:,j])**2)
    return matrix

def stdeuclidean_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    std = np.nanvar(occurences_list.T,axis=0,ddof=1)
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
             matrix[i][j]= matrix[j][i] = np.sqrt(np.sum(np.where(std!=0,((occurences_list[:,i] - occurences_list[:,j])**2/std),0)))
    return matrix

def weighted_euclidean_distance(seqs_number, occurences_list,weights_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
             matrix[i][j]= matrix [j][i] = np.sum(((occurences_list[:,i] - occurences_list[:,j])**2) * weights_list[:,0])
    return matrix

def angle_cos(i,j,occurences_list):
    return np.sum(occurences_list[:,i] * occurences_list[:,j])/(np.sqrt(np.sum(occurences_list[:,i]**2)) * np.sqrt(np.sum(occurences_list[:,j]**2)))

def angle_cos_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
             matrix[i][j] = matrix[j][i] = (1 - angle_cos(i,j,occurences_list)).round(10) + 0.0
    return matrix

def evolutionary_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
              matrix[i][j] = matrix[j][i] =(-np.log((1+angle_cos(i,j,occurences_list))/2)).round(10) + 0.0
    return matrix

def kullback_leibler(seqs_number, frequencies_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
            log2 = np.log2((frequencies_list[:,i]+1)/(frequencies_list[:,j]+1))
            matrix[i][j] = (np.sum(log2 * (frequencies_list[:,i]+1))).round(10) + 0.0
    return matrix

def pearson(seqs_number, frequencies_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
             m = np.vstack((frequencies_list[:,i], frequencies_list[:,j]))
             matrix[i][j] = matrix[j][i] = (np.corrcoef(m[:, np.apply_along_axis(np.count_nonzero, 0, m) >= 1])[0,1]).round(10) + 0.0
    return matrix

def stdeuclidean_distance_diffrest_resolutions(lenght,maximum_lenght,input_file,seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for l in range(lenght,maximum_lenght+1):
        teiresias_run(l, input_file)
        occurences_list = teiresias_patterns(seqs_number)
        matrix+=stdeuclidean_distance(seqs_number, occurences_list)
    return matrix

def weighted_euclidean_distance_diffrest_resolutions(lenght,maximum_lenght,input_file,seqs_number):
    matrix = np.zeros([seqs_number, seqs_number])
    for l in range(lenght,maximum_lenght+1):
        teiresias_run(l, input_file)
        occurences_list = teiresias_patterns(seqs_number)
        weights_list = calculate_weights(seqs_number)
        matrix+=weighted_euclidean_distance(seqs_number, occurences_list,weights_list)
    return matrix

def main():
    arguments = parse()
    seqs_number = teiresias_run(arguments.lenght, arguments.input.name)
    occurences_list = teiresias_patterns(seqs_number)    #można zastanowic sie nad poprawa
    weights_list = calculate_weights(seqs_number)        #można zastanowic sie nad poprawa
    frequencies_list = (calculate_frequencies(seqs_number, occurences_list)).T #można zastanowic sie nad poprawa
    seqs_identifiers = identify_headers(arguments.input.name)
    euclidean_dist = euclidean_distance(seqs_number, occurences_list)
    print "\nRaw euclidean distance\n", np.sqrt(euclidean_dist)
    print "\nSquared euclidean distance\n", euclidean_dist
    print "\nStandarized euclidean distance\n", stdeuclidean_distance(seqs_number, occurences_list)
    print "\nWeighted euclidean distance\n", weighted_euclidean_distance(seqs_number, occurences_list,weights_list)
    print "\nAngle metrics\n", angle_cos_distance(seqs_number, occurences_list)         
    print "\nEvolutionary distance\n", evolutionary_distance(seqs_number, occurences_list) 
    print "\nKullback–Leibler divergence\n", kullback_leibler(seqs_number, frequencies_list)
    print "\nPearson product-moment correlation coefficient\n", pearson(seqs_number, frequencies_list)
    print "\nPearson product-moment correlation coefficient - normalized and scaled\n", (pearson(seqs_number, frequencies_list) + 1)/2
    if arguments.maximum_lenght is not None:
        print "\nStandarized euclidean distance from "+str(arguments.lenght)+" to "+str(arguments.maximum_lenght)+"-tuples\n",stdeuclidean_distance_different_resolutions(arguments.lenght,arguments.maximum_lenght,arguments.input.name,seqs_number)
        print "\nWeighted euclidean distance from "+str(arguments.lenght)+" to "+str(arguments.maximum_lenght)+"-tuples\n",weighted_euclidean_distance_different_resolutions(arguments.lenght,arguments.maximum_lenght,arguments.input.name,seqs_number)

if __name__ == '__main__':
    main()
