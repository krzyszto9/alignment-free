#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import subprocess
import numpy as np
import math
import os.path

np.set_printoptions(threshold=np.nan) #potem do usuniecia - wyświetlanie całych macierzy
np.seterr(divide='ignore',invalid='ignore')

def parse():
    parser = argparse.ArgumentParser(description = 'Methods based on word (oligomer) frequency', add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', help='input file', required=True, type=argparse.FileType('r'))
    required.add_argument('--lenght', '-l', help='word (oligomer) lenght', required=True, type=int)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help="show this help message and exit")
    try:
        return parser.parse_args()
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
    subprocess.Popen("rm output.txt", stdout=subprocess.PIPE, stderr=None, shell=True).wait()
    return np.vstack(occurences_list)

def calculate_frequencies(seqs_number, occurences_list):
    frequences_list =[]
    for i in range(0,seqs_number):
        frequences_list.append(occurences_list[:,i]/np.sum(occurences_list[:,i]))
    return np.vstack(frequences_list)

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
             matrix[i][j]= matrix[j][i] = np.sum(np.where(std!=0,((occurences_list[:,i] - occurences_list[:,j])**2/std),0))
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

arguments = parse()
seqs_number = teiresias_run(arguments.lenght, arguments.input.name)
occurences_list = teiresias_patterns(seqs_number)    #można zastanowic sie nad poprawa
frequencies_list = (calculate_frequencies(seqs_number, occurences_list)).T #można zastanowic sie nad poprawa
#Dystans euklidesowy    			DZIAŁA DOBRZE - wszystkie trzy
euclidean_dist = euclidean_distance(seqs_number, occurences_list)
print "Raw euclidean distance matrix\n", np.sqrt(euclidean_dist)
print "Squared euclidean distance matrix\n", euclidean_dist
print "Standarized euclidean distance matrix\n", np.sqrt(stdeuclidean_distance(seqs_number, occurences_list))
#Cosinus kata pomiedzy sekwencjami 		DZIAŁA DOBRZE
print "Angle metrics\n", angle_cos_distance(seqs_number, occurences_list)         
#Dystans ewolucyjny 				DZIAŁA DOBRZE
print "Evolutionary distance\n", evolutionary_distance(seqs_number, occurences_list) 
#Dywergencja Kullbacka-Leiblera 		DZIAŁA DOBRZE
print "Kullback–Leibler divergence\n", kullback_leibler(seqs_number, frequencies_list)
#Wspolczynnik korelacji liniowej Pearsona
print "Pearson product-moment correlation coefficient\n", pearson(seqs_number, frequencies_list)
