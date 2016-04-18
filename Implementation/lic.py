#!/usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import subprocess
import numpy as np
import math

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
    return int(p.stdout.readline())

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

def euclidean_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(i,seqs_number):
             matrix[i][j]= matrix [j][i] = np.sum((occurences_list[:,i] - occurences_list[:,j])**2)
    return matrix

def angle_cos(i,j,occurences_list):
    return np.sum(occurences_list[:,i] * occurences_list[:,j])/(math.sqrt(np.sum(occurences_list[:,i]**2)) * math.sqrt(np.sum(occurences_list[:,j]**2)))

def angle_cos_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
             matrix[i][j] = math.fabs(1 - angle_cos(i,j,occurences_list))  
    return matrix

def evolutionary_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
              matrix[i][j] = -math.log((1+angle_cos(i,j,occurences_list))/2) + 0
    return matrix

def calculate_frequencies(seqs_number, occurences_list):
    frequences_list =[]
    for rows in occurences_list:
        frequences_list.append(rows/ np.sum(rows))
    return np.vstack(frequences_list)

'''def pearson(seqs_number, frequencies_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
             (x_2, y_2, x, y, xy) = (0,0,0,0,0)
             l=len(frequencies_list)
             for dicts in frequencies_list:
                 value_x = 0 if (i not in dicts) else dicts[i] 
                 value_y = 0 if (j not in dicts) else dicts[j]
                 if value_x == 0 and value_y == 0: l-=1
                 x_2+= value_x**2
                 y_2+= value_y**2
                 x+= value_x
                 y+= value_y
                 xy+= value_x*value_y
             result = math.sqrt((l*x_2-x**2)*(l*y_2-y**2))
             matrix[i][j] = float(0) if result==0 else (l * xy - (x * y))/ math.sqrt((l*x_2-x**2)*(l*y_2-y**2))
    return matrix'''

'''def pearson(seqs_number, frequencies_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
             l=4324
             result = math.sqrt((l*np.sum(frequencies_list[:,i]**2)-np.sum(frequencies_list[:,i])**2)*(l*np.sum(frequencies_list[:,j]**2)-np.sum(frequencies_list[:,j])**2))
             matrix[i][j] = float(0) if result==0 else (l * np.sum(frequencies_list[:,i]*frequencies_list[:,j]) - (np.sum(frequencies_list[:,i]) * np.sum(frequencies_list[:,j])))/ result
    return matrix'''

arguments = parse()
seqs_number = teiresias_run(arguments.lenght, arguments.input.name)
occurences_list = teiresias_patterns(seqs_number)    #można zastanowic sie nad poprawa, tak samo frequencies list
np.set_printoptions(suppress=True)
#Dystans euklidesowy
#print "Euclidean distance matrix\n", euclidean_distance(seqs_number, occurences_list).round(4)
#Cosinus kata pomiedzy sekwencjami
#print "Angle metrics\n", angle_cos_distance(seqs_number, occurences_list).round(4)
#dystans ewolucyjny
#print "Evolutionary distance\n", evolutionary_distance(seqs_number, occurences_list).round(4)

frequencies_list = calculate_frequencies(seqs_number, occurences_list)
#wspolczynnik korelacji liniowej Pearsona
print "Pearson product-moment correlation coefficient\n", pearson(seqs_number, frequencies_list).round(4)
#Dywergencja Kullbacka-Leiblera
#print "Kullback–Leibler divergence\n", kullback_leibler(seqs_number, frequencies_list).round(4)


