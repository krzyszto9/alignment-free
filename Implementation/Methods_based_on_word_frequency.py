#!/usr/bin/python
import argparse
import subprocess
import numpy as np
import math

'''
seq_number - liczba "obrabianych sekwencji"
occurences_list - lista zawierajaca slowniki dla kazdego wzorca, key - numer sekwencji, d[key] - liczba wystapien danego wzorca wsekwencji

chodzi o to, ze w outpucie, szukamy lini zaczynajacych sie od liczby
potem szukamy w jakich sekwencjach wystepuje dany wzorzec i nastepnie dodajemy zwiekszamy liczbe wystapien danego wzorca w slowniku - jeden slownik, jednen motyw 
'''


def parse():
    parser = argparse.ArgumentParser(description = 'Methods based on word (oligomer) frequency', add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', help='input file', required=True, type=argparse.FileType('r'))
    required.add_argument('--lenght', '-l', help='word (oligomer) lenght', required=True, type=int)
    required.add_argument('--method', '-m', help='used method', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-h', action='help', help="show this help message and exit")
    #optional.add_argument('--std', '-s', help="maximum lenght of word (oligomer) length, must be greater than -l/--lenght", type=int)
    try:
        pars = parser.parse_args()
        #if(pars.std and pars.lenght>=pars.std):
            #parser.error( "maximum lenght of word (oligomer) length, must be greater than -l/--lenght")
        return pars
    except IOError, msg:
        parser.error(str(msg))

def teiresias_run(lenght,input_file):
    cmd = 'teiresias -l{0} -w{0} -k1 -i{1} -p -s'.format(lenght,input_file)
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None,shell=True).wait()
    cmd = "grep -c '^>' {0}".format(input_file)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
    return int(p.stdout.readline())

def teiresias_patterns():
    occurences_list =[]
    for line in open('output.txt'):
        if line[0].isdigit():
            line = line.rstrip('\n').split()
            d={}
            for i in range(3, len(line),2):
                if int(line[i]) not in d:
                    d[int(line[i])]=0
                d[int(line[i])]+=1
            occurences_list.append(d)
    subprocess.Popen("rm output.txt", stdout=subprocess.PIPE, stderr=None, shell=True).wait()
    return occurences_list

def identify_headers(input_file):
    list_of_ids=[]
    for line in open(input_file):
        if line.startswith('>'):
            line = line.replace('>','').split()
            list_of_ids.append(line[0])
    return list_of_ids

def euclidean_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
             euclidean_dist=0
             for dicts in occurences_list:
                 key1 = 0 if (i not in dicts) else dicts[i] 
                 key2 = 0 if (j not in dicts) else dicts[j] 
                 euclidean_dist+=(key1-key2)**2
             matrix[i][j] = euclidean_dist
    return matrix

def angle_cos(i,j,occurences_list):
    euclidean = 0
    sum_key1 = 0 
    sum_key2 = 0
    for dicts in occurences_list:
        key1 = 0 if (i not in dicts) else dicts[i] 
        key2 = 0 if (j not in dicts) else dicts[j] 
        euclidean+=(key1*key2)
        sum_key1+=key1**2
        sum_key2+=key2**2
    return euclidean/(math.sqrt(sum_key1) * math.sqrt(sum_key2))

def angle_cos_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
             #value = (1 - self.__angle_cos(seqnum1, seqnum2)) / 2 <-decaf+py
             # w innych zrodlach samo 1- https://en.wikipedia.org/wiki/Cosine_similarity
             matrix[i][j] = round(1 - angle_cos(i,j,occurences_list),4) + 0
    return matrix

def evolutionary_distance(seqs_number, occurences_list):
    matrix = np.zeros([seqs_number, seqs_number])
    for i in range(0,seqs_number):
        for j in range(0,seqs_number):
              matrix[i][j] = round(-math.log((1+angle_cos(i,j,occurences_list))/2),4) + 0 
    return matrix

arguments = parse()
seqs_number = teiresias_run(arguments.lenght, arguments.input.name)
occurences_list = teiresias_patterns()
seqs_identifiers = identify_headers(arguments.input.name)
print "Euclidean distance matrix\n", euclidean_distance(seqs_number, occurences_list)
print "Angle metrics\n", angle_cos_distance(seqs_number, occurences_list)
print "Evolutionary distance\n", evolutionary_distance(seqs_number, occurences_list)
#else:
    #matrix = euclidean(seqs_number, occurences_list)
    #for i in range(arguments.lenght+1, arguments.std+1):
        #teiresias_run(i, arguments.input.name)
        #occurences_list = teiresias_patterns()
        #matrix = matrix + euclidean(seqs_number, occurences_list)
    #print "Macierz dystansu euklidesowego, st\n"
