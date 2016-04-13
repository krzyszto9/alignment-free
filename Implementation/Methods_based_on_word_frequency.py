#!/usr/bin/python
import argparse
import subprocess
import time

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
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-help', action='help', help="show this help message and exit")
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

def euclidean(seqs_number, occurences_list,seqs_identifiers):
    for i in range(0,seqs_number):
        print seqs_identifiers[i], "\t",
        for j in range(0,seqs_number):
             euclidean_dist=0
             for dicts in occurences_list:
                 key1 = 0 if (i not in dicts) else dicts[i] 
                 key2 = 0 if (j not in dicts) else dicts[j] 
                 euclidean_dist+=(key1-key2)**2
             print euclidean_dist, "\t",
        print ""

arguments = parse()
seqs_number = teiresias_run(arguments.lenght, arguments.input.name)
occurences_list = teiresias_patterns()
seqs_identifiers = identify_headers(arguments.input.name)
euclidean(seqs_number, occurences_list, seqs_identifiers)

