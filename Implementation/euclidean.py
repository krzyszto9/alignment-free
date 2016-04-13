#!/usr/bin/python
import argparse
import subprocess

def teiresias_run(lenght,input_file):
    cmd = 'teiresias -l{0} -w{0} -k1 -i{1} -p -s'.format(lenght,input_file)
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None, shell=True)


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


arguments = parse()
teiresias_run(arguments.lenght, arguments.input.name)

