#!/usr/bin/python
import argparse

def parse():
    parser = argparse.ArgumentParser(description = 'Methods based on word (oligomer) frequency', add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', '-i', help='input file', required=True, type=file, nargs=1)
    required.add_argument('--lenght', '-l', help='word (oligomer) lenght', required=True, type=int, nargs=1)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--help', '-help', action='help', help="show this help message and exit")
    return parser.parse_args()

arguments = parse()

