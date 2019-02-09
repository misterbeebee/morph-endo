import os
import sys
import datetime

# http://stutzbachenterprises.com/blist/blist.html
from blist import blist

import build
from util import *
from arrow import Arrow
from dna import *
import test


def decode_pattern(pattern_str):
    return Arrow(prefix_str=pattern_str).consume_pattern()


PRE_TEST_PREFIX = 'IIPIFFCPICICIICPIICIPPPICIIC'
SELF_TEST_PREFIX = 'IIPIFFCPICFPPICIICCIICIPPPFIIC'

import zipfile


def load_endo_dna():
    zip_ref = zipfile.ZipFile(os.getcwd() + "/../data/endo.dna.zip", 'r')
    dna_str = zip_ref.read("endo.dna").decode("utf-8")
    zip_ref.close()
    return dna_str


def remove_comments_and_space(string):
    # remove comments and then linebreaks:
    noncomments = ''.join([line.split('#')[0] for line in string.split('\n')])
    return ''.join([x for x in noncomments if not x.isspace()])


def read_dna_file(dna_dir, dna_file_basename):
    dna_file = open(file=(dna_dir + dna_file_basename), mode='r')
    dna = remove_comments_and_space(dna_file.read())
    dna_file.close()
    print('read DNA file done:', dna_file_basename)
    return dna


def run_prefix_from_file(filename):
    dna = read_dna_file(
        dna_dir=(os.getcwd() + '/../data/prefix/'), dna_file_basename=filename)
    arrow = Arrow(
        prefix_name=filename,
        prefix_str=dna,
        endo_dna_str=load_endo_dna(),
        archive_consumption=True)
    arrow.execute()
    build.buildRNA(arrow.rna, arrow.prefix_name)


def main():
    log(FORCED, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    test.test()

    prefix_name = 'pre_test'
    prefix_str = PRE_TEST_PREFIX

    # prefix_name = 'self_test'
    # prefix_str = SELF_TEST_PREFIX

    # prefix_name = 'original_endo'
    # prefix_str = ''

    run_prefix_from_file('guide-0000100.dna')


main()
