import os
import sys
# http://stutzbachenterprises.com/blist/blist.html
from blist import blist
import build

from util import *
from arrow import Arrow


def test():
    print(Pattern([
        PBase('C'),
        POpen,
        PSkip(10),
        PSearch('ICFPPCFI'),
        PClose,
    ]))

    icfp = DNA('ICFP')
    assert_eq(icfp[0:2], DNA('IC'))
    assert_eq(icfp[2:0], empty())
    assert_eq(icfp[2:2], empty())
    assert_eq(icfp[2:3], DNA('F'))
    assert_eq(icfp[2], F)
    assert_eq(icfp[2:6], DNA('FP'))
    assert_eq(icfp[2:None], DNA('FP'))
    assert_eq(icfp[6], '')  # Spec is poorly typed for this case.

    dna = DNA('ABC')
    dna.pop(1)
    assert_eq(dna, 'BC')

    cases = {'CIIC': 'I', 'IIPIPICPIICICIIF': '(<2>)P'}
    for prefix in cases:
        a = Arrow(prefix_dna=prefix)
        p = a.consume_pattern()
        assert_eq(
            str(p), cases[prefix],
            '\n\tdna_in=%s\n\tdna_out=%s ' % (prefix, str(a.dna)))

    cases = {'IPPPICIIC': '[0^0]P'}
    for prefix in cases:
        a = Arrow(prefix_dna=prefix)
        p = a.consume_template()
        assert_eq(
            str(p), cases[prefix],
            '\n\tdna_in=%s\n\tdna_out=%s ' % (prefix, str(a.dna)))


PRE_TEST_PREFIX = 'IIPIFFCPICICIICPIICIPPPICIIC'
SELF_TEST_PREFIX = 'IIPIFFCPICFPPICIICCIICIPPPFIIC'

import zipfile


def load_endo_dna():
    zip_ref = zipfile.ZipFile(os.getcwd() + "/../data/endo.dna.zip", 'r')
    dna_str = zip_ref.read("endo.dna").decode("utf-8")
    zip_ref.close()
    return dna_str


def main():
    # test()
    print("hi")
    build.testColors()
    #prefix_name = 'pre_test'
    #prefix_str = PRE_TEST_PREFIX

    prefix_name = 'self_test'
    prefix_str = SELF_TEST_PREFIX

    # prefix_name = 'original_endo'
    # prefix_str = ''
    arrow = Arrow(
        prefix_name=prefix_name,
        prefix_str=prefix_str,
        endo_dna_str=load_endo_dna())
    arrow.execute()
    build.buildRNA(arrow.rna, arrow.prefix_name)


main()
