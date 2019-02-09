import os
import sys
import datetime

import build
from util import *
from arrow import Arrow
from dna import *


def test():

    print(Pattern([
        PBase('C'),
        POpen,
        PSkip(10),
        PSearch(DNA('ICFPPCFI')),
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

    mp3_clue = 'IIPIFFCPICPICIICIPCCCPIICIPPPFFCFFFFIIC'
    prefix = mp3_clue
    a = Arrow(prefix_str=prefix)
    p = a.consume_pattern()
    t = a.consume_template()
    assert_eq(
        str(p), '(<IFPFP>)<7>',
        '\n\tdna_in=%s\n\tdna_out=%s ' % (prefix, str(a.dna)))
    assert_eq(
        str(t), '[0^0]CCICCCC',
        '\n\tdna_in=%s\n\tdna_out=%s ' % (prefix, str(a.dna)))

    cases = {'CIIC': 'I', 'IIPIPICPIICICIIF': '(<2>)P'}
    for prefix in cases:
        a = Arrow(prefix_str=prefix)
        p = a.consume_pattern()
        assert_eq(
            str(p), cases[prefix],
            '\n\tdna_in=%s\n\tdna_out=%s ' % (prefix, str(a.dna)))

    cases = {'IPPPICIIC': '[0^0]P'}
    for prefix in cases:
        a = Arrow(prefix_str=prefix)
        t = a.consume_template()
        assert_eq(
            str(t), cases[prefix],
            '\n\tdna_in=%s\n\tdna_out=%s ' % (prefix, str(a.dna)))

    build.testColors()
