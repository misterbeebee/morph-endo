import blist
from util import *


class DNA():
    "Like BList, but has DNA semantics for index/slice"

    def __init__(self, sequence):
        if isinstance(sequence, blist.blist):
            self.blist = sequence
        elif isinstance(sequence, str):
            self.blist = blist.blist(sequence)
        else:
            raise TypeError

    def __eq__(self, sequence):
        "Flexible equality that implicitly converts 'str' and 'blist' types"
        if isinstance(sequence, DNA):
            return self.blist == sequence.blist
        if isinstance(sequence, blist.blist):
            return self.blist == sequence
        elif isinstance(sequence, str):
            return self.blist == blist.blist(sequence)
        else:
            raise TypeError

    def raw_str(self):
        return ''.join(self.blist)

    def head_str(self, prefix_len=10):
        # max length to print
        stop = min(len(self), prefix_len)
        string = ''.join(self.blist[0:stop])
        return "%s%s %8s" % (string, '' if len(self) <= prefix_len else '...',
                             len(self))

    def __str__(self):
        prefix_len = 10  # max length to print
        stop = min(len(self), prefix_len)
        string = ''.join(self.blist[0:stop])
        return "%s%s %8s" % (string, '' if len(self) <= prefix_len else '...',
                             len(self))

    def __add__(self, other):
        return DNA(self.blist + other.blist)

    def append(self, base):
        self.blist.append(base)

    def extend(self, dna):
        self.blist += dna.blist

    def __len__(self):
        return len(self.blist)

    # bounds-safe item/slice (clamps out-of-range index to end)
    # Warning: Does not suppport negative indexes, which are not in DNA spec.
    # Warning: Does not suppport slice step, which are not in DNA spec.
    # Warning: DNA spec's bounds-safety is incompatible with for...in loop
    #    (including functional methods like 'join')  (which expects IndexError)
    def __getitem__(self, key):
        if isinstance(key, slice):
            assert key.start >= 0
            assert key.stop is None or key.stop >= 0
            assert key.step is None
            stop = key.stop if key.stop is not None else len(self)
            return DNA(self.blist[max(0, key.start):(min(stop, len(self)))])
        elif isinstance(key, int):
            assert key >= 0
            if key < len(self):
                return self.blist[key]
            else:
                return ''

    def insert(self, *args):
        self.blist.insert(*args)

    def pop(self, n):
        "Remove and return n elements from head of sequence"
        stop = min(n, len(self))
        ret = DNA(self.blist[0:stop])
        del self.blist[0:stop]
        return ret

    def index(self, sub_dna, begin):
        "Returns index of sub_dna in self, or else None"
        log(DEBUG, "index(): Searching for %s" % sub_dna)
        # Unwrap to avoid excess allocations of DNA object.
        sub_dna_b = sub_dna.blist
        dna_b = self.blist
        found_after_begin = ''.join(dna_b[begin:]).find(''.join(sub_dna_b))
        return None if found_after_begin == -1 else begin + found_after_begin


class Base:
    I = 'I'
    C = 'C'
    F = 'F'
    P = 'P'


I = Base.I
C = Base.C
F = Base.F
P = Base.P

# Base :: char: I C F P
# DNA :: blist[Base]

empty = lambda: DNA('')
# FIXME: add enforced freeze
const_empty = empty()


class PItem:
    pass


class PBase(PItem):
    def __init__(self, base):
        self.base = base

    def __str__(self):
        return self.base


PI = PBase(I)
PC = PBase(C)
PF = PBase(F)
PP = PBase(P)


class PSkip(PItem):
    length = None

    def __init__(self, length):
        self.length = length

    def __str__(self):
        return '<%s>' % self.length


class PSearch(PItem):
    def __init__(self, dna):
        self.dna_literal = (dna if isinstance(dna, DNA) else DNA(dna))

    def __str__(self):
        return '<%s>' % self.dna_literal.raw_str()


class PGroup(PItem):
    def __init__(self, str):
        self.str = str

    def __str__(self):
        return self.str


POpen = PGroup('(')
PClose = PGroup(')')


def p_capture(pitems):
    return [POpen] + pitems + [PClose]


class TProtectedReference(PItem):
    "Writes a matched capture group, quoted 0 or more times"

    def __init__(self, protection_level, reference_id):
        self.protection_level = protection_level
        self.reference_id = reference_id

    def __str__(self):
        return "[%s^%s]" % (self.reference_id, self.protection_level)


class TReferenceLength(PItem):
    "Writes the length of a matched capture-group, as a literal number"

    def __init__(self, reference_id):
        self.reference_id = reference_id

    def __str__(self):
        return "|%s|" % self.reference_id


def protect(level, dna):
    "Protect (quote/escape) dna, level times."
    if level == 0: return dna
    else: return protect(level - 1, quote(dna))


QUOTE_BASE = {
    'I': 'C',
    'C': 'F',
    'F': 'P',
    'P': 'IC',
}

UNQUOTE_BASE = {
    'C': 'I',
    'F': 'C',
    'P': 'F',
    'IC': 'P',
}


def quote(dna):
    return DNA(''.join([QUOTE_BASE[b] for b in dna]))


def parse_nat(nat_str):
    "Convert a non-P-terminated binary-IC string to int."
    n = 0
    for power in range(0, len(nat_str)):
        h = nat_str[power]
        if h in ['I', 'F']:
            # 0 bit in this position
            pass
        elif h == 'C':
            n += 1 << power
    return n


def asnat(n):
    rev_dna = []
    while n > 0:
        if n % 2 == 0:
            rev_dna.append('I')
        else:
            rev_dna.append('C')
        n = n // 2
    rev_dna.append('P')
    return DNA(''.join(reversed(rev_dna)))


# Data model hack!!
# Pattern is also used for Template, since some are both!
class Pattern:
    def __init__(self, pitems=None):
        self.pitems = (pitems or [])[:]  # defensive copy from input

    def __eq__(self, other):
        return self.pitems == other.pitems

    def __str__(self):
        return ''.join([str(pitem) for pitem in self.pitems])

    def append(self, pitem):
        self.pitems.append(pitem)
