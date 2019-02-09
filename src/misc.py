#  Code from offline to merge into the main project.

import io
import struct
import itertools
import math

# low-order bits of color are fiddled to make it easier to spot differences
# when printed as text bytes.
red = (250, 0, 0)
green = (0, 251, 0)
yellow = (253, 253, 0)
blue = (0, 0, 252)

white = (255, 255, 255)

# Known DNA marker prefixes:
#   IFP
#   IFPICFPP
# The very end of Endo is two markers with a 'P' (0) betwen them:
#  IFPICFPPCIFF P IFPICFPPCIFP
#  IFPICFPPCIFF  "red zone" start
#  IFPICFPPCIFP "red zone" end End of endo

# (<7509409><IFPICFPPCIFF>)(<1>)IFPICFPPCIFP)
# This copies some data into the "red zone"
# That data is later removedf with:

# In analysis below: X+ means a long block of X, unspecified size
# Endo is approx 2743x2743 (~7.5Mbases)
# with several distinct regions (in order):
# [0-13500] 1350 blocks of:
#   IIIPxP, where x is  elem of {IIIII, CCCCC, FFFFF, CCIFF, FFICC}, most often x
#       but first block is IIIPIPIIPC  (bases 6,9,10 are "special")
#   TODO: decode this as a string ins a 5-bit alphabset.941
# IIIII looks like "background")
# [  13500:  13580] 110 bases of approximately the above pattern, "degraded"
# [  13580:  13712] 102 bases of (too complex to describe, but not random)
# [  13712:  13762] 50 bases: (Cx24 IC Cx24 IC)   (== quoted(24-bit 0) x2
# [  13708:  14860] 1152 bases: (Ix8 P)x143   (8-bit naturals, all = 0)
# [  14860:  14935] 75 bases: IFPCFFP  Ix23 P  IFPFP(PI) Ix7 P IFPFP(IP) Ix23 P   (23bit 0 and 7 bit 0 naturals)
# [  14935:  27735] 12800 bases: I on C background. When viewed as ((16x16) bitmap)x50:
#                  picture of: "    ACHTUNG! PORTABLE NETWORK GRAPHICS FOLLOWS    "
# [  27735:  27736] 1 base: P  (end of previous message)
# [  27736:  29520] 1784 (223*8) bases: binary data  [png.dna]
# [  29520:  29522] 2 bases: PP
# [  29522: 223634] 194112 (3033*8*8) bases: binary data [png2.dna]
# [ 223634: 223635] 1 base: P
# [ 223635: 223642] 7 bases: marker IFPFCCF
# [ 223642: 223666] 14 bases:   CCIICCIIIIIIIIIIIIIIIIIP: binary:51
# [ 223666: 223675] 9 bases:   CCCCCCCCP: binary:127
# [ 223675: 224820] 1145 (143*8 + 1) bases: I*1143 C P
# [ 224820: 224869] 49 bases: CPFPPCFCIIPFICCFCFIFCCCFCIIIFIFCCICCPCPFFFCCCFPIC
# [ 224869: 224889] 20 bases: markers:
#                    IFPPFI
#                    IFPIPP    occurs at 224875,
#                             957390 (9 more)  984182    (specifically: 959075, 962071, 962557, 967824, 968975, 971655, 973163, 974571, 980168)
#                            2196089, (10-20 more) 2918470;
#                            4398959, (5-10 more) 4910991,
#                            5719651
#                    IFPCIFIC  occurs at 224881 and 2199375 (almost 2M bases later)
#
#                    Note: 1M to 2.19M, 3M-4.39M, 4.4M-5.7M  have none of these markers
# [ 224889: 224914] 25 bases: Ix24 P  24-bit binary 0.
# [ 224914: 224916] 2 bases: FF
# [ 224916: 234132] 9216 bases: 384 (128*3) 23-bit Binary numbers.
# [ 234132: 295182] 61050 bases: 13180 opcodes with header PIIC (or IIC, and ends with P ?)
#                                 Many codes are pairs of PPICPCCPCCP PxN C
#                                 Most codes are many binary numbers
#                                 After 8 opcodes is an extra insertion of base 'C'. Perhaps an error to delete.

# [295182: 304398] 9216 bases: 384 (128*3) 23-bit Binary numbers, again.
# [305398: 367588] 62190 bases: IIC opcodes again
# [367588: 376804] 9216 bases: 384 (128*3) 23-bit Binary numbers, again.
#  376ish: 500ish: IIC opcodes
#  500ish : 513ish 24bit 0
#         : 520ish ICP mixed, then ICFP mixed
#          : 530ish 24bit numbers
# a couple more times
# 8 bit numbers

LSB_FIRST_IN_BYTE = "LSB_FIRST_IN_BYTE"
LSB_LAST_IN_BYTE = "LSB_LAST_IN_BYTE"
PAD_END = "PAD_END"
PAD_START = "PAD_START"


def dna_bin_to_char_list(dna_str, least_significant_bit, pad_at):
    """Converts a long binary string to a list of 8-bit chars.
    If length is not a multiple of 0, pads at `pad_at`.
       pad_at=PAD_END may be weird if least_significant bit is LAST_IN_BYTE!
    """
    ret = []
    padding = ['I'] * ((8 - (len(dna_str) % 8)) % 8)
    dna_list = list(dna_str)
    if pad_at == PAD_END:
        dna_list.extend(padding)
    elif pad_at == START:
        dna_list = padding + dna_list

    else:
        raise Exception

    if least_significant_bit == LSB_LAST_IN_BYTE:
        # Move LSB to front of bytes (side_effect: reverse order of chars)
        dna_list.reverse()

    # Core algo assumes LSB is first in byte
    for char_pos in range(0, len(dna_list) // 8):
        n = 0
        ret.append(parse_nat(dna_list[char_pos * 8:(char_pos * 8) + 8]))
    if least_significant_bit == LSB_LAST_IN_BYTE:
        # un-reverse order of chars from side_effect aboce
        ret.reverse()

    return ret


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


def remove_whitespace(string):
    return ''.join([x for x in string if not x.isspace()])


#  IC (near-random mix)
#  I with short blocks of C, and single Ps at semi-regular frequency
#  CP with some I and very little F (genome errors?)
colormap = {
    'I': red,
    'C': green,
    'F': yellow,
    'P': blue,
}


# create a rectangle off area equal to or slightly larger than n
# return dimensions, and padding (excess area of rectangle)
def rectangle(n, rows=None, cols=None):
    rows = rows if rows else int(math.ceil(n / cols if cols else math.sqrt(n)))
    cols = cols if cols else int(math.ceil(n / rows if rows else math.sqrt(n)))
    padding = (cols * rows) - n
    print("rectangle: rows:%s cols:%s padding:%s" % (rows, cols, padding))
    return (rows, cols, padding)


def read_dna_file(dna_dir, dna_file_basename):
    dna_file = open(file=(dna_dir + dna_file_basename), mode='r')
    dna = remove_whitespace(dna_file.read())
    dna_file.close()
    print('read DNA file done')
    return dna


def write_sparkline_from_file(dna_file_basename,
                              dna_dir,
                              print_binary=True,
                              print_text=False,
                              rows=None,
                              cols=None):
    return sparkline(
        read_dna_file(dna_file_basename=dna_file_basename, dna_dir=dna_dir),
        print_binary=print_binary,
        print_text=print_text,
        dna_file_basename=dna_file_basename,
        output_dir=dna_dir,
        rows=rows,
        cols=cols)


def sparkline(dna,
              dna_file_basename,
              output_dir,
              print_binary=True,
              print_text=False,
              rows=None,
              cols=None):
    """
     Print_text    P3 format is expensove, for debugging.
 """

    (rows, cols, padding) = rectangle(len(dna), rows=rows, cols=cols)

    out_filename = file = ('%s/%s.%sx%s.p6binary.ppm' %
                           (output_dir, dna_file_basename, cols, rows))
    fileio = io.FileIO(out_filename, mode='w')
    pixel_tuples = [colormap[x] for x in list(dna)]

    if print_binary:
        fileio.write(bytes('P6 %s %s 255 ' % (cols, rows), 'ascii'))
        int_subpixels_iter = itertools.chain.from_iterable(pixel_tuples)
        print('built int_subpixels_iter')
        for rgb in pixel_tuples:
            fileio.write(struct.pack('BBB', *rgb))
        fileio.write(b''.join([struct.pack('B', 0)] * 3 * padding))
        fileio.close()
        print("wrote", out_filename)

    if print_text:
        text_out = open(
            file=('%s/%s.%sx%s.p3text.ppm' % (output_dir, dna_file_basename,
                                              cols, rows)),
            mode='w')
        print('P3 %s %s 255 ' % (cols, rows), file=text_out)
        print(
            '  '.join([' %3s %3s %3s' % rgb for rgb in pixel_tuples]),
            file=text_out)
        print('  ', '  '.join(['   0   0   0'] * padding), file=text_out)
        text_out.close()


def read_dna_to_binary_file(dna_file_basename, header_bytes,
                            least_significant_bit, pad_at, dna_dir):
    dna = read_dna_file(dna_dir=dna_dir, dna_file_basename=dna_file_basename)
    byte_list = dna_bin_to_char_list(
        dna_str=dna,
        least_significant_bit=least_significant_bit,
        pad_at=pad_at)
    out_filename = ('%s/%s.%s.%s.binary' % (dna_dir, dna_file_basename,
                                            least_significant_bit, pad_at))
    fileio = io.FileIO(file=out_filename, mode='w')
    fileio.write(header_bytes)
    for char_byte in byte_list:
        fileio.write(struct.pack('B', char_byte))
    fileio.close()
    print("wrote", out_filename)


def main():
    endo_dna_base = 'endo.dna'
    dna_dir = '/Users/michaelroger/Desktop/personal/endo/dna/'
    for cols in [
            # 1000,
            # 1536,
            1144
    ]:
        write_sparkline_from_file(
            dna_file_basename='endo.dna', dna_dir=dna_dir, cols=cols)
# write_sparkline_from_file(dna_file_basename='hidden_message.14935-27735.dna', cols=16, dna_dir=dna_dir)
    for binary_dna_file in [
            "mp3.29522-223634.dna", "sample_binary_file.01space.dna",
            "png.27736-29520.dna"
    ]:
        pass


#    read_dna_to_binary_file(dna_file_basename=binary_dna_file,
#                   least_significant_bit=LSB_FIRST_IN_BYTE, pad_at=PAD_END, header_bytes=b'   ', dna_dir=dna_dir)

main()
