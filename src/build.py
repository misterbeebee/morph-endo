from util import *
import os
from PIL import Image
import numpy
import blist


def isRGB(color):
    return isinstance(color, tuple)


def isTransparency(color):
    return isinstance(color, int)


# must be in clockwise order, for %4 math.
N = 0
E = 1
S = 2
W = 3

black = (0, 0, 0)
red = (255, 0, 0)
green = (0, 255, 0)
yellow = (255, 255, 0)
blue = (0, 0, 255)
magenta = (255, 0, 255)
cyan = (0, 255, 255)
white = (255, 255, 255)

transparent = 0
opaque = 255


def transparentBitmap():
    return [[(black, transparent) for y in range(0, 600)]
            for x in range(0, 600)]


def addColor(color):
    raise Error


def build(rna):
    Builder().build(rna)


def addColorsTo(colors, bktRGB, bktA):
    for c in colors:
        addColorTo(c, bktRGB, bktA)


def addColorTo(color, bktRGB, bktA):
    if isRGB(color):
        bktRGB.insert(0, color)
    elif isTransparency(color):
        bktA.insert(0, color)


def currentPixel(bktRGB, bktA):
    (r, g, b) = averageRGB(bktRGB, 0)
    a = average(bktA, 255)
    return (tuple(map(lambda c: (c * a) // 255, (r, g, b))), a)


def averageRGB(values, default):
    return (average([r for (r, g, b) in values], default),
            average([g for (r, g, b) in values], default),
            average([b for (r, g, b) in values], default))


def average(values, default):
    if values == []:
        return default
    else:
        return sum(values) // len(values)


def testColor(colors, pixel):
    bktRGB = []
    bktA = []
    addColorsTo(colors, bktRGB, bktA)
    assert_eq(currentPixel(bktRGB, bktA), pixel)


def testColors():
    b = black
    r = red
    m = magenta
    w = white
    y = yellow
    c = cyan
    t = transparent
    o = opaque
    testColor([t, o, o], ((0, 0, 0), 170))
    testColor([b, y, c], ((85, 170, 85), 255))
    testColor([y, t, o], ((127, 127, 0), 127))
    testColor([b] * 18 + [r] * 7 + [m] * 39 + [w] * 10 + [o] * 3 + [t] * 1,
              ((143, 25, 125), 191))


def move(pos, d):
    (x, y) = pos
    if d == N: return (x, (y - 1) % 600)
    elif d == E: return ((x + 1) % 600, y)
    elif d == S: return (x, (y + 1) % 600)
    elif d == W: return ((x - 1) % 600, y)
    else: raise Exception


def turnCounterClockwise(d):
    return (d - 1) % 4


def turnClockwise(d):
    return (d + 1) % 4


def getPixel(bitmap, pos):
    (x, y) = pos
    # indices reverse to make array in row-major order
    return (bitmap)[y][x]


def setPixel(bitmap, pos, bktRGB, bktA):
    pixel = currentPixel(bktRGB, bktA)
    setPixelFaster(bitmap, pos, pixel)


def setPixelFaster(bitmap, pos, pixel):
    (x, y) = pos
    # log(INFO, lambda: "setPixel(): bitmap[%s][%s] = %s " % (x, y, pixel))
    # indices reverse to make array in row-major order
    bitmap[y][x] = pixel


def line(pos0, pos1, bitmap, bktRGB, bktA):
    # log(INFO, lambda:
    #    "line(%s, %s) -> %s" % (pos0, pos1, currentPixel(bktRGB, bktA)))
    (x0, y0) = pos0
    (x1, y1) = pos1
    deltax = x1 - x0
    deltay = y1 - y0
    d = max(abs(deltax), abs(deltay))
    c = 1 if (deltax * deltay <= 0) else 0  # why??
    x = x0 * d + (d - c) // 2
    y = y0 * d + (d - c) // 2
    pixel = currentPixel(bktRGB, bktA)
    for i in range(0, d):
        #log("INFO", "setPixelFaster(%s, %s)" % (((x // d), (y // d)), pixel))
        setPixelFaster(bitmap, ((x // d), (y // d)), pixel)

        x = x + deltax
        y = y + deltay
    setPixelFaster(bitmap, (x1, y1), pixel)


def tryfill(position, bitmap, bktRGB, bktA):
    new = currentPixel(bktRGB, bktA)
    old = getPixel(bitmap, position)
    log(INFO, lambda: "tryfill(%s) %s->%s" % (position, old, new))
    if new != old: fill(position, old, bitmap, bktRGB, bktA)


def fill(pos, initial, bitmap, bktRGB, bktA):
    log(INFO, lambda: "fill(%s, %s)" % (pos, initial))
    queue = blist.sortedset([pos])
    # Change all adjacent pixels of same initial color to new color
    queued = 1
    filled = 0
    iters = 0
    skipped = 0
    pixel = currentPixel(bktRGB, bktA)
    setPixelFaster(bitmap, pos, pixel)
    while queue:
        if iters % 40000 == 0:
            log(
                INFO, "fill() iters=%s skipped=%s len(queue)=%s q=%s pixel=%s"
                % (iters, skipped, len(queue), queue[0],
                   bitmap[queue[0][0]][queue[0][1]]))
        iters += 1
        (x, y) = queue.pop()
        filled += 1
        if x > 0:
            next = (x - 1, y)
            if getPixel(bitmap, next) == initial:
                setPixelFaster(bitmap, next, pixel)
                queue.add(next)
                queued += 1
            else:
                skipped += 1
        if x < 599:
            next = (x + 1, y)
            if getPixel(bitmap, next) == initial:
                setPixelFaster(bitmap, next, pixel)
                queue.add(next)
                queued += 1
            else:
                skipped += 1
        if y > 0:
            next = (x, y - 1)
            if getPixel(bitmap, next) == initial:
                setPixelFaster(bitmap, next, pixel)
                queue.add(next)
                queued += 1
            else:
                skipped += 1
        if y < 599:
            next = (x, y + 1)
            if getPixel(bitmap, next) == initial:
                setPixelFaster(bitmap, next, pixel)
                queue.add(next)
                queued += 1
            else:
                skipped += 1
    log(
        INFO, "fill() queued %s and filled %s pixels in %s iterations" %
        (queued, filled, iters))


addCount = 0


def addBitmap(bitmaps, b):
    global addCount
    log(INFO, lambda: "addBitmap()")
    draw([bitmaps[0]], "pre-add.%s" % addCount)
    addCount += 1
    if len(bitmaps) < 10:
        bitmaps.insert(0, b)
    else:
        log(WARNING, "Ignoring bitmap over limit (10)")


popCount = 0


# opposite of python's  "pop"!
def pop_bitmap(bitmaps):
    global popCount
    bitmap = bitmaps[0]
    draw([bitmap], "popped.%s" % popCount)
    popCount += 1
    del bitmaps[0]


def compose(bitmaps):
    log(INFO, lambda: "compose()")
    # Merge front 2 bitmaps
    if len(bitmaps) >= 2:
        for x in range(0, 600):
            for y in range(0, 600):
                ((r0, g0, b0), a0) = (bitmaps[0])[x][y]
                ((r1, g1, b1), a1) = (bitmaps[1])[x][y]
                pixel = ((r0 + ((r1 * (255 - a0)) // 255),
                          g0 + ((g1 * (255 - a0)) // 255),
                          b0 + ((b1 * (255 - a0)) // 255)),
                         a0 + ((a1 * (255 - a0)) // 255))
                bitmaps[1][x][y] = pixel
        pop_bitmap(bitmaps)

    else:
        log(WARNING, "Tried to compose single bitmap")


def clip(bitmaps):
    log(INFO, lambda: "clip()")
    # Apply bitmap1's alpha channel as a filter to next bitmap
    if len(bitmaps) >= 2:
        for x in range(0, 600):
            for y in range(0, 600):
                ((r0, g0, b0), a0) = (bitmaps[0])[x][y]
                ((r1, g1, b1), a1) = (bitmaps[1])[x][y]
                pixel = (((r1 * a0) // 255, (g1 * a0) // 255,
                          (b1 * a0) // 255), (a1 * a0) // 255)
                bitmaps[1][x][y] = pixel
                log(
                    INFO,
                    lambda: "clip():  bitmaps[0][%s][%s] = %s " % (x, y, pixel)
                )
        pop_bitmap(bitmaps)
    else:
        log(WARNING, "Tried to clip single bitmap")


def buildRNA(rna, prefix_name):
    bktRGB = []  # [color]
    bktA = []  # [color]
    position = (0, 0)
    mark = (0, 0)
    direction = E
    bitmaps = [transparentBitmap()]  # [bitmap]
    addColor = lambda c: addColorTo(c, bktRGB, bktA)

    for r in rna:
        if (r) == 'PIPIIIC': addColor(black)
        elif r == 'PIPIIIP': addColor(red)
        elif r == 'PIPIICC': addColor(green)
        elif r == 'PIPIICF': addColor(yellow)
        elif r == 'PIPIICP': addColor(blue)
        elif r == 'PIPIIFC': addColor(magenta)
        elif r == 'PIPIIFF': addColor(cyan)
        elif r == 'PIPIIPC': addColor(white)
        elif r == 'PIPIIPF': addColor(transparent)
        elif r == 'PIPIIPP': addColor(opaque)

        elif r == 'PIIPICP': (bktRGB, bktA) = ([], [])

        elif r == 'PIIIIIP': position = move(position, direction)
        elif r == 'PCCCCCP': direction = turnCounterClockwise(direction)
        elif r == 'PFFFFFP': direction = turnClockwise(direction)

        elif r == 'PCCIFFP': mark = position
        elif r == 'PFFICCP': line(position, mark, bitmaps[0], bktRGB, bktA)
        elif r == 'PIIPIIP': tryfill(position, bitmaps[0], bktRGB, bktA)
        elif r == 'PCCPFFP': addBitmap(bitmaps, transparentBitmap())
        elif r == 'PFFPCCP': compose(bitmaps)
        elif r == 'PFFICCF': clip(bitmaps)
        else:
            log(SIDE_CHANNEL, "build(): junk RNA: %s" % r)
    draw(bitmaps, prefix_name)  # all alpha values are set to 255!


def draw(bitmaps, prefix_name):
    for i in range(len(bitmaps)):
        bitmap = bitmaps[i]
        outfile_base_name = (os.getcwd() + "/../output/%s.%s") % (prefix_name,
                                                                  i)
        ppm = False
        if ppm:
            outfile_name = outfile_base_name + ".ppm.txt"
            outfile = open(outfile_name, 'w')
            # http://netpbm.sourceforge.net/doc/ppm.html
            # P3 is RAW ASCII PPM -- human readable but least efficient.
            print("P3 600 600 255", end='\n', file=outfile)
            for x in range(0, 600):
                print(
                    *[
                        "%3s %3s %3s" % (r, g, b)
                        for ((r, g, b), a) in bitmap[x]
                    ],
                    sep='   ',
                    end='\n',
                    file=outfile)
            outfile.close()
        # FIXME: need to flip x<->y axes?
        bitmap_array = numpy.array(
            [[[r, g, b] for ((r, g, b), a) in col] for col in bitmap],
            order='C',
            dtype=numpy.uint8)
        log(FORCED,
            "draw(): drawing endo to %s" % (outfile_base_name + ".png"))
        Image.fromarray(bitmap_array, 'RGB').save(outfile_base_name + ".png",
                                                  'png')
        # raise Exception  # quit for debug
