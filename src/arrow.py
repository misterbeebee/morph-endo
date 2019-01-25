from dna import *
from util import *
import traceback
import os


class Arrow:
    # safety throttle FIXME
    MAX_ITERATIONS = 5000
    # Sentinel for finish() -- invalid DNA that (intentionally) terminates execute()
    FINISH_NOW = None

    def __init__(self,
                 prefix_str='',
                 prefix_name=None,
                 prefix_file=None,
                 endo_dna_str=''):
        if prefix_file is not None:
            raise ValueError("prefix_file Not implemented yet")
            # prefix_dna = ...

        endo_dna = DNA(endo_dna_str)
        prefix_dna = DNA(prefix_str)
        self.prefix_name = (prefix_name
                            or prefix_dna_str[0:min(32, len(prefix_str))])
        self.dna = prefix_dna.__add__(endo_dna)
        self.rna = []
        self.iteration = 0

        self.log(INFO, "prefix dna: %s" % str(prefix_dna))
        self.log(INFO, "endo dna: %s" % str(endo_dna))

    def log(self, level, message):
        if ((self.iteration <= 21 and level not in [VERBOSE, SIDE_CHANNEL])
                or level not in [VERBOSE, SIDE_CHANNEL, DEBUG]):
            print(level[0], message)

    def execute(self):
        try:
            while True:
                self.log(
                    INFO if (self.iteration < 100
                             or self.iteration % 100 == 0) else DEBUG,
                    "execute(): iteration=%s  dna=%s" % (self.iteration,
                                                         self.dna.head_str()))
                if (self.iteration > Arrow.MAX_ITERATIONS):
                    print("Max iterations reached! Exiting early.")
                    break
                self.iteration += 1
                p = self.consume_pattern()
                t = self.consume_template()
                self.log(
                    VERBOSE, "execute(): Pattern=%s Template=%s   DNA=%s" %
                    (p, t, self.dna.head_str()))

                self.matchreplace(pat=p, template=t)
        except FinishNow as e:
            print("Finish!")
            self.write()
            return
        except Exception as e:
            self.log(ERROR, "execute() Aborting after error. Dumping RNA")
            self.write()
            traceback.print_exc()

    def consume_pattern(self):
        # WARNING: The spec calls for silently truncating out-of-range slices,
        # but this method will raise exception on then.
        p = Pattern()
        level = 0
        while True:
            # dna_head_str
            h = self.dna[0:3].raw_str()
            # self.log(DEBUG, 'pattern: dna[0:3]: %s' % h)
            if h[0] == 'C':
                self.dna.pop(1)
                p.append(PI)
            elif h[0] == 'F':
                self.dna.pop(1)
                p.append(PC)
            elif h[0] == 'P':
                self.dna.pop(1)
                p.append(PF)
                # Note the use of 'I' as "escape" character
            elif h[0:2] == 'IC':
                self.dna.pop(2)
                p.append(PP)
            elif h[0:2] == 'IP':
                self.dna.pop(2)
                p.append(PSkip(self.consume_nat()))
            elif h[0:2] == 'IF':
                # three bases consumed!
                self.log(SIDE_CHANNEL,
                         "consts pattern consumed junk base: %s" % self.dna[2])
                self.dna.pop(3)
                p.append(PSearch(self.consume_consts()))
            elif h[0:3] == 'IIP':
                self.dna.pop(3)
                level += 1
                p.append(POpen)
            elif h[0:3] in ['IIC', 'IIF']:  # PClose; also EOF
                self.dna.pop(3)
                if level == 0:
                    # FIXME: print more context info, in case this needs debugging.
                    self.log(
                        WARNING if len(p.pitems) == 0 else VERBOSE,
                        "consumed pattern: %s" % p)

                    return p
                else:
                    level = level - 1
                    p.append(PClose)
            elif h[0:3] == 'III':
                self.dna.pop(3)
                self.rna.append(self.dna[0:7].raw_str())
                self.dna.pop(7)
            else:
                self.log(
                    INFO, "consume_pattern() encountered non-decodable DNA."
                    " Finishing execution. (in-spec.)")
                raise FinishNow

    def consume_nat(self):
        if self.dna == empty(): raise FinishNow
        n = 0
        power = 0
        for power in range(0, 1 << 62):
            if len(self.dna) == 0: raise FinishNow
            h = self.dna[0]
            self.dna.pop(1)
            # self.log(VERBOSE, "consume_nat: %s" % h)
            if h == 'P':
                # self.log(DEBUG, "consume_nat: %s" % n)
                return n
            elif h in ['I', 'F']:
                # 0 bit in this position
                pass
            elif h == 'C':
                n += 1 << power
        else:
            self.log(ERROR, "nat() value > 1<<62")
            raise FinishNow

    def consume_consts(self):
        # WARNING: This could raise IndexOutOfBounds -- which would be out-of-spec code bug.
        h1 = self.dna[0]
        if h1 == 'C':
            self.dna.pop(1)
            s = self.consume_consts()
            s.insert(0, I)
            return s
        if h1 == 'F':
            self.dna.pop(1)
            s = self.consume_consts()
            s.insert(0, C)
            return s
        if h1 == 'P':
            self.dna.pop(1)
            s = self.consume_consts()
            s.insert(0, F)
            return s
        h2 = self.dna[0:2].raw_str()
        if h2 == 'IC':
            self.dna.pop(2)
            s = self.consume_consts()
            s.insert(0, P)
            return s
        return empty()

    def consume_template(self):
        # WARNING: The spec calls for silently truncating out-of-range slices,
        # but this method will raise exception on then.
        p = Pattern()  # reusing datatype!
        while True:
            # dna_head_str
            h = self.dna[0:3].raw_str()
            if (len(h) < 1): raise FinishNow
            if h[0] == 'C':
                self.dna.pop(1)
                p.append(PI)
            elif h[0] == 'F':
                self.dna.pop(1)
                p.append(PC)
            elif h[0] == 'P':
                self.dna.pop(1)
                p.append(PF)
            elif (len(h) < 2):
                raise FinishNow
            elif h[0:2] == 'IC':
                # Note the use of 'I' as "escape" character
                self.dna.pop(2)
                p.append(PP)
            elif h[0:2] in ['IP', 'IF']:
                self.dna.pop(2)
                p.append(
                    TProtectedReference(
                        protection_level=self.consume_nat(),
                        reference_id=self.consume_nat()))
            elif (len(h) < 3):
                raise FinishNow
            elif h[0:3] == 'IIP':
                self.dna.pop(3)
                p.append(TReferenceLength(reference_id=self.consume_nat()))

            elif h[0:3] in ['IIC', 'IIF']:
                self.dna.pop(3)
                self.log(VERBOSE, "consumed template: %s " % p)
                return p
            elif h[0:3] == 'III':
                self.dna.pop(3)
                self.rna.append(self.dna[0:7].raw_str())
                self.dna.pop(7)
            else:
                self.log(
                    INFO,
                    "consume_template() encountered non-decodable DNA %s."
                    " Finishing execution. (in-spec.)" % h)
                raise FinishNow

    def matchreplace(self, pat, template):
        self.log(
            DEBUG, "matchreplace(): DNA=%s Pattern=%s     Template=%s" %
            (self.dna.head_str(), pat, template))
        i = 0  # current index in DNA, position after DNA matched so far
        env = []  # environment: [DNA] List of DNA substrings in capture groups
        capture_positions = [
        ]  # capture group positions: [int] Each element is position in DNA where captured group starts.
        for pitem in pat.pitems:
            self.log(VERBOSE, "matchreplace(): Pattern item: %s" % pitem)
            if pitem in [PI, PC, PF, PP]:
                if self.dna[i] == pitem.base:
                    i += 1
                else:
                    # FIXME: this happens a lot. spammy, but could be a bug too.
                    self.log(
                        DEBUG,
                        "Match failed at dna[%i], pattern base=%s; dna[%i:]=%s"
                        % (i, pitem.base, i, self.dna[i:].head_str()))
                    return
            if isinstance(pitem, PSkip):
                i = i + pitem.length
                if i > len(self.dna):
                    self.log(
                        WARNING,
                        "matchreplace(): Match failed dna[%i], Skip=%s too long:"
                        % (i - pitem.length, pitem.length))
                    return
            if isinstance(pitem, PSearch):
                match_start = self.dna.index(pitem.dna_literal, begin=i)
                if match_start is None:
                    self.log(
                        WARNING, "matchreplace(): failed search for %s" %
                        pitem.dna_literal)
                    return
                else:
                    i = match_start + len(pitem.dna_literal)
                    # Advance to character after end of match.
                    self.log(
                        DEBUG,
                        "matchreplace(): Search succeeded, ending at i=%s" % i)
            if pitem == POpen: capture_positions = [i] + capture_positions
            if pitem == PClose:
                env.append(self.dna[capture_positions[0]:i])
                capture_positions = capture_positions[1:]
        # Remove matched DNA.
        self.dna = self.dna[i:]
        # Insert replacement DNA
        self.log(
            DEBUG, "matchreplace(): Match (consumed) length=%s   new DNA=%s" %
            (str(i), self.dna.head_str()))

        self.replace(template=template, env=env, pattern=pat)
        return

    def replace(self, template, env, pattern):
        self.log(
            DEBUG, "replace(): template: %s   env=%s" %
            (template, [str(e) for e in env]))
        r = empty()
        for t in template.pitems:
            if isinstance(t, PBase):
                r.append(t.base)
            elif isinstance(t, TProtectedReference):
                protected = protect(t.protection_level, env[t.reference_id])
                r.extend(protected)
            elif isinstance(t, TReferenceLength):
                if t.reference_id >= len(env):
                    # FIXME: save popped dna from previous iterations, to
                    # analyze how we got here.
                    self.log(
                        ERROR,
                        "env[ref_id] out of range pattern=%s template=%s   env=%s  ref_id=%s "
                        % (pattern, template, env, t.reference_id))
                e = env[t.reference_id]
                r.extend(asnat(len(e)))
        self.log(DEBUG, "replace(): Replacement: %s" % str(r))
        self.dna = r + self.dna
        return

    def write(self):
        outfile_name = (os.getcwd() + "/../output/%s.rna") % self.prefix_name
        self.log(FORCED, "Writing RNA to %s" % outfile_name)
        outfile = open(outfile_name, 'w')
        rna = self.rna
        print(*rna, sep='\n', end='\n', file=outfile)
        outfile.close()
