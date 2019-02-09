import dna
from dna import *
from util import *
import traceback
import os


class Arrow:
    # safety throttle FIXME
    MAX_ITERATIONS = 5000
    # Sentinel for finish() -- invalid DNA that (intentionally) terminates execute()
    FINISH_NOW = None
    MAX_MATCHREPLACE_CONSECUTIVE_FAILURES = 10

    def __init__(
            self,
            prefix_str='',
            prefix_name=None,
            prefix_file=None,
            endo_dna_str='',
            # Keep copy of consumed DNA for debugging
            archive_consumption=True):
        if prefix_file is not None:
            raise ValueError("prefix_file Not implemented yet")
            # prefix_dna = ...

        endo_dna = DNA(endo_dna_str)
        prefix_dna = DNA(prefix_str)
        self.prefix_name = (prefix_name
                            or prefix_str[0:min(32, len(prefix_str))])
        self.dna = prefix_dna.__add__(endo_dna)
        # transient remembrance of recently consumed DNA
        # reset frequently
        self.consumed = empty()
        self.archive_consumption = archive_consumption

        self.matchreplace_consecutive_failures = 0
        self.rna = []
        self.iteration = 0

        self.log(INFO, lambda: "prefix DNA: %s" % str(prefix_dna))
        self.log(INFO, lambda: "endo's DNA: %s" % str(endo_dna))

    def log(self, level, message):
        if ((self.iteration <= 40 and level not in [VERBOSE, SIDE_CHANNEL])
                or level not in [VERBOSE, SIDE_CHANNEL, DEBUG]):
            print(level[0], message() if callable(message) else message)

    def execute(self):
        try:
            while True:
                self.log(
                    INFO if (self.iteration < 100
                             or self.iteration % 100 == 0) else DEBUG,
                    lambda: "execute():\n\titeration:\t%s\n\tDNA:       \t%s" % (self.iteration,
                                                         self.dna.head_str()))
                if (self.iteration > Arrow.MAX_ITERATIONS):
                    log(ERROR, "Max iterations reached! Exiting early.")
                    break
                assert_eq(len(self.consumed), 0)
                p = self.consume_pattern()
                # FIXME: print more context info, in case this needs debugging.
                self.log(
                    WARNING if
                    ((len(p.pitems) == 0) or
                     (len([i for i in p.pitems if i == POpen]) != len(
                         [i for i in p.pitems if i == PClose]))) else VERBOSE,
                    lambda: "pattern():\n\tpattern:\t%s\n\tfrom:\t%s" % (p, str(self.consumed))
                )
                self.consumed = empty()

                t = self.consume_template()
                self.log(
                    VERBOSE,
                    lambda: "\ttemplate:\t%s\n\tfrom:\t%s " % (t, self.consumed)
                )
                self.consumed = empty()
                self.log(
                    VERBOSE,
                    lambda: "execute():\n\tPattern:\t%s\n\tTemplate:\t%s" % (p, t)
                )

                replacement = self.matchreplace(pat=p, template=t)
                if replacement == const_empty:
                    self.matchreplace_consecutive_failures += 1

                    self.log(
                        "WARNING",
                        lambda: ("Iteration: %s   matchreplace() failed\t\n"
                                 'Pattern: %s \tTemplate:%s') % (self.iteration, p, t)
                    )
                    if self.matchreplace_consecutive_failures > Arrow.MAX_MATCHREPLACE_CONSECUTIVE_FAILURES:
                        log(
                            "ERROR",
                            lambda: "Too many matchreplace failures; giving up"
                        )
                        raise FinishNow

                else:
                    self.matchreplace_consecutive_failures = 0
                if (len(replacement) - len(self.consumed)) > 100000:
                    self.log(
                        FORCED,
                        lambda: ("matchreplace(): DNA grew a lot!!!"
                                 "\n\tPattern:  \t%s"
                                 "\n\tTemplate:  \t%s"
                                 "\n\tConsumed: \t%s"
                                 "\n\tAdded:    \t%s"
                                 "\n\tnew DNA:  \t%s" % (p, t, self.consumed, replacement, self.dna))
                    )
                self.log(
                    DEBUG,
                    lambda: ("matchreplace():\n\tconsumed:\t%s\n\tadded:    \t%s" % (self.consumed, replacement))
                )
                self.consumed = empty()
                self.iteration += 1
        except FinishNow as e:
            print("Finish!")
            self.write()
            return
        except Exception as e:
            self.log(ERROR, "execute() Aborting after error. Dumping RNA")
            self.write()
            traceback.print_exc()

    def consume(self, n):
        # Archiving consumed DNA is expensive, so we can disable it
        if self.archive_consumption:
            self.consumed.extend(self.dna.pop(n))
        else:
            self.dna.pop(n)

    def consume_pattern(self):
        # WARNING: The spec calls for silently truncating out-of-range slices,
        # but this method will raise exception on then.
        p = Pattern()
        level = 0
        # consume n bases from DNA, and record them in consumed
        while True:
            # dna_head_str
            h = self.dna[0:3].raw_str()
            # self.log(DEBUG, lambda: 'pattern: dna[0:3]: %s' % h)
            if len(h) < 1:
                log(WARNING, "end of DNA during pattern!")
                return p
            if h[0] == 'C':
                self.consume(1)
                p.append(PI)
            elif h[0] == 'F':
                self.consume(1)
                p.append(PC)
            elif h[0] == 'P':
                self.consume(1)
                p.append(PF)
                # Note the use of 'I' as "escape" character
            elif h[0:2] == 'IC':
                self.consume(2)
                p.append(PP)
            elif h[0:2] == 'IP':
                self.consume(2)
                p.append(PSkip(self.consume_nat()))
            elif h[0:2] == 'IF':
                # three bases consumed!
                self.log(SIDE_CHANNEL,
                         "consts pattern consumed junk base: %s" % self.dna[2])
                self.consume(3)
                p.append(PSearch(self.consume_consts()))
            elif h[0:3] == 'IIP':
                self.consume(3)
                level += 1
                p.append(POpen)
            elif h[0:3] in ['IIC', 'IIF']:  # PClose; also EOF
                self.consume(3)
                if level == 0:
                    return p
                else:
                    level = level - 1
                    p.append(PClose)
            elif h[0:3] == 'III':
                self.consume(3)
                self.rna.append(self.dna[0:7].raw_str())
                self.consume(7)
            else:
                self.log(
                    INFO, "consume_pattern() encountered non-decodable DNA."
                    " Finishing execution. (in-spec.)")
                raise FinishNow

    def consume_nat(self):
        if self.dna == empty(): raise FinishNow
        i = self.dna.blist.index('P')
        nat = self.dna[0:i].raw_str()
        if len(nat) > 30:
            self.log(INFO, lambda: "consume_nat: %s" % (len(nat)))
        self.consume(i + 1)
        return parse_nat(nat)

    def consume_consts(self):
        # WARNING: This could raise IndexOutOfBounds -- which would be out-of-spec code bug.
        s = []
        while True:
            h1 = self.dna[0]
            if h1 in ['C', 'F', 'P']:
                unquoted = dna.UNQUOTE_BASE[h1]
                self.consume(1)
                s.append(unquoted)
            else:
                # WARNING: This could raise IndexOutOfBounds -- which would be out-of-spec code bug.
                h2 = self.dna[0:2].raw_str()
                if h2 == 'IC':
                    unquoted = dna.UNQUOTE_BASE[h2]
                    self.consume(2)
                    s.append(unquoted)
                else:
                    break
        return DNA(''.join(s))

    def consume_template(self):
        # WARNING: The spec calls for silently truncating out-of-range slices,
        # but this method will raise exception on then.
        p = Pattern()  # reusing datatype!
        while True:
            # dna_head_str
            h = self.dna[0:3].raw_str()
            if (len(h) < 1): raise FinishNow
            if h[0] == 'C':
                self.consume(1)
                p.append(PI)
            elif h[0] == 'F':
                self.consume(1)
                p.append(PC)
            elif h[0] == 'P':
                self.consume(1)
                p.append(PF)
            elif (len(h) < 2):
                raise FinishNow
            elif h[0:2] == 'IC':
                # Note the use of 'I' as "escape" character
                self.consume(2)
                p.append(PP)
            elif h[0:2] in ['IP', 'IF']:
                self.consume(2)
                p.append(
                    TProtectedReference(
                        protection_level=self.consume_nat(),
                        reference_id=self.consume_nat()))
            elif (len(h) < 3):
                raise FinishNow
            elif h[0:3] == 'IIP':
                self.consume(3)
                p.append(TReferenceLength(reference_id=self.consume_nat()))

            elif h[0:3] in ['IIC', 'IIF']:
                self.consume(3)
                return p
            elif h[0:3] == 'III':
                self.consume(3)
                self.rna.append(self.dna[0:7].raw_str())
                self.consume(7)
            else:
                self.log(
                    INFO,
                    "consume_template() encountered non-decodable DNA %s."
                    " Finishing execution. (in-spec.)" % h)
                raise FinishNow

    def matchreplace(self, pat, template):
        "Returns replacement."
        self.log(
            INFO if (self.iteration <= 40) else DEBUG,
            "matchreplace():\n\tDNA:    \t%s\n\tPattern:\t%s\n\tTemplate:\t%s"
            % (self.dna.head_str(), pat, template))
        i = 0  # current index in DNA, position after DNA matched so far
        # environment: [DNA] List of DNA substrings in capture groups, ordered
        # by position of pClose
        env = []
        # unclosed capture group positions: [int] Each element is position in
        # DNA where yet-unclosed captured group starts.
        capture_positions = []
        for pitem in pat.pitems:
            self.log(VERBOSE,
                     lambda: "matchreplace(): Pattern item: %s" % pitem)
            if pitem in [PI, PC, PF, PP]:
                if self.dna[i] == pitem.base:
                    i += 1
                else:
                    # FIXME: this happens a lot. spammy, but could be a bug too.
                    self.log(
                        DEBUG,
                        "Match failed at:\n\tpattern base=%s\n\tdna[%i:]=%s" %
                        (pitem.base, i, self.dna[i:].head_str()))
                    return const_empty
            if isinstance(pitem, PSkip):
                i = i + pitem.length
                if i > len(self.dna):
                    self.log(WARNING,
                             ("matchreplace(): Match failed at\n\tdna[%i]"
                              "\n\tSkip=%s too long") % (i - pitem.length,
                                                         pitem.length))
                    return const_empty
            if isinstance(pitem, PSearch):
                match_start = self.dna.index(pitem.dna_literal, begin=i)
                if match_start is None:
                    self.log(
                        WARNING, "matchreplace(): failed search for %s" %
                        pitem.dna_literal)
                    return const_empty
                else:
                    i = match_start + len(pitem.dna_literal)
                    # Advance to character after end of match.
                    self.log(
                        DEBUG,
                        "matchreplace(): Search succeeded, ending at i=%s dna=%s..."
                        % (i, self.dna[match_start:match_start +
                                       len(pitem.dna_literal)]))
            if pitem == POpen: capture_positions.insert(0, i)
            if pitem == PClose:
                env.append(self.dna[capture_positions[0]:i])
                capture_positions = capture_positions[1:]
        # Remove matched DNA.
        self.consume(i)
        # Insert replacement DNA
        return self.prepend_replacement(
            template=template, env=env, pattern=pat)

    def prepend_replacement(self, template, env, pattern):
        "Returns reference to prefix added to DNA"
        self.log(
            DEBUG, "replace():\n\ttemplate:\t%s\n\tenv:       \t%s" %
            (template, '\n\t            \t'.join([str(e) for e in env])))
        references_used = blist.sortedset()
        for t in template.pitems:
            if isinstance(t, TProtectedReference) or isinstance(
                    t, TReferenceLength):
                references_used.add(t.reference_id)
        if references_used != blist.sortedset(range(0, len(env))):
            self.log(
                WARNING, "references_used %s != range(0,len(env)=%s)" %
                (references_used, len(env)))
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
                    self.log(ERROR,
                             ("env[ref_id] out of range\n\tpattern:\t%s"
                              "\n\ttemplate:\t%s\n\tenv:\t%s\n\tref_id:\t%s") %
                             (pattern, template, env, t.reference_id))
                e = env[t.reference_id]
                r.extend(asnat(len(e)))
        self.dna = r + self.dna
        return r

    def write(self):
        outfile_name = (os.getcwd() + "/../output/%s.rna") % self.prefix_name
        self.log(FORCED, "Writing RNA to %s" % outfile_name)
        outfile = open(outfile_name, 'w')
        rna = self.rna
        print(*rna, sep='\n', end='\n', file=outfile)
        outfile.close()
