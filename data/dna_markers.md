# DNA markers

Example:

*  Pattern: `(<7509409><IFPICFPPCIFF>)(<1>)(IFPICFPPCIFP)`
*  Template: `[0^0][1^0][1^0][1^0][1^0][2^0]`

DNA Layout implied by Pattern:

1. pattern
1. template
1. FIXED_ZONE: 7509409 bases
1. variable zone
1. `IFPICFPPCIFF` VARIABLE\_ZONE_END marker
1. special base
1. `IFPICFPPCIFP` TEMP\_DATA_START marker
1. trailing temp data that is discarded this iteration


I assume that pattern like this is placed there (with accompanying template)
with the knowledge that the template ends right before a known fixed chunk
of DNA that never changes size (but might be spot-patched).


The template:

1. `[0^0]` Preserves the FIXED_ZONE and variable zone up to the `IFPICFPPCIFF` marker.
1. `[1^0][1^0][1^0][1^0]` replicates a snippet (to 4x).
1. `[2^0]` preserves the tail data up to the "EOF" marker, discarding temp data

Note that both markers have form `IFPICFPPCIF?` aka `IFP ICFPPC IF?`,
which decodes to pattern: `<PCFFI><` (is truncated at start of const-string)


### Example: Patch in-place:

Search for a marker (starts with `IFP` as above), and preserve it (`[0^0]`)
and match a base and replace it:

*  Pattern:        `(<IFPP>)F`
*  Template:       `[0^0]P`

changes `F` to `P`


### Example: Insert a constant binary-integer

* Pattern:        `<0>(<4536645>(<800>)<2971964>)`
* Template:       `[0^0][1^0]IIIIIIIIIIIIIIIIIIIIIIIPIIIIIIIIIIIIIIIIIIIIIIIP`


FIXME...