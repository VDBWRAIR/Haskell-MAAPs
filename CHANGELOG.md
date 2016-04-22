0.01
====
Basic functionality.


0.02
====
Codons which translated to the same AA multiple times were reported as non-synonymous. This is because the code checked the length of the AAs rather than the uniqueness. Fixed by adding a call to `nub`.
Dropped support for Frameshifts; they were invalidly called when a sequence had a length not divisible by 3. In reality, need to account for gaps within a sequence.
Fixed amiguous Nucleotide table, which was just wrong
