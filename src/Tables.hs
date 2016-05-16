module Tables 
  (codonTable, degens) where

import qualified Data.Map.Strict as H

{- | table from nucleotide to the bases it represents if its ambiguities are "expanded."
 bases which are not ambiguous expand to themselves. e.g.
>>> H.lookup 'A' degens 
Just "A"
-}
degens = H.fromList  [('A', "A"), ('C', "C"), ('G', "G"), ('T', "T"), ('R', "AG"), ('Y',"CT"), ('S', "GC"), ('W', "AT")]

type CodonTable = H.Map String String

codonTable :: CodonTable
{- | table from codons to Amino Acids.
>>> H.lookup "TTT" codonTable
Just "Phenylalanine"
-}
codonTable = foldr f (H.empty) $ zip codons aas
  where
    f :: ([String], String) -> CodonTable -> CodonTable
    f (ks, v) m = foldr (\k m' -> H.insert k v m') m ks
    
codons =  [
    ["GCT", "GCC", "GCA","GCG"],
    ["TGT", "TGC"],
    ["GAT", "GAC"],
    ["GAA", "GAG"],
    ["TTT", "TTC"],
    ["GGT", "GGC","GGA","GGG"],
    ["CAT", "CAC"],
    ["ATT", "ATC", "ATA"],
    ["AAA", "AAG"],
    ["TTA", "TTG", "CTT","CTC","CTA","CTG"],
    ["ATG"],
    ["AAT", "AAC"],
    ["CCT", "CCC", "CCA","CCG"],
    ["CAA", "CAG"],
    ["CGT", "CGC", "CGA","CGG","AGA","AGG"],
    ["TCT", "TCC", "TCA","TCG","AGT","AGC"],
    ["ACT", "ACC", "ACA","ACG"],
    ["GTT", "GTC", "GTA","GTG"],
    ["TGG"],
    ["TAT", "TAC"]]
aas = ["Alanine",
    "Cysteine",
    "Aspartic acid",
    "Glutamic acid",
    "Phenylalanine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Lysine",
    "Leucine",
    "Methionine, Start",
    "Asparagine",
    "Proline",
    "Glutamine",
    "Arginine",
    "Serine",
    "Threonine",
    "Valine",
    "Tryptophan",
    "Tyrosine"]
