module Lib
    ( someFunc
    ) where

import qualified Data.HashMap.Strict as H
import Data.List (unfoldr, splitAt)
--import Tables (codonTable, degens)


someFunc = do
  print $ expand  "ATR" -- "Isoleucine", -- "Methionine Start",
  print $ expand  "ATC"  -- returns its normal AA (synonymous, without degen)
  print $ expand  "zzz"  -- Nothing, not in `degen` list
  print $ expand  "ATRYCSA"  -- Nothing, not divisible by 3
  --print (length aas, length codons) -- these need to be equal

expand xs =  sequence $ map expandTriple $ triples xs
  where triples xs = takeWhile (not . null) $ unfoldr (Just . (splitAt 3)) xs
{- Post-processing
 - * note stop codons
 - check frame shifts (length check)
-}
-- this should also report position of degenerate NT (could be done after return) -- do this by just keeping the triple as you go. AA position is just its index
-- Should check if synonymous -- do this by checking size of AA list
-- report inserts, Ns
--
-- in order to do this need AAs, triples (codons), index
expandTriple xs = do
    degens' <- lookups allBases xs
    let perms = sequence degens'
    aas' <- lookups codonTable perms
    return aas'
 where lookups m xs' = sequence $ map (`H.lookup` m) xs'
    

--newtype Codon = Codon String
--type CodonTable = H.HashMap Codon AA
--type CodonTable = H.HashMap Char String
--codonTable :: CodonTable
codonTable = H.fromList $ zip (concat codons) aas'
  where
    aas' = concat $ zipWith replicate lengths aas
    lengths = map length codons

allBases :: H.HashMap Char String
allBases = H.fromList (ambig ++ nonAmbig ++ otherBases) 

nonAmbig = [('A', "A"), ('C', "C"), ('G', "G"), ('T', "T")]

otherBases = [('N', "N"), ('-', "-")] 

ambig =  [ ('R', "AG"), ('Y',"CT"),
           ('S', "GC"), ('W', "AT"),
           ('K', "TG"), ('M', "CA"),
           ('V', "ACG"), ('H', "ACT"),
           ('B', "CGT")]  -- Doesn't include `N`

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

--codonTable :: CodonTable
--codonTable = H.fromList $ zip  codons aas
--  where 
--    codons = map Codon $ sequence $ replicate 3 "ACGT" 
--    aas = map char2AA "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
--    char2AA = lookup $ zip (map (head . show) [R..V]) [R..V]

--table = H.fromList $ zip (replicateM 3 "ACGT") aa
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
data AA =  R |  N |  D |  C |  E |  Q |  G |  H |  I |  L |  K |  M |  F |  P |  S |  T |  W |  Y |  V
  deriving (Show, Bounded, Eq)
