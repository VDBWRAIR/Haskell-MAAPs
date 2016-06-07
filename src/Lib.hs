{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
module Lib where
import qualified Data.HashMap.Strict as H
import Data.Hashable 
import Data.Maybe (fromMaybe)
import Data.List (unfoldr, splitAt, findIndices, intersperse, intercalate, intersect, nub)
import GHC.Generics
import Data.Char (ord, toUpper)
import Data.Csv hiding (lookup)
import Data.Text (Text, pack)
import Control.Monad (forM_, zipWithM_)
import qualified Data.ByteString.Lazy as B
import qualified Data.ByteString.Lazy.Char8 as C
--import Bio.Sequence.Fasta (readFasta, toStr, seqdata, seqid)
import Bio.Sequence.Fasta 
import Bio.Core.Sequence
import Options.Generic 
import Types
import Data.FixedList hiding (length, head)

{-
X TODO: Fix formatting somehow
X TODO: Report fasta IDs somehow
TODO: More tests better
TODO: CI
X TODO: Convert maybes to Either with error descriptions
X TODO: Row should not be stringly typed

 Post-processing
 * note stop codons
 check frame shifts (length check)


 this should also report position of degenerate NT (could be done after return) -- do this by just keeping the triple as you go. AA position is just its index
 Should check if synonymous -- do this by checking size of AA list
 report inserts, Ns

 in order to do this need AAs, triples (codons), index
 case on info to get Degen; case on Degen to get Report
 frame shift should be handled seperately
 (but a frame shift is likely to be noticed sooner, and needs to be handled explicitly
 before / within `expand`, because we will get a triple of size 1 or 2.
 Not sure how to handle frame shifts within polypeptide, because there are no stop codons to go off of.
 | NT pos | NTs | AAs   | AAposition | Type |
 | 9      | AKK | R/Z/Z |    3       | Non-synonymous |
 | 4      | ANA |  *    |    *       | N     |
 | 4      | A-A |  *    |    *       | Gap |
 | 4      | ATA |  A    |    2       | Stop Codon|
 | 4      | ATA |  A    |    2       | Stop Codon|
 how handle a *possible* stop codon?
 should codons with multiple degens output NT degen indices?
-}

type Error = String 
{- For every record in the fasta file found in opts, print to stdout the result CSV
followed by newline. -}
run :: Options -> IO ()
run opts = do
  recs <- readFasta fastaFile
  C.putStrLn header' 
  forM_ recs printOne
  where
    fastaFile = unHelpful $ fasta opts -- unwrap the help-annotated type
    printOne x = either putStr (`forM_` C.putStrLn) $ process opts x
    header' = B.intercalate "\t" $ toList fields

{- | remove "Degens" that are synonymous or normal.
Also removes the stop codons that are at the end--as these are legitimate stops and shouldn't be marked.
>>> filterDegens [NormalCodon, WithN (Codon "TNT") 2, StopCodon Z 3 (Codon "ATC") []]
[WithN (Codon "TNT") 2]
-}
filterDegens :: Options -> [Degen] -> [Degen]
filterDegens opts xs = filter (\x -> (not $ isSynonymous x) || keepSyn ) $ filter (not . isNormal) $ dropStopCodon $ xs
  where
    keepSyn = unHelpful $ syn opts
    isSynonymous (Synonymous' _ _ _ _) = True
    isSynonymous _                    = False
    isNormal NormalCodon = True
    isNormal _           = False
    dropStopCodon [] = [] -- this pattern suggests that there is no stop codon!
    dropStopCodon ((StopCodon _ _ _ _ ):[]) = []
    dropStopCodon (x:xs)                    = x : dropStopCodon xs
    
{- pull out the ids for each sequence by removing content after the seperator.
create the CSV list for each degen and encode. -}
process :: Options -> Sequence -> (Either Error [B.ByteString])
process opts s@(Seq header' (SeqData {unSD=seq}) _ ) = do
  xs <- filterDegens opts <$> getDegens (map toUpper seq')
  let rows = map (record . toList . (fieldList id')) xs
  return $ [(encodeWith outOptions rows)]
  where
    id'  =  Id $  takeWhile (not . (`elem` [' ', '_', '-', '/'])) $ C.unpack $ unSL $ header'
    seq' =  C.unpack seq
    outOptions = defaultEncodeOptions {encDelimiter = fromIntegral (ord '\t')}
    getDegens :: String -> Either Error [Degen] 
    getDegens s = do
      (cds, aas) <- unzip <$> expand s
      return $ zipWith3 toDegen cds (nub aas) [1..]
    
toList :: FieldList a -> [a]  
toList = foldr (:) []

-- | convert the codon and AA information into the Degen type, representing wether it's an insert, contains an N, etc.
-- >>> toDegen (Codon "ATC") [F] 1
-- NormalCodon
-- >>> toDegen (Codon "AAA") [Z] 3
-- StopCodon Z 3 (Codon "AAA") []
-- >>> toDegen (Codon "TNT") [] 5
-- WithN (Codon "TNT") 5
toDegen :: Codon -> [AA] -> Index -> Degen
toDegen cdn@(Codon nts) aas i
  | '-' `elem` nts = Insert cdn i
  | 'N' `elem` nts = WithN  cdn i
  | otherwise = case aas of
       ([])  ->   error $ "Invalid state; uncaught bad trasnlation at index " ++ show i ++ "codon " ++ nts
       (Z:[])  -> StopCodon     Z   i cdn  ntIdxs
       (aa:[]) | (not $ doIntersect nts ambigNts) -> NormalCodon 
       (aa:[])    -> Synonymous' aa  i cdn  ntIdxs
       aas   ->   NonSynonymous  aas i cdn  ntIdxs
  where
    ntIdxs = map (\n -> ((i - 1) * 3) + 1 + n) $ findIndices (`elem` ambigNts) nts
    doIntersect x y = not $ null $ intersect x y


-- | Expand a string of arbitrary size into a the Codons and amino acids the string represents after expanding
-- degenerate nucleotides and seeeing what is coded for. Splits the string up into triples and maps `expandTriple` over it.
-- >>> expand "SAGTAA" 
-- Right [(Codon "SAG",[E,Q]),(Codon "TAA",[Z])]
--
-- Expects the string to divide evenly into triplets, and for every base to be a valid nucleotide or gap (represtented as '-')
-- >>> expand "AGTCC"
-- Left "Permutation AA lookup in codon position 2, CC not found."
-- >>> expand "CCzCC"
-- Left "Bases in codon position 1, CCz not found." 
--
-- codons with gaps come back as codons which code to nothing.
-- >>> expand "CC-ATG"
-- Right [(Codon "CC-",[]),(Codon "ATG",[M])] 
expand :: String -> Either Error [(Codon, [AA])] 
expand xs = (zip codons) <$> expandeds
  where
    codons = map Codon $ takeWhile (not . null) $ unfoldr (Just . (splitAt 3)) xs
    expandeds = sequence $ zipWith tryExpand [1..] codons
    tryExpand i cdn@(Codon x) = if (not $ null $ intersect x "-N") then Right [] else expandTriple i  cdn 

toEither :: b -> Maybe a -> Either b a
toEither b a = maybe (Left b) Right a

-- Given the index of the first NT and the codon, expand the degenerate bases and return all possible
-- amino acids the codon could code to, or an error message.
expandTriple :: Int -> Codon -> Either Error [AA]
expandTriple i (Codon xs) = do
    degens' <- toEither ("Bases in codon position " ++ show i ++ ", " ++ xs ++ " not found.") $ lookups allBases xs
    -- sequence here is used as a sort of ordered permutation,
    -- e.g. sequence ["A","ATG"] == ["AA","AT","AG"] 
    let perms = map Codon $ sequence degens'
    aas' <- toEither ("Permutation AA lookup in codon position " ++ show i ++ ", " ++ xs ++ " not found.") $ lookups codonTable perms
    return aas'
 where
  lookups m xs' = sequence $ map (`H.lookup` m) xs' 
  allBases :: H.HashMap Char String
  allBases = H.fromList (ambig ++ nonAmbig ++ otherBases)
    where
      nonAmbig   = zipString "ACTG"
      otherBases = zipString "N-"
      ambigExp   = map snd ambigTable
      ambig      = zip ambigNts ambigExp
      zipString :: String -> [(Char, String)]
      zipString xs = zip xs $ map (:[]) xs
    
{- | table from codons to Amino Acids.
>>> H.lookup (Codon "TTT") codonTable
Just F
-}
codonTable :: H.HashMap Codon AA
codonTable = H.fromList $ zip codons aas
  where
    -- TODO: why this works and tests
    codons = map Codon $ sequence $ replicate 3 "ACGT" 
    aas = map char2AA ("KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF" :: String)
    char2AA x = fromMaybe (error ("Bad AA " ++ show x) ) $ lookup x (zip (map (head . show) [K ..] ) [K ..])

-- (length aas, length codons) -- these need to be equal 

fields :: FieldList B.ByteString
fields = fromFoldable' $  ["ID", "Codon", "NTPos", "AA", "AAPos", "RowType"]

-- Convert degen information into a list of "Fields" (wrapped text) for the CSV decoder to consume
fieldList :: Id -> Degen -> FieldList Field
fieldList (Id id') x = case x of
  (Insert                (Codon nts)  idx) -> tf' id' :. tf' nts :. tf idx       :. "-"        :. "-"     :. tf Is_Gap :. Nil
  (WithN                 (Codon nts)  idx) -> tf' id' :. tf' nts :. tf idx       :. "-"        :. "-"     :. tf Has_N :. Nil
  (StopCodon      aa aaI (Codon nts)  ntI) -> tf' id' :. tf' nts :. "-"          :. aaf aa      :. tf  aaI :. tf Stop_Codon :. Nil
  (Synonymous'     aa aaI (Codon nts)  ntI) -> tf' id' :. tf' nts :. jf "," ntI :. aaf aa      :. tf aaI  :. tf Synonymous :. Nil
  (NonSynonymous aas aaI (Codon nts)  ntI) -> tf' id' :. tf' nts :. jf "," ntI :. aajf "/" aas :. tf aaI  :. tf Non_Synonymous :. Nil
  NormalCodon -> error "NormalCodon shouldn't be output"
  where 
      aaf = toField . pack . aaShow 
      tf' = toField
      tf a = toField $ pack (show a)
      jf c xs = toField $ pack (intercalate c $ map show xs) 
      aajf c xs = toField $ pack (intercalate c $ map aaShow xs) 

ambigNts = map fst ambigTable

{- | table from nucleotide to the bases it represents if its ambiguities are "expanded."
>>> lookup 'R' ambigTable
Just "AG"
-}
ambigTable = [('R', "AG"),
              ('Y', "CT"),
              ('S', "GC"),
              ('W', "AT"),
              ('K', "GT"),
              ('M', "AC"),
              ('B', "CGT"),
              ('D', "AGT"),
              ('H', "ACT"),
              ('V', "ACG")]
--ambigNts = ['M' ,  'R' , 'W' , 'S' ,  'Y' , 'K' , 'V' , 'H' , 'D' , 'B'] 
--ambigExp = ["AG", "GC", "TG", "ACG", "CGT", "CT", "AT", "CA", "ACT"]
