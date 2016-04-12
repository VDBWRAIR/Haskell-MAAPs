module Lib
    ( someFunc,
      expand,
      AA(..),
      Codon) where

import Data.Hashable 
import qualified Data.HashMap.Strict as H
import Data.Maybe (fromMaybe)
import Data.List (unfoldr, splitAt, findIndices, intersperse, intercalate)
{- Post-processing
 - * note stop codons
 - check frame shifts (length check)
-}

-- this should also report position of degenerate NT (could be done after return) -- do this by just keeping the triple as you go. AA position is just its index
-- Should check if synonymous -- do this by checking size of AA list
-- report inserts, Ns
--
-- in order to do this need AAs, triples (codons), index
-- case on info to get Degen; case on Degen to get Report
-- frame shift should be handled seperately
-- (but a frame shift is likely to be noticed sooner, and needs to be handled explicitly
-- before / within `expand`, because we will get a triple of size 1 or 2.
-- Not sure how to handle frame shifts within polypeptide, because there are no stop codons to go off of.
-- | NT pos | NTs | AAs   | AAposition | Type |
-- | 9      | AKK | R/Z/Z |    3       | Non-synonymous |
-- | 4      | ANA |  *    |    *       | N     |
-- | 4      | A-A |  *    |    *       | Gap |
-- | 4      | ATA |  A    |    2       | Stop Codon|
-- | 4      | ATA |  A    |    2       | Stop Codon|


type Index = Int

data AA = K | N | T | R | S | I | M | Q | H | P | L | E | D | A | G | V | Z | Y | C | W | F
  deriving (Show, Eq, Enum)

newtype Codon = Codon String
  deriving (Eq, Show)
instance Hashable Codon where
  hashWithSalt salt (Codon s) = hashWithSalt salt s
-- how handle a *possible* stop codon?
-- should codons with multiple degens output NT degen indices?
data Degen  = Insert Codon Index
            | WithN Codon Index
            | StopCodon AA Index Codon [Index]
            | Synonymous AA Index Codon [Index]
            | NonSynonymous [AA] Index Codon [Index]


type CodonTable = H.HashMap Codon AA
data Row = Row { codon :: Codon
                ,aa :: Maybe [AA]
                ,ntPos :: Maybe [Index]
                ,aaPos :: Index
                ,rowType :: RowType} -- etc.

data RowType = InsertT | WithNT | StopCodonT | SynonymousT | NonSynonymousT
  deriving (Show, Eq)

toRow :: Degen -> Row
toRow  (Insert cdn idx)                 =  Row cdn Nothing Nothing idx InsertT
toRow  (WithN  cdn idx)                 =  Row cdn Nothing Nothing idx WithNT
toRow  (StopCodon      aa aaI cdn  ntI) =  Row cdn (Just [aa]) (Just ntI) aaI StopCodonT
toRow  (Synonymous     aa aaI cdn  ntI) =  Row cdn (Just [aa]) (Just ntI) aaI SynonymousT
toRow  (NonSynonymous aas aaI cdn  ntI) =  Row cdn (Just aas)  (Just ntI) aaI NonSynonymousT

header = intercalate "/" ["NTs", "NT_pos ","AAs", "AAposition","Type"]

showRow :: Row -> String
showRow Row{codon=(Codon nts),aa=aas,ntPos=ntI,aaPos=aaI,rowType=rowT} =
  intercalate ['\t'] cols
    where
      cols = [nts, ntI', aas', aaI', rowT']
      aas'  = intercalate "/" $ fromMaybe ["-"] $ (map show) <$> aas
      aaI'  = show aaI
      ntI'  = defaultStr $ show <$> ntI
      rowT' = show rowT
      defaultStr x = fromMaybe "-" x

toDegen :: Codon -> [AA] -> Index -> Degen
toDegen cdn@(Codon nts) aas i
  | '-' `elem` nts = Insert cdn i
  | 'N' `elem` nts = WithN  cdn i
  | otherwise = case aas of
       (Z:[])  -> StopCodon     Z   i cdn  ntIdxs
       (aa:[]) -> Synonymous    aa  i cdn  ntIdxs
       aas   ->   NonSynonymous aas i cdn  ntIdxs
  where
    ntIdxs = map (+ (i * 3) ) $ findIndices (`elem` ambigNts) nts

someFunc = do
  print $ expand  "ATR" -- "Isoleucine", -- "Methionine Start",
  print $ expand  "ATC"  -- returns its normal AA (synonymous, without degen)
  print $ expand  "zzz"  -- Nothing, not in `degen` list
  print $ expand  "ATRYCSA"  -- Nothing, not divisible by 3
  --print (length aas, length codons) -- these need to be equal
  
getDegens :: String -> Maybe [Degen]
getDegens s = do
  (cds, aas) <- unzip <$> expand s
  return $ zipWith3 toDegen cds aas [1..]
  
  

expand :: String -> Maybe [(Codon, [AA])]
expand xs = (zip codons) <$> expandeds
  where
    codons = map Codon $ takeWhile (not . null) $ unfoldr (Just . (splitAt 3)) xs
    expandeds = sequence $ map expandTriple $ codons 
   
expandTriple :: Codon -> Maybe [AA]
expandTriple (Codon xs) = do
    degens' <- lookups allBases xs
    let perms = map Codon $ sequence degens'
    aas' <- lookups codonTable perms
    return aas'
 where lookups m xs' = sequence $ map (`H.lookup` m) xs' 

allBases :: H.HashMap Char String
allBases = H.fromList (ambig ++ nonAmbig ++ otherBases) 
  where
    nonAmbig   = zipString "ACTG"
    otherBases = zipString "N-"
    ambigExp   = ["AG", "GC", "TG", "ACG", "CGT", "CT", "AT", "CA", "ACT"]
    ambig      = zip ambigNts ambigExp
    zipString :: String -> [(Char, String)]
    zipString xs = zip xs $ map (:[]) xs
codonTable :: CodonTable
codonTable = H.fromList $ zip codons aas
  where 
    codons = map Codon $ sequence $ replicate 3 "ACGT" 
    aas = map char2AA "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF"
    char2AA x = fromMaybe (error ("Bad AA " ++ show x) ) $ lookup x (zip (map (head . show) [K ..] ) [K ..])



--ambig =  [ ('R', "AG"), ('Y', "CT"),
--           ('S', "GC"), ('W', "AT"),
--           ('K', "TG"), ('M', "CA"),
--           ('V', "ACG"),('H', "ACT"),
--           ('B', "CGT")]  -- Doesn't include `N`

ambigNts = ['R', 'S', 'K', 'V', 'B', 'Y', 'W', 'M', 'H']
