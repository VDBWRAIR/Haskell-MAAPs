{-# LANGUAGE DeriveGeneric #-}
module Lib
    ( someFunc,
      expand,
      AA(..),
      Codon) where
import Data.Hashable 
import qualified Data.HashMap.Strict as H
import Data.Maybe (fromMaybe)
import Data.List (unfoldr, splitAt, findIndices, intersperse, intercalate)
import Data.Csv hiding (lookup)
import GHC.Generics
import Data.Char (ord)
import Data.Text (Text, pack)
import qualified Data.ByteString.Lazy as B
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
-- how handle a *possible* stop codon?
-- should codons with multiple degens output NT degen indices?

type Index = Int

data AA = K | N | T | R | S | I | M | Q | H | P | L | E | D | A | G | V | Z | Y | C | W | F
  deriving (Show, Eq, Enum, Generic)

newtype Codon = Codon String
  deriving (Eq, Show)
instance Hashable Codon where
  hashWithSalt salt (Codon s) = hashWithSalt salt s
  
data Degen  = Insert Codon Index
            | WithN Codon Index
            | StopCodon AA Index Codon [Index]
            | Synonymous AA Index Codon [Index]
            | NonSynonymous [AA] Index Codon [Index]


type CodonTable = H.HashMap Codon AA
data Row = Row { codon ::   !Text
                ,aa ::      !Text
                ,ntPos ::   !Text
                ,aaPos ::   !Index
                ,rowType :: !Text} -- etc. 
 deriving Generic

instance FromNamedRecord Row
instance ToNamedRecord   Row
instance DefaultOrdered  Row

data RowType = InsertT | WithNT | StopCodonT | SynonymousT | NonSynonymousT
  deriving (Show, Eq)
toRow :: Degen -> Row
toRow  (Insert                (Codon nts)  idx)  =  Row (pack nts) (pack "-") (pack "-") idx $ text InsertT
toRow  (WithN                 (Codon nts)  idx)  =  Row (pack nts) (pack "-") (pack "-") idx $ text WithNT
toRow  (StopCodon      aa aaI (Codon nts)  ntI) =   Row (pack nts) (text aa)  (join ntI) aaI $ text StopCodonT
toRow  (Synonymous     aa aaI (Codon nts)  ntI) =   Row (pack nts) (text aa)  (join ntI) aaI $ text SynonymousT
toRow  (NonSynonymous aas aaI (Codon nts)  ntI) =   Row (pack nts) (join aas) (join ntI) aaI $ text NonSynonymousT

text a = pack (show a)
join xs = pack (intercalate "/" $ map show xs)
header =  ["NTs", "NT_pos ","AAs", "AAposition","Type"]

process s = do
  xs <- getDegens s
  return $ encodeDefaultOrderedByNameWith outOptions $ map toRow xs
  where outOptions = defaultEncodeOptions {encDelimiter = fromIntegral (ord '\t')}

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
  B.putStrLn $ fromMaybe (error "Error!") $ process "ATR"
  print $ expand  "ATC"  -- returns its normal AA (synonymous, without degen)
  print $ expand  "zzz"  -- Nothing, not in `degen` list
  print $ expand  "ATRYCSA"  -- Nothing, not divisible by 3
  --print 
  
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
    aas = map char2AA ("KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF" :: String)
    char2AA x = fromMaybe (error ("Bad AA " ++ show x) ) $ lookup x (zip (map (head . show) [K ..] ) [K ..])

ambigNts = ['R', 'S', 'K', 'V', 'B', 'Y', 'W', 'M', 'H']
-- (length aas, length codons) -- these need to be equal 
