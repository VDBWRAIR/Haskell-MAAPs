{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE DataKinds         #-}  -- these two needed for help output
{-# LANGUAGE TypeOperators     #-}
module Types where
import Data.Hashable 
import qualified Data.HashMap.Strict as H
import Data.Csv hiding (lookup)
import Data.Text (Text, pack)
import GHC.Generics
import Data.List (intercalate)
import Options.Generic
import Data.FixedList 
import qualified Data.ByteString.Lazy as B
import Data.ByteString.Internal as I

type Index = Int
newtype CodonIndex = CodonIndex Int -- make this part of codon?
newtype AAIndex = AAIndex Int -- make this part of codon?

-- type FilePath = String
data Seperator = Tab | Comma
  deriving (Show, Eq, Generic)
-- , geneSource :: Maybe GenBankSource <?> "Genbank ID, Genbank file, or csv file"
      
deriving instance Read Seperator
instance ParseField Seperator
      
deriving instance Read RowType
instance ParseField RowType
data Options = Options {   rowFilter :: [RowType] <?> "What show -- Insert, etc.."
                         , sep :: First Seperator <?> "Comma or Tab output"
                         , fasta :: FilePath <?> "Input aligned fasta file"
                         , align :: Bool <?> "Align the sequence first?"}
             deriving (Generic, Show)

instance ParseRecord Options

data AA = K | N | T | R | S | I | M | Q | H | P | L | E | D | A | G | V | Z | Y | C | W | F
  deriving (Show, Eq, Enum, Generic)

newtype Codon = Codon String
  deriving (Eq, Show)
instance Hashable Codon where
  hashWithSalt salt (Codon s) = hashWithSalt salt s
  
type CodonTable = H.HashMap Codon AA

data RowType = InsertT | WithNT | StopCodonT | SynonymousT | NonSynonymousT | FrameShiftT
  deriving (Show, Eq, Generic)

data Degen  = Insert Codon Index
            | WithN Codon Index
            | StopCodon AA Index Codon [Index]
            | Synonymous AA Index Codon [Index]
            | NonSynonymous [AA] Index Codon [Index]
            | FrameShift Index -- Codon index or AA Index? Should make newtypes
            | NormalCodon 
              
-- *** Exception: Data.Csv.Encoding.namedRecordToRecord: header contains name "RowType" which is not present in the named record

type FieldList = FixedList6
--instance ToRecord Degen where
--  toRecord = record . toList . fieldList
toList :: FieldList a -> [a]  
toList = foldr (:) []
newtype Id = Id String 
fields :: FieldList B.ByteString
fields = fromFoldable' $  ["ID", "Codon", "NTPos", "AA", "AAPos", "RowType"]

fieldList :: Id -> Degen -> FieldList Field
fieldList (Id id') x = case x of
  (Insert                (Codon nts)  idx) -> tf' id' :. tf' nts :. tf idx       :. "-"        :. "-"     :. tf InsertT :. Nil
  (WithN                 (Codon nts)  idx) -> tf' id' :. tf' nts :. tf idx       :. "-"        :. "-"     :. tf WithNT :. Nil
  (FrameShift idx)                         -> "-"     :. "-"     :. tf idx       :. "-"        :. "-"     :. tf FrameShiftT :. Nil
  (StopCodon      aa aaI (Codon nts)  ntI) -> tf' id' :. tf' nts :. jf " :." ntI :. tf aa      :. tf  aaI :. tf StopCodonT :. Nil
  (Synonymous     aa aaI (Codon nts)  ntI) -> tf' id' :. tf' nts :. jf " :." ntI :. tf aa      :. tf aaI  :. tf SynonymousT :. Nil
  (NonSynonymous aas aaI (Codon nts)  ntI) -> tf' id' :. tf' nts :. jf " :." ntI :. jf "/" aas :. tf aaI  :. tf NonSynonymousT :. Nil
  NormalCodon -> error "NormalCodon shouldn't be output"
  where
      tf' = toField
      tf a = toField $ pack (show a)
      jf c xs = toField $ pack (intercalate c $ map show xs)
