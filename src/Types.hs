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


-- type FilePath = String
-- , geneSource :: Maybe GenBankSource <?> "Genbank ID, Genbank file, or csv file"
-- commandline stuff
data Seperator = Tab | Comma
  deriving (Show, Eq, Generic)
      
deriving instance Read Seperator
instance ParseField Seperator
      
deriving instance Read RowType
instance ParseField RowType

data Options = Options {   rowFilter :: [RowType] <?> "What show -- Insert, etc.."
                         , sep :: First Seperator <?> "Comma or Tab output"
                         , fasta :: FilePath <?> "Input aligned fasta file"
                         , align :: Bool <?> "Align the sequence first?"
                         , syn   :: Bool <?> "Display synonymous AAs?"}
             deriving (Generic, Show)

instance ParseRecord Options

-- Primary types

type Index = Int
newtype Id = Id String

data AA = K | N | T | R | S | I | M | Q | H | P | L | E | D | A | G | V | Z | Y | C | W | F
  deriving (Show, Eq, Enum, Generic)

aaShow Z = "!"
aaShow x = show x

newtype Codon = Codon String
  deriving (Eq, Show)

instance Hashable Codon where
  hashWithSalt salt (Codon s) = hashWithSalt salt s
  
type CodonTable = H.HashMap Codon AA

data RowType = Is_Gap | Has_N | Stop_Codon | Synonymous | Non_Synonymous
  deriving (Show, Eq, Generic)

data Degen  = Insert Codon Index
            | WithN Codon Index
            | StopCodon AA Index Codon [Index]
            | Synonymous' AA Index Codon [Index]
            | NonSynonymous [AA] Index Codon [Index] -- Codon index or AA Index? Should make newtypes
            | NormalCodon
  deriving (Show, Eq)

-- type for the csv row list              
type FieldList = FixedList6

--instance ToRecord Degen where
--  toRecord = record . toList . fieldList


--data GenBankSource = GBFile FilePath | CSVFile FilePath | GBID String
--  deriving (Show, Generic, ParseRecord)

-- unused newtype CodonIndex = CodonIndex Int -- make this part of codon?
-- unused newtype AAIndex = AAIndex Int -- make this part of codon?

-- *** Exception: Data.Csv.Encoding.namedRecordToRecord: header contains name "RowType" which is not present in the named record
