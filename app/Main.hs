{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE StandaloneDeriving     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}  -- these two needed for help output
{-# LANGUAGE TypeOperators     #-}
module Main where

import Lib  (process, RowType(..))
import Options.Generic
import Options.Applicative.Types
-- type FilePath = String
data Seperator = Tab | Comma
  deriving (Show, Eq, Generic)
-- , geneSource :: Maybe GenBankSource <?> "Genbank ID, Genbank file, or csv file"
      
deriving instance Read Seperator
instance ParseField Seperator
      
deriving instance Read RowType
instance ParseField RowType
data Options = Options {   filter :: [RowType] <?> "What show -- Insert, etc.."
                         , sep :: First Seperator <?> "Comma or Tab output"
                         , fasta :: FilePath <?> "Input aligned fasta file"
                         , align :: Bool <?> "Align the sequence first?"}
             deriving (Generic, Show)

instance ParseRecord Options

main :: IO ()
main = do
  x <- getRecord "Running Program"
  print (x :: Options)

  
--data GenBankSource = GBFile FilePath | CSVFile FilePath | GBID String
--  deriving (Show, Generic, ParseRecord)

--instance ParseField GenBankSource where
--      parseField m = fmap GenBankSource (parseField m)
      
  --process (x :: Options)

--instance Read Seperator
--instance ParseField Seperator where
--     --parseField m s = fmap (\x -> if (x == "tab") then Tab else Comma) (parseField m s)
--     parseField (Just m) (Just s) = NilP $ Just (if (m == "tab") then Tab else Comma)
