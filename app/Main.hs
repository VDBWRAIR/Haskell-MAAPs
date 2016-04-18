{-# LANGUAGE OverloadedStrings #-}
module Main where 
import Lib  (process, run)
import Types (Options, RowType(..))
import Options.Generic (getRecord)

main :: IO ()
main = do
  x <- getRecord "Running Program"
  print (x :: Options)
  run   (x :: Options)

  
--data GenBankSource = GBFile FilePath | CSVFile FilePath | GBID String
--  deriving (Show, Generic, ParseRecord)

--instance ParseField GenBankSource where
--      parseField m = fmap GenBankSource (parseField m)
      
  --process (x :: Options)

--instance Read Seperator
--instance ParseField Seperator where
--     --parseField m s = fmap (\x -> if (x == "tab") then Tab else Comma) (parseField m s)
--     parseField (Just m) (Just s) = NilP $ Just (if (m == "tab") then Tab else Comma)
