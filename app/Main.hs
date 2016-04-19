{-# LANGUAGE OverloadedStrings #-}
module Main where 
import Lib  (codonTable, process, run)
import Types (Options, RowType(..))
import Options.Generic (getRecord)

main :: IO ()
main = do
  x <- getRecord "Running Program"
  run   (x :: Options)

