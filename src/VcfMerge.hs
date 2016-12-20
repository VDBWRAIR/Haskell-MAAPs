{-# LANGUAGE OverloadedStrings #-}

import Safe (atMay, headMay, lastMay, readMay)
import Data.Maybe (fromMaybe)
import qualified Data.Text as T
import Data.Text (Text())
import Data.List (find, intercalate)
import Control.Monad (when)

pos :: Text -> Int
pos s = fromMaybe (-1) (toInt =<< (nth 1 $ T.splitOn "\t" s))

toInt x = readMay $ T.unpack x :: Maybe Int
nth = flip atMay
readNuc x = readMay $ T.unpack x :: Maybe Nucleotide
data Nucleotide = A | T | C | G
  deriving (Show, Eq, Read)

data Stats = Stats { dp :: Int, alts :: [Nucleotide], pacs :: [Int] }
  deriving (Show, Eq)
--
formatStats Stats {dp=dp', alts=alts', pacs=pacs'} = s where
  s = (show dp') ++ "\t" ++ (intercalate "," $ zipWith (\x y -> (show x) ++ "==" ++ (show y)) alts' pacs')
--
--
s = "Den4/KDH0146A/Thailand/Unk/Den4_1\t6\t.\tA\tG,C,T\t.\t.\tDP=2998;RC=2973;RAQ=37;PRC=99;AC=1,10,14;AAQ=37,37,37;PAC=3,2,1;CBD=2973;CB=A;HPOLY"
---- get DP; then depending on number of alts, expect and get that number of PAC (note PAC may be zero)
-- | hey
-- >>> parseInfo s
-- Just (Stats {dp = 2998, alts = [G,C,T], pacs = [3,2,1]}) 
parseInfo :: Text -> Maybe Stats
parseInfo line = do
  info <- lastMay $ T.splitOn "\t" line 
  let fields = T.splitOn ";" info
  dp'  <- toInt =<< fieldValue "DP" fields
  rawAlts' <- T.splitOn "," <$> nth 4 cols
  alts' <- sequence $ map readNuc rawAlts' -- note mapM = sequence . map
  return alts'
  pac' <-  fieldValue "PAC" fields
  pacs' <- sequence $ map toInt $ T.splitOn "," pac' 
  if (length alts' == length pacs') then
    return $ Stats { dp = dp', alts = alts', pacs = pacs'}
    else Nothing
  where 
    cols   = T.splitOn "\t" line
    fieldValue pre xs = nth 1 =<< T.splitOn "=" <$> find (T.isPrefixOf pre) xs

-- output as TSV, then join with MAAPs output on CHROM and POS

main = putStrLn $ show $ parseInfo s
