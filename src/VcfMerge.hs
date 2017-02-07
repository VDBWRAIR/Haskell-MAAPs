{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE DataKinds         #-}  -- last three for Options
{-# LANGUAGE TypeOperators     #-}
{-# LANGUAGE DeriveGeneric #-}

import Safe (atMay, headMay, lastMay, readMay, tailMay)
import Data.Maybe (catMaybes, fromMaybe)
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.Text (Text())
import Data.List (find, groupBy, group)
import Control.Monad (when)
import Options.Generic

--pos :: Text -> Int
--pos s = fromMaybe (-1) (toInt =<< (nth 1 $ T.splitOn "\t" s))
data Options = Options {   maaps :: FilePath <?> "Tab delimited MAAPs output"
                         , vcf :: FilePath <?> "VCF ouptut of ngsmapper." }
  deriving (Generic, Show)
instance ParseRecord Options

toInt x = readMay $ T.unpack x :: Maybe Int
nth = flip atMay
readNuc x = readMay $ T.unpack x :: Maybe Nucleotide
data Nucleotide = A | T | C | G | (*)
  deriving (Show, Eq, Read)


data Stats = Stats { dp :: Int, alts :: [Nucleotide], pacs :: [Int], pos :: Int, ref :: Text }
  deriving (Show, Eq)
--
formatStats Stats {dp=dp', alts=alts', pacs=pacs', pos=_, ref=ref'} = s where
  s = T.intercalate "\t" [(tshow dp'), showAlts alts' pacs', ref'] 
  showAlts :: [Nucleotide] -> [Int] -> Text
  showAlts [] [] = "."
  showAlts xs ys = T.intercalate "," $ zipWith (\x y -> (tshow x) `T.append` "==" `T.append` (tshow y)) xs ys
 
tshow :: (Show a) => a -> Text
tshow = T.pack . show
--
--T
s = "Den4/KDH0146A/Thailand/Unk/Den4_1\t6\t.\tA\tG,C,T\t.\t.\tDP=2998;RC=2973;RAQ=37;PRC=99;AC=1,10,14;AAQ=37,37,37;PAC=3,2,1;CBD=2973;CB=A;HPOLY"

byRef xs = groupBy (\x y -> (headMay x) == (headMay y)) xs

readMaaps = filter (\x -> length x > 1) . map T.words . tail . T.lines 
-- note that ngs_mapper vcfs have an entry for every position, even if there are no ALTs there. so suck!
--readVCF = catMaybes . map parseInfo . dropWhile (T.isPrefixOf "#") . T.lines
readVCF :: Text -> [Stats]
readVCF x = fromMaybe (error "Failed to parse vcf file.") $ traverse parseInfo $ dropWhile (T.isPrefixOf "#") $ T.lines x
main = do
  opts <- (getRecord "Starting maaps" :: IO Options)
  let vcf' = unHelpful $ vcf opts
  let tsv  = unHelpful $ maaps opts
  stats <- readVCF   <$> T.readFile vcf'
  maaps <- readMaaps <$> T.readFile tsv
  maapsHeader <- head <$> T.lines <$> T.readFile tsv
  let header = maapsHeader `T.append` "\tDP\tAlts\tRef"
  putStrLn $ T.unpack header
  putStrLn $ process stats maaps

process ms ss = T.unpack $ T.unlines $ map (uncurry output) $ joinMaapsStats ms ss


--output :: Maybe ([Text], Stats) -> String
output :: [Text] -> Maybe Stats -> Text
output mr (Just stats) = (T.intercalate "\t" mr) `T.append` "\t" `T.append` (formatStats stats)
output mr Nothing      = T.intercalate "\t" (mr ++ ["-", "-"])

joinMaapsStats :: [Stats] -> [[Text]] -> [([Text], Maybe Stats)]
joinMaapsStats vcf maaps = map (\m -> (m, find (matcher m) vcf)) maaps
  where matcher = if (length $ group $ fmap ref vcf) == 1 then matchPos else matchBoth
--  f m = do
--    r <- find (match m) tsvs
--    return (r, m)
--matchPos x y = (toInt =<< nth 2 x) == (Just $ pos y)
matchPos :: [Text] -> Stats -> Bool
matchPos x y = (read $ T.unpack  (x !! 2)) == (pos y)
matchBoth :: [Text] -> Stats -> Bool
matchBoth x y = ((ref'  x) == (ref y)) && matchPos x y
--ref' x = T.tail <$> snd <$> T.breakOn "_" <$> headMay x
ref' x = T.tail $ snd $ T.breakOn "_" $ head x

   
--groupBy (\x y -> ( (pos y)) 
---- get DP; then depending on number of alts, expect and get that number of PAC (note PAC may be zero)
-- | hey
-- >>> parseInfo s
-- Just (Stats {dp = 2998, alts = [G,C,T], pacs = [3,2,1]}) 
parseInfo :: Text -> Maybe Stats
parseInfo line = do
  info <- lastMay cols
  let fields = T.splitOn ";" info
  dp'  <- toInt =<< fieldValue "DP" fields
  rawAlts' <- T.splitOn "," <$> nth 4 cols
  let altsMaybe' = mapM readNuc rawAlts' -- note mapM = sequence . map
  let pacMaybe' = mapM toInt =<< T.splitOn "," <$> fieldValue "PAC" fields
  let pacs' = fromMaybe [] pacMaybe'  :: [Int]
  let alts' = fromMaybe [] altsMaybe'
  --pacs' <- sequence $ map toInt $ T.splitOn "," pac' 
  pos'  <- toInt =<< nth 1 cols
  ref'  <- nth 0 $ T.splitOn "\t" line
  if (length alts' == length pacs') then
    return $ Stats { dp = dp', alts = alts', pacs = pacs', pos = pos', ref = ref'}
    else Nothing
  where 
    cols   = T.splitOn "\t" line
    fieldValue pre xs = nth 1 =<< T.splitOn "=" <$> find (T.isPrefixOf pre) xs

-- output as TSV, then join with MAAPs output on CHROM and POS

--main = putStrLn $ show $ parseInfo s
