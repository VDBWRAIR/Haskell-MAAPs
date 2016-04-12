import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck
import Lib (expand, AA(..), Codon)

main :: IO ()
main = defaultMain tests

tests = [testGroup "Sorting Group 1" [
                testProperty "example ATR" prop_example_expand
           ]

        ]

--prop_frame_shift_reported 
--prop_N_reported
--prop_insert_reported
--prop_synonymous_is_synonymous
--prop_nonsynonymous_is_nonsynonymous
prop1 b = b == False
  where types = (b :: Bool)

       
prop_example_expand = expand "ATR" == Just [(Codon "ATR", [I,M])]
prop_lower_case_equivalent_upper x = expand x == expand $ toUpper x
  where types = (x :: String)
prop_bad_nts_fail = expand "zzz" == Nothing
bad_char_count_fails x = expand x == Nothing
  where types = (x :: String) suchThat (((length x) `mod` 3) /= 0)
