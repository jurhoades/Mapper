# Imports
import os
import sys
import unittest

#os.chdir("..")
from mapper import Mapper

class TestMapperMethod(unittest.TestCase):

    def setUp(self):
        self.mapper = Mapper("data/refFlatMm10.txt")

    
    def test_normal_case(self):
        cds_pos, aa_pos = self.mapper.map(101153495, "NM_146145")
        self.assertEqual(cds_pos, 3462)
        self.assertEqual(aa_pos, 1154)


    def test_coord_in_intron(self):
        cds_pos, aa_pos = self.mapper.map(101153494, "NM_146145")
        self.assertEqual(cds_pos, None)
        self.assertEqual(aa_pos, None)


    def test_refseq_id_not_in_file(self):
        cds_pos, aa_pos = self.mapper.map(101153495, "NM_899287")
        self.assertEqual(cds_pos, None)
        self.assertEqual(aa_pos, None)



if __name__ == "__main__":
    unittest.main()
