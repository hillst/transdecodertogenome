#!/usr/bin/env python
import unittest
import string
from transdecoder_to_genome import Exon
from transdecoder_to_genome import Transcript
from transdecoder_to_genome import build_coordinates

class GFFLineTest(unittest.TestCase):
    
    def setUp(self):
        ret = unittest.TestCase.setUp(self)
        exon1 ="TCONS_00001425  .       gene    1439    1564    .       -       .       ID=g.11422;Name=ORF%20g.11422%20m.11422%20type%3A5prime_partial%20len%3A42%20%28-%29"
        exon2 ="TCONS_00001425  .       gene    1043    1231    .       +       .       ID=g.11419;Name=ORF%20g.11419%20m.11419%20type%3Acomplete%20len%3A63%20%28%2B%29"
        exon3 ="TCONS_00001425  .       gene    1       1062    .       +       .       ID=g.11414;Name=ORF%20g.11414%20m.11414%20type%3A5prime_partial%20len%3A354%20%28%2B%29"
        self.actual_exon = Exon("Spolyrrhiza9509S001     Cufflinks       exon    26722   26764   .       -       .       Parent=TCONS_00001425")
        self.actual_trans = Transcript(exon1)
        actual_exon1 = "Spolyrrhiza9509S001     Cufflinks       exon    24497   24656   .       -       .       Parent=TCONS_00001425"
        actual_exon2 = "Spolyrrhiza9509S001     Cufflinks       exon    24809   24938   .       -       .       Parent=TCONS_00001425"
        actual_exon3 = "Spolyrrhiza9509S001     Cufflinks       exon    26237   26599   .       -       .       Parent=TCONS_00001425"
        self.actual_exon1 = Exon(actual_exon1)
        self.actual_exon2 = Exon(actual_exon2)
        self.actual_exon3 = Exon(actual_exon3)
        self.exon1 = Transcript(exon1)
        self.exon2 = Transcript(exon2)
        self.exon3 = Transcript(exon3)
        return ret
    """
    Tests comparison on RNAseq format
    """
    def testTranscriptComparison(self):
        self.assertGreater(self.exon2, self.exon1)
        self.assertGreater(self.exon3, self.exon1)
        self.assertGreater(self.exon2, self.exon3)
        self.assertLess(self.exon1, self.exon2)
        self.assertLess(self.exon1, self.exon3)
        self.assertLess(self.exon3, self.exon2)
        self.assertEqual(self.exon1, self.exon1)
        
    def testBaseGetDicFormat(self):
        actual = self.actual_exon.getDicFormat()
        expected = {"name": "TCONS_00001425", "start": 26722, "stop": 26764}
        self.assertEqual(actual, expected)
        
    def testBaseConstructor(self):
        self.assertEqual(self.exon1.start, 1439)
        self.assertEqual(self.exon1.stop, 1564)
        self.assertEqual(self.exon1.getLabel(), "gene")
        self.assertEqual(self.exon1.parent, "TCONS_00001425")
        self.assertEqual(len(self.exon1), 125)
        self.assertEqual(str(self.exon1), "TCONS_00001425  .       gene    1439    1564    .       -       .       ID=g.11422;Name=ORF%20g.11422%20m.11422%20type%3A5prime_partial%20len%3A42%20%28-%29")
        self.assertEqual(self.exon1.partial, "partial")
        self.assertEqual(self.exon1.strand, "-")
        
    def testBadPositions(self):
        tester = "TCONS_00001425  .       gene    2000    1564    .       -       .       ID=g.11422;Name=ORF%20g.11422%20m.11422%20type%3A5prime_partial%20len%3A42%20%28-%29"
        try:
            Exon(tester)
            self.fail("No exception thrown")
        except:
            pass
        tester = "TCONS_00001425  .       gene    1564    1564    .       -       .       ID=g.11422;Name=ORF%20g.11422%20m.11422%20type%3A5prime_partial%20len%3A42%20%28-%29"

        try:
            Exon(tester)
            self.fail("No exception thrown")
        except:
            pass
    def testSetLenth(self):
        self.exon1.setLength(2222)
        self.assertEqual(len(self.exon1), 2222)
    
    def testSwitchStrand(self):
        self.exon1.switchStrand()
        self.assertEquals(self.exon1.strand, "+")
        
    def testGetExonDicFormat(self):
        temp = self.actual_trans
        temp.addExon(self.actual_exon1)
        temp.addExon(self.actual_exon2)
        temp.addExon(self.actual_exon3)
        actual = temp.getExonDicList()
        expected = [{"start": 24497, "stop":24656 ,"name": "TCONS_00001425"}, {"start":24809 , "stop": 24938 ,"name": "TCONS_00001425"}, {"start": 26237, "stop": 26599 ,"name": "TCONS_00001425"}]
        self.assertEqual(actual, expected)

"""
    The following tests show that the sign showed in the gff file are irrelevant, the coordinates
    remain the same regardless.
    
    That is, if a negative sign is shown in the transdecoder.gff3 file, the output coordinates work, and to get the data you must
    slice the subsquence from exons.fasta, then reverse transcribe it.
    
    This will output the sequence in the transdecoder.cds file.
"""        
class TransDecoderStrandednessTest(unittest.TestCase):

    def setUp(self):
        ret = unittest.TestCase.setUp(self)
        self.plusplusorf, self.plusminusorf, self.minusminusorf, self.minusplusorf = "", "", "", ""
        self.plusplusgen, self.plusminusgen, self.minusminusgen, self.minusplusgen = "", "", "", ""
        self.trans = string.maketrans("ATCGR", "TAGCY")
        with open("test/plusplus.gm", "r") as gm:
            for line in gm: self.plusplusgen += line.strip()
        with open("test/plusplus.orf", "r") as orf:
            for line in orf: self.plusplusorf += line.strip()
        with open("test/plusminus.gm", "r") as gm:
            for line in gm: self.plusminusgen += line.strip()
        with open("test/plusminus.orf", "r") as orf:
            for line in orf: self.plusminusorf += line.strip() 
        with open("test/minusplus.gm", "r") as gm:
            for line in gm: self.minusplusgen += line.strip()
        with open("test/minusplus.orf", "r") as orf:
            for line in orf: self.minusplusorf += line.strip() 
        with open("test/minusminus.gm", "r") as gm:
            for line in gm: self.minusminusgen += line.strip()
        with open("test/minusminus.orf", "r") as orf:
            for line in orf: self.minusminusorf += line.strip()         


        return ret
    def testCorrectPlusPlus(self):
        """
            Plus Plus works exactly as expected. No reverse transcription required. ATG appears immediately
            at final[522:3012], and the sequence in the transdecoder.cds match the subsequence.
        """
        gene = {'start': 74641, 'stop': 78761, 'name':'plusplus', 
                'exons': ({ 'start': 74641, 'stop': 76014 } , 
                          { 'start': 76111, 'stop': 76277 } , 
                          { 'start': 76368, 'stop': 76465 } ,
                          { 'start': 76558, 'stop': 76914 } ,
                          { 'start': 77029, 'stop': 77445 } ,
                          { 'start': 77628, 'stop': 78761 }
                         )}
        transcript = {'name':'plusplus', 'start': 523, 'stop':3012}
        cds = build_coordinates(gene, transcript) #??? what should this do... 
        # should equal start-1:stop
        self.assertEqual(self.plusplusgen[522:3012], self.plusplusorf)
    def testPlusMinus(self):
        # + - CASE FOR PARTIAL
        """
            final[878:1190][::-1].translate(transtab) where final is the exons.fasta for something that was originally +.
            and the gff3 file said it is negative works correctly. that is, our coordinates provided by the gff3 file are
            correct for + in org - in the new. This means to get the sequence it needs to be reverse transcribed, but the
            coordinate themselves are correct.
            
            TCONS_00000005
        """
        gene = {'start': 67878, 'stop': 74308, 'name':'plusmin', 
                'exons': ({ 'start': 67878, 'stop': 68085 } , 
                          { 'start': 68882, 'stop': 68943 } , 
                          { 'start': 69034, 'stop': 69141 } ,
                          { 'start': 69683, 'stop': 69762 } ,
                          { 'start': 71555, 'stop': 71717 } ,
                          { 'start': 72187, 'stop': 72261 } ,
                          { 'start': 72370, 'stop': 72458 } ,
                          { 'start': 73151, 'stop': 73270 } ,
                          { 'start': 73918, 'stop': 73962 } ,
                          { 'start': 74066, 'stop': 74308 } 
                         )}
        transcript = {'name':'plusmin', 'start': 879, 'stop': 1190}
        #reverse
        self.assertEqual(self.plusminusgen[878:1190][::-1].translate(self.trans), self.plusminusorf)
    
    def testMinusPlus(self):
        # - + case
        """
            Minus plus immediately finds the ATG at the start of the genome sequence, this is as expected,
            as the gff3 provided by transdecoder said +, and the genome sequence in exons.fasta and the 
            transdecoder coding sequence agree.
            
            found at final[0:957]
            TCONS_00001426
        """
        gene = {'start': 26167, 'stop': 27852, 'name':'minplus', 
                'exons': ({ 'start': 26167, 'stop': 26599 } , 
                          { 'start': 26722, 'stop': 26764 } , 
                          { 'start': 27215, 'stop': 27494 } ,
                          { 'start': 27650, 'stop': 27852 } 
                         )}
         
        transcript = {'name':'minplus', 'start': 1, 'stop': 957}
        
        #forward
        self.assertEqual(self.minusplusgen[0:957], self.minusplusorf)
    
    def testMinusMinus(self):
        # - - case
        """
            final[1003:1327][::-1].translate(trans) where final is the exons.fasta for somethinga that was originally -.
            The gff3 file from transdecoder also said -, which means to reverse transcribe the fasta file to get what was
            in the coding sequence. These are the coordinates found in the gff3, so we know they are correct.
            TCONS_00015061
        """
        gene = {'start': 5731958, 'stop': 5744344, 'name':'minmin', 
                'exons': ({ 'start': 5731958, 'stop': 5732148 } , 
                          { 'start': 5732445, 'stop': 5733573 } , 
                          { 'start': 5733673, 'stop': 5733929 } ,
                          { 'start': 5735453, 'stop': 5736066 } ,
                          { 'start': 5736111, 'stop': 5736553 } ,
                          { 'start': 5736642, 'stop': 5736953 } ,
                          { 'start': 5737039, 'stop': 5738345 } ,
                          { 'start': 5738598, 'stop': 5738961}  ,
                          { 'start': 5739248, 'stop': 5739641 } ,
                          { 'start': 5739724, 'stop': 5740829 } ,
                          { 'start': 5741250, 'stop': 5743397 } ,
                          { 'start': 5743651, 'stop': 5744344 } 
                         )}
        transcript = {'name':'minmin', 'start': 1004, 'stop':1327}
        #reverse
        self.assertEqual(self.minusminusgen[1003:1327][::-1].translate(self.trans),self.minusminusorf)
class BuildCoordinateTest(unittest.TestCase):
    
    """
        previous bugs found:
        If the transcript coordinate was longer than the physical sequence it would break, 
        
        off by 2 in 3prime partial and on the last exon.
        
        CDS would skil the 2nd to last place and put it at the end
        
        Transcript names were incorrect (off by one or so)
    """
    def setUp(self):
        ret = unittest.TestCase.setUp(self)
        #101 - 1000 length of 900
        self.gene = { "start": 101, "stop": 1000, "name": "off1end", "exons" : 
                                    [{ "start": 101, "stop": 105}, #1-5
                                     { "start": 111, "stop": 200}, #6-95
                                     { "start": 225, "stop": 275}, #96-146
                                     { "start": 301, "stop": 700}, #147-546
                                     { "start": 801, "stop": 900}, #547-646
                                     { "start": 999, "stop": 1000}]} #647
                                                                    
        return ret
    
    def testPerfectMatch(self):
        trans = {"start": 1, "stop": 647, "name":"off1end"}
        cds = build_coordinates(self.gene, trans)
        self.assertEqual(self.gene['exons'], cds, "Full gene does not match.") #should map perfectly
        
    def testTooLong(self):
        trans = {"start": 1, "stop": 648, "name":"off1end"}
        cds = build_coordinates(self.gene, trans)
        self.assertEqual(self.gene['exons'], cds, "Did not correct for an input that goes past maximum length") 

    def testTooShort(self):
        trans = {"start": 1, "stop": 646, "name":"off1end"}
        cds = build_coordinates(self.gene, trans)
        expect =   [{ "start": 101, "stop": 105},
                                     { "start": 111, "stop": 200},
                                     { "start": 225, "stop": 275},
                                     { "start": 301, "stop": 700},
                                     { "start": 801, "stop": 900}
                                     ]
        self.assertEqual(expect, cds, "included last region") #should map off by one
    
    def testMiddle(self):
        trans = {"start": 147, "stop": 550, "name":"off1end"}
        cds = build_coordinates(self.gene, trans)
        expect =   [{ "start": 301, "stop": 700},
                    { "start": 801, "stop": 804}
                    ]
        self.assertEqual(expect, cds, "included last region") 
    
    def testSingleExon(self):
        trans = {"start":6, "stop": 48, "name":"off1end"}
        cds = build_coordinates(self.gene, trans)
        expect = [{"start": 111, "stop":153}]
        self.assertEqual(expect, cds)
    
    def testMismatchNames(self):
        trans = {"start":1, "stop": 648, "name":"not a name"}
        try:
            cds = build_coordinates(self.gene,trans)
            self.fail("expected exception for name mismatch.")
        except:
            pass
    
    def testNegativeNumbers(self):
        trans = {"start": -200, "stop": 400, "name": "off1end"}
        try:
            cds = build_coordinates(self.gene,trans)
            self.fail("expected exception for negative number.")
        except:
            pass
        trans = {"start": 200, "stop": -400, "name": "off1end"}
        try:
            cds = build_coordinates(self.gene, trans)
            self.fail("expected exception for negative number")
        except:
            pass
        trans = {"start": -200, "stop": -100, "name": "off1end"}
        try:
            cds = build_coordinates(self.gene, trans)
            self.fail("expected exception for negative number")
        except:
            pass
                                
        
if __name__ == "__main__":
    unittest.main()
