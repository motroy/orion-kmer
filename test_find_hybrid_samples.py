import unittest
from find_hybrid_samples import classify_platform

class TestFindHybridSamples(unittest.TestCase):
    def test_classify_platform(self):
        self.assertEqual(classify_platform('Illumina MiSeq'), 'SHORT')
        self.assertEqual(classify_platform('MinION'), 'LONG')
        self.assertEqual(classify_platform('GridION'), 'LONG')
        self.assertEqual(classify_platform('PacBio RS II'), 'LONG')
        self.assertEqual(classify_platform('NextSeq 500'), 'SHORT')
        self.assertEqual(classify_platform('DNBSEQ-T7'), 'SHORT')
        self.assertEqual(classify_platform('Ion Torrent PGM'), 'SHORT')
        self.assertEqual(classify_platform('Unknown'), 'OTHER')
        self.assertEqual(classify_platform(None), 'OTHER')
        self.assertEqual(classify_platform(123), 'OTHER')

if __name__ == '__main__':
    unittest.main()
