#!/usr/bin/env python3

"""tests.test_compare.py: tests for exfi.compare"""


from unittest import main, TestCase

from exfi.compare import \
    bedtools_intersect, \
    classify, \
    compute_stats_per_exon, \
    compute_true_bases, \
    compute_pred_bases, \
    compute_true_positive_bases, \
    compute_false_positive_bases, \
    compute_false_negative_bases, \
    compute_stats_per_base


from tests.compare import \
    BED3_EMPTY_FN, BED3_TRUE_FN, BED3_PRED_FN, \
    BED3_EMPTY, BED3_TRUE, BED3_PRED, \
    TP_DF, FP_DF, FN_DF, \
    CLASSIFICATION, \
    STATS_PER_EXON, STATS_PER_BASE



class TestBedtoolsIntersect(TestCase):
    """Tests for exfi.compare.bedtools_intersect"""

    def test_empty(self):
        """exfi.compare.bedtools_intersect: empty case"""
        observed = bedtools_intersect(BED3_EMPTY_FN, BED3_EMPTY_FN, [])
        print('Observed:', observed.values, 'Expected:', BED3_EMPTY.values, sep='\n')
        self.assertEqual(observed.values.tolist(), BED3_EMPTY.values.tolist())

    def test_tp(self):
        """exfi.compare.bedtools_intersect: true positive case"""
        observed = bedtools_intersect(
            bed1_fn=BED3_PRED_FN,
            bed2_fn=BED3_TRUE_FN,
            additional_flags=['-f', '0.99', '-r', '-wo']
        )
        print('Observed:', observed.values.tolist(), 'Expected:', TP_DF, sep='\n')
        self.assertTrue(observed.equals(TP_DF))

    def test_fp(self):
        """exfi.compare.bedtools_intersect: false positive case"""
        observed = bedtools_intersect(
            bed1_fn=BED3_PRED_FN,
            bed2_fn=BED3_TRUE_FN,
            additional_flags=['-f', '0.99', '-r', '-v']
        )
        print('Observed:', observed, 'Expected:', FP_DF, sep='\n')
        self.assertTrue(observed.equals(FP_DF))

    def test_fn(self):
        """exfi.compare.bedtools_intersect: false negative case"""
        observed = bedtools_intersect(
            bed1_fn=BED3_TRUE_FN,
            bed2_fn=BED3_PRED_FN,
            additional_flags=['-f', '0.99', '-r', '-v']
        )
        print('Observed:', observed, 'Expected:', FN_DF, sep='\n')
        self.assertTrue(observed.equals(FN_DF))



class TestClassify(TestCase):
    """Tests for exfi.compare.classify"""

    def test_complex(self):
        """exfi.compare.classify: Test some exons"""
        observed = classify(
            bed3_true=BED3_TRUE, bed3_pred=BED3_PRED, fraction=0.99
        )
        print(
            'Observed:', observed['true_positives'],
            'Expected:', CLASSIFICATION['true_positives'],
            sep='\n'
        )
        self.assertTrue(
            observed['true_positives'].equals(CLASSIFICATION['true_positives'])
        )
        print(
            'Observed:', observed['false_positives'],
            'Expected:', CLASSIFICATION['false_positives'],
            sep='\n'
        )
        self.assertTrue(
            observed['false_positives'].equals(CLASSIFICATION['false_positives'])
        )
        print(
            'Observed:', observed['false_negatives'],
            'Expected:', CLASSIFICATION['false_negatives'],
            sep='\n'
        )
        self.assertTrue(
            observed['false_negatives'].equals(CLASSIFICATION['false_negatives'])
        )



class TestComputeStatsPerExon(TestCase):
    """Tests for exfi.compare.compute_stats_per_exon"""

    def test_complex(self):
        """exfi.compare.compute_stats_per_exon: some exons"""
        observed = compute_stats_per_exon(CLASSIFICATION)
        print(observed.values.tolist())
        print(STATS_PER_EXON.values.tolist())
        self.assertTrue(observed.equals(STATS_PER_EXON))



class TestComputeTrueBases(TestCase):
    """Tests for exfi.compare.compute_true_bases"""

    def test_complex(self):
        """exfi.compare.compute_true_bases: some exons"""
        observed = compute_true_bases(CLASSIFICATION)
        print(observed)
        self.assertEqual(observed, 5080.0)



class TestComputePredBases(TestCase):
    """Tests for exfi.compare.compute_pred_bases"""

    def test_complex(self):
        """exfi.compare.compute_pred_bases: some exons"""
        observed = compute_pred_bases(CLASSIFICATION)
        print(observed)
        self.assertEqual(observed, 5090.0)



class TestComputeTruePositiveBases(TestCase):
    """Tests for exfi.compare.compute_true_positive_bases"""

    def test_complex(self):
        """exfi.compare.compute_true_positive_bases: some exons"""
        observed = compute_true_positive_bases(CLASSIFICATION)
        self.assertEqual(observed, 4159.0)



class TestComputeFalsePositiveBases(TestCase):
    """Tests for exfi.compare.compute_false_positive_bases"""

    def test_complex(self):
        """exfi.compare.compute_false_positive_bases: some exons"""
        observed = compute_false_positive_bases(CLASSIFICATION)
        self.assertEqual(observed, 931)



class TestComputeFalseNegativeBases(TestCase):
    """Tests for exfi.compare.compute_false_negative_bases"""

    def test_complex(self):
        """exfi.compare.compute_false_negative_bases: some exons"""
        observed = compute_false_negative_bases(CLASSIFICATION)
        self.assertEqual(observed, 911)


class TestComputeStatsPerBase(TestCase):
    """Tests for exfi.compare.compute_stats_per_base"""

    def test_complex(self):
        """exfi.compare.compute_stats_per_base: some exons"""
        observed = compute_stats_per_base(CLASSIFICATION)
        print(observed.values.tolist())
        print(STATS_PER_BASE.values.tolist())
        self.assertTrue(observed.equals(STATS_PER_BASE))



if __name__ == '__main__':
    main()
