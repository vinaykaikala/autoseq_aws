import unittest

from autoseq.util.vcfutils import vt_split_and_leftaln, fix_ambiguous_cl, remove_dup_cl


class TestVcfUtils(unittest.TestCase):
    def test_vt_split_and_leftaln(self):
        """
        test that -n parameter is passed when allow_ref_mismatches is True
        and vice versa when it's False
        """
        reference_sequence = "ref.fasta"
        cmd = vt_split_and_leftaln(reference_sequence, allow_ref_mismatches=True)
        self.assertIn(' -n ', cmd)

        cmd = vt_split_and_leftaln(reference_sequence, allow_ref_mismatches=False)
        self.assertNotIn(' -n ', cmd)

    def test_fix_ambiguous_cl(self):
        """
        test that the passed integer get's put in the right place and that the complete string is correct
        """
        cmd = fix_ambiguous_cl(10)
        self.assertEqual(cmd,
                         r"""awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $10) } {print}'""")

    def test_remove_dup_cl(self):
        cmd = remove_dup_cl()
        self.assertEqual(cmd,
                         r""" awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}'""")