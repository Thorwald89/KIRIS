"""
Unit Tests per KIRIS PRO - Motore di Validazione e Calcolo Deterministico (ISO 13485)
"""

import sys
import os
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.kir_analysis import (
    validate_and_extract_alleles,
    map_hla_to_ligands,
    calculate_b_content_logic,
    assess_donor_education,
    calculate_vectors_quantitative,
    calculate_donor_recommendation,
    generate_radar_chart_svg,
    analyze_single_pair
)

class TestKirisImmunogenetics(unittest.TestCase):

    def test_allele_validation_and_normalization(self):
        # Test formato standard
        self.assertEqual(validate_and_extract_alleles("A*03:01, A*11:01", "HLA-A"), ["A*03:01", "A*11:01"])
        # Test normalizzazione automatica (manca prefisso locus)
        self.assertEqual(validate_and_extract_alleles("03:01, 11:01", "HLA-A"), ["A*03:01", "A*11:01"])
        # Test normalizzazione automatica (manca asterisco)
        self.assertEqual(validate_and_extract_alleles("B08:01, B44:02", "HLA-B"), ["B*08:01", "B*44:02"])
        # Test errore di locus
        with self.assertRaises(ValueError):
            validate_and_extract_alleles("B*08:01", "HLA-A")

    def test_c16_and_c12_epitopes(self):
        # C*16:01 e C*16:04 devono essere C1; C*16:02 deve essere C2
        ligs_c1 = map_hla_to_ligands(["C*16:01", "C*16:04"], [], [])
        self.assertIn("C1", ligs_c1)
        self.assertNotIn("C2", ligs_c1)

        ligs_c2 = map_hla_to_ligands(["C*16:02"], [], [])
        self.assertIn("C2", ligs_c2)

        # C*12:02 C1 vs C*12:04 C2
        ligs_1204 = map_hla_to_ligands(["C*12:04"], [], [])
        self.assertIn("C2", ligs_1204)

    def test_bw4_residue_classification(self):
        # B*13:01 -> L80 (Anergica)
        ligs_b13 = map_hla_to_ligands([], ["B*13:01"], [])
        self.assertIn("Bw4-L80", ligs_b13)

        # B*51:01 -> I80 (Alta affinità)
        ligs_b51 = map_hla_to_ligands([], ["B*51:01"], [])
        self.assertIn("Bw4-I80", ligs_b51)

        # B*44:02 -> T80 (Media affinità)
        ligs_b44 = map_hla_to_ligands([], ["B*44:02"], [])
        self.assertIn("Bw4-T80", ligs_b44)

        # B*08:01 -> Bw6
        ligs_b08 = map_hla_to_ligands([], ["B*08:01"], [])
        self.assertIn("Bw6", ligs_b08)

    def test_kir2ds1_education_stratification(self):
        """
        Verifica conformità biologica per KIR2DS1:
        - Donatore C1/C1 -> Hyperexcitable
        - Donatore C1/C2 -> Responsivo
        - Donatore C2/C2 -> Anergico/Tollerizzato
        """
        kir_set = {"KIR2DS1", "KIR3DL3", "KIR2DL4", "KIR3DL2", "KIR3DP1", "KIR2DP1"}

        # Caso 1: C1/C1
        d_ligs_c1 = {"C1"}
        edu_c1 = assess_donor_education(kir_set, d_ligs_c1, ["B*08:01"], ["C*07:01", "C*08:01"])
        row_c1 = next(r for r in edu_c1 if r['KIR'] == 'KIR2DS1')
        self.assertIn("Hyperexcitable", row_c1['Edu'])
        self.assertIn("C1/C1", row_c1['Lig'])

        # Caso 2: C1/C2
        d_ligs_c1c2 = {"C1", "C2"}
        edu_c1c2 = assess_donor_education(kir_set, d_ligs_c1c2, ["B*08:01"], ["C*07:01", "C*02:02"])
        row_c1c2 = next(r for r in edu_c1c2 if r['KIR'] == 'KIR2DS1')
        self.assertIn("Responsivo", row_c1c2['Edu'])
        self.assertIn("C1/C2", row_c1c2['Lig'])

        # Caso 3: C2/C2
        d_ligs_c2 = {"C2"}
        edu_c2 = assess_donor_education(kir_set, d_ligs_c2, ["B*08:01"], ["C*04:01", "C*02:02"])
        row_c2 = next(r for r in edu_c2 if r['KIR'] == 'KIR2DS1')
        self.assertIn("Anergico", row_c2['Edu'])
        self.assertIn("C2/C2", row_c2['Lig'])

    def test_kir2ds1_gvl_activation_on_recipient_target(self):
        """
        Verifica che KIR2DS1 dia boost GvL SOLO se il ricevente esprime C2.
        """
        kir_set = {"KIR2DS1", "KIR3DL3", "KIR2DL4", "KIR3DL2", "KIR3DP1", "KIR2DP1"}
        don_b = ["B*08:01", "B*08:01"]
        don_c = ["C*07:01", "C*08:01"]  # Donatore C1/C1 -> Hyperexcitable
        d_ligs = {"C1", "Bw6"}
        edu_res = assess_donor_education(kir_set, d_ligs, don_b, don_c)
        k2ds1_edu = next(r['Edu'] for r in edu_res if r['KIR'] == 'KIR2DS1')

        # Ricevente 1: C2-positivo (es. C*02:01)
        r_ligs_c2_pos = {"C2", "Bw6"}
        _, _, _, _, _, _, gvl_score_c2, gvl_status_c2, _ = calculate_vectors_quantitative(
            edu_res, d_ligs, r_ligs_c2_pos, "B/X", 1, k2ds1_edu
        )
        self.assertIn("KIR2DS1 Ipereccitabile vs C2", gvl_status_c2)
        self.assertGreaterEqual(gvl_score_c2, 3)

        # Ricevente 2: C2-negativo (C1/C1, es. C*07:01, C*08:01)
        r_ligs_c2_neg = {"C1", "Bw6"}
        _, _, _, _, _, _, gvl_score_c1, gvl_status_c1, _ = calculate_vectors_quantitative(
            edu_res, d_ligs, r_ligs_c2_neg, "B/X", 1, k2ds1_edu
        )
        self.assertNotIn("KIR2DS1 Ipereccitabile vs C2 (+1.5)", gvl_status_c1)
        self.assertIn("Ricevente C2-negativo", gvl_status_c1)

    def test_b_content_and_cooley_score(self):
        # Aplotipo A/A
        haplo_a, score_a, cen_a, tel_a = calculate_b_content_logic({
            "KIR3DL3", "KIR2DL3", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR2DS4", "KIR3DL2"
        })
        self.assertEqual(haplo_a, "A/A")
        self.assertEqual(score_a, 0)

        # Aplotipo B/X con marker Centromerico e Telomerico (es. 2DS2 e 2DS1)
        haplo_b, score_b, _, _ = calculate_b_content_logic({
            "KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DP1", "KIR2DL1", "KIR3DP1", "KIR2DL4", "KIR2DS1", "KIR3DL2"
        })
        self.assertEqual(haplo_b, "B/X")
        self.assertGreaterEqual(score_b, 2)

    def test_donor_recommendation_and_radar_svg(self):
        # Test Donatore Ottimale
        rec_opt = calculate_donor_recommendation(gvh_score=4, hvg_score=1, gvl_score=5, haplo="B/X", b_score=3)
        self.assertEqual(rec_opt["level"], "PRIORITÀ ELEVATA")

        # Test Donatore Subottimale (HvG alto)
        rec_sub = calculate_donor_recommendation(gvh_score=1, hvg_score=4, gvl_score=2, haplo="A/A", b_score=0)
        self.assertEqual(rec_sub["level"], "SUBOTTIMALE")

        # Test generazione SVG
        svg = generate_radar_chart_svg(3, 1, 4)
        self.assertIn("<svg", svg)
        self.assertIn("polygon", svg)
        self.assertIn("GvL", svg)

if __name__ == '__main__':
    unittest.main()
