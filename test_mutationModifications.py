import unittest
import mutationModifications 
import pandas as pd
class mutationModification(unittest.TestCase):

    def test_delecion(self):
        row = {"referenceAllele": "A", "variantAllele": "-"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Deleción")

    def test_insercion(self):
        row = {"referenceAllele": "-", "variantAllele": "A"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Inserción")

    def test_transicion(self):
        row = {"referenceAllele": "A", "variantAllele": "G"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Transición")
        row = {"referenceAllele": "C", "variantAllele": "T"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Transición")

    def test_transversion(self):
        row = {"referenceAllele": "A", "variantAllele": "T"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Transversión")
        row = {"referenceAllele": "C", "variantAllele": "G"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Transversión")

    def test_complex_transicion(self):
        row = {"referenceAllele": "AG", "variantAllele": "GA"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Transición")

    def test_complex_transversion(self):
        row = {"referenceAllele": "AG", "variantAllele": "CT"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Transversión")

    def test_general_delecion(self):
        row = {"referenceAllele": "AGT", "variantAllele": "A"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Deleción")

    def test_general_insercion(self):
        row = {"referenceAllele": "A", "variantAllele": "AGT"}
        self.assertEqual(mutationModifications.clasificar_mutacion(row), "Inserción")

    def test_generate_mutated_peptides_none_sequence(self):
        self.assertEqual(mutationModifications.generate_mutated_peptides(None, "A1B"), [])

    def test_generate_mutated_peptides_basic(self):
        sequence = "ABCDEFGHIABCDEFGHIABCDEFGHI"
        mutation = "I9B"
        expected = ["BABCDEFGH", "HBABCDEFG", "GHBABCDEF", "FGHBABCDE", "EFGHBABCD", "DEFGHBABC", "CDEFGHBAB", "BCDEFGHBA", "ABCDEFGHB"]
        self.assertEqual(mutationModifications.generate_mutated_peptides(sequence, mutation), expected)
    
    def test_classify_binding_strong(self):
        row = {"presentation_percentile": 0.3}
        self.assertEqual(mutationModifications.classify_binding(row), "SB")

    def test_classify_binding_weak(self):
        row = {"presentation_percentile": 1.5}
        self.assertEqual(mutationModifications.classify_binding(row), "WB")

    def test_classify_binding_na(self):
        row = {"presentation_percentile": 2.5}
        self.assertEqual(mutationModifications.classify_binding(row), "N/A")

    def test_classify_binding_custom_thresholds(self):
        row = {"presentation_percentile": 1.0}
        self.assertEqual(mutationModifications.classify_binding(row, presentation_percentile_hard=1.0, presentation_percentile_soft=3.0), "SB")
        row = {"presentation_percentile": 2.5}
        self.assertEqual(mutationModifications.classify_binding(row, presentation_percentile_hard=1.0, presentation_percentile_soft=3.0), "WB")
        row = {"presentation_percentile": 3.5}
        self.assertEqual(mutationModifications.classify_binding(row, presentation_percentile_hard=1.0, presentation_percentile_soft=3.0), "N/A")

    def test_calcularNeoantigenosPaciente_basic(self):
        clinical_df = pd.DataFrame({
            'Patient ID': ['P001', 'P002', 'P003']
        })
        predictions_df = pd.DataFrame({
            'patientId': ['P001', 'P001', 'P002', 'P003', 'P003', 'P003'],
            'Binding_Classification': ['SB', 'WB', 'SB', 'SB', 'WB', 'WB']
        })
        expected_df = pd.DataFrame({
            'Patient ID': ['P001', 'P002', 'P003'],
            'Neoantigen_SB_Count': [1, 1, 1],
            'Neoantigen_WB_Count': [1, 0, 2]
        })
        result_df = mutationModifications.calcularNeoantigenosPaciente(clinical_df, predictions_df)
        pd.testing.assert_frame_equal(result_df, expected_df)

    def test_calcularNeoantigenosPaciente_extra_patients(self):
        clinical_df = pd.DataFrame({
            'Patient ID': ['P001', 'P002']
        })
        predictions_df = pd.DataFrame({
            'patientId': ['P001', 'P001', 'P002', 'P003', 'P003', 'P003'],
            'Binding_Classification': ['SB', 'WB', 'SB', 'SB', 'WB', 'WB']
        })
        expected_df = pd.DataFrame({
            'Patient ID': ['P001', 'P002'],
            'Neoantigen_SB_Count': [1, 1],
            'Neoantigen_WB_Count': [1, 0]
        })
        result_df = mutationModifications.calcularNeoantigenosPaciente(clinical_df, predictions_df)
        pd.testing.assert_frame_equal(result_df, expected_df)

    def test_calcularNeoantigenosPaciente_missing_patients(self):
        clinical_df = pd.DataFrame({
            'Patient ID': ['P001', 'P002', 'P003', 'P004']
        })
        predictions_df = pd.DataFrame({
            'patientId': ['P001', 'P001', 'P002', 'P003', 'P003', 'P003'],
            'Binding_Classification': ['SB', 'WB', 'SB', 'SB', 'WB', 'WB']
        })
        expected_df = pd.DataFrame({
            'Patient ID': ['P001', 'P002', 'P003', 'P004'],
            'Neoantigen_SB_Count': [1.0, 1.0, 1.0, 0.0],
            'Neoantigen_WB_Count': [1.0, 0.0, 2.0, 0.0]
        })
        result_df = mutationModifications.calcularNeoantigenosPaciente(clinical_df, predictions_df)
        pd.testing.assert_frame_equal(result_df, expected_df)


if __name__ == '__main__':
    unittest.main()