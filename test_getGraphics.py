import unittest
import pandas as pd
import os
import getGraphics

class TestGetGraphics(unittest.TestCase):

    def setUp(self):
        # Crear un DataFrame de muestra para pruebas
        data = {
            'Diagnosis Age': [25, 35, 45, 55, 65],
            'Tumor Size': [2.5, 3.0, 3.5, 4.0, 4.5],
            'Patient Status': ['Alive', 'Deceased', 'Alive', 'Deceased', 'Alive']
        }
        self.clinical_df = pd.DataFrame(data)
        self.attribute = 'Patient Status'
        
        # Crear directorio para guardar figuras si no existe
        if not os.path.exists('Figuras 2'):
            os.makedirs('Figuras 2')

        # Crear directorio para guardar figuras si no existe
        if not os.path.exists('Figuras'):
            os.makedirs('Figuras')

        # Crear un DataFrame de muestra para pruebas de mutaciones
        mutations_data = {
            'Clasificación': ['Tipo1', 'Tipo2', 'Tipo1', 'Tipo3', 'Tipo2', 'Tipo1']
        }
        self.mutations_df = pd.DataFrame(mutations_data)
        self.mutations_file_path = 'resultados/mutations.csv'
        self.mutations_df.to_csv(self.mutations_file_path, index=False)

    def test_boxplotsVariablesNumericasPorAtributo(self):
        # Llamar a la función para generar boxplots
        getGraphics.boxplotsVariablesNumericasPorAtributo(self.clinical_df, self.attribute)
        
        # Verificar si los archivos fueron creados
        for column in self.clinical_df.select_dtypes(include=['number']).columns:
            file_path = f'Figuras 2/Boxplot of {column} by {self.attribute}.png'
            self.assertTrue(os.path.exists(file_path))
            # Limpiar los archivos generados
            if os.path.exists(file_path):
                os.remove(file_path)

    def test_mutacionesTipo(self):
        # Llamar a la función para generar el gráfico de mutaciones
        getGraphics.mutacionesTipo()
        
        # Verificar si el archivo fue creado
        file_path = 'Figuras/número_mutaciones_por_tipo.png'
        self.assertTrue(os.path.exists(file_path))
        # Limpiar los archivos generados
        if os.path.exists(file_path):
            os.remove(file_path)

    def test_neoantigenosFuertesVsDebiles(self):
        # Crear un DataFrame de muestra para pruebas de neoantígenos
        neoantigen_data = {
            'Binding_Classification': ['SB', 'WB', 'SB', 'WB', 'SB', 'WB', 'SB']
        }
        self.neoantigen_df = pd.DataFrame(neoantigen_data)
        self.neoantigen_file_path = 'resultados/unique_predictions.csv'
        self.neoantigen_df.to_csv(self.neoantigen_file_path, index=False)

        # Llamar a la función para generar el gráfico de neoantígenos
        getGraphics.neoantigenosFuertesVsDebiles()
        
        # Verificar si el archivo fue creado
        file_path = 'Figuras/número_antígenos_fuertes_vs_debiles.png'
        self.assertTrue(os.path.exists(file_path))
        # Limpiar los archivos generados
        if os.path.exists(file_path):
            os.remove(file_path)

        # Eliminar el archivo de neoantígenos de prueba
        if os.path.exists(self.neoantigen_file_path):
            os.remove(self.neoantigen_file_path)

    def test_neoantigenosPorMutacion(self):
        # Crear un DataFrame de muestra para pruebas de mutaciones a tratar
        mutations_to_be_treated_data = {
            'patientId': [1, 2, 3, 4, 5],
            'Gene': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
            'Clasificación': ['Tipo1', 'Tipo2', 'Tipo1', 'Tipo3', 'Tipo2']
        }
        self.mutations_to_be_treated_df = pd.DataFrame(mutations_to_be_treated_data)
        self.mutations_to_be_treated_file_path = 'resultados/mutationsToBeTreated.csv'
        self.mutations_to_be_treated_df.to_csv(self.mutations_to_be_treated_file_path, index=False)

        # Crear un DataFrame de muestra para pruebas de predicciones únicas
        unique_predictions_data = {
            'patientId': [1, 2, 3, 4, 5],
            'gen': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
            'Binding_Classification': ['SB', 'WB', 'SB', 'WB', 'SB']
        }
        self.unique_predictions_df = pd.DataFrame(unique_predictions_data)
        self.unique_predictions_file_path = 'resultados/unique_predictions.csv'
        self.unique_predictions_df.to_csv(self.unique_predictions_file_path, index=False)

        # Llamar a la función para generar el gráfico de neoantígenos por mutación
        getGraphics.neoantigenosPorMutacion()
        
        # Verificar si el archivo fue creado
        file_path = 'Figuras/neoantígenos_por_tipo_mutación.png'
        self.assertTrue(os.path.exists(file_path))
        # Limpiar los archivos generados
        if os.path.exists(file_path):
            os.remove(file_path)

        # Eliminar los archivos de prueba
        if os.path.exists(self.mutations_to_be_treated_file_path):
            os.remove(self.mutations_to_be_treated_file_path)
        if os.path.exists(self.unique_predictions_file_path):
            os.remove(self.unique_predictions_file_path)
        if os.path.exists('resultados/combinaciónmutaciones.csv'):
            os.remove('resultados/combinaciónmutaciones.csv')

    def test_matrizCorrelacion(self):
        # Llamar a la función para generar el heatmap de la matriz de correlación
        getGraphics.matrizCorrelacion(self.clinical_df)
        
        # Verificar si el archivo fue creado
        file_path = 'Figuras 2/matriz_correlacion.png'
        self.assertTrue(os.path.exists(file_path))
        # Limpiar los archivos generados
        if os.path.exists(file_path):
            os.remove(file_path)

    def test_createIntervalosEdad(self):
        # Llamar a la función para crear intervalos de edad
        result_df = getGraphics.createIntervalosEdad(self.clinical_df, num_intervals=3)
        
        # Verificar si la columna 'Age Interval' fue creada
        self.assertIn('Age Interval', result_df.columns)
        
        # Verificar si los intervalos de edad son correctos
        expected_intervals = ['(0, 21]', '(21, 42]', '(42, 63]']
        actual_intervals = result_df['Age Interval'].cat.categories.astype(str).tolist()
        self.assertEqual(expected_intervals, actual_intervals)
        
        # Verificar si los valores de 'Age Interval' son correctos
        expected_values = ['(21, 42]', '(21, 42]', '(42, 63]', '(42, 63]', 'nan']
        actual_values = result_df['Age Interval'].astype(str).tolist()
        self.assertEqual(expected_values, actual_values)

    def test_boxplotConjuntoNeoantigenoPorAtributo(self):
        # Crear un DataFrame de muestra para pruebas de neoantígenos
        neoantigen_data = {
            'Patient Status': ['Alive', 'Deceased', 'Alive', 'Deceased', 'Alive', 'Deceased'],
            'Neoantigen_SB_Count': [10, 20, 30, 40, 50, 60],
            'Neoantigen_WB_Count': [5, 15, 25, 35, 45, 55]
        }
        self.neoantigen_df = pd.DataFrame(neoantigen_data)
        self.attribute = 'Patient Status'

        # Llamar a la función para generar el boxplot
        getGraphics.boxplotConjuntoNeoantigenoPorAtributo(self.neoantigen_df, self.attribute)
        
        # Verificar si el archivo fue creado
        file_path = f'Figuras 2/Boxplot del Número de Neoantígenos SB y WB por {self.attribute}.png'
        self.assertTrue(os.path.exists(file_path))
        # Limpiar los archivos generados
        if os.path.exists(file_path):
            os.remove(file_path)



    def tearDown(self):
        # Eliminar el directorio si está vacío
        if os.path.exists('Figuras 2') and not os.listdir('Figuras 2'):
            os.rmdir('Figuras 2')

        # Eliminar el directorio si está vacío
        if os.path.exists('Figuras') and not os.listdir('Figuras'):
            os.rmdir('Figuras')

        # Eliminar el archivo de mutaciones de prueba
        if os.path.exists(self.mutations_file_path):
            os.remove(self.mutations_file_path)


if __name__ == '__main__':
    unittest.main()