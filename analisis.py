import getGraphics
import getDuplicates
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, kstest


##### Preprocesamiento de los datos  #########


# Leer el archivo de datos clínicos
clinical_df = pd.read_csv('es_dfarber_broad_2014_clinical_data_with_neoantigens.csv')

# Filtrar solo las muestras que son tumores
clinical_df = clinical_df[clinical_df['Sample Class'] == 'Tumor']

# Obtener los datos clínicos de los pacientes únicos. Se quitan las muestras de los pacientes
clinical_df = clinical_df.drop_duplicates(subset=['Patient ID'])

# Eliminar columna númeor de muestras
clinical_df = clinical_df.drop(columns=['Number of Samples Per Patient'])

# Guardar en .csv el archivo con el resumen de los datos clínicos
clinical_df.describe().to_csv('resultados/es_dfarber_broad_2014_clinical_data_summary.csv')

# Crear los intervalos de edad basados en la edad máxima
clinical_df = getGraphics.createIntervalosEdad(clinical_df)

#print(clinical_df["Age Interval"])

##### Ejecución de pruebas estadísticas para comprobar la normalidad de las variables numéricas  #########

df = clinical_df.select_dtypes(include=['number']) # Seleccionar solo las columnas numéricas
df = df.dropna() # Elimina filas con valores nulos
for column in df.columns:
    print(f"\nPruebas para la columna: {column}") 
    series = df[column] 
    # Prueba de Shapiro-Wilk 
    stat, p = shapiro(series) 
    print('Prueba de Shapiro-Wilk: Estadístico=%.3f, p=%.3f' % (stat, p)) 
    # Prueba de Kolmogorov-Smirnov 
    stat, p = kstest(series, 'norm') 
    print('Prueba de Kolmogorov-Smirnov: Estadístico=%.3f, p=%.3f' % (stat, p))

##### Generación de gráficos  #########

# Obtener matriz de correlación de variables numéricas
getGraphics.matrizCorrelacion(clinical_df)

# Crear boxplots de neoantígenos por atributo
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Sex')
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Age Interval')
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Overall Survival Status')

box_pairs=[
            (("White/Europe", "Neoantigen_SB_Count"), ("White/Latin America", "Neoantigen_SB_Count")),
            (("White/Europe", "Neoantigen_WB_Count"), ("White/Latin America", "Neoantigen_WB_Count"))
            ]
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Ethnicity Category', box_pairs)

# Obtener gráfico de número de mutaciones por tipo de mutación
getGraphics.mutacionesTipo()

# Obtener gráfico de número de neoantígenos por tipo de neoantígeno
getGraphics.neoantigenosFuertesVsDebiles()

# Obtener gráfico de número de neoantígenos por tipo de mutación
getGraphics.neoantigenosPorMutacion()

##### Buscar neoantígenos combinados  #########

getDuplicates.genesConNeoantigenosRepetidos()
getDuplicates.neoantigenosRepetidosPorPeptido()
getDuplicates.pacientesConNeoantigenosRepetidos('peptide')
getDuplicates.pacientesConNeoantigenosRepetidos('gen')
getDuplicates.genesConNeoantigenosPorPetido()


######### Generación de gráficos auxiliares  #########


# Crear boxplots para cada columna numérica en función de los intervalos de edad
for column in clinical_df.select_dtypes(include=['number']).columns:
    if column != 'Diagnosis Age':  # Eliminar la columna de edad
        #column = clinical_df.dropna(subset=[column])
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='Age Interval', y=column, data=clinical_df)
        plt.title(f'Boxplot of {column} by Age Interval')
        plt.xlabel('Age Interval')
        plt.ylabel(column)
        plt.xticks(rotation=45)
        plt.savefig(f'Figuras 2/Boxplot of {column} by Age Interval.png', bbox_inches='tight')

# Crear boxplots para cada columna numérica en función del sexo
getGraphics.boxplotsVariablesNumericasPorAtributo(clinical_df, 'Sex')

# Crear boxplots para cada columna numérica en función de la etnia
getGraphics.boxplotsVariablesNumericasPorAtributo(clinical_df, 'Ethnicity Category')

# Crear boxplots para cada columna numérica en función de la etnia
getGraphics.boxplotsVariablesNumericasPorAtributo(clinical_df, 'Overall Survival Status')