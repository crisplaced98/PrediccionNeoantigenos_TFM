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

# Eliminar columna número de muestras
clinical_df = clinical_df.drop(columns=['Number of Samples Per Patient'])

# Guardar en .csv el archivo con el resumen de los datos clínicos
clinical_df.describe().to_csv('resultados/es_dfarber_broad_2014_clinical_data_summary.csv')

# Crear los intervalos de edad basados en la edad máxima
clinical_df = getGraphics.createIntervalosEdad(clinical_df)


##### Ejecución de pruebas estadísticas para comprobar la normalidad de las variables numéricas  #########

df = clinical_df.select_dtypes(include=['number']) # Seleccionar solo las columnas numéricas
df = df.dropna() # Elimina filas con valores nulos

custom_titles = {
    'Diagnosis Age': 'Edad de diagnóstico',
    'Mutation Count': 'Número de mutaciones',
    'TMB (nonsynonymous)': 'TMB (no sinónimas)',
    'Neoantigen_SB_Count': 'Número de neoantígenos fuertes',
    'Neoantigen_WB_Count': 'Número de neoantígenos débiles',
}

# Generación de gráficos Q-Q para cada variable numérica
for column in df.columns:
    series = df[column] 
    # crear el gráfico qqplot
    title = custom_titles.get(column)
    getGraphics.graficoqq(series, title)

for idx, column in enumerate(df.columns, 1):
    print(f"\nPruebas para la columna: {column}") 
    series = df[column] 
    # Prueba de Shapiro-Wilk 
    stat, p = shapiro(series) 
    print(f'Prueba de Shapiro-Wilk: Estadístico={stat}, p={p}') 
    # Prueba de Kolmogorov-Smirnov 
    stat, p = kstest(series, 'norm') 
    print(f'Prueba de Kolmogorov-Smirnov: Estadístico={stat}, p={p}')

    # Crear subgráfico para cada columna
    plt.subplot(2, len(df.columns)//2 + 1, idx)  # 2 filas, columnas según el número de variables
    sns.boxplot(y=column, data=clinical_df)
    title = custom_titles.get(column)
    plt.title(f'Boxplot {title}')
    plt.ylabel(title)

plt.tight_layout()  # Ajusta el espaciado entre subgráficos
plt.savefig('Figuras/boxplots_columnas_numericas.png')  # Guarda la figura como imagen
plt.close()

##### Generación de gráficos  #########

# Obtener matriz de correlación de variables numéricas
getGraphics.matrizCorrelacion(clinical_df)

# Crear boxplots de neoantígenos por atributo
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Sex', title='Boxplot de neoantígenos por sexo', traduccion_atributo='Sexo')
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Age Interval', title='Boxplot de neoantígenos por intervalo de edad', traduccion_atributo='Intervalo de edad')
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Overall Survival Status', title='Boxplot de neoantígenos por estado de supervivencia', traduccion_atributo='Estado de supervivencia')

box_pairs=[
            (("White/Europe", "SB"), ("White/Latin America", "SB")),
            (("White/Europe", "WB"), ("White/Latin America", "WB"))
            ]
getGraphics.boxplotConjuntoNeoantigenoPorAtributo(clinical_df, 'Ethnicity Category', box_pairs=box_pairs, title='Boxplot de neoantígenos por etnia', traduccion_atributo='Etnia')

# Obtener gráfico de número de mutaciones por tipo de mutación
getGraphics.mutacionesTipo()

# Obtener gráfico de número de mutaciones por tipo de mutación general
getGraphics.mutacionesTipoGeneral()

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
