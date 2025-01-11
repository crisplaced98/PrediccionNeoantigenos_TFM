import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from itertools import combinations


def mutacionesTipo():
    """
    Lee un archivo CSV de mutaciones, cuenta el número de mutaciones por tipo,
    y genera un gráfico de barras que muestra el número de mutaciones por tipo.
    El archivo CSV debe tener una columna llamada "Clasificación" que contiene
    los tipos de mutaciones.
    El gráfico generado se guarda en un archivo PNG en la carpeta "Figuras".
    Returns:
        None
    """
    # Leer el archivo de mutaciones
    df = pd.read_csv("resultados/mutations.csv")

    # Contar el número de mutaciones por tipo
    mutations_count = df["Clasificación"].value_counts().reset_index()
    mutations_count.columns = ["Tipo de Mutación", "Número de Mutaciones"]

    # Crear el gráfico de barras
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Tipo de Mutación", y="Número de Mutaciones", data=mutations_count, palette="viridis")

    # Añadir títulos y etiquetas
    plt.title("Número de Mutaciones por Tipo de Mutación")
    plt.xlabel("Tipo de Mutación")
    plt.ylabel("Número de Mutaciones")

    # Rotar las etiquetas del eje x para mayor legibilidad
    plt.xticks(rotation=45)

    # Guardar el gráfico en un archivo 
    plt.savefig(f'Figuras/número_mutaciones_por_tipo.png', bbox_inches='tight')

def neoantigenosFuertesVsDebiles():
    """
    Genera un gráfico de barras comparando el número de neoantígenos fuertes (SB) y débiles (WB).
    Esta función lee un archivo CSV que contiene predicciones de neoantígenos, cuenta el número de neoantígenos fuertes y débiles,
    y crea un gráfico de barras para visualizar la comparación. El gráfico se guarda como un archivo PNG.
    El archivo CSV debe tener una columna llamada "Binding_Classification" que contiene la clasificación de neoantígenos
    como fuertes (SB) o débiles (WB).
    El gráfico resultante se guarda en 'Figuras/número_antígenos_fuertes_vs_debiles.png'.
    Returns:
        None
    """
    predictions_df = pd.read_csv("resultados/unique_predictions.csv")
    # Contar el número de neoantígenos fuertes (SB) y débiles (WB)
    neoantigen_counts = predictions_df["Binding_Classification"].value_counts().reset_index()
    neoantigen_counts.columns = ["Clasificación", "Número de Neoantígenos"]

    # Crear el gráfico de barras
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Clasificación", y="Número de Neoantígenos", data=neoantigen_counts, palette="viridis")

    # Añadir títulos y etiquetas
    plt.title("Número de Neoantígenos Fuertes (SB) vs Débiles (WB)")
    plt.xlabel("Clasificación")
    plt.ylabel("Número de Neoantígenos")

    plt.savefig('Figuras/número_antígenos_fuertes_vs_debiles.png', bbox_inches='tight')

def neoantigenosPorMutacion():
    """
    Genera un gráfico de barras que muestra el número de neoantígenos débiles (WB) y fuertes (SB) para cada tipo de mutación.
    Esta función realiza los siguientes pasos:
        1. Lee los datos de mutaciones y predicciones desde archivos CSV.
        2. Combina los datos de mutaciones con los datos de predicciones basándose en 'patientId' y 'Gene'.
        3. Cuenta el número de neoantígenos débiles y fuertes para cada tipo de mutación.
        4. Crea y guarda un gráfico de barras visualizando el número de neoantígenos por tipo de mutación y clasificación de unión.
        5. El gráfico resultante se guarda como 'neoantígenos_por_tipo_mutación.png' en el directorio 'Figuras'.
    Archivos CSV de entrada:
        'resultados/mutationsToBeTreated.csv': Contiene datos de mutaciones a tratar.
        'resultados/unique_predictions.csv': Contiene datos de predicciones únicas.
    Archivo CSV de salida:
        'resultados/combinaciónmutaciones.csv': Contiene los datos combinados de mutaciones y predicciones.
    Archivo de salida del gráfico:
        'Figuras/neoantígenos_por_tipo_mutación.png': El gráfico de barras generado.
    Returns:
        None
    """
    # Leer las tablas
    mutations_to_be_treated_df = pd.read_csv('resultados/mutationsToBeTreated.csv')
    predictions_df = pd.read_csv('resultados/unique_predictions.csv')

    # Unir las mutaciones tratadas con las predicciones usando 'patientId' y 'Gene'
    mutations_combined = mutations_to_be_treated_df.merge(predictions_df, left_on=['patientId', 'Gene'], right_on=['patientId', 'gen'])
    mutations_combined.to_csv("resultados/combinaciónmutaciones.csv", index=False)

    # Contar el número de neoantígenos débiles (WB) y fuertes (SB) para cada tipo de mutación
    neoantigen_counts = mutations_combined.groupby(['Clasificación', 'Binding_Classification']).size().reset_index(name='Número de Neoantígenos')

    # Crear el gráfico de barras
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Clasificación', y='Número de Neoantígenos', hue='Binding_Classification', data=neoantigen_counts, palette='viridis')

    # Añadir títulos y etiquetas
    plt.title("Número de Neoantígenos Débiles (WB) vs Fuertes (SB) para Transiciones y Transversiones")
    plt.xlabel("Tipo de Mutación")
    plt.ylabel("Número de Neoantígenos")
    plt.legend(title='Clasificación de Unión')

    plt.savefig('Figuras/neoantígenos_por_tipo_mutación.png', bbox_inches='tight')

def matrizCorrelacion(clinical_df):
    """
    Genera y guarda un heatmap de la matriz de correlación de las variables numéricas en un DataFrame clínico.
    El gráfico resultante se guarda como 'matriz_correlacion.png' en el directorio 'Figuras 2'.
    Args:
        clinical_df (pandas.DataFrame): DataFrame que contiene los datos clínicos.
    Returns:
        None
    """
    # Seleccionar solo las columnas numéricas 
    numerical_df = clinical_df.select_dtypes(include=['number']) 
    # Calcular la matriz de correlación 

    correlation_matrix = numerical_df.corr()

    # Configurar el tamaño del gráfico
    plt.figure(figsize=(12, 10))

    # Crear un heatmap usando seaborn
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0)

    # Añadir títulos
    plt.title('Matriz de Correlación de Variables Clínicas')

    # Guardar el gráfico
    plt.savefig('Figuras 2/matriz_correlacion.png', bbox_inches='tight')

def createIntervalosEdad(clinical_df, num_intervals=5):
    """
    Crear intervalos de edad para un DataFrame clínico.
    Esta función toma un DataFrame que contiene datos clínicos y crea intervalos de edad
    basados en la columna 'Diagnosis Age'. Los intervalos se crean dividiendo el 
    rango de edades en un número especificado de intervalos.
    Parámetros:
        clinical_df (pd.DataFrame): Un DataFrame que contiene datos clínicos con una columna 'Diagnosis Age'.
        num_intervals (int, opcional): El número de intervalos a crear. El valor predeterminado es 5.
    Retorna:
        pd.DataFrame: El DataFrame original con una columna adicional 'Age Interval' que contiene los intervalos de edad.
    """
    edad_maxima = clinical_df['Diagnosis Age'].max()
    bins = pd.interval_range(start=0, end=edad_maxima, freq=edad_maxima//num_intervals)
    labels = [f'{int(interval.left)}-{int(interval.right)}' for interval in bins]

    clinical_df['Age Interval'] = pd.cut(clinical_df['Diagnosis Age'], bins=bins, labels=labels, right=False)

    return clinical_df

def boxplotsVariablesNumericasPorAtributo(clinical_df, atributo):
    """
    Genera y guarda boxplots para variables numéricas en el DataFrame dado agrupadas por un atributo especificado.
    Parámetros:
        clinical_df (pandas.DataFrame): El DataFrame que contiene los datos clínicos.
        atributo (str): El atributo por el cual agrupar las variables numéricas para los boxplots.
    Retorna:
        None
    """
    for column in clinical_df.select_dtypes(include=['number']).columns:
        plt.figure(figsize=(10, 6))
        sns.boxplot(x=atributo, y=column, data=clinical_df)
        plt.title(f'Boxplot of {column} by {atributo}')
        plt.xlabel(atributo)
        plt.ylabel(column)
        plt.savefig(f'Figuras 2/Boxplot of {column} by {atributo}.png', bbox_inches='tight')

def boxplotConjuntoNeoantigenoPorAtributo(clinical_df, atributo, box_pairs=None):
    """
    Genera un boxplot comparando los conteos de dos tipos de neoantígenos (SB y WB) a través de diferentes categorías de un atributo dado.
    Parámetros:
        clinical_df (pd.DataFrame): DataFrame que contiene datos clínicos, incluyendo conteos de neoantígenos y el atributo a comparar.
        atributo (str): El atributo en el DataFrame por el cual agrupar y comparar.
        box_pairs (list of tuples, opcional): Pares específicos de cajas a comparar para la anotación estadística. Si es None, se comparan todos los pares posibles.
    Retorna:
        None
    La función realiza los siguientes pasos:
        1. Convierte el DataFrame a formato largo adecuado para seaborn.
        2. Identifica categorías únicas del atributo dado.
        3. Crea pares de cajas para comparación si no se proporcionan.
        4. Genera un boxplot con seaborn.
        5. Añade anotaciones estadísticas al gráfico.
        6. Guarda el gráfico como un archivo PNG.
    """
    # Crear un DataFrame en formato largo (long format) para seaborn
    long_df = pd.melt(clinical_df, id_vars=[atributo], value_vars=['Neoantigen_SB_Count', 'Neoantigen_WB_Count'],
                    var_name='Neoantigen_Type', value_name='Count')
    
    # Obtener las categorías únicas del atributo para hacer las comparaciones
    categoria_atributo = long_df[atributo].dropna().unique()
    neoantigen_types = ['Neoantigen_SB_Count', 'Neoantigen_WB_Count']

    # Crear pares de cajas para comparar
    if box_pairs is None:
        box_pairs = [((valor1, type_), (valor2, type_)) for valor1, valor2 in combinations(categoria_atributo, 2) for type_ in neoantigen_types]
    
    # Crear un boxplot de las dos variables combinadas por atributo
    plt.figure(figsize=(12, 8))
    ax = sns.boxplot(x=atributo, y='Count', hue='Neoantigen_Type', data=long_df)

    # Añadir las pruebas estadísticas, para mostar las significancias
    add_stat_annotation(ax, data=long_df, x=atributo, y='Count', hue='Neoantigen_Type',
                    box_pairs=box_pairs,
                    test='Mann-Whitney', text_format='star', loc='inside')
    plt.title(f'Boxplot del Número de Neoantígenos SB y WB por {atributo}')
    plt.xlabel(atributo)
    plt.ylabel('Número de Neoantígenos')
    plt.legend(title='Tipo de Neoantígeno')

    plt.savefig(f'Figuras 2/Boxplot del Número de Neoantígenos SB y WB por {atributo}.png', bbox_inches='tight')