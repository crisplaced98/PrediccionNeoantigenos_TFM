import pandas as pd


def neoantigenosRepetidosPorPeptido():
    """
    Identifica y guarda neoantígenos duplicados basados en secuencias de péptidos.
    Esta función lee un archivo CSV que contiene predicciones de neoantígenos, filtra los datos para incluir solo 
    las filas con información sobre la clasificación de unión, identifica neoantígenos duplicados basados en la 
    columna 'peptide', y guarda los neoantígenos duplicados en un nuevo archivo CSV.
    Retorna:
        pandas.DataFrame: Un DataFrame que contiene los neoantígenos duplicados basados en secuencias de péptidos.
                            El archivo contiene los péptidos que que han sido considerados neoantígenos más de una vez.
    """
    # Leer la tabla de predicciones
    predictions_df = pd.read_csv('resultados/unique_predictions.csv')

    # Filtrar los datos para incluir solo las filas con información sobre la fuerza del neoantígeno
    neoantigen_data = predictions_df[predictions_df['Binding_Classification'].notna()]

    # Identificar los neoantígenos duplicados en la columna 'péptido'
    duplicated_neoantigens = neoantigen_data[neoantigen_data.duplicated(subset=['peptide'], keep=False)]

    # Crear DataFrame con los péptidos mutados
    duplicated_neoantigens_df = pd.DataFrame(duplicated_neoantigens)

    # Guardar los péptidos mutados en un archivo .csv
    duplicated_neoantigens_df.to_csv("resultados neoantigenos combinados/neoantigenosRepetidosPorPeptido.csv", index=False, sep=",")

    return duplicated_neoantigens

def genesConNeoantigenosRepetidos():
    """
    Identifica y guarda genes con neoantígenos repetidos.
    Esta función lee un archivo CSV que contiene datos de predicción, filtra los datos para incluir solo 
    las filas con información sobre la fuerza del neoantígeno, identifica neoantígenos duplicados basados 
    en la columna 'gen', y guarda el DataFrame resultante con los neoantígenos duplicados en un nuevo archivo CSV.
    El archivo CSV de entrada debe estar ubicado en 'resultados/unique_predictions.csv' y el archivo CSV de salida 
    se guardará en 'resultados neoantigenos combinados/genesConNeoantigenosRepetidos.csv'.
    Retorna:
        None
        El archivo contiene los genes que han sido considerados neoantígenos más de una vez.
    """
    # Leer la tabla de predicciones
    predictions_df = pd.read_csv('resultados/unique_predictions.csv')

    # Filtrar los datos para incluir solo las filas con información sobre la fuerza del neoantígeno
    neoantigen_data = predictions_df[predictions_df['Binding_Classification'].notna()]

    # Identificar los neoantígenos duplicados en la columna 'gen'
    duplicated_neoantigens = neoantigen_data[neoantigen_data.duplicated(subset=['gen'], keep=False)]

    # Crear DataFrame con los péptidos mutados
    duplicated_neoantigens_df = pd.DataFrame(duplicated_neoantigens)

    # Guardar los péptidos mutados en un archivo .csv
    duplicated_neoantigens_df.to_csv("resultados neoantigenos combinados/genesConNeoantigenosRepetidos.csv", index=False, sep=",")


def pacientesConNeoantigenosRepetidos(atributo):
    """
    Identifica y guarda los neoantígenos repetidos entre diferentes pacientes basándose en un atributo específico.
    Args:
        atributo (str): El nombre de la columna en la que se basará la identificación de neoantígenos repetidos.
    Returns:
        None: La función guarda un archivo .csv con los neoantígenos repetidos entre diferentes pacientes.
    Proceso:
        1. Lee la tabla de predicciones desde 'resultados/unique_predictions.csv'.
        2. Filtra los datos para incluir solo las filas con información sobre la fuerza del neoantígeno.
        3. Identifica los neoantígenos duplicados en la columna especificada entre diferentes pacientes.
        4. Agrupa los neoantígenos duplicados para ver cuáles están repetidos entre pacientes.
        5. Filtra los neoantígenos repetidos por más de un paciente.
        6. Crea un DataFrame con los péptidos mutados.
        7. Guarda los péptidos mutados en un archivo .csv en 'resultados neoantigenos combinados/'.
    """
    # Leer la tabla de predicciones
    predictions_df = pd.read_csv('resultados/unique_predictions.csv')

    # Filtrar los datos para incluir solo las filas con información sobre la fuerza del neoantígeno
    neoantigen_data = predictions_df[predictions_df['Binding_Classification'].notna()]

    # Identificar los neoantígenos duplicados en la columna 'Binding_Classification' entre diferentes pacientes
    duplicated_neoantigens = neoantigen_data[neoantigen_data.duplicated(subset=[atributo], keep=False)]

    # Agrupar por 'Binding_Classification' para ver cuáles neoantígenos están repetidos entre pacientes
    grouped_duplicated_neoantigens = duplicated_neoantigens.groupby(atributo)['patientId'].unique().reset_index()

    # Filtrar los neoantígenos repetidos por más de un paciente
    multiple_patients_neoantigens = grouped_duplicated_neoantigens[grouped_duplicated_neoantigens['patientId'].apply(lambda x: len(x) > 1)]

    # Crear DataFrame con los péptidos mutados
    multiple_patients_neoantigens_df = pd.DataFrame(multiple_patients_neoantigens)

    # Guardar los péptidos mutados en un archivo .csv
    multiple_patients_neoantigens_df.to_csv(f"resultados neoantigenos combinados/pacientesConNeoantigenosRepetidosPor_{atributo}.csv", index=False, sep=",")


def genesConNeoantigenosPorPetido():
    """
    Identifica y guarda neoantígenos que están duplicados entre diferentes pacientes.
    Esta función realiza los siguientes pasos:
        1. Lee un archivo CSV que contiene datos de predicción.
        2. Filtra los datos para incluir solo las filas con información sobre la fuerza del neoantígeno.
        3. Identifica neoantígenos duplicados en la columna 'Binding_Classification' entre diferentes pacientes.
        4. Agrupa los neoantígenos duplicados por 'gen' para ver cuáles neoantígenos están repetidos entre pacientes.
        5. Filtra los neoantígenos agrupados para incluir solo aquellos repetidos en más de un paciente.
        6. Crea un DataFrame con los péptidos mutados.
        7. Guarda el DataFrame en un archivo CSV.
    Returns:
        None: La función guarda un archivo .csv con los péptidos considerados neoantígenos agrupados por gen.
    El archivo CSV de entrada debe estar ubicado en 'resultados/unique_predictions.csv'.
    El archivo CSV de salida se guardará en 'resultados neoantigenos combinados/genesConNeoantigenosPorPetido.csv'.
    """
    # Leer la tabla de predicciones
    predictions_df = pd.read_csv('resultados/unique_predictions.csv')

    # Filtrar los datos para incluir solo las filas con información sobre la fuerza del neoantígeno
    neoantigen_data = predictions_df[predictions_df['Binding_Classification'].notna()]

    # Identificar los neoantígenos duplicados en la columna 'Binding_Classification' entre diferentes pacientes
    duplicated_neoantigens = neoantigen_data[neoantigen_data.duplicated(subset='gen', keep=False)]

    # Agrupar por 'Binding_Classification' para ver cuáles neoantígenos están repetidos entre pacientes
    grouped_duplicated_neoantigens = duplicated_neoantigens.groupby('gen')['peptide'].unique().reset_index()

    # Filtrar los neoantígenos repetidos por más de un paciente
    multiple_patients_neoantigens = grouped_duplicated_neoantigens[grouped_duplicated_neoantigens['peptide'].apply(lambda x: len(x) > 1)]

    # Crear DataFrame con los péptidos mutados
    multiple_patients_neoantigens_df = pd.DataFrame(multiple_patients_neoantigens)

    # Guardar los péptidos mutados en un archivo .csv
    multiple_patients_neoantigens_df.to_csv("resultados neoantigenos combinados/genesConNeoantigenosPorPetido.csv", index=False, sep=",")


