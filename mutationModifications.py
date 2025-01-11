# Description: Este script contiene funciones para modificar y clasificar mutaciones genéticas.

def create_mutations_dict(muts):
    """
    Convierte una lista de objetos de mutación en una lista de diccionarios con atributos específicos.
    Estructurar los datos correctamente con los datos que son relevantes.
    Args:
        muts (list): Una lista de objetos de mutación. Se espera que cada objeto de mutación tenga los siguientes atributos:
            - chr (str): Identificador del cromosoma.
            - startPosition (int): Posición inicial de la mutación.
            - endPosition (int): Posición final de la mutación.
            - referenceAllele (str): Alelo de referencia.
            - variantAllele (str): Alelo variante.
            - variantType (str): Tipo de variante.
            - gene.hugoGeneSymbol (str): Símbolo HUGO del gen.
            - proteinChange (str): Notación del cambio de proteína.
            - patientId (str): Identificador del paciente.
            - sampleId (str): Identificador de la muestra.
            - tumorAltCount (int): Conteo de alelos alternativos en el tumor.
            - tumorRefCount (int): Conteo de alelos de referencia en el tumor.
            - mutationType (str): Tipo de mutación.
            - molecularProfileId (str): Identificador del perfil molecular.
            - studyId (str): Identificador del estudio.
    Returns:
        list: Una lista de diccionarios, cada uno conteniendo los atributos de un objeto de mutación.
    """

    mutation_dicts = [
        {
            "chr": mutation.chr,
            "startPosition": mutation.startPosition,
            "endPosition": mutation.endPosition,
            "referenceAllele": mutation.referenceAllele,
            "variantAllele": mutation.variantAllele,
            "variantType": mutation.variantType,
            "Gene": mutation.gene.hugoGeneSymbol,
            "Protein Change": mutation.proteinChange,
            "patientId": mutation.patientId,
            "sampleId": mutation.sampleId,
            "tumorAltCount": mutation.tumorAltCount,
            "tumorRefCount": mutation.tumorRefCount,
            "Mutation Type": mutation.mutationType,
            "molecularProfileId": mutation.molecularProfileId,
            "studyId": mutation.studyId
        }
        for mutation in muts
    ]
    return mutation_dicts

# Definir función para clasificar mutaciones
def clasificar_mutacion(row):
    """
    Clasifica el tipo de mutación genética basado en los alelos de referencia y alternativos.
    Args:
        row (dict): Un diccionario que contiene las claves "referenceAllele" y "variantAllele",
                    representando los nucleótidos de referencia y alternativos respectivamente.
    Returns:
        str: Una cadena que describe el tipo de mutación. Puede ser uno de los siguientes valores:
             - "Deleción": Si el alelo alternativo es "-" y si la longitud del alelo de referencia es mayor que la del alelo alternativo
             - "Inserción": Si el alelo de referencia es "-" y si la longitud del alelo de referencia es menor que la del alelo alternativo
             - "Transición": Si ambos alelos tienen longitud 1 o tienen la misma longitud mayor a 1 y son del mismo tipo (purina a purina o pirimidina a pirimidina)
             - "Transversión": Si ambos alelos tienen longitud 1 o  tienen la misma longitud mayor a 1 y son de diferente tipo (purina a pirimidina o viceversa)
             - "Otro": Para cualquier otro tipo de mutación no cubierto por las reglas anteriores
    """
    ref = row["referenceAllele"]  # Nucleótido(s) de referencia
    alt = row["variantAllele"]  # Nucleótido(s) alternativo
    
    # Deleción: alelo alternativo es "-"
    if alt == "-":
        return "Deleción"
    
    # Inserción: alelo de referencia es "-"
    if ref == "-":
        return "Inserción"
    
    # Mutaciones puntuales (longitud de 1 en ambos alelos)
    if len(ref) == 1 and len(alt) == 1:
        if (ref in "AG" and alt in "AG") or (ref in "CT" and alt in "CT"):
            return "Transición"
        else:
            return "Transversión"
    
    # Mutaciones complejas: misma longitud mayor a 1
    if len(ref) == len(alt) and len(ref) > 1:
        transversion_detected = any((r in "AG" and a in "CT") or (r in "CT" and a in "AG") for r, a in zip(ref, alt))
        if transversion_detected:
            return "Transversión"
        else:
            return "Transición"
    
    # Inserciones y deleciones más generales
    if len(ref) > len(alt):
        return "Deleción"
    if len(ref) < len(alt):
        return "Inserción"
    
    # Caso general para cualquier otro tipo de mutación
    return "Otro"



# 
def generate_mutated_peptides(sequence, mutation, length=9):
    """
    Genera péptidos con la mutación en cada posición a partir de una secuencia dada de longitud length.
    Parámetros:
    sequence (str): La secuencia original del péptido.
    mutation (str): La mutación en el formato 'A1B', donde 'A' es el aminoácido original,
                    '1' es la posición (índice basado en 1), y 'B' es el aminoácido mutado.
    length (int, optional): La longitud de los péptidos a generar. El valor predeterminado es 9.
    Retorna:
    list: Una lista de péptidos mutados de la longitud especificada.
            Si la secuencia es None, devuelve una lista vacía.
    """
    if sequence is None:
        return []  # Devolver lista vacía si no hay secuencia
    
    position = int(mutation[1:-1]) - 1  # Ajustar el índice
    mutated_peptides = []

    # Generar péptidos con la mutación en cada posición del péptido
    for i in range(length):
        start = max(0, position - i)
        end = min(len(sequence), position - i + length)
        segment = sequence[start:end]

        # Aplicar la mutación en la posición específica dentro del segmento
        if i < len(segment):
            segment = segment[:i] + mutation[-1] + segment[i+1:]

        # Asegurarse de que el péptido tenga la longitud especificada
        if len(segment) == length:
            mutated_peptides.append(segment)

    return mutated_peptides


def generate_peptides(row):
    """
    Genera una lista de péptidos mutados a partir de una fila de datos.
    Args:
        row (dict): Un diccionario que contiene las siguientes claves:
            - "Gene" (str): El nombre del gen.
            - "Protein Change" (str): La mutación en la proteína.
            - "Protein_Sequence" (str): La secuencia base de la proteína.
            - "patientId" (str): El ID del paciente.
            - "sampleId" (str): El ID de la muestra.
    Returns:
        list: Una lista de diccionarios, cada uno conteniendo:
            - "peptido" (str): La secuencia del péptido mutado.
            - "gen" (str): El nombre del gen.
            - "patientId" (str): El ID del paciente.
            - "sampleId" (str): El ID de la muestra.
    """
    gene = row["Gene"]
    mutation = row["Protein Change"]
    base_sequence = row["Protein_Sequence"]
    patient_id = row["patientId"]
    sample_id = row["sampleId"]
    peptides = generate_mutated_peptides(base_sequence, mutation)
    return [{"peptido": peptide, "gen": gene, "patientId": patient_id, "sampleId": sample_id} for peptide in peptides]

# Clasificar las predicciones en WB y SB
def classify_binding(row, presentation_percentile_hard = 0.5, presentation_percentile_soft = 2):
    """
    Clasifica la fuerza de unión del neoantígeno contra MHC basada en la afinidad y el percentil de presentación.
    Parámetros:
        row (dict): Un diccionario que contiene las claves 'affinity' y 'presentation_percentile' con sus respectivos valores.
        presentation_percentile_hard (float, opcional): El umbral estricto para el percentil de presentación. Por defecto es 0.5.
        presentation_percentile_soft (float, opcional): El umbral suave para el percentil de presentación. Por defecto es 2.
    Retorna:
        str: Una cadena que indica la fuerza de unión:
            - "SB" para unión fuerte si el percentil de presentación es menor que presentation_percentile_hard.
            - "WB" para unión débil si el percentil de presentación está entre presentation_percentile_hard y presentation_percentile_soft.
            - "N/A" si ninguna de las condiciones anteriores se cumple.
    """

    if row["presentation_percentile"] <= presentation_percentile_hard:
        return "SB" # Unión fuerte
    elif (row["presentation_percentile"] < presentation_percentile_soft and row["presentation_percentile"] >= presentation_percentile_hard):
        return "WB" # Unión débil
    else:
        return "N/A"  # No aplica

def calcularNeoantigenosPaciente(clinical_df, predictions_df):
    # Contar el número de neoantígenos SB y WB por paciente
    neoantigen_counts = predictions_df.groupby('patientId')['Binding_Classification'].value_counts().unstack(fill_value=0)

    # Renombrar las columnas para mayor claridad
    neoantigen_counts.rename(columns={'SB': 'Neoantigen_SB_Count', 'WB': 'Neoantigen_WB_Count'}, inplace=True)

    # Añadir columnas para los pacientes que no tienen neoantígenos SB o WB
    clinical_df = clinical_df.merge(neoantigen_counts, left_on='Patient ID', right_index=True, how='left')

    # Rellenar los valores NaN con 0 para los pacientes sin neoantígenos SB o WB
    clinical_df['Neoantigen_SB_Count'] = clinical_df['Neoantigen_SB_Count'].fillna(0)
    clinical_df['Neoantigen_WB_Count'] = clinical_df['Neoantigen_WB_Count'].fillna(0)

    return clinical_df

   
