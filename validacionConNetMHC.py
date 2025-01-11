import pandas as pd
import glob
import os
from venn import venn
import matplotlib.pyplot as plt

######### Crear los archivos para pasar a NetMHC #########


# Leer el archivo CSV
df = pd.read_csv('resultados/unique_predictions.csv')

# Seleccionar solo la columna 'peptide'
df_peptide = df[['peptide']]

# Dividir el DataFrame en partes de 5000 registros
chunk_size = 5000
for i in range(0, len(df_peptide), chunk_size):
    chunk = df_peptide[i:i + chunk_size]
    chunk.to_csv(f'archivosNetMHC/unique_predictions_part_{i//chunk_size + 1}.csv', header=False, index=False)


########## Unificar predicciones de netMHC ##########

# Definir el directorio donde se encuentran los archivos CSV
carpeta = 'archivosGeneradosNetMHC'

# Buscar todos los archivos CSV en la carpeta
archivos_csv = glob.glob(os.path.join(carpeta, '*.csv'))

# Leer y concatenar todos los archivos CSV
lista_df = []
for archivo in archivos_csv:
    # Leer el archivo CSV
    df = pd.read_csv(archivo, sep=';')
    # Añadir el DataFrame a la lista
    lista_df.append(df)

# Concatenar todos los DataFrames
df_concatenado = pd.concat(lista_df, ignore_index=True)

# Guardar el archivo resultante
df_concatenado.to_csv('predicciones_netMHC.csv', index=False)


########## Clasificar las predicciones de netMHC ##########

def classify_binding(row, presentation_percentile_hard = 0.5, presentation_percentile_soft = 2):
    """
    Clasifica la fuerza de unión del neoantígeno contra MHC basada en la afinidad y el percentil de presentación.
    Parámetros:
        row (dict): Un diccionario que contiene las claves 'affinity' y 'presentation_percentile' con sus respectivos valores.
        presentation_percentile_hard (float, opcional): El umbral estricto para el percentil de presentación. Por defecto es 0.5.
        presentation_percentile_soft (float, opcional): El umbral suave para el percentil de presentación. Por defecto es 2.
    Retorna:
        str: Una cadena que indica la fuerza de unión:
            - "SB" para unión fuerte si la afinidad es menor que affinity_threshold_hard y el percentil de presentación es menor que presentation_percentile_hard.
            - "WB" para unión débil si la afinidad está entre affinity_threshold_hard y affinity_threshold_soft, y el percentil de presentación está entre presentation_percentile_hard y presentation_percentile_soft.
            - "N/A" si ninguna de las condiciones anteriores se cumple.
    """
    if row["Rank"] <= presentation_percentile_hard:
        return "SB" # Unión fuerte
    elif (row["Rank"] < presentation_percentile_soft and row["Rank"] >= presentation_percentile_hard):
        return "WB" # Unión débil
    else:
        return "N/A"  # No aplica
    
# Clasificar las predicciones en SB (Strong Binding) y WB (Weak Binding)
df_concatenado["Binding_Classification"] = df_concatenado.apply(classify_binding, axis=1)

df_concatenado.to_csv('predicciones_netMHC.csv', index=False)


########## Diagrama de Venn ##########


# Leer los archivos CSV
netMHC = pd.read_csv('predicciones_netMHC.csv')
MHCFlurry = pd.read_csv('resultados/unique_predictions.csv')

# Crear conjuntos
sb_peptides_net = set(netMHC[netMHC['Binding_Classification'] == 'SB']['Peptide'])
wb_peptides_net = set(netMHC[netMHC['Binding_Classification'] == 'WB']['Peptide'])

sb_peptides_flurry = set(MHCFlurry[MHCFlurry['Binding_Classification'] == 'SB']['peptide'])
wb_peptides_flurry = set(MHCFlurry[MHCFlurry['Binding_Classification'] == 'WB']['peptide'])

# Crear el diccionario para la librería venn
data = {
    'SB (netMHC)': sb_peptides_net,
    'WB (netMHC)': wb_peptides_net,
    'SB (MHCFlurry)': sb_peptides_flurry,
    'WB (MHCFlurry)': wb_peptides_flurry,
}

# Generar el diagrama de Venn
venn(data)
plt.title('Diagrama de Venn para comparar las predicciones de netMHC y MHCFlurry')
plt.savefig(f'Figuras 2/Diagrama de Venn.png', bbox_inches='tight')





