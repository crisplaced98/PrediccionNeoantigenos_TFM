# **************************************************************************** #

import pandas as pd
from mhcflurry import Class1PresentationPredictor

import getInformation
import mutationModifications
import uniProtCache


######### Obtener información clínica y mutaciones  #########


# Leer el archivo de datos clínicos
clinical_df = pd.read_csv('es_dfarber_broad_2014_clinical_data.tsv', sep='\t')

# Filtrar solo las muestras que son tumores
tumor_samples = clinical_df[clinical_df['Sample Class'] == 'Tumor']['Patient ID'].tolist()
'''
# Obtener las mutaciones del estudio
muts = getInformation.get_mutations_cBioPortal("es_dfarber_broad_2014")

# Estructurar los datos correctamente con los datos que son relevantes
mutation_dicts = mutationModifications.create_mutations_dict(muts)

# Crear un DataFrame de las mutaciones obtenidas
df = pd.DataFrame(mutation_dicts)


######### Filtrar datos  #########


# Filtrar mutaciones solo para los tumores
df = df[df['patientId'].isin(tumor_samples)]

# Clasificar las mutaciones y añadir una nueva columna al DataFrame
df["Clasificación"] = df.apply(mutationModifications.clasificar_mutacion, axis=1)

# Guardar las mutaciones que se han descargado en un archivo CSV
df.to_csv("resultados/mutations.csv", index=False)

# Filtrar solo las mutaciones de tipo 'Missense_Mutation' 
df = df[df["Mutation Type"] == "Missense_Mutation"]

# Guardar las mutaciones que se van a tratar en un archivo CSV
df.to_csv("resultados/mutationsToBeTreated.csv", index=False)
print("Datos clínicos cargados y filtrados")


######### Obtener secuencias de UniProt  #########


# Crear una instancia de UniProtCache
cache = uniProtCache.UniProtCache()

# Aplicar el método cached_get_uniprot_info a la columna "Gene"
uniprot_info = df["Gene"].apply(cache.cached_get_uniprot_info)
df["UniProt_ID"], df["Protein_Sequence"] = zip(*uniprot_info)

# Crear DataFrame con los péptidos mutados
uniprot_info_df = pd.DataFrame(uniprot_info)

# Guardar los péptidos mutados en un archivo .csv
uniprot_info_df.to_csv("resultados/uniprot_info_df.csv", index=False, sep=",")

print("Información de UniProt obtenida y añadida al DataFrame")


######### Generar péptidos mutados #########


# Generar todas las secuencias mutadas 
mutated_peptides = df.apply(mutationModifications.generate_peptides, axis=1).explode().tolist()
print("Péptidos mutados generados")

# Filtrar los datos válidos (solo diccionarios)
valid_peptides = [item for item in mutated_peptides if isinstance(item, dict)]

# Verificar si hay elementos inválidos en la lista
invalid_peptides = [item for item in mutated_peptides if not isinstance(item, dict)]
if invalid_peptides:
    print("Se encontraron elementos inválidos en la lista y fueron ignorados:", invalid_peptides)

# Crear DataFrame con los péptidos mutados
mutated_peptides_df = pd.DataFrame(valid_peptides)

# Guardar los péptidos mutados en un archivo .csv
mutated_peptides_df.to_csv("resultados/mutated_peptides.csv", index=False, sep=",")
print("Péptidos mutados guardados en mutated_peptides.csv")'''

######### Predecir neoanígenos y clasificarlos #########

mutated_peptides_df = pd.read_csv('resultados/mutated_peptides.csv')

# Cargar MHCflurry predictor
predictor = Class1PresentationPredictor.load()

# Predecir la afinidad de unión usando MHCflurry
predictions = predictor.predict(peptides=mutated_peptides_df["peptido"].tolist(), alleles=["HLA-A*02:01"])
print("Predicciones de afinidad de unión realizadas")

# Convertir las predicciones en un DataFrame 
predictions_df = pd.DataFrame(predictions) 

# Agregar la información del gen, el ID del paciente y la muestra a las predicciones 
predictions_df["gen"] = mutated_peptides_df["gen"] 
predictions_df["patientId"] = mutated_peptides_df["patientId"]
predictions_df["sampleId"] = mutated_peptides_df["sampleId"]

# Guardar las predicciones en un archivo .csv 
predictions_df.to_csv("resultados/predictions.csv", index=False, sep=",")
predictions_df = pd.read_csv('resultados/predictions.csv')

# Clasificar las predicciones en SB (Strong Binding) y WB (Weak Binding)
predictions_df["Binding_Classification"] = predictions_df.apply(mutationModifications.classify_binding, axis=1)

# Guardar las clasificaciones de las predicciones en un archivo .csv 
predictions_df.to_csv("resultados/predictions.csv", index=False, sep=",")


######### Tratamiento de los datos para unificar y separar en archivos  #########


# Eliminar filas duplicadas basadas en 'peptido', 'gen' y 'patientId'
unique_predictions = predictions_df.drop_duplicates(subset=['peptide', 'gen', 'patientId'])

# Guardar las predicciones únicas en un archivo .csv
unique_predictions.to_csv("resultados/unique_predictions.csv", index=False, sep=",")
print("Predicciones únicas guardadas en unique_predictions.csv")

unique_predictions = pd.read_csv('resultados/unique_predictions.csv')
# Filtrar y guardar los péptidos con alta probabilidad de presentación
strong_binding_peptides = predictions_df[predictions_df["Binding_Classification"] == "SB"]
strong_binding_peptides.to_csv("resultados/strong_binding_peptides.csv", index=False)
print("Predicciones con alta probabilidad de presentación guardadas en strong_binding_peptides.csv")

# Filtrar y guardar los péptidos con alta afinidad
weak_binding_peptides = predictions_df[predictions_df["Binding_Classification"] == "WB"]
weak_binding_peptides.to_csv("resultados/weak_binding_peptides.csv", index=False)
print("Predicciones de alta afinidad guardadas en weak_binding_peptides.csv")

# Calcular el número de neoantígenos por paciente y actualizar el DataFrame clínico
clinical_df = mutationModifications.calcularNeoantigenosPaciente(clinical_df, unique_predictions)

# Guardar el DataFrame actualizado en un nuevo archivo CSV
clinical_df.to_csv('es_dfarber_broad_2014_clinical_data_with_neoantigens.csv', index=False)

print("Archivo actualizado con los contajes de neoantígenos SB y WB guardado como 'es_dfarber_broad_2014_clinical_data_with_neoantigens.csv'")


