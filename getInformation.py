import requests
from bravado.client import SwaggerClient

### Función para obtener las mutaciones del estudio
def get_mutations_cBioPortal(studyId):
    """
    Obtener mutaciones de cBioPortal para un ID de estudio dado.
    Esta función establece una conexión con la API de cBioPortal y recupera
    datos de mutaciones para el ID de estudio especificado. Las mutaciones se
    obtienen del perfil molecular del estudio e incluyen información detallada
    de los genes.
    Args:
        studyId (str): El ID del estudio para el cual se desea recuperar datos de mutaciones.
    Returns:
        list: Una lista de mutaciones para el estudio especificado, incluyendo
              información detallada de los genes.
    """
    # Establecer la conexión con cBioPortal a su API
    cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/v2/api-docs',
                                        config={"validate_requests":False,"validate_responses":False,"validate_swagger_spec": False})

    # Obtener las mutaciones por estudio 
    muts = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
        molecularProfileId=f"{studyId}_mutations", # obtiene las mutaciones del perfil molecular del estudio 
        sampleListId=f"{studyId}_all", # obtiene todas las muestras
        projection="DETAILED" # obtiene la información de los genes
    ).result()

    return muts

# Función para buscar el identificador de UniProt a partir del nombre del gen y el organismo
def get_uniprot_id(gene_name, taxonomy_id="9606"):
    """
    Recuperar el ID de UniProt para un nombre de gen y un ID de taxonomía dados.
    Args:
        gene_name (str): El nombre del gen a buscar.
        taxonomy_id (str, opcional): El ID de taxonomía del organismo. Por defecto es "9606" (Homo sapiens).
    Returns:
        str o None: El primer ID de UniProt encontrado para el nombre del gen y el ID de taxonomía dados, o None si no se encuentra ningún ID o si ocurre un error.
    Raises:
        None: Esta función no lanza ninguna excepción, pero imprime mensajes de error en la consola.
    """
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene_name}+AND+taxonomy_id:{taxonomy_id}&format=tsv&fields=accession"
    response = requests.get(url)
    if response.status_code == 200:
        uniprot_ids = response.text.split('\n')[1:-1]  # Eliminar las primeras y últimas líneas vacías
        if uniprot_ids:
            return uniprot_ids[0]  # Devolver el primer identificador encontrado
        else:
            print(f"No se encontraron identificadores de UniProt para el gen {gene_name} y el organismo {taxonomy_id}.")
            return None
    else:
        print(f"Error al buscar el identificador de UniProt: {response.status_code}")
        return None

# Función para obtener la secuencia proteica de UniProt
def get_protein_sequence(uniprot_id):
    """
    Obtiene la secuencia proteica de la base de datos UniProt para un ID de UniProt dado.
    Args:
        uniprot_id (str): El ID de UniProt de la proteína.
    Returns:
        str: La secuencia proteica si la solicitud es exitosa.
        None: Si hay un error al obtener la secuencia.
    Raises:
        requests.exceptions.RequestException: Si hay un problema con la solicitud de red.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        # Extraer la secuencia del archivo FASTA
        lines = fasta_data.split('\n')
        sequence = ''.join(line for line in lines if not line.startswith('>'))
        return sequence
    else:
        print("Error al obtener la secuencia:", response.status_code)
        return None


