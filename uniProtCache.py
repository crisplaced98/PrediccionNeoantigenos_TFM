import getInformation

class UniProtCache:
    def __init__(self):
        self.uniprot_cache = {}
        self.sequence_cache = {}

    def cached_get_uniprot_info(self, gene):
        """
        Recupera información de UniProt para un gen dado, utilizando una caché para minimizar llamadas redundantes a la API.
        Args:
            gene (str): El símbolo del gen para el cual se desea recuperar información de UniProt.
        Returns:
            tuple: Una tupla que contiene el ID de UniProt y la secuencia de la proteína. Si no se encuentra el ID de UniProt,
               el segundo elemento de la tupla será None.
        """
        if gene not in self.uniprot_cache:
            self.uniprot_cache[gene] = getInformation.get_uniprot_id(gene)
        uniprot_id = self.uniprot_cache[gene]
        
        if uniprot_id:
            if uniprot_id not in self.sequence_cache:
                    self.sequence_cache[uniprot_id] = getInformation.get_protein_sequence(uniprot_id)
            return uniprot_id, self.sequence_cache[uniprot_id]
        else:
            return uniprot_id, None
