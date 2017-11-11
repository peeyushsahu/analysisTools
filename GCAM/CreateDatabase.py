__author__ = 'sahu'



class BuildGcamDatabase:
    '''
    This will create occurrence database for all the genes.
    '''
    def __init__(self, synonym_df, pubmed_db):
        self.synonym_df = synonym_df
        self.pubmed_db = pubmed_db
        self.gene_synonym_dict = None
        self.gene_2_pmid_dict = None

    def combine_gene_synonym(self):
        '''
        Combing synonyms for genes.
        ex. {'A1BG': ['A1B', 'ABG', 'GAB', 'HYST2477']}
        :return:
        '''
        gene_synonym_dict = {}
        synonym_df = self.synonym_df
        for ind, row in synonym_df.iterrows():
            gene = row['Symbol']
            synonym = row['Synonyms'].split('|')
            if len(synonym) == 1 and synonym[0] == '-':
                gene_synonym_dict[gene] = [gene]
            else:
                synonym.append(gene)
                gene_synonym_dict[gene] = synonym
        self.gene_synonym_dict = gene_synonym_dict
        print(gene_synonym_dict)

    def get_pmids(self, gene):
        '''
        This function fetches pmids for provided genes.
        :return:
        '''
        from Bio import Entrez
        Entrez.email = "peeyush215@gmail.com"
        data = Entrez.esearch(db="pubmed", retmax=15000, term=gene)
        res = Entrez.read(data)
        PMID = res["IdList"]
        #print 'Pubmed ids for '+gene+':', len(PMID)
        #print('Gene:', gene, ' Lenght of pmid:', len(PMID))
        return PMID

    def build_gene_2_pmid(self):
        '''
        Creating list of PMIDS associated with gene and its synonym.
        ex. {'A1BG': ['29103060', '28652196', '27132531', .....]}
        :return:
        '''
        gene_2_pmid_dict = {}
        gene_synonym_dict = self.gene_synonym_dict
        for gene, synonym in gene_synonym_dict.items():
            pmid_list = []
            for i in synonym:
                ids = self.get_pmids(i)
                pmid_list = pmid_list + ids
            print("PMIDS for gene", gene, ":", len(pmid_list))
            gene_2_pmid_dict[gene] = list(set(pmid_list))
        self.gene_2_pmid_dict = gene_2_pmid_dict
        print(gene_2_pmid_dict)

    def get_pmid_pos(self, pubmed_db):
        celltype_list = []
        for pmid in self.pmids:
            if int(pmid) in pubmed_db.index:
                celltype_list.append(pubmed_db.loc[int(pmid)][0])
        self.cellinpmid = celltype_list
