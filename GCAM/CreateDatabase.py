__author__ = 'sahu'
import GCAM.FilesFolders as read
import pandas as pd
import time
import sys


class BuildGcamDatabase:
    '''
    This will create occurrence database for all the genes.
    '''
    def __init__(self, resource_path, synonym_df):
        self.synonym_df = synonym_df
        self.pubmed_db = read.read_annotation_database(resource_path)
        self.cell_subtract_df = read.read_cell_subtractdf(resource_path)
        self.cell_name_df = read.celltype_DB(resource_path)
        self.cell_synonym = read.cell_synonym(resource_path)
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
                #print('pmids for snonym:', i, ':', len(ids))
                pmid_list = pmid_list + ids
            #print("PMIDS for gene", gene, ":", len(pmid_list))
            gene_2_pmid_dict[gene] = list(set(pmid_list))
        self.gene_2_pmid_dict = gene_2_pmid_dict
        #print(gene_2_pmid_dict)

    def get_pmid_pos(self, pmid_list):
        celltype_list = []
        pubmed_db = self.pubmed_db
        for pmid in pmid_list:
            if int(pmid) in pubmed_db.index:
                cells = pubmed_db.loc[int(pmid)][0]
                cells = cells.split('%%')
                celltype_list = celltype_list + cells
                #print(cells)
        return celltype_list

    def get_occurrence(self):
        '''
        Calculate celltype occurrence for each gene.
        :param genes_dict:
        :param cellDB:
        :return:
        '''
        gene2cell_occur = {}
        for gene, pmids in self.gene_2_pmid_dict.items():
            celloccu = dict.fromkeys(self.cell_name_df, 0)
            celltype = self.get_pmid_pos(pmids)
            for cells in self.cell_name_df:
                for found in celltype:
                    if cells in found.lower():
                        occu = found.lower().count(cells)
                        celloccu[cells] += occu
            gene2cell_occur[gene] = celloccu
        #print(gene2cell_occur)
        gene2cell_occur = self.subtract_cellnamepeat(gene2cell_occur)
        print("########################################################################")
        #print(gene2cell_occur)
        gene2cell_occur = self.joincellsynonym(gene2cell_occur)
        #return {self.gene: celloccu}

    def subtract_cellnamepeat(self, gene2cell_occur):
        '''
        This will subtract count from ex. t lymphocyte that belongs to cd4 t cells
        :return:
        '''
        subtract_df = self.cell_subtract_df
        for gene, occur in gene2cell_occur.items():
            for k, v in subtract_df.iterrows():
                for cell in v['subs'].split(','):
                    occur[v['celltype']] = occur[v['celltype']] - occur[cell.lower()]
            #print (celloccu)
        return gene2cell_occur

    def joincellsynonym(self, gene2cell_occur):
        '''
        Join multiple cell synonym to one.
        :param celloccu:
        :param cellSyn:
        :return:
        '''
        cellSyn = self.cell_synonym
        #print (celloccu)
        for gene, celldict in gene2cell_occur.items():
            for k, v in cellSyn.iterrows():
                index = v['cell'].lower()
                #print(gene)
                #print(celldict)
                #print(v['synonyms'].split(','))
                for cell in v['synonyms'].split(','):
                    indexsyn = cell.lower()
                    ## adding synonym
                    #print(indexsyn)
                    celldict[index] = celldict[index] + celldict[indexsyn]
                    celldict.pop(indexsyn)
        print('#################################################################')
        print(gene2cell_occur)
        return gene2cell_occur


class BuildGcamDatabaseInOneGo:
    '''
    This will create occurrence database for all the genes.
    '''
    def __init__(self, resource_path, synonym_df):
        self.synonym_df = synonym_df
        self.pubmed_db = read.read_annotation_database(resource_path)
        self.cell_subtract_df = read.read_cell_subtractdf(resource_path)
        self.cell_name_df = read.celltype_DB(resource_path)
        self.cell_synonym = read.cell_synonym(resource_path)
        self.gene_synonym_dict = None
        self.gene_2_pmid_dict = None

    def combine_gene_synonym(self):
        '''
        Combing synonyms for genes.
        ex. {'A1BG': ['A1B', 'ABG', 'GAB', 'HYST2477']}
        :return:
        '''
        gene_synonym_dict = {}
        cellOccurDict = {}
        genes_left = []
        synonym_df = self.synonym_df
        index = 0

        for ind, row in synonym_df.iterrows():
            index = ind
            gene_synonym_dict = {}
            gene = row['Symbol']
            #print(gene)
            sys.stdout.write("\r%d %s" % (ind, gene))
            sys.stdout.flush()
            synonym = row['Synonyms'].split('|')
            gene_synonym_dict[gene] = [gene]

            if len(synonym) == 1 and synonym[0] == '-':
                gene_synonym_dict[gene] = [gene]
            else:
                synonym.append(gene)
                gene_synonym_dict[gene] = synonym

            # Associate gene names to pmids
            pmids=None
            try:
                pmids = self.build_gene_2_pmid(gene_synonym_dict)

            except Exception as e:
                print(e)
                time.sleep(5)
                genes_left.append(index)
                pass
            # Calculate celtype occurrence from associated abstract
            cellOccur = self.get_occurrence(pmids)
            cellOccurDict.update(cellOccur)

        #self.gene_synonym_dict = gene_synonym_dict
        #print(gene_synonym_dict)
        celloccudf = pd.DataFrame.from_dict(cellOccurDict, orient='index')
        celloccudf.to_csv('/ps/imt/resources/celltype_occu_widout_synonym_DB', sep='\t', header=True)
        return cellOccurDict, genes_left

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

    def build_gene_2_pmid(self, gene_synonym_dict):
        '''
        Creating list of PMIDS associated with gene and its synonym.
        ex. {'A1BG': ['29103060', '28652196', '27132531', .....]}
        :return:
        '''
        gene_2_pmid_dict = {}
        for gene, synonym in gene_synonym_dict.items():
            pmid_list = []
            for i in synonym:
                ids = self.get_pmids(i)
                #print('pmids for snonym:', i, ':', len(ids))
                pmid_list = pmid_list + ids
            #print("PMIDS for gene", gene, ":", len(pmid_list))
            gene_2_pmid_dict[gene] = list(set(pmid_list))
        return gene_2_pmid_dict
        #print(gene_2_pmid_dict)

    def get_pmid_pos(self, pmid_list):
        '''
        Find PMID in pubmed_db and find annotated celltype
        :param pmid_list:
        :return:
        '''
        celltype_list = []
        pubmed_db = self.pubmed_db
        for pmid in pmid_list:
            if int(pmid) in pubmed_db.index:
                cells = pubmed_db.loc[int(pmid)][0]
                cells = cells.split('%%')
                celltype_list = celltype_list + cells
                #print(cells)
        return celltype_list

    def get_occurrence(self, gene_2_pmid_dict):
        '''
        Calculate celltype occurrence for each gene.
        :param genes_dict:
        :param cellDB:
        :return:
        '''
        gene2cell_occur = {}
        for gene, pmids in gene_2_pmid_dict.items():
            celloccu = dict.fromkeys(self.cell_name_df, 0)
            celltype = self.get_pmid_pos(pmids)
            for cells in self.cell_name_df:
                for found in celltype:
                    if cells in found.lower():
                        occu = found.lower().count(cells)
                        celloccu[cells] += occu
            gene2cell_occur[gene] = celloccu
        #print(gene2cell_occur)
        gene2cell_occur = self.subtract_cellnamepeat(gene2cell_occur)
        #print("########################################################################")
        #print(gene2cell_occur)
        gene2cell_occur = self.joincellsynonym(gene2cell_occur)
        return gene2cell_occur

    def subtract_cellnamepeat(self, gene2cell_occur):
        '''
        This will subtract count from ex. t lymphocyte that belongs to cd4 t cells
        :return:
        '''
        subtract_df = self.cell_subtract_df
        for gene, occur in gene2cell_occur.items():
            for k, v in subtract_df.iterrows():
                for cell in v['subs'].split(','):
                    occur[v['celltype']] = occur[v['celltype']] - occur[cell.lower()]
            #print (celloccu)
        return gene2cell_occur

    def joincellsynonym(self, gene2cell_occur):
        '''
        Join multiple cell synonym to one.
        :param celloccu:
        :param cellSyn:
        :return:
        '''
        cellSyn = self.cell_synonym
        #print (celloccu)
        for gene, celldict in gene2cell_occur.items():
            for k, v in cellSyn.iterrows():
                index = v['cell'].lower()
                #print(gene)
                #print(celldict)
                #print(v['synonyms'].split(','))
                for cell in v['synonyms'].split(','):
                    indexsyn = cell.lower()
                    ## adding synonym
                    #print(indexsyn)
                    celldict[index] = celldict[index] + celldict[indexsyn]
                    celldict.pop(indexsyn)
        #print('#################################################################')
        #print(gene2cell_occur)
        return gene2cell_occur


def background_prob_4_celltype_occurence_per_gene():
    '''
    Calculate background probability for celltype occurrrence in all genes.
    This will be used in calculating significance for significance of celltype in all gene.
    :return:
    '''
    occudf = pd.read_csv('/home/peeyush/Desktop/gcam_test_data/resources/celltype_occu_widout_synonym_DB', header=0, sep='\t')
    print(occudf.head())
    occudf.index = occudf['Unnamed: 0']
    occudf = occudf.drop('Unnamed: 0', axis=1)
    empty_gene = occudf.sum(axis=1)
    non_empty_gene = empty_gene.index[empty_gene!=0].tolist()
    filtered_occur_df = occudf[occudf.index.isin(non_empty_gene)]
    celltype_background_prob = {}

    # populating background prob dict
    for ind, col in filtered_occur_df.iteritems():
        #print(ind)
        #print(len(col[col>0]))
        celltype_background_prob[ind] = {'found_in':len(col[col>0]),
                                         'total_gene':len(col),
                                         'background_prob':len(col[col>0])/len(col)}
    pd.DataFrame(celltype_background_prob).T.to_csv('/home/peeyush/Desktop/gcam_test_data/resources/celltype_occu_binom_prob',
                                                    header=True, sep='\t')


def create_db_4_synonym_2_symbol_check():
    gene_synonym = pd.read_csv('/home/peeyush/Desktop/gcam_test_data/resources/Homo_sapiens.gene_info', header=0,
                               index_col=None, sep='\t')
    gene_synonym_2_symbol = []
    for ind, col in gene_synonym.iterrows():
        genesymbol = col['Symbol']
        geneid = col['GeneID']
        gene_synonym_2_symbol.append({'synonym':genesymbol, 'geneID':geneid, 'symbol':genesymbol})
        synonyms = col['Synonyms'].split('|')
        if len(synonyms) >= 1 and synonyms[0] != '-':
            for synonym in synonyms:
                gene_synonym_2_symbol.append({'synonym':synonym, 'geneID':geneid, 'symbol':genesymbol})
    pd.DataFrame(gene_synonym_2_symbol).to_csv('/home/peeyush/Desktop/gcam_test_data/resources/map_synonym_2_symbol_DB',
                                                    header=True, sep='\t')