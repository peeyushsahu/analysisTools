__author__ = 'peeyush'
from pandas import read_csv
import os, sys, stat
import time


def create_folders(path):
    print('Output directory created: ' + os.path.join(path, 'GCAM_output_'+str(time.strftime("%d_%m_%Y"))+'_'+str(time.strftime("%H_%M_%S"))))
    npath = os.path.join(path, 'GCAM_output_'+str(time.strftime("%d_%m_%Y"))+'_'+str(time.strftime("%H_%M_%S")))
    if not os.path.exists(npath):
        os.makedirs(npath)
        #print('Changed permission')
        os.chmod(npath, 0o777)
        #os.chmod(npath, mode=777)
    return npath


def zipdir(dirPath=None, zipFilePath=None, includeDirInZip=False):
    import zipfile
    """Create a zip archive from a directory.

    Note that this function is designed to put files in the zip archive with
    either no parent directory or just one parent directory, so it will trim any
    leading directories in the filesystem paths and not include them inside the
    zip archive paths. This is generally the case when you want to just take a
    directory and make it into a zip file that can be extracted in different
    locations.

    Keyword arguments:

    dirPath -- string path to the directory to archive. This is the only
    required argument. It can be absolute or relative, but only one or zero
    leading directories will be included in the zip archive.

    zipFilePath -- string path to the output zip file. This can be an absolute
    or relative path. If the zip file already exists, it will be updated. If
    not, it will be created. If you want to replace it from scratch, delete it
    prior to calling this function. (default is computed as dirPath + ".zip")

    includeDirInZip -- boolean indicating whether the top level directory should
    be included in the archive or omitted. (default True)

"""
    if not zipFilePath:
        zipFilePath = dirPath + ".zip"
    if not os.path.isdir(dirPath):
        raise OSError("dirPath argument must point to a directory. "
            "'%s' does not." % dirPath)
    parentDir, dirToZip = os.path.split(dirPath)
    #Little nested function to prepare the proper archive path
    def trimPath(path):
        archivePath = path.replace(parentDir, "", 1)
        if parentDir:
            archivePath = archivePath.replace(os.path.sep, "", 1)
        if not includeDirInZip:
            archivePath = archivePath.replace(dirToZip + os.path.sep, "", 1)
        return os.path.normcase(archivePath)

    outFile = zipfile.ZipFile(zipFilePath, "w",
        compression=zipfile.ZIP_DEFLATED)
    for (archiveDirPath, dirNames, fileNames) in os.walk(dirPath):
        for fileName in fileNames:
            filePath = os.path.join(archiveDirPath, fileName)
            outFile.write(filePath, trimPath(filePath))
        #Make sure we get empty directories as well
        if not fileNames and not dirNames:
            zipInfo = zipfile.ZipInfo(trimPath(archiveDirPath) + "/")
            outFile.writestr(zipInfo, "")
    outFile.close()


def read_annotation_database(path):
    '''
    Reading annotation DB as pandas dataframe
    :param path:
    :return:
    '''
    try:
        annoDB = read_csv(path + os.path.sep + 'pmid_celltype_index_final.txt', header=0, sep="\t")
    except Exception as e:
        print(e)
        raise ValueError("annotation db doesnot exist.")
    annoDB = annoDB.set_index(['pmid'])
    return annoDB


def celltype_DB(path):
    '''
    Import celltype database.
    :param path:
    :return:
    try:
        cellDB = read_csv(path + os.path.sep + 'cellTypes.csv', header=None, sep=',')
    except:
        raise ValueError("celltype db does not exist.")
    cellDB.columns = ['celltype']
    cellDB['celltype'] = cellDB['celltype'].str.lower()
    #print cellDB
    return cellDB
    '''
    celltypeList = []
    try:
        with open(path + os.path.sep + 'cellTypes.txt') as file:
            for celltype in file:
                celltype = celltype.rstrip()
                if len(celltype) > 0:
                    celltypeList.append(celltype.lower())
    except Exception as e:
        print(e)
        raise ValueError("celltype db does not exist.")
    celltypeList = list(set(celltypeList))
    #print (celltypeList)
    return celltypeList


def cell2tissue_DB(path):
    '''
    Import celltype database.
    :param path:
    :return:
    '''
    try:
        cell2tissue = read_csv(os.path.join(path, 'cell2tissue.txt'), header=0, sep='\t', index_col=0)
    except:
        raise ValueError("cell2tissue db does not exist.")
    #print cellDB
    return cell2tissue


def cell_synonym(path):
    '''
    Import cell synonym database.
    :param path:
    :return:
    '''
    try:
        cellSyn = read_csv(path + os.path.sep + 'cell_type_synonyms_python.csv', header=0, sep=',')
    except:
        raise ValueError("Synonym db does not exist.")
    return cellSyn

def get_genes(path):
    '''
    Read input gene list.
    :param path:
    :return:
    '''
    geneList = []
    try:
        with open(path) as file:
            for gene in file:
                gene = gene.strip()
                if len(gene) > 0:
                    geneList.append(gene.lower())
    except IOError:
        print ("Error: Query file does not appear to exist.")
        sys.exit(1)
    f_genelist = list(set(geneList))
    print ('Size of user provided gene list:' + str(len(geneList)))
    print ('No. of genes after removing duplicates:' + str(len(f_genelist)))
    if len(f_genelist) < 5:
        raise IOError("Minimum gene requirement for the analysis is 5.")
    return f_genelist

def key_celltypes(path):
    '''
    Read input gene list.
    :param path:
    :return:
    '''
    celltypeList = []
    path = path + os.path.sep + 'key_cellTypes.csv'
    try:
        with open(path) as file:
            for cells in file:
                cells = cells.strip()
                celltypeList.append(cells.lower())
    except IOError:
        print ("Error: Key celltype db does not exist.")
        sys.exit(1)
    f_celltypes = list(set(celltypeList))
    print ('Size of key gene list:' + str(len(celltypeList)))
    #print 'No. of genes after removing duplicates:', len(f_genelist)
    return f_celltypes

def gene_synonym(path, organism):
    '''
    Reads synonym file for genes
    :param path:
    :param organism:
    :return:
    '''
    try:
        if organism == 'human':
            geneSyn = read_csv(path + os.path.sep + 'Human_synonym.txt', header=None, sep='\t')
        elif organism == 'mouse':
            geneSyn = read_csv(path + os.path.sep + 'Mouse_synonym.txt', header=None, sep='\t')
    except:
        raise ValueError("Synonym db does not exist.")
    geneSyn.columns = ['gene', 'synonym']
    geneSyn['gene'] = geneSyn['gene'].str.lower()
    geneSyn = geneSyn.set_index(geneSyn['gene'])
    return geneSyn


def read_expression_file(path):
    '''
    Reads expression data for analysis
    :param path:
    :return:
    '''
    try:
        # print ('exprs path:'+ path)
        expressiondf = read_csv(path, header=0, sep="\t", index_col=0)
        # print(expressiondf.head(1))
        expressiondf.index = expressiondf.index.to_series().str.lower()
        expressiondf = expressiondf[expressiondf.index.to_series().str.len() >= 3] ##Removing gene names with only 2 letters.
    except IOError:
        print ("Error: Expression File does not appear to exist.")
        sys.exit(0)
    return expressiondf


def read_previous_occurrence_table(resource_path):
    '''
    This will read the occurrence database for already analysed genes to save time.
    :param resource_path:
    :return:
    '''
    try:
        print ('reading previously analysed genes')
        gene_occu_db = read_csv(os.path.join(resource_path, 'gene_occu_db.tsv'), header=0, sep="\t", index_col=0)
    except:
        print ("Warning: Creating Gene occurrence db, analysis will take longer.")
        return None, False
    return gene_occu_db, True


def read_occurrence_table(resource_path, synonym=False):
    '''
    This will read the occurrence database for already analysed genes to save time.
    :param resource_path:
    :return:
    '''
    if synonym:
        print('reading celltype occurrence table with synonym...')
        gene_occu_db = read_csv(os.path.join(resource_path, 'celltype_occu_with_synonym_DB'), header=0, sep="\t", index_col=0)
    else:
        print('reading celltype occurrence table...')
        gene_occu_db = read_csv(os.path.join(resource_path, 'celltype_occu_widout_synonym_DB'), header=0, sep="\t", index_col=0)
    return gene_occu_db


def read_binom_prob_occu(resource_path):
    '''
    This background probability is for occurrence for a celltype in all celltype occurrence.
    :param resource_path:
    :return:
    '''
    try:
        #print ('binomial probability calculation...')
        cell_binom = read_csv(os.path.join(resource_path, 'cellType_binom_prob.txt'), header=0, sep="\t")
    except:
        raise ValueError('Binomial probabilities not found')
    return cell_binom


def read_binom_prob_occu_genes(resource_path):
    '''
    This background probability is for occurrence for celltype in all genes .
    :param resource_path:
    :return:
    '''
    try:
        #print ('binomial probability calculation...')
        cell_binom = read_csv(os.path.join(resource_path, 'celltype_occu_binom_prob'), header=0, sep="\t")
    except:
        raise ValueError('Binomial probabilities not found')
    return cell_binom


def read_genename_2_genesymbol_map(resource_path):
    '''

    :return:
    '''
    try:
        #print ('binomial probability calculation...')
        name2symbol = read_csv(os.path.join(resource_path, 'map_synonym_2_symbol_DB'), header=0, sep="\t", index_col=0)
        print(name2symbol.head())
    except:
        raise ValueError('Gene name symbol DB: map_synonym_2_symbol_DB not found')
    return name2symbol


def read_pheno_data(path):
    '''
    Reads pheno-data for analysis
    :param path:
    :return:
    '''
    try:
        #print ('pheno path:'+ path)
        phenodf = read_csv(path, header=0, sep="\t")
        if not 'phenotype' in phenodf.columns or not 'sample' in phenodf.columns:
            print ('Error phenotype: please name columns as phenotype and sample')
            sys.exit(0)
    except IOError:
        print ("Error: pheno File does not appear to exist.")
        sys.exit(0)
    return phenodf

def read_cell_subtractdf(path):
    '''
    Reads pheno-data for analysis
    :param path:
    :return:
    '''
    try:
        cellsub_df = read_csv(os.path.join(path, 'celltype_subtract.txt'), header=0, sep="\t")
        #print cellsub_df
    except IOError:
        print ("Error: cell_subtract db does not appear to exist.")
        sys.exit(0)
    return cellsub_df