import timeit, os
from GCAM import Occurrence
from GCAM import FilesFolders
from GCAM import SignificanceTesting
from GCAM import ExpressionClustering
from GCAM import Previous_genecheck
from GCAM import plots
import pandas as pd
import logging

__author__ = 'peeyush'


def gcam_analysis(args, outpath, resource_path):
    '''
    Main GCAM function.
    :param args:
    :param resource_path:
    :return: cell occurrence dataframe
    '''
    outdir = FilesFolders.create_folders(outpath)
    logging.basicConfig(filename=os.path.join(outdir, 'GCAM.log'), level=logging.INFO)
    logging.info('Started')
    resource_path = resource_path
    subcommand = args['subcommand']
    tstart = timeit.default_timer()

    # Reading require databases
    logging.info('Reading required resource dbs..')
    warnings(args)
    write_parameter(args, outdir)

    if subcommand == 'genebased':
        try:
            genenames = args['genelist']
        except:
            print("Reading genenames from:", args['path'])
            genenames = FilesFolders.get_genes(args['path'])
        gene_based(args, resource_path, genenames, outdir)

    if subcommand == 'exprbased':
        expressiondf = args['exppath']
        expressiondf.index = expressiondf['SYMBOL']
        expressiondf = expressiondf.drop('SYMBOL', axis=1)
        print(expressiondf)
        pheno_data = args['phenopath']

        if args['controlsample'] is not None and args['controlsample'] not in list(pheno_data['phenotype']):
            logging.error(args['controlsample']+": control sample name not found in phenotype list")
            raise KeyError(args['controlsample']+": control sample name not found in phenotype list")
        genenames, newexprdf = expr_based(outdir, expressiondf, pheno_data, args)

        significance_Df = gene_based(args, resource_path, genenames, outdir)
        plotdf = ExpressionClustering.exprdf4plot(significance_Df, newexprdf, pheno_data, args,
                                                  control=args['controlsample'], path=outdir, clusterSize=int(args['celltypeClusterSize']))

    tstop = timeit.default_timer()
    print ('Total time elapsed: ' + str(tstop - tstart) + ' sec')
    logging.info('Total time elapsed: ' + str(tstop - tstart) + ' sec')
    logging.info('Finished')
    FilesFolders.zipdir(outdir) # write zip dir of output folder for download
    return outdir


def write_parameter(args, outdir):
    '''
    Write parameters in a txt file.
    :param args:
    :param subcommand:
    :param outdir:
    :return:
    '''
    if args['subcommand'] == 'exprbased':
        parameter = open(os.path.join(outdir, 'parameter.txt'), 'w')
        parameter.write("Analysis type: "+args['subcommand']+"\n")
        parameter.write("Som grid size: "+str(args['som_gridsize'])+"\n")
        parameter.write("Minimum no of genes for cell fraction analysis: "+str(args['celltypeClusterSize'])+"\n")
        parameter.write("Regression method used: nuSVR \n")
        if args['meanAsControl'] == 'sample':
            parameter.write("Control sample name: "+args['controlsample']+"\n")
        else:
            parameter.write("Consider mean of all samples as reference")
        parameter.close()


def warnings(args):
    '''
    Gives warning based on paramenter collection.
    :return:
    '''
    subcommand=args['subcommand']
    if subcommand == 'exprbased':
        userCelltype = args['selectCelltypes']
        if userCelltype is not None:
            userCelltype = [x.strip() for x in userCelltype.split(',')]
            if len(userCelltype) > 26:
                raise ValueError("Selected cell-types for analysis should be less than 26 and more than 0.")


def gene_based(args, resource_path, genenames, outdir):
    synonym = args['synonym']
    organism = args['org']
    genenames = genenames
    primarygene = genenames
    binom_prob = FilesFolders.read_binom_prob(resource_path)
    pd.DataFrame(genenames, columns=['genenames']).to_csv(outdir + os.path.sep + 'input_gene_list.txt', sep='\t', encoding='utf-8', index=False)

    ocstart = timeit.default_timer()
    if synonym:
        geneSyn = FilesFolders.gene_synonym(resource_path, organism)
        genenames = Occurrence.gene2synonym(genenames, geneSyn)
        print ('Gene count after synonym:' + str(len(genenames)))
    cellOccu = Previous_genecheck.occurrence_df(genenames, resource_path) # subquery is deprecated
    #cellOccu = Previous_genecheck.joincellsynonym(occuDF, resource_path)
    if synonym:
        cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)

    # Reduced celltypes
    if args['key_celltype_list']:
        key_celltypes = FilesFolders.key_celltypes(resource_path)
        cellOccu = cellOccu[key_celltypes]
        #print('-------------------------------------')
    #print(cellOccu.head())

    ocstop = timeit.default_timer()
    logging.info("TC in occurrence analysis:"+str(ocstop - ocstart)+'sec')

    # Scale df for heatmap and do further analysis
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu, binom_prob, resource_path, outdir, args)
    significanceDF.fisher_occurrence_test()
    write_result(significanceDF, outdir, args)
    return significanceDF


def expr_based(outdir, expressiondf, pheno_data, args):
    # Expression analysis of celltype
    return ExpressionClustering.SOMclustering(expressiondf, pheno_data, outdir, float(args['som_foldifference']), iteration=int(args['somiter']))


def write_result(significanceDF, outdir, args):
    '''
    Print all the output for genebased analysis.
    :param significanceDF:
    :return:
    '''
    cellgenedf = significanceDF.cellgenedf  # [significanceDF.cellgenedf['p-val'] < 0.05]
    cellgenedf.sort_values('p-val', ascending=True)
    if len(cellgenedf)>0:
        filtereddf = filter_df(cellgenedf)
        #filtereddf.to_excel(os.path.join(outdir, 'GCAM_sigenes.xlsx'), index=False)
        filtereddf.to_csv(os.path.join(outdir, 'GCAM_sigenes.tsv'), sep='\t', index=False)
    else: print('No significant genes for celltype')

    significanceDF.heatmapdf_create()
    significanceDF.plot_heatmap(outdir)
    significanceDF.heatmapdf.to_csv(os.path.join(outdir, 'GCAM_heatmap_df.tsv'), sep='\t')
    significanceDF.occurrencedf.to_csv(os.path.join(outdir, 'GCAM_occurence_df.tsv'), sep='\t')

    sigCelltypedf = significanceDF.sigCelltypedf
    if sigCelltypedf is not None:
        #sigCelltypedf = sigCelltypedf[sigCelltypedf['binom_pval'] < 0.05]
        significanceDF.data4radarplot()

        sigCelltypedf = sigCelltypedf.sort_values('genecluster', ascending=False)
        sigCelltypedf.index = range(len(sigCelltypedf))
        sigCelltypedf = sigCelltypedf[:20]
        sigCelltypedf = sigCelltypedf.sort_values('genecluster', ascending=True)
        plots.plot_celltypesignificance(outdir, sigCelltypedf)
        #sigCelltypedf.to_excel(os.path.join(outdir, 'GCAM_sigCelltypes.xlsx'), index=False)
        sigCelltypedf.to_csv(os.path.join(outdir, 'GCAM_sigCelltypes.tsv'), sep='\t', index=False)
    else: print('No significant celltypes')


def filter_df(df):
    '''
    Filter sigenes df to print only 3 most significant gene-celltype entries.
    :param df:
    :return:
    '''
    new_df = pd.DataFrame(columns=df.columns)
    df_gr = df.groupby('celltype')
    for k, gr in df_gr:
        new_gr = gr.sort_values('p-val', ascending=True)
        for i, r in new_gr.iterrows():
            if r['p-val'] <= 0.05:
                new_df = new_df.append(r)
    return new_df