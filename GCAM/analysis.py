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
        genenames = args['genelist']
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
    cellSyn = FilesFolders.cell_synonym(resource_path)
    binom_prob = FilesFolders.read_binom_prob(resource_path)
    pd.DataFrame(genenames, columns=['genenames']).to_csv(outdir + os.path.sep + 'input_gene_list.txt', sep='\t', encoding='utf-8', index=False)

    ocstart = timeit.default_timer()
    if synonym:
        geneSyn = FilesFolders.gene_synonym(resource_path, organism)
        genenames = Occurrence.gene2synonym(genenames, geneSyn)
        print ('Gene count after synonym:' + str(len(genenames)))
    occuDF = Previous_genecheck.occurrence_df(genenames, resource_path) # subquery is deprecated
    cellOccu = Occurrence.joincellsynonym(occuDF, cellSyn)
    if synonym:
        cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)

    # Reduced celltypes
    if args['key_celltype_list']:
        key_celltypes = FilesFolders.key_celltypes(resource_path)
        cellOccu = cellOccu[cellOccu['celltype'].isin(key_celltypes)]

    # print ('size of new df', len(cellOccu))
    cellOccu = cellOccu.set_index(cellOccu['celltype'])
    cellOccu = cellOccu.drop(['celltype'], axis=1)
    ocstop = timeit.default_timer()
    logging.info("TC in occurrence analysis:"+str(ocstop - ocstart)+'sec')

    # Scale df for heatmap and do further analysis
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu, binom_prob, resource_path, outdir)
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

    significanceDF.heatmapdf_create(thres=(20,20))
    significanceDF.plot_heatmap(outdir)
    sigCelltypedf = significanceDF.sigCelltypedf
    if sigCelltypedf is not None:
        sigCelltypedf = significanceDF.sigCelltypedf[significanceDF.sigCelltypedf['FDR'] < 0.05]
        significanceDF.data4radarplot()
        sigCelltypedf.sort_values('genecluster', ascending=True)
        plots.plot_celltypesignificance(outdir, sigCelltypedf, args)
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
    df_gr = df.groupby('gene')
    for k, gr in df_gr:
        new_gr = gr.sort_values('p-val', ascending=True)
        rows = 1
        for i, r in new_gr.iterrows():
            if rows <= 3:
                new_df = new_df.append(r)
                rows += 1
    return new_df