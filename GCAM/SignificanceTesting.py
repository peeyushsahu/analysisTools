from GCAM.plots import HiearchicalHeatmap
__author__ = 'peeyush'
import pandas as pd
import scipy.stats as stats
import os
import logging, timeit
from GCAM import FilesFolders as read
from GCAM import plots
import seaborn as sns
import statsmodels.sandbox.stats.multicomp as statsmodels


class SignificanceObject:
    '''
    This Object will hold dataframes used and created in analysis
    '''
    def __init__(self, occurrencedf, binom_prob, resourcepath, outdir, args, heatmapdf=None):
        self.outdir = outdir
        self.resourcepath = resourcepath
        self.occurrencedf = occurrencedf
        self.binom_prob = binom_prob
        self.binom_prob_occu = read.read_binom_prob_occu_genes(resourcepath)
        self.heatmapdf = heatmapdf
        self.args = args
        self.filheatmapdf = None
        self.pvaldf = None
        self.adjpvaldf = None
        self.cellgenedf = None
        self.sigCelltypedf = None
        self.binom_pval_df = None

    def heatmapdf_create(self):
        '''
        This function will generate df for HeatMap by scaling the occurrencedf
        '''
        occurrencedf = self.occurrencedf
        scaled_df = scale_dataframe(occurrencedf)
        scaled_df.columns = occurrencedf.columns
        scaled_df = scaled_df.set_index(occurrencedf.index)
        self.heatmapdf = scaled_df
        self.filter_heatmapdf()

    def filter_heatmapdf(self):
        '''
        This method filters rows and columns with sum < 1 in HeatMapdf
        '''
        df = self.heatmapdf
        filheatmapdf = df.loc[df.sum(1) > 10, df.sum(0) > 100]
        #filheatmapdf.index = range(len(filheatmapdf))
        #print(filheatmapdf.head())
        self.filheatmapdf = filheatmapdf

    def plot_heatmap(self, path):
        '''
        This method will plot the HeatMap dataframe. Using package HclustHeatMap.
        :param HiearchicalHeatmap:
        :param df:
        :param path:
        :return:
        '''
        #hclustHeatmap = HiearchicalHeatmap()
        #hclustHeatmap.frame = self.filheatmapdf
        #hclustHeatmap.path = os.path.sep.join([path, 'GCAM_heatMap.svg'])
        #fig, axm, axcb, cb = hclustHeatmap.plot()
        dataframe = self.filheatmapdf
        path = os.path.sep.join([path, 'GCAM_heatMap.svg'])
        sns.set(context="talk")
        cmap = sns.diverging_palette(250, 10, as_cmap=True)
        if len(dataframe) > 30:
            ax = sns.clustermap(dataframe, cmap=cmap, yticklabels=False)
            ax.ax_heatmap.set_ylabel('Genes')
        else:
            ax = sns.clustermap(dataframe, cmap=cmap)
            for ticks in ax.ax_heatmap.get_yticklabels():
                ticks.set_rotation(0)
        ax.savefig(path)

    def filter_occuDf(self):
        '''
        Filter occurrence df for removing gene wid less than 5 celltype tags
        :return:
        '''
        #print('########### Shape of occu df #############')
        #print(self.occurrencedf.shape)
        occuDf = self.occurrencedf
        Columns=[]
        for k, v in occuDf.iteritems():
            if v.sum() < 5:
                Columns.append(k)
        #print(len(Columns))
        self.occurrencedf = occuDf.drop(Columns, axis=1)
        '''
        rows=[]
        for k, r in occuDf.iterrows():
            if r.sum() < 5:
                rows.append(k)
        print(len(rows))
        self.occurrencedf = occuDf.drop(rows, axis=0)
        '''
        #print(self.occurrencedf.shape)

    def pergene_celltype_occurrence_test(self):
        '''
        This method will calculate significance of celltypes per gene using their occurrence.
        Statistical test used is Fisher Exact Test
        '''
        fsstart = timeit.default_timer()
        self.filter_occuDf()
        occu_df = self.occurrencedf
        binom_prob = self.binom_prob
        if self.args['key_celltype_list']:
            key_celltypes = read.key_celltypes(self.resourcepath)
            binom_prob = binom_prob[binom_prob['celltype'].isin(key_celltypes)]
            binom_prob.index = range(len(binom_prob))
            self.binom_prob = binom_prob

        pvaldf = pd.DataFrame()
        binom_pvaldf = pd.DataFrame()
        enrichmentdf = pd.DataFrame()
        #print(occu_df.head())
        for k, v in occu_df.iterrows():
            key = v.keys()
            rowsum = v.sum()
            for i in range(0, v.shape[0]):
                value = v[i]
                if not value == 0:
                    #print(k, v, rowsum, value, i)
                    enrichment = float(value)/occu_df.iloc[:, i].sum()
                    celltype = key[i]
                    index = binom_prob[binom_prob['celltype'] == celltype].index.tolist()
                    #background_prob = binom_prob.iloc[index[0], 1] / sum(binom_prob['occurrence'])
                    background_prob = binom_prob.loc[index[0], 'background_prob']
                    #print(binom_prob.iloc[index[0], 1], sum(binom_prob['occurrence']), background_prob)
                    #print(celltype, index, value, rsum+value)
                    #print(binom_prob.head())
                    b_pval = stats.binom_test(value, rowsum, background_prob, alternative='two-sided')
                    binom_pvaldf.loc[k, key[i]] = b_pval
                    enrichmentdf.loc[k, key[i]] = enrichment
                else:
                    binom_pvaldf.loc[k, key[i]] = 1
                    enrichmentdf.loc[k, key[i]] = 0
        self.binom_pval_df = binom_pvaldf
        self.pvaldf = pvaldf
        fsstop = timeit.default_timer()
        logging.info('TC in sig occur test:'+str(fsstop-fsstart)+'sec')
        self.celltype_overrepresntation_list(enrichmentdf)

    def celltype_overrepresntation_list(self, enrichmentdf):
        '''
        This method will save the result of significance in one DF.
        '''
        significance = 1
        column = ['celltype', 'gene', 'enrichment', 'binom_pval', 'FDR']
        cellgenedf = pd.DataFrame()
        #print(self.binom_pval_df.head())
        for gene, celltype in self.binom_pval_df.iterrows():
            for cell, pval in celltype.iteritems():
                if pval < significance:
                    cellgenedf = cellgenedf.append(
                        pd.Series([cell, gene, enrichmentdf.loc[gene, cell], pval, 0]), ignore_index=True)
        #print cellgenedf.head(10)
        cellgenedf.columns = column
        cellgenedf = cellgenedf.sort_values(['celltype', 'binom_pval'], ascending=[True, True])
        cellgenedf.index = range(len(cellgenedf))

        pvals = cellgenedf['binom_pval'].values
        corr_pvals = statsmodels.multipletests(pvals=pvals, alpha=0.05, method='fdr_bh')
        #print(pvals)
        #print(corr_pvals)
        cellgenedf['FDR'] = 0
        cellgenedf['FDR'] = corr_pvals[1]
        '''
        for ind, row in cellgenedf.iterrows():
            fdr = (row['Binom p-val'] * len(cellgenedf)) / (ind + 1)
            cellgenedf.iloc[ind, 4] = fdr
        '''
        print('cellgenedf shape:', cellgenedf.shape)
        #cellgenedf = self.filter_df(cellgenedf)
        self.cellgenedf = cellgenedf
        print('cellgenedf shape after:', cellgenedf.shape)
        #self.filter_cellgenedf()  # Filter single cell multigene enrihment
        self.overall_significant_celltypes()

    def overall_significant_celltypes(self):
        '''
        This method will test the combined significance of celltype in the data and help predicts
        its association with user given data.
        '''
        sigcelltype = pd.DataFrame()
        #print self.cellgenedf
        cellgroup = self.cellgenedf.groupby(self.cellgenedf['celltype'])
        for celltype, val in cellgroup:
            #print celltype
            if len(val[val['FDR'] <= 0.05]) > 1:
                #print val
                a = len(val[val['FDR'] <= 0.05])
                sigcelltype = sigcelltype.append(pd.Series([celltype, a]), ignore_index=True)
        sigcelltype.columns = ['celltype', 'genecluster']
        #print('Sig cell type\n', sigcelltype)
        self.sigCelltypedf = sigcelltype
        self.binom_significant_celltypes()
        self.hypergeometric_significant_celltypes()

    def hypergeometric_significant_celltypes(self):
        '''
        hypergeometric test for significance of celltype enrichment.
        '''
        print('Testing celltype enrichment....')
        sigcelltype = self.sigCelltypedf
        cellgroup = self.cellgenedf.groupby(self.cellgenedf['celltype'])
        totalgenes = self.occurrencedf.shape[0]

        allsiggenes = self.cellgenedf
        allsiggenes = allsiggenes[allsiggenes['FDR'] <= 0.05]
        allsiggenes = len(set(allsiggenes['gene']))

        sigcelltype.loc[:, 'hyper_pval'] = 1
        col = sigcelltype.columns.get_loc('hyper_pval')
        for index, row in sigcelltype.iterrows():
            #print(row['celltype'])
            #print(row['genecluster'], totalgenes, len(cellgroup.get_group(row['celltype'])), allsiggenes)
            ## stats.hypergeom.sf(x, M, n, N)
            hyper_pval = stats.hypergeom.sf(row['genecluster']-1, totalgenes, allsiggenes, len(cellgroup.get_group(row['celltype'])))
            #print(hyper_pval)
            sigcelltype.iloc[index, col] = hyper_pval

        sigcelltype.loc[:, 'hyper_FDR'] = 1
        #ind_fdr = sigcelltype.columns.get_loc('hyper_FDR')
        sigcelltype = sigcelltype.sort_values('hyper_pval', ascending=True)
        sigcelltype.index = range(len(sigcelltype))

        pvals = sigcelltype['hyper_pval'].values
        corr_pvals = statsmodels.multipletests(pvals=pvals, alpha=0.05, method='fdr_bh')
        sigcelltype['hyper_FDR'] = corr_pvals[1]
        self.sigCelltypedf = sigcelltype

    def binom_significant_celltypes(self):
        '''
        Binomial test for significance of celltype enrichment.
        '''
        print('Testing celltype enrichment....')
        sigcelltype = self.sigCelltypedf
        cellgroup = self.cellgenedf.groupby(self.cellgenedf['celltype'])

        binom_prob_occu = self.binom_prob_occu

        sigcelltype.loc[:, 'binom_pval'] = 1
        col = sigcelltype.columns.get_loc('binom_pval')
        for index, row in sigcelltype.iterrows():
            #print(row['celltype'])
            #print(row['genecluster'], totalgenes, len(cellgroup.get_group(row['celltype'])), allsiggenes)

            bprob_ind = binom_prob_occu[binom_prob_occu['celltype'] == row['celltype']].index[0]
            #print(bprob_ind)
            background_prob = binom_prob_occu.loc[bprob_ind, 'background_prob']
            #print(background_prob)
            binom_pval = stats.binom_test(row['genecluster']-1, len(cellgroup.get_group(row['celltype'])), background_prob, alternative='two-sided')
            sigcelltype.iloc[index, col] = binom_pval

        sigcelltype.loc[:, 'binom_FDR'] = 1
        sigcelltype = sigcelltype.sort_values('binom_pval', ascending=True)
        sigcelltype.index = range(len(sigcelltype))

        pvals = sigcelltype['binom_pval'].values
        corr_pvals = statsmodels.multipletests(pvals=pvals, alpha=0.05, method='fdr_bh')
        #print(pvals)
        #print(corr_pvals)
        sigcelltype['binom_FDR'] = corr_pvals[1]
        self.sigCelltypedf = sigcelltype

    def data4radarplot(self):
        '''
        Preparing data for radar plot
        '''
        import numpy as np
        tissue = []
        genes = []
        sigCelltypedf = self.sigCelltypedf
        cell2tissue = read.cell2tissue_DB(self.resourcepath)
        #print(cell2tissue.head(3))
        cell2tissue['genes'] = 0
        for k, v in sigCelltypedf.iterrows():
            if v['celltype'] in list(cell2tissue.index):
                cell2tissue.loc[v['celltype'], 'genes'] = v['genecluster']

        norm_cell2gene = cell2tissue['tissue'].value_counts()
        #print(norm_cell2gene)
        df_cell_gr = cell2tissue.groupby('tissue')
        for cell in df_cell_gr.groups:
            #print(cell)
            df = df_cell_gr.get_group(cell)
            tissue.append(cell)
            genes.append((sum(df['genes'])*1.)/norm_cell2gene.loc[cell])
        gene2plot = np.divide(genes, sum(genes))
        #print(gene2plot, tissue, self.outdir)
        plots.plot_radar(gene2plot, tissue, self.outdir)


def scale(val, src, dst):
    '''
    This returns scaled value
    :param val: value to be scaled
    :param src: min and max of values to be scaled
    :param dst: range to scale
    :return:
    '''
    return (((val - src[0]) * (dst[1] - dst[0])) / (src[1] - src[0])) + dst[0]


def scale_dataframe(df):
    '''
    This will scale a Pandas.DataFrame with every row min-max
    :param df: DataFrame
    :return: Scaled DataFrame
    '''
    new_max = 100
    new_min = 0
    list_of_rows = []
    #print 'Process: scaling of dataframe'
    for r, v in df.iterrows():
        #print v.sum()
        rows = []
        for val in v:
            #print val
            old_min = min(v)
            old_max = max(v)
            if not isinstance(val, str):
                if val == 0:
                    rows.append(0)
                else:
                    sc_val = scale(val, (old_min, old_max), (new_min, new_max))
                    rows.append(sc_val)
                #rows.append(float(val)/v.sum())
        list_of_rows.append(rows)
    return pd.DataFrame(data=list_of_rows)
