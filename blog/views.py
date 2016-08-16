from django.shortcuts import render, redirect
from django.http import Http404
from .models import Post
from .forms import GenebasedForm, ExprbasedForm
import GCAM.analysis as ganalysis
import pandas as pd


def posts(request):
    posts = Post.objects.all()
    return render(request, 'blog/posts.html', {'posts': posts})

# view for genebased analysis


def genebased(request):
    parameter = {'subcommand': 'genebased'}
    form = GenebasedForm
    if request.method == 'POST':
        form = form(request.POST, request.FILES)

        if form.is_valid():
            gene_text_field = request.POST.get('gene_text_field', '')
            gene_list = list(gene_text_field.split('\r\n'))
            print(len(gene_list))
            # if upload file does not contain any gene names
            try:
                if len(gene_list) < 2:
                    print('Here')
                    file = request.FILES['file']
                    gene_list = []
                    for gene in file:
                        print(gene.decode('utf8').strip('\n'))
                        gene_list.append(gene.decode('utf8').strip('\n'))
            except:
                print('Give cell')

            only_key_celltypes = form.cleaned_data['only_key_celltypes']
            gene_cluster = request.POST.get('gene_cluster')
            #if gene_cluster == '': gene_cluster = 20
            synonym = form.cleaned_data['synonym']#request.POST.get('synonym', False)
            #if synonym:
            species = request.POST.get('species', '')
            parameter.update({'org': species})
            print(gene_list, only_key_celltypes, gene_cluster, synonym, species)
            parameter.update({'genelist': list(set(gene_list)),'synonym': synonym,
                              'key_celltype_list': only_key_celltypes,
                              'celltypeClusterSize': int(gene_cluster)})
            # This will take some time to  redirect
            #time.sleep(10)
            '''
            Here we will connect gcam
            '''
            result_path = ganalysis.gcam_analysis(parameter, outpath='/home/peeyush/Desktop/gcam_test_data',
                                                  resource_path='/home/peeyush/Desktop/gcam_test_data/resources')
            result_folder = result_path.split('/')[-1].strip()

            return redirect('genebased_results', path=result_folder)

    return render(request, 'blog/genebased.html', {'form': form})


def genebased_res(request, path):
    print(path)
    return render(request, 'blog/genebased_results.html', {'rfolder': path})

# view for exprbased analysis


def exprbased(request):
    parameter = {'subcommand': 'exprbased'}
    form = ExprbasedForm
    if request.method == 'POST':
        form = form(request.POST, request.FILES)

        if form.is_valid():
            exprDf = read_file_as_dataframe(request.FILES['exprfile'])
            print(exprDf)

            phenoDf = read_file_as_dataframe(request.FILES['phenofile'])
            print(phenoDf)
            if not inspect_exprbased_dataframes(exprDf, phenoDf):
                raise Http404("Column name mismatch b/w expression data and pheno data.")

            control = request.POST.get('select_control', '')
            parameter.update({'org': control})
            control_sample_name = None
            if control == 'sample':
                control_sample_name = request.POST.get('control_sample_name')
                print(control_sample_name)
                if not inspect_control_samplename_dataframe(control_sample_name, phenoDf):
                    raise Http404("Control sample name not found in phenodata column phenotype.")

            only_key_celltypes = form.cleaned_data['only_key_celltypes']
            fold_change = request.POST.get('fold_change')
            gene_cluster = request.POST.get('gene_cluster')
            synonym = form.cleaned_data['synonym']
            species = request.POST.get('species', '')
            parameter.update({'org': species})
            print(control, fold_change, only_key_celltypes, gene_cluster, synonym, species)

            parameter.update({'exppath': exprDf, 'phenopath': phenoDf, 'som_foldifference': fold_change,
                              'meanAsControl': control, 'controlsample': control_sample_name, 'synonym': synonym,
                              'key_celltype_list': only_key_celltypes, 'celltypeClusterSize': int(gene_cluster),
                              'selectCelltypes': None, 'som_gridsize': 10, 'somiter': 1000})

            # This will take some time to  redirect
            #time.sleep(10)

            result_path = ganalysis.gcam_analysis(parameter, outpath='/home/peeyush/Desktop/gcam_test_data',
                                                  resource_path='/home/peeyush/Desktop/gcam_test_data/resources')
            result_folder = result_path.split('/')[-1].strip()
            print(result_folder)
            return redirect('exprbased_results', path=result_folder)  # path=result_folder

    return render(request, 'blog/exprbased.html', {'form': form})


def exprbased_res(request, path):
    print(path)
    return render(request, 'blog/exprbased_results.html', {'rfolder': path})


def read_file_as_dataframe(fileobj):
    '''
    Read uploaded file as data frame.
    '''
    exprList = []
    for row in fileobj:
        nrow = row.decode('utf8').strip('\r\n').split('\t')
        exprList.append(nrow)
    headers = exprList.pop(0)
    dataframe = pd.DataFrame(exprList, columns=headers)
    return dataframe


def inspect_exprbased_dataframes(exprDf, phenoDf):
    '''
    Inspect if all sample names of phenodf are present in exprDf column name.
    '''
    return set(list(phenoDf['sample'])).issubset(list(exprDf.columns))


def inspect_control_samplename_dataframe(control_sample_name, phenoDf):
    '''
    Inspect if control sample name is present in phenodata.
    :param control_sample_name:
    :param phenoDf:
    :return:
    '''
    print(control_sample_name in list(phenoDf['phenotype']))
    return control_sample_name in list(phenoDf['phenotype'])


def contact(request):
    return render(request, 'blog/contact.html', {})