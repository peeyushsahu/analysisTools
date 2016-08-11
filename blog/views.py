from django.shortcuts import render, redirect
from django.utils import timezone
from .models import Post
from .forms import GenebasedForm
import time
import GCAM.analysis as ganalysis


def posts(request):
    posts = Post.objects.all()
    return render(request, 'blog/posts.html', {'posts': posts})


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
                        gene_list.append(gene)
            except:
                print('Give cell')

            only_key_celltypes = form.cleaned_data['only_key_celltypes']
            gene_cluster = request.POST.get('gene_cluster', '')
            synonym = form.cleaned_data['synonym']#request.POST.get('synonym', False)
            #if synonym:
            species = request.POST.get('species', '')
            parameter.update({'org': species})
            print(set(gene_list), only_key_celltypes, gene_cluster, synonym, species)
            parameter.update({'synonym': synonym,
                              'key_celltype_list': only_key_celltypes,
                              'celltypeClusterSize': int(gene_cluster)})
            # This will take some time to  redirect
            #time.sleep(10)
            '''
            Here we will connect gcam
            '''
            result_path = ganalysis.gcam_analysis(parameter, outpath='/home/peeyush/Desktop/gcam_test_data',
                                                  resource_path='/home/peeyush/PycharmProjects/GCAM_python/GCAM/resources', genelist=list(set(gene_list)))
            result_folder = result_path.split('/')[-1].strip()

            return redirect('genebased_results', path=result_folder)

    return render(request, 'blog/genebased.html', {'form': form})


def genebased_res(request, path):
    print(path)
    return render(request, 'blog/genebased_results.html', {'rfolder': path})


def contact(request):
    return render(request, 'blog/contact.html', {})