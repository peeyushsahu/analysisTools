from django.shortcuts import render, redirect
from django.utils import timezone
from .models import Post
from .forms import genebasedform
import time
import GCAM.analysis as ganalysis


def posts(request):
    posts = Post.objects.all()
    return render(request, 'blog/posts.html', {'posts': posts})


def genebased(request):
    parameter = {'subcommand': 'genebased'}
    form = genebasedform
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
            only_key_celltypes = request.POST.get('only_key_celltypes', True)
            gene_cluster = request.POST.get('gene_cluster', '')
            synonym = request.POST.get('synonym', False)
            if synonym:
                species = request.POST.get('species', '')
            else:
                species = None
            print(gene_list, only_key_celltypes, gene_cluster, synonym, species)
            parameter.update({'key_celltype_list':only_key_celltypes, 'synonym':synonym,
                              'celltypeClusterSize':int(gene_cluster), 'org':species})
            # This will take some time to  redirect
            time.sleep(10)
            ganalysis.gcam_analysis('', '', parameter, gene_list)
            '''
            Here we will connect gcam
            '''
            return redirect('genebased_results', path='GCAM86482')

    return render(request, 'blog/genebased.html', {'form': form})


def genebased_res(request, path):
    print(path)
    return render(request, 'blog/genebased_results.html', {'rfolder': path})


def contact(request):
    return render(request, 'blog/contact.html', {})