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
            if synonym:
                species = request.POST.get('species', '')
                parameter.update({'org': species})
            else:
                species = None
            print(gene_list, only_key_celltypes, gene_cluster, synonym, species)
            parameter.update({'synonym':synonym})
            # This will take some time to  redirect
            #time.sleep(10)
            ganalysis.gcam_analysis(parameter, outpath='', resource_path='/home/sahu/PycharmProjects/GCAM_python/GCAM/resources', genelist=gene_list)
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