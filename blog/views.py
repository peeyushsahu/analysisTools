from django.shortcuts import render
from django.utils import timezone
from .models import Post
from .forms import genebasedform


def posts(request):
    posts = Post.objects.all()
    return render(request, 'blog/posts.html', {'posts': posts})


def genebased(request):
    form = genebasedform()
    return render(request, 'blog/genebased.html', {'form': form})