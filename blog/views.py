from django.shortcuts import render
from django.utils import timezone
from .models import Post


def posts(request):
    posts = Post.objects.all()
    return render(request, 'blog/posts.html', {'posts': posts})