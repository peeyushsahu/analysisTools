__author__ = 'sahu'

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.posts, name='posts'),
    url(r'^blog/genebased/$', views.genebased, name='genebased'),
]