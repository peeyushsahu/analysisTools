__author__ = 'sahu'

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^gcam/$', views.posts, name='posts'),
    url(r'^gcam/genebased/$', views.genebased, name='genebased'),
    url(r'^gcam/contact/$', views.contact, name='contact'),
    url(r'^gcam/genebased_results/(?P<path>GCAM[0-9]+)/$', views.genebased_res, name='genebased_results'),
]