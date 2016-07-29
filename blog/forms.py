__author__ = 'peeyush'

from django import forms


class genebasedform(forms.Form):
    gene_text_field = forms.CharField(label='Paste gene list', widget=forms.Textarea)
    # gene_upload
    only_key_celltypes = forms.BooleanField(label='Association for only key cell types', required=False)
    gene_cluster = forms.CharField(label='Minimum gene in cluster', max_length=5, strip=True)
    synonym = forms.BooleanField(label='Consider gene synonym', required=True)
    options = (('hg','Human'), ('mm','Mouse'))
    species = forms.ChoiceField(options)

