__author__ = 'peeyush'

from django import forms


class genebasedform(forms.Form):
    gene_text_field = forms.CharField(label='Paste gene list', widget=forms.Textarea)
    # gene_upload
    only_key_celltypes = forms.BooleanField(label='Association for only key cell types', required=False)
    synonym = forms.BooleanField(label='Consider gene synonym', required=True)
    options = (('Human','hs'), ('Mouse','mm'))
    species = forms.ChoiceField(options)

