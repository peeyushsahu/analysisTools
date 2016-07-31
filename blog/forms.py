__author__ = 'peeyush'

from django import forms


class genebasedform(forms.Form):
    gene_text_field = forms.CharField(widget=forms.Textarea(attrs={'rows': 10}), required=False)
    file = forms.FileField(required=False)
    # gene_upload
    only_key_celltypes = forms.BooleanField(required=False)
    gene_cluster = forms.CharField(max_length=3, strip=True, required=False)
    synonym = forms.BooleanField(required=False)
    options = (('hg','Human'), ('mm','Mouse'))
    species = forms.ChoiceField(options, required=False)

