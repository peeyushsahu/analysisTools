__author__ = 'peeyush'

from django import forms


class GenebasedForm(forms.Form):
    gene_text_field = forms.CharField(widget=forms.Textarea(attrs={'rows': 10}), required=False)
    file = forms.FileField(required=False)  # gene_upload
    only_key_celltypes = forms.BooleanField(required=False)
    gene_cluster = forms.CharField(max_length=3, strip=True, required=False)
    synonym = forms.BooleanField(required=False)
    options = (('human', 'Human'), ('mouse', 'Mouse'))
    species = forms.ChoiceField(choices=options, required=False)

    def __init__(self, data=None, *args, **kwargs):
        super(GenebasedForm, self).__init__(data, *args, **kwargs)
        # if 'gene_text_field' is empty and 'file' is not given raise error
        if data and data.get('gene_text_field', '') == '':
            self.fields['file'].required = True
            self.fields['file'].error_messages = {'required': 'Please upload file or paste gene names!!!!'}
