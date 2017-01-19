from django import forms
from ... import performancetest as pt

class NWChemOutputField(forms.Field):
    def clean(self, value):
        if not value:
            return
        file_exts = ('.out',)
        from os.path import splitext
        if splitext(value)[1].lower() not in file_exts:
            raise forms.ValidationError("Only following Output types accepted: %s"
                                    % ", ".join([f.strip('.') for f in file_exts]))
        return value
    
class NWChemInputField(forms.Field):
    def clean(self, value):
        if not value:
            return
        file_exts = ('.nw',)
        from os.path import splitext
        if splitext(value)[1].lower() not in file_exts:
            raise forms.ValidationError("Only following Input types accepted: %s"
                                    % ", ".join([f.strip('.') for f in file_exts]))
        return value
    
class NWChemOtherField(forms.Field):
    def clean(self, value):
        if not value:
            return
        file_exts = ('.db',)
        from os.path import splitext
        if splitext(value)[1].lower() not in file_exts:
            raise forms.ValidationError("Only following Database types accepted: %s"
                                    % ", ".join([f.strip('.') for f in file_exts]))
        return value

class NWC_UploadFileForm(forms.Form):
    output_file = NWChemOutputField(label="NWChem Output File", widget=forms.FileInput())
    input_file = NWChemInputField(label="NWChem Input File", widget=forms.FileInput())
    database_file = NWChemOtherField(label="NWChem Database File", widget=forms.FileInput())

class Octopus_UploadFileForm(forms.Form):
    pass
    #output_file = NWChemOputputField(label="Octopus Output File", widget=forms.FileInput())
    #input_file = NWChemInputField(label="Octopus Input File", widget=forms.FileInput())
    #database_file = NWChemOtherField(label="Octopus Database File", widget=forms.FileInput())
