from .. import models 
from django.forms import ModelForm
from django import forms
#from django.forms.widgets import *
#from django.forms.extras.widgets import *


class CalculationGroupsForm(ModelForm):
    class Meta:
        model = models.CalculationGroups
        exclude = ("calc_group_id", )
        
class ElementsForm(ModelForm): 
    class Meta:
        model = models.Elements

class GeometriesForm(ModelForm):
    class Meta:
        model = models.Geometries
        exclude = ("geom_id")
        
class MolecularSpeciesForm(ModelForm):
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    class Meta:
        model = models.MolecularSpecies
        exclude = ("species_id", )


        
class PartialMolecularSpeciesForm(ModelForm):
    class Meta:
        model = models.MolecularSpecies
        exclude = ("species_id", "time_stamp", "qual_index", )
        
class TheoryLevelsForm(ModelForm):
    description = forms.CharField()
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    xc_description = forms.CharField(required=False,  widget=forms.TextInput(attrs={'size':'40'}))
    bibliographies = forms.ModelMultipleChoiceField(queryset=models.Bibliography.alive_objects.all(), required = False)
    class Meta:
        model = models.TheoryLevels
        exclude = ("thlevel_id", )
        
class ChemistryCodesForm(ModelForm):
    description = forms.CharField(required=False)
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    bibliographies = forms.ModelMultipleChoiceField(queryset=models.Bibliography.alive_objects.all(), required = False)
    class Meta:
        model = models.ChemistryCodes
        exclude = ("code_id", )
        

class CalculationsForm(ModelForm):
    input_md5 = forms.CharField(required = False,  widget=forms.TextInput(attrs={'size':'40'}))
    output_md5 = forms.CharField(required = False,  widget=forms.TextInput(attrs={'size':'40'}))
    other_output_md5 = forms.CharField(required = False,  widget=forms.TextInput(attrs={'size':'40'}))
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    class Meta:
        model = models.Calculations
        exclude = ("calc_id", "code",)
        
        
class BasisSetsForm(ModelForm):
    name = forms.CharField(required=False) #REMOVE required=False
    description = forms.CharField(required=False)
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    bibliographies = forms.ModelMultipleChoiceField(queryset=models.Bibliography.alive_objects.all(), required = False)
    class Meta:
        model = models.BasisSets
        exclude = ("basisset_id", )
        
class TasksForm(ModelForm):
    description = forms.CharField()
    comments = forms.CharField(required=False)
    class Meta:
        model = models.Tasks
        exclude = ("task_id", "thlevel", "calc", )
        
class ElectonicStatesForm(ModelForm):
    bibliographies = forms.ModelMultipleChoiceField(queryset=models.Bibliography.alive_objects.all(), required = False)
    class Meta:
        model = models.ElectronicStates
        exclude = ("state_id", "species", "geom", "task", )
        
class DipoleMomentsForm(ModelForm):
    bibliographies = forms.ModelMultipleChoiceField(queryset=models.Bibliography.alive_objects.all(), required = False)
    class Meta:
        model = models.DipoleMoments
        exclude = ("dip_id", "state", "task", )
        
class ElementSpeciesForm(ModelForm):
    class Meta:
        model = models.ElementSpecies
        exclude = ("elem_species_id", "element", "species",)
        
class VibrationalAnalysesArmonicForm(ModelForm):
    polimode = forms.CharField(required=False)
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    bibliographies = forms.ModelMultipleChoiceField(queryset=models.Bibliography.alive_objects.all(), required = False)
    class Meta:
        model = models.VibrationalAnalysesArmonic
        exclude = ("vibalnalysis_id", "state","task",)

class RotationalConstantsForm(ModelForm):
    time_stamp = forms.DateTimeField(required=False)
    comments = forms.CharField(required=False)
    class Meta:
        model = models.RotationalConstants
        exclude = ("rot_id","state",)

class TabulatedVibrationsForm(ModelForm):
    sym_type = forms.CharField(required=False)
    eigenvectors = forms.CharField(required=False)
    class Meta:
        model = models.TabulatedVibrations
        exclude = ("vib_id", "vibanalysis_armonic",)
