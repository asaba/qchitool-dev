'''
Created on 19/dic/2011

@author: asaba
'''

from .. import models
from .. import performancetest as pt

from django.core.exceptions import ObjectDoesNotExist

def count_calculation_by_md5(md5):
    calculations = models.Calculations.alive_objects.defer("input", "output", "other_output", "comments").filter(output_md5 = md5)
    #calculations = models.Calculations.objects.filter(output_md5 = md5)
    matching_calculus = False
    first_calculus_found = None
    tasks_list = None
    if len(calculations)> 0:
        matching_calculus = True
        first_calculus_found = calculations[0].pk
        pt.p("start search tasks")
        tasks_list =  models.Tasks.alive_objects.filter(calc = first_calculus_found)
        pt.p("end search tasks")
    else:
        pass
    return matching_calculus, first_calculus_found, tasks_list

def return_molecule_specie_by_inchikey(inchikey):
    
    ms = models.MolecularSpecies.alive_objects.filter(inchikey=inchikey).order_by("-qual_index")
    if len(ms) > 0:
        return ms[0]
    else:
        return None
    
class DummyObjects():
    def __init__(self):
        self.geometryclass = models.Geometryclasses.objects.filter(geometry_md5 = "DUMMY", geometry = "DUMMY", comments = "DUMMY")
        if len(self.geometryclass) == 0:
            self.geometryclass = models.Geometryclasses(geometry_md5 = "DUMMY", geometry = "DUMMY", comments = "DUMMY")
            self.geometryclass.save()
        else:
            self.geometryclass = self.geometryclass[0]
        self.xcclasses = models.Xcclasses.objects.filter(name = "DUMMY", description = "DUMMY", comments = "DUMMY")
        if len(self.xcclasses) == 0:
            self.xcclasses = models.Xcclasses(name = "DUMMY", description = "DUMMY", comments = "DUMMY")
            self.xcclasses.save()
        else:
            self.xcclasses = self.xcclasses[0]
        
dummyobject = DummyObjects()