from .. import models
from .. import performancetest as pt

from modelforms import TheoryLevelsForm, BasisSetsForm
from modelforms import TasksForm, ElementsForm, MolecularSpeciesForm
from modelforms import GeometriesForm, DipoleMomentsForm,ElectonicStatesForm
from modelforms import VibrationalAnalysesHarmonicForm, RotationalConstantsForm, TabulatedVibrationsForm

from ..parser.tools import parser_tools

def BuildCalculationForm():
    pass



def BuildTasksForm(selected_tasks):
    result = {}
    if len(selected_tasks) > 0:
        result["Tasks"]= []
        result["TheoryLevels"] = []
        #result["BasisSets"] = []
        for i in range(len(selected_tasks)):
            all_theorylevel = models.TheoryLevels.alive_objects.filter(name = selected_tasks[i][1].theorylevel.name, xc_description = selected_tasks[i][1].theorylevel.xc_description).order_by("-qual_index")
            if len(all_theorylevel) > 0:  
                #TheoryLevel Found in DB
                selected_tasks[i][1].theorylevel = all_theorylevel[0]
            else:
                #complete new chemistrycode_object
                selected_tasks[i][1].theorylevel.qual_index = 0
            
            #if len(all_basis_sets) > 0:
            #    #basis set found in DB
            #    selected_tasks[i].basissets = all_basis_sets[0]
            #else:
            #    #complete new Basis Set
            #    selected_tasks[i].basissets.qual_index = 0
            selected_tasks[i][1].task.qual_index = 0 
            result["Tasks"].append(TasksForm(instance = selected_tasks[i][1].task, prefix = "Tasks" + str(i)))
            result["TheoryLevels"].append(TheoryLevelsForm(instance = selected_tasks[i][1].theorylevel, prefix = "TheoryLevels" + str(i)))
            #result["BasisSets"].append(BasisSetsForm(instance = selected_tasks[i].basissets, prefix = "BasisSets" + str(i)))
            
    return result

def BuildMolecularSpeciesForm(molecularspecies, index, forms):
    all_MolecularSpecies_list = models.MolecularSpecies.alive_objects.filter(inchikey=molecularspecies.inchikey).order_by("-qual_index")
                    
    if len(all_MolecularSpecies_list) > 0:
        #Molecular Species Found in DB
        molecularspecies = all_MolecularSpecies_list[0]
    else:
        #Complete new Molecular species Info
        molecularspecies.qual_index = 0
    if not ("MolecularSpecies" in forms):
        forms["MolecularSpecies"] = []
        
    forms["MolecularSpecies"].append(MolecularSpeciesForm(instance = molecularspecies,
                                                     prefix = "MolecularSpecies" + str(index)))      




def BuildMolecularSpeciesandRelatedForm(molecularspecies, elements, electronicstates, dipolemoments, index, forms):
    
    BuildMolecularSpeciesForm(molecularspecies, index, forms)
        
    if elements:
        #CHECK if it is an isotope
        if not ("Elements" in forms):
            forms["Elements"] = []
        forms["Elements"].append([])
        
        if not ("BasisSets" in forms):
            forms["BasisSets"] = []
        forms["BasisSets"].append([])
        #database_elements = models.Elements.objects.all()
        #database_basissets_alive = models.BasisSets.alive_objects.all()
        element_index = 0
        for element in elements:
            element.update_elements_by_database()
            #filtered_elements = database_elements.filter(atomic_number = element[0].atomic_number, 
            #                                                   atomic_mass = round(element[0].atomic_mass))
            #if len(filtered_elements)>0:
            #    #Element from DB
            #    element[0] = filtered_elements[0]
            #filtered_basissets = database_basissets_alive.filter(name = element[1].name).order_by("-qual_index")    
            #if len(filtered_basissets) > 0:
            #    #Basis Set from DB
            #    element[1] = filtered_basissets[0]

            forms["Elements"][-1].append(ElementsForm(instance = element.element, prefix = "Elements" + str(index) + "-" + str(element_index)))
            forms["BasisSets"][-1].append(BasisSetsForm(instance = element.basisset, prefix = "BasisSets"+ str(index) + "-" + str(element_index)))
            element_index += 1

    BuildElectronicStateForm(electronicstates, index, forms)
    
    BuildDipoleMomentForm(dipolemoments, index, forms)

def BuildChemistryCode():
    pass

def BuildElectronicStateForm(electronicstates, index, forms):
    if not ("ElectronicStates"  in forms):
        forms["ElectronicStates"] = []
    
    forms["ElectronicStates"].append(ElectonicStatesForm(instance = electronicstates, 
                                                     prefix = "ElectronicStates" + str(index)))    

def BuildDipoleMomentForm(dipolemoments, index, forms):
    if not ("DipoleMoments" in forms):
        forms["DipoleMoments"] = []
    if dipolemoments:
        forms["DipoleMoments"].append(DipoleMomentsForm(instance = dipolemoments, 
                                                     prefix = "DipoleMoments" + str(index)))
    else:
        forms["DipoleMoments"].append(None)

def BuildVibrationalAnalysisHarmonic(vibrationalanalisysarmonic, rotationalconstants, tabulatedvibrations, index, forms):
    #this function build all object related to vibrational analisys
    
    if not ("VibrationalAnalysesHarmonic" in forms):
        forms["VibrationalAnalysesHarmonic"] = []
    
    if vibrationalanalisysarmonic:
        forms["VibrationalAnalysesHarmonic"].append(VibrationalAnalysesHarmonicForm(instance = vibrationalanalisysarmonic,
                                                     prefix = "VibrationalAnalysesHarmonic" + str(index)))
    else:
        forms["VibrationalAnalysesHarmonic"].append(None)
        
    if not ("RotationalConstants" in forms):
        forms["RotationalConstants"] = []    
    if rotationalconstants:
        forms["RotationalConstants"].append(RotationalConstantsForm(instance = rotationalconstants, 
                                                     prefix = "RotationalConstants" + str(index)))
    else:
        forms["RotationalConstants"].append(None)
        
    if not ("TabulatedVibrations" in forms):
        forms["TabulatedVibrations"] = []
        
    forms["TabulatedVibrations"].append([])    
    if tabulatedvibrations:
        for vibration_index in range(len(tabulatedvibrations)):
            forms["TabulatedVibrations"][-1].append(TabulatedVibrationsForm(instance = tabulatedvibrations[vibration_index], prefix = "Vibrations" + str(index) + "-" + str(vibration_index)))
    else:
        forms["TabulatedVibrations"][-1].append(None)

def BuildGeometriesForm(geometries, converged, index, forms):
    all_geometry_list = models.Geometries.alive_objects.filter(geometry_md5 = geometries.geometry_md5).order_by("-qual_index")
    if len(all_geometry_list)>0:
        #Geometry found in Database
        geometries = all_geometry_list[0]
    else:
        #Complete new geometry Info
        if converged <> {}:
            pass
            #symmetry_section = parser_tools.returnsection("SymmetryInformation", converged["section"], "First")
            ##for Geometry ASK to user
            ##calculation_object["Geometries"].sym_group = #calculate
            ##for total_energy
            #NwchemDftModule_section = parser_tools.returnsection("NwchemDftModule", converged["section"], "First")
            ##for spin multiplicity
            #GeneralInformation_section = parser_tools.returnsection("GeneralInformation", NwchemDftModule_section, "First")
    if not ("Geometries" in forms):
        forms["Geometries"] = []
        
    forms["Geometries"].append(GeometriesForm(instance = geometries,
                                         prefix = "Geometries" + str(index)))

    