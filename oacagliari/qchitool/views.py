'''
Created on 26/ago/2011
@author: asaba
'''
# Create your views here.

import os

from django.contrib.sessions.models import Session
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response, get_object_or_404
from django.core.context_processors import csrf
from django.db import transaction
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout

from models import *
from modelforms import modelforms, build_form
from modelforms.nwchem import forms
from tools import qchitool_tools, modeltools
from parser import impnwc, commonsobjects
from parser.nwchem import commons, check_import_object
from parser.tools import parser_tools
from version import version
from tools.repairtool.nwchem import bugfix 
from tools import geometry_tools as gt
from tools.functionals import functional_tools as ft
from tools.ionizationenergies import ionizationenergies_tools as iet
import performancetest as pt

def global_index(request, PATH_FOR_MOLECULE_FILES):
    if not os.path.exists(PATH_FOR_MOLECULE_FILES):
	try:	    
    	    os.makedirs(PATH_FOR_MOLECULE_FILES)
	except:
	    pass
    currente_session = Session.objects.filter(session_key=request.session.session_key) 
    if currente_session:
        currente_session.delete()
    logout(request)
    t = loader.get_template("index.html")
    c = Context({"version": version, })
    #output = ",".join([e.name for e in all_elements_list])
    return HttpResponse(t.render(c))

def explore_index(request, PATH_FOR_MOLECULE_FILES):
    if not os.path.exists(PATH_FOR_MOLECULE_FILES):
        os.makedirs(PATH_FOR_MOLECULE_FILES)
    currente_session = Session.objects.filter(session_key=request.session.session_key) 
    if currente_session:
        currente_session.delete()
    logout(request)
    #molecules = MolecularSpecies.alive_objects.get_alive_query_set(timecheck = datetime.datetime(year = 2012, month = 6, day = 10)).all().order_by("name")
    molecules = MolecularSpecies.alive_objects.all().order_by("name")
    molecules.query.group_by = ["name", "formula"]
    t = loader.get_template("explore/index.html")
    c = Context({"molecules": molecules,
                 "version": version, })
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
def nwc_import_index(request):
    post_data = request.POST.copy()
    expert_form = forms.NWC_UploadFileForm(post_data)  
    #delete session info
    if expert_form.is_valid():
            t = loader.get_template("nwc_import/index.html")
            c = Context({"expert_form": expert_form, })
            c.update(csrf(request)) 
            return HttpResponse(t.render(c))
    else:
            return HttpResponse("This Form isn't valid")
        
        
@login_required(login_url='/accounts/login/')
def nwc_import_file_info(request, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    post_data.update(request.FILES)
    enable_jmol = False
    if "output_file" in post_data:
        output_file = post_data["output_file"]
    else:
        output_file = None
    if "input_file" in post_data:
        input_file = post_data["input_file"]
    else:
        input_file = None
    if "database_file" in post_data:
        database_file = post_data["database_file"]
    else:
        database_file = None
    if "jmol" in post_data:
        enable_jmol = post_data["jmol"]
    else:
        enable_jmol = False
        
    md5 = qchitool_tools.returnmd5(request.session.session_key)
    session_path = PATH_FOR_MOLECULE_FILES + str(md5) + "/"
    if not os.path.exists(session_path):
        os.makedirs(session_path)
    uploaded_output_file = qchitool_tools.handle_uploaded_file(output_file, session_path)
    uploaded_input_file = qchitool_tools.handle_uploaded_file(input_file, session_path)
    uploaded_database_file = qchitool_tools.handle_uploaded_file(database_file, session_path)
    #new_xyz_file = newfilename = session_path + output_file.name


    request.session.update(check_import_object.return_parsed_object_for_session(str(output_file), uploaded_output_file, uploaded_input_file, uploaded_database_file, session_path, is_dbfile_humanreadable = False))

    t = loader.get_template("nwc_import/file_info_detail.html")
    c = Context({"tasks": request.session["tasks"],
                 "frl": commons.fields_repr_by_list, 
                 "frd": commons.fields_repr_by_dict,
                 "parsing_result": None, #request.session["parsing_result"],
                 "calculation_object": request.session["calculation_object"],
                 "calculation_found_tasks_list" : request.session["calculation_found_tasks_list"],
                 "calculation_by_md5": request.session["calculation_found"],
                 "URL_FOR_MOLECULE_FILES": URL_FOR_MOLECULE_FILES,
                 "enable_jmol": enable_jmol,
             })
    c.update(csrf(request))
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
def nwc_import_file_atoms_info(request):
    post_data = request.POST.copy()
    post_data.update(request.FILES)
    database_atoms = []
    
    for atom in request.session["calculation_object"]["Elements"]:
        list_of_database_elements = Elements.objects.filter(atomic_number=atom.atomic_number, 
                                                            atomic_mass = atom.atomic_mass)
        if len(list_of_database_elements) > 0:
            #Element in DB
            database_atoms.append([atom, list_of_database_elements[0]])
        else:
            database_atoms.append([atom, None])
    t = loader.get_template("nwc_import/molecule_elements_info.html")
    c = Context({
                     "atoms": database_atoms,
                 })
    c.update(csrf(request))
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
def nwc_import_file_task_info(request):
    #Get selected tasks and show the tasks' information that will be saved
    post_data = request.POST.copy()
    post_data.update(request.FILES)
    
    if post_data["selected_action"] == "Save":
        selected_tasks = post_data.getlist("task")
        results = associate_task_operation(selected_tasks, request)

        t = loader.get_template("nwc_import/file_info_task_detail.html")
        c = Context({
                 "results": results,
                 "frl": commons.fields_repr_by_list, 
                 "frd": commons.fields_repr_by_dict,
             })
    else:
            t = loader.get_template("nwc_import/index.html")
            c = Context({
                 })

    c.update(csrf(request))
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
def nwc_post_selected_task_to_save(request):
    post_data = request.POST.copy()
    post_data.update(request.FILES)
    calculation_object = request.session["calculation_object"] 

    if post_data["selected_action"] == "Save":
        selected_tasks = request.session["selected_tasks"]

    
        forms, tasks_forms_obj = build_tasks_forms(calculation_object, selected_tasks)
        request.session["calculation_object"] = calculation_object 
        
        t = loader.get_template("nwc_import/calculus_information_to_save.html")
        c = Context({
                 "forms": forms,
                 "tasks": tasks_forms_obj,
             })
        c.update(csrf(request))
        return HttpResponse(t.render(c))   
        
@login_required(login_url='/accounts/login/')
@transaction.commit_on_success
def nwc_calculus_information_to_save(request):

    forms = {}
    if request.method == 'POST': # If the calculation_form has been submitted...
        calculation_object = request.session["calculation_object"] 
        selected_tasks = request.session["selected_tasks"]
        try:
            calculations_form = modelforms.CalculationsForm(request.POST,  prefix = "Calculations") # A calculation_form bound to the POST data
            chemistrycodes_form = modelforms.ChemistryCodesForm(request.POST, prefix = "ChemistryCodes")
            tasks_form = []
            theorylevels_form = []
            geometries_form = []
            molecularspecies_form = []
            elements_form = []
            basissets_form = []
            electronicstates_form = []
            dipolemoments_form = []
            tabualtedvibrations_form = []
            rotationalconstants_form = []
            vibrationalanalisysarmonic_form = []
            electronicstates_by_geom = []
            for i in range(len(selected_tasks)):
                tasks_form.append(modelforms.TasksForm(request.POST, prefix = "Tasks" + str(i)))
                theorylevels_form.append(modelforms.TheoryLevelsForm(request.POST, instance = TheoryLevels(), prefix = "TheoryLevels" + str(i)))
                geometries_form.append(modelforms.GeometriesForm(request.POST, prefix = "Geometries" + str(i)))

                molecularspecies_form.append(modelforms.MolecularSpeciesForm(request.POST, prefix = "MolecularSpecies" + str(i)))

                elements_form.append([])
                basissets_form.append([])
                for l in range(len(calculation_object["Elements"][i])):
                    elements_form[-1].append(modelforms.ElementsForm(request.POST, prefix = "Elements" + str(i) + "-" + str(l)))
                    basissets_form[-1].append(modelforms.BasisSetsForm(request.POST, prefix = "BasisSets" + str(i) + "-" + str(l)))
                
                electronicstates_form.append(modelforms.ElectonicStatesForm(request.POST, prefix = "ElectronicStates" + str(i)))
                dipolemoments_form.append(modelforms.DipoleMomentsForm(request.POST, prefix = "DipoleMoments" + str(i)))
                
                tabualtedvibrations_form.append([])
                if calculation_object["VibrationalAnalysesHarmonic"][i]:
                    for l in range(len(calculation_object["TabulatedVibrations"][i])):
                        tabualtedvibrations_form[-1].append(modelforms.TabulatedVibrationsForm(request.POST, prefix = "Vibrations" + str(i) + "-" + str(l)))

                rotationalconstants_form.append(modelforms.RotationalConstantsForm(request.POST, prefix = "RotationalConstants" + str(i)))
                vibrationalanalisysarmonic_form.append(modelforms.VibrationalAnalysesHarmonicForm(request.POST, prefix = "VibrationalAnalysesHarmonic" + str(i)))
        except ValueError:  
            return HttpResponse("Error" + str(ValueError))
        
        
        #BEGIN Transaction
        if calculation_object["ChemistryCodes"].pk:
            #ChemistryCode present in DB
            
            c = chemistrycodes_form.save(commit = False)
            if (calculation_object["ChemistryCodes"].name == 
                c.name) and (calculation_object["ChemistryCodes"].version == 
                             c.version) and (calculation_object["ChemistryCodes"].description == 
                                             c.description) and (calculation_object["ChemistryCodes"].comments == c.comments):
                
                #do not save
                pass
            else:
                #update information
                calculation_object["ChemistryCodes"].name = c.name
                calculation_object["ChemistryCodes"].version = c.version
                calculation_object["ChemistryCodes"].description = c.description
                calculation_object["ChemistryCodes"].comments = c.comments
                calculation_object["ChemistryCodes"].now(True)
                calculation_object["ChemistryCodes"].save()
        else:
            #new chemistry code
            calculation_object["ChemistryCodes"] = chemistrycodes_form.save(commit = False)
            calculation_object["ChemistryCodes"].now(True)
            calculation_object["ChemistryCodes"].save()
            
        if calculation_object["Calculations"].pk:
            #Calculation Present on DB
            c = calculations_form.save(commit = False)
            if (calculation_object["Calculations"].input_md5 == 
                c.input_md5) and (calculation_object["Calculations"].output_md5 == 
                             c.output_md5) and (calculation_object["Calculations"].other__output_md5 == 
                                             c.other__output_md5) and (calculation_object["Calculations"].qual_index == c.qual_index):
                
                #do not save
                pass
            else:
                #update information
                calculation_object["Calculations"].input = c.input
                calculation_object["Calculations"].output = c.output
                calculation_object["Calculations"].other_output = c.other_output
                calculation_object["Calculations"].qual_index = c.qual_index
                calculation_object["Calculations"].comments = c.comments
                calculation_object["Calculations"].calculate_md5()
                calculation_object["Calculations"].now(True)
                calculation_object["Calculations"].save()
        else:
            calculation_object["Calculations"] = calculations_form.save(commit = False)
            calculation_object["Calculations"].code = calculation_object["ChemistryCodes"]
            calculation_object["Calculations"].now()
            calculation_object["Calculations"].save()               
        
        for i in range(len(selected_tasks)):
            if selected_tasks[i][1].theorylevel.pk:
                #Theory Level present in DB
                t = theorylevels_form[i].save(commit = False)
                selected_tasks[i][1].theorylevel.qual_index = t.qual_index
                selected_tasks[i][1].theorylevel.now()
                selected_tasks[i][1].theorylevel.save()
            else:
                selected_tasks[i][1].theorylevel = theorylevels_form[i].save(commit = False)
                selected_tasks[i][1].theorylevel.now()
                selected_tasks[i][1].theorylevel.save()
                        
            selected_tasks[i][1].task = tasks_form[i].save(commit = False)
            selected_tasks[i][1].task.calc = calculation_object["Calculations"]
            selected_tasks[i][1].task.thlevel = selected_tasks[i][1].theorylevel
            selected_tasks[i][1].task.now()
            selected_tasks[i][1].task.save()                
    
            if calculation_object["Geometries"][i].pk:
                g = geometries_form[i].save(commit = False)
                if (calculation_object["Geometries"][i].geometry == 
                    g.geometry) and (calculation_object["Geometries"][i].geometry_md5 ==
                    g.geometry_md5) and (calculation_object["Geometries"][i].sym_group == 
                    g.sym_group) and (calculation_object["Geometries"][i].sym_elements == 
                    g.sym_elements) and (calculation_object["Geometries"][i].comments == 
                    g.comments) and (calculation_object["Geometries"][i].geometryclass_id == g.geometryclass_id) :
                    #do not save
                    pass
                else:
                    #update information
                    calculation_object["Geometries"][i].geometry = g.geometry 
                    calculation_object["Geometries"][i].geometry_md5 = g.geometry_md5
                    calculation_object["Geometries"][i].sym_group = g.sym_group
                    calculation_object["Geometries"][i].sym_elements = g.sym_elements
                    calculation_object["Geometries"][i].comments = g.comments
                    calculation_object["Geometries"][i].now(True)
                    calculation_object["Geometries"][i].save()
            else:
                calculation_object["Geometries"][i] = geometries_form[i].save(commit = False)
                calculation_object["Geometries"][i].now()
                calculation_object["Geometries"][i].save()
            
            #Check Geometry Class
            
            
            if calculation_object["MolecularSpecies"][i].pk:
                #add new molecular species if quality index is up
                c = molecularspecies_form[i].save(commit = False)
                if (calculation_object["MolecularSpecies"][i].name == 
                    c.name) and (calculation_object["MolecularSpecies"][i].formula == 
                    c.formula) and (calculation_object["MolecularSpecies"][i].inchi ==
                    c.inchi) and (calculation_object["MolecularSpecies"][i].inchikey ==
                    c.inchikey) and (calculation_object["MolecularSpecies"][i].aromatic_cycles == 
                    c.aromatic_cycles) and (calculation_object["MolecularSpecies"][i].charge ==
                    c.charge) and (calculation_object["MolecularSpecies"][i].comments ==
                    c.comments) and (calculation_object["MolecularSpecies"][i].isotopologue_of == c.isotopologue_of):
                    pass
                else:
                    #update information
                    calculation_object["MolecularSpecies"][i].name = c.name
                    calculation_object["MolecularSpecies"][i].formula = c.formula
                    calculation_object["MolecularSpecies"][i].inchi = c.inchi
                    calculation_object["MolecularSpecies"][i].inchikey = c.inchikey
                    calculation_object["MolecularSpecies"][i].aromatic_cycles = c.aromatic_cycles
                    calculation_object["MolecularSpecies"][i].charge = c.charge
                    calculation_object["MolecularSpecies"][i].comments = c.comments
                    calculation_object["MolecularSpecies"][i].isotopologue_of = c.isotopologue_of
                    calculation_object["MolecularSpecies"][i].now()
                    calculation_object["MolecularSpecies"][i].save()
                
                #create all elementspeciesbasisset record by task, molecular specie and elements (with basisset)
                create_elementspeciesbasisset(task = selected_tasks[i][1].task, elementspecies = create_elementspecies(molecularspecie = calculation_object["MolecularSpecies"][i], elements = calculation_object["Elements"][i]), elements = calculation_object["Elements"][i])

 
            else:
                calculation_object["MolecularSpecies"][i] = molecularspecies_form[i].save(commit = False)
                calculation_object["MolecularSpecies"][i].now()
                calculation_object["MolecularSpecies"][i].save()
                elementspecies = []
                for j in range(len(calculation_object["Elements"][i])):
                    if calculation_object["Elements"][i][j].element.pk:
                        #Element from DB
                        e = elements_form[i][j].save(commit = False)
                        if (calculation_object["Elements"][i][j].element.name == 
                            e.name) and (calculation_object["Elements"][i][j].element.symbol == 
                            e.symbol) and (calculation_object["Elements"][i][j].element.element_group == 
                            e.element_group) and (calculation_object["Elements"][i][j].element.standard_atomic_weight ==
                            e.standard_atomic_weight) and (calculation_object["Elements"][i][j].element.atomic_number == 
                            e.atomic_number) and (calculation_object["Elements"][i][j].element.atomic_mass == 
                            e.atomic_mass):
                            pass
                        else:
                            #update
                            calculation_object["Elements"][i][j].element.name = e.name
                            calculation_object["Elements"][i][j].element.symbol = e.symbol
                            #calculation_object["Elements"][i].atomic_number = e.atomic_number
                            calculation_object["Elements"][i][j].element.element_group = e.element_group
                            calculation_object["Elements"][i][j].element.standard_atomic_weight = e.standard_atomic_weight
                            calculation_object["Elements"][i][j].element.save() #update
                    else:
                        #new element
                        calculation_object["Elements"][i][j].element = elements_form[i][j].save(commit = False)
                        calculation_object["Elements"][i][j].element.save() 
                        #Check other elements 
                        for element in calculation_object["Elements"][i]:
                            if not element.element.pk:
                                if (element.element.atomic_mass == calculation_object["Elements"][i][j].element.atomic_mass) and (element.element.atomic_mass == calculation_object["Elements"][i][j].element.atomic_mass):
                                    element.element.pk = calculation_object["Elements"][i][j].element.pk
                                
                                
                        
                    if calculation_object["Elements"][i][j].basisset.pk:
                        #BasisSet from DB
                        b = basissets_form[i][j].save(commit = False)
                        if (calculation_object["Elements"][i][j].basisset.name == 
                            b.name) and (calculation_object["Elements"][i][j].basisset.description == 
                            b.description) and (calculation_object["Elements"][i][j].basisset.comments ==
                            b.comments) and (calculation_object["Elements"][i][j].basisset.qual_index == 
                            b.qual_index):
                            pass
                        else:
                            calculation_object["Elements"][i][j].basisset.name = b.name
                            calculation_object["Elements"][i][j].basisset.description = b.description
                            calculation_object["Elements"][i][j].basisset.comments = b.comments
                            calculation_object["Elements"][i][j].basisset.qual_index = b.qual_index
                            calculation_object["Elements"][i][j].basisset.bibliographies = b.bibliographies
                            calculation_object["Elements"][i][j].basisset.now()
                            calculation_object["Elements"][i][j].basisset.save()
                    else:
                        #New Basis set
                        calculation_object["Elements"][i][j].basisset = basissets_form[i][j].save(commit = False)
                        calculation_object["Elements"][i][j].basisset.now()
                        calculation_object["Elements"][i][j].basisset.save()
                            #Check other basis sets 
                        for element in calculation_object["Elements"][i]:
                            if not element.basisset.pk:
                                if (element.basisset.name == calculation_object["Elements"][i][j].basisset.name) and (element.basisset.description == calculation_object["Elements"][i][j].basisset.description):
                                    element.basisset.pk = calculation_object["Elements"][i][j].basisset.pk
                                    #element[1].bibliographies = calculation_object["Elements"][i][j][1].bibliographies
                                    #element[1].fields["bibliographies"] = calculation_object["Elements"][i][j][1].bibliographies
                                    #write values in fields
                                    
                    #check the presence of elementspecies
                    #elementspecies_list = ElementSpecies.objects.filter(element = calculation_object["Elements"][i][j][0], species = calculation_object["MolecularSpecies"][i])
                    #if len(elementspecies_list) > 0:
                    #    #elementspecies present
                    #    elementspecies = elementspecies_list[0]
                    #else:
                    #    elementspecies = ElementSpecies(element = calculation_object["Elements"][i][j][0], 
                    #                                species = calculation_object["MolecularSpecies"][i])
                    #    elementspecies.save()
                    elementspecies = ElementSpecies(element = calculation_object["Elements"][i][j].element, 
                                                    species = calculation_object["MolecularSpecies"][i])
                    elementspecies.save()
                    elementspeciesbasisset = ElementSpeciesBasisSet(basisset = calculation_object["Elements"][i][j].basisset,
                                                                    task = selected_tasks[i][1].task,
                                                                    elementspecies = elementspecies)
                    elementspeciesbasisset.save()


            calculation_object["ElectronicStates"][i] = electronicstates_form[i].save(commit = False)
            calculation_object["ElectronicStates"][i].task = selected_tasks[i][1].task
            calculation_object["ElectronicStates"][i].species = calculation_object["MolecularSpecies"][i]
            calculation_object["ElectronicStates"][i].geom = calculation_object["Geometries"][i]
            calculation_object["ElectronicStates"][i].now()
            calculation_object["ElectronicStates"][i].save()
            
            #calculate energy ionization
            iet.build_ionizationenergies(calculation_object["ElectronicStates"][i])
            #check for is_minimun in other electronic states
            electronicstates_by_geom.append(gt.electronicstatebygeometryclass(calculation_object["ElectronicStates"][i], None, None, None))
            
            
            if calculation_object["DipoleMoments"][i]:
                if calculation_object["DipoleMoments"][i].pk:
                    pass
                else:
                    calculation_object["DipoleMoments"][i] = dipolemoments_form[i].save(commit = False)
                    calculation_object["DipoleMoments"][i].state = calculation_object["ElectronicStates"][i]
                    calculation_object["DipoleMoments"][i].task = selected_tasks[i][1].task
                    calculation_object["DipoleMoments"][i].now()
                    calculation_object["DipoleMoments"][i].save()
                
            if calculation_object["RotationalConstants"][i]:
                if calculation_object["RotationalConstants"][i].pk:
                    pass
                else:
                    calculation_object["RotationalConstants"][i] = rotationalconstants_form[i].save(commit = False)
                    calculation_object["RotationalConstants"][i].state = calculation_object["ElectronicStates"][i]
                    calculation_object["RotationalConstants"][i].now()
                    calculation_object["RotationalConstants"][i].save()
                
            if calculation_object["VibrationalAnalysesHarmonic"][i]:
                if calculation_object["VibrationalAnalysesHarmonic"][i].pk:
                    pass
                else:
                    calculation_object["VibrationalAnalysesHarmonic"][i] = vibrationalanalisysarmonic_form[i].save(commit = False)
                    calculation_object["VibrationalAnalysesHarmonic"][i].state = calculation_object["ElectronicStates"][i]
                    calculation_object["VibrationalAnalysesHarmonic"][i].task = selected_tasks[i][1].task
                    calculation_object["VibrationalAnalysesHarmonic"][i].now()
                    calculation_object["VibrationalAnalysesHarmonic"][i].save()
                    
                if calculation_object["TabulatedVibrations"][i]:
                    for j in range(len(calculation_object["TabulatedVibrations"][i])):
                        calculation_object["TabulatedVibrations"][i][j] = tabualtedvibrations_form[i][j].save(commit = False)
                        calculation_object["TabulatedVibrations"][i][j].vibanalysis_armonic = calculation_object["VibrationalAnalysesHarmonic"][i]
                        calculation_object["TabulatedVibrations"][i][j].save()
                    
        request.session["calculation_object"]  = calculation_object
        request.session["electronicstates"] = []
        index = 0
        for i in range(len(electronicstates_by_geom)):
            eslist = electronicstates_by_geom[i]
            for es in eslist:
                request.session["electronicstates"].append(es)
                build_form.BuildElectronicStateForm(es, index, forms)
                index += 1
            
        t = loader.get_template("nwc_import/similar_electronic_states_to_save.html")
        c = Context({"forms": forms["ElectronicStates"],
             })
        c.update(csrf(request))
        return HttpResponse(t.render(c))     

    else:
        pass
        #calculation_form = ContactForm() # An unbound calculation_form
    
        #return render_to_response('contact.html', {'calculation_form': calculation_form, })
        return HttpResponse("Request Method isn't POST")
    
@login_required(login_url='/accounts/login/')
@transaction.commit_on_success
def nwc_electronic_states_similar(request):
    index = 0
    if request.method == 'POST': # If the calculation_form has been submitted...
        for es in request.session["electronicstates"]:
            electronicstates_form = modelforms.ElectonicStatesForm(request.POST, prefix = "ElectronicStates" + str(index))

            es = electronicstates_form.save(commit = False)
            es.now()
            es.save()


        t = loader.get_template("nwc_import/calculus_saved.html")
        c = Context({
             })
        c.update(csrf(request))
        return HttpResponse(t.render(c))    



@login_required(login_url='/accounts/login/')
def octopus_import_index(request):
    post_data = request.POST.copy()
    expert_form = forms.Octopus_UploadFileForm(post_data)  
    #delete session info
    if expert_form.is_valid():
            t = loader.get_template("octopus_import/index.html")
            c = Context({"expert_form": expert_form, })
            c.update(csrf(request)) 
            return HttpResponse(t.render(c))
    else:
            return HttpResponse("This Form isn't valid")
        
        
@login_required(login_url='/accounts/login/')
def octopus_import_file_info(request, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    post_data.update(request.FILES)
    enable_jmol = False
    if "output_file" in post_data:
        output_file = post_data["output_file"]
    else:
        output_file = None
    if "input_file" in post_data:
        input_file = post_data["input_file"]
    else:
        input_file = None
    if "database_file" in post_data:
        database_file = post_data["database_file"]
    else:
        database_file = None
    if "jmol" in post_data:
        enable_jmol = post_data["jmol"]
    else:
        enable_jmol = False
        
    md5 = qchitool_tools.returnmd5(request.session.session_key)
    session_path = PATH_FOR_MOLECULE_FILES + str(md5) + "/"
    if not os.path.exists(session_path):
        os.makedirs(session_path)
    uploaded_output_file = qchitool_tools.handle_uploaded_file(output_file, session_path)
    uploaded_input_file = qchitool_tools.handle_uploaded_file(input_file, session_path)
    uploaded_database_file = qchitool_tools.handle_uploaded_file(database_file, session_path)
    #new_xyz_file = newfilename = session_path + output_file.name

    parsing_result = impnwc.parsing(uploaded_output_file, True, "NwChemOutput") 
    calculation_object = check_import_object.build_calculation_object(uploaded_output_file, 
                                                                      uploaded_input_file, 
                                                                      uploaded_database_file, 
                                                                      parsing_result, 
                                                                      session_path)
    
    molecule = None
    #last_molecular_species = None
    calculation_found = False
    geometry = []
    tasks_description = []
    relative_file_path = []
    if "Calculations" in calculation_object:
        calculation_by_md5 = Calculations.alive_objects.filter(output_md5 = calculation_object["Calculations"].output_md5)
        if len(calculation_by_md5) > 0:
            calculation_found = True
            
    if  "Tasks" in calculation_object:
        if len(calculation_object["Tasks"]) > 0:
            for i in range(len(calculation_object["Tasks"])):
                    
                relative_file_path = re.sub(r".*media", "", calculation_object["molecule_info"][i].geometries["sdf"]["filename"])
                #Delete molecule object because it'snt serializable in session
                calculation_object["molecule_info"][i].moleculeobject = None
                #Search molecular species in DB
                all_MolecularSpecies_list = MolecularSpecies.alive_objects.filter(inchikey=calculation_object["MolecularSpecies"][i].inchikey).order_by("-qual_index")
                
                if len(all_MolecularSpecies_list) > 0:
                    #Almost one MolecularSpecies identified by inchikey is present
                    calculation_object["MolecularSpecies"][i] = all_MolecularSpecies_list[0]
                    
                tasks_description.append(commonsobjects.TaskCode(calculation_object["TheoryLevels"][i], 
                                                                 calculation_object["Tasks"][i],
                                                                 calculation_object["molecule_info"][i],
                                                                 calculation_object["ElectronicStates"][i],
                                                                 calculation_object["RotationalConstants"][i],
                                                                 calculation_object["VibrationalAnalysesHarmonic"][i],
                                                                 calculation_object["DipoleMoments"][i],
                                                                 calculation_object["MolecularSpecies"][i],
                                                                 calculation_object["Geometries"][i],
                                                                 relative_file_path,
                                                                 calculation_object["Elements"][i]))
            request.session["tasks"] = tasks_description

    request.session["calculation_object"] = calculation_object
    
    t = loader.get_template("nwc_import/file_info_detail.html")
    c = Context({"tasks": request.session["tasks"],
                 "frl": commons.fields_repr_by_list, 
                 "frd": commons.fields_repr_by_dict,
                 "parsing_result": parsing_result,
                 "calculation_object": calculation_object,
                 "calculation_by_md5": calculation_found,
                 "URL_FOR_MOLECULE_FILES": URL_FOR_MOLECULE_FILES,
                 "enable_jmol": enable_jmol,
             })
    c.update(csrf(request))
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
def utility_index(request):
    #post_data = request.POST.copy()
    #delete session info

    t = loader.get_template("utility/bugfixutility.html")
    c = Context()
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
def utility_geom(request, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    md5 = qchitool_tools.returnmd5(request.session.session_key)
    session_path = PATH_FOR_MOLECULE_FILES + str(md5) + "/"
    bugfix.re_run_parsing(session_path) 
    t = loader.get_template("utility/bugfixutility.html")
    c = Context()
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
#@transaction.commit_on_success
def utility_geom_classes(request, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    md5 = qchitool_tools.returnmd5(request.session.session_key)
    session_path = PATH_FOR_MOLECULE_FILES + str(md5) + "/"
    if not os.path.exists(session_path):
        os.makedirs(session_path)
    outputresult, classes = gt.buildequivalenceclasses2(session_path, limit = 0)
    #elimina tutti record della tebella delle classi di geometria
    Geometries.alive_objects.all().update(geometryclass = modeltools.dummyobject.geometryclass)
    Geometryclasses.objects.exclude(pk = modeltools.dummyobject.geometryclass.geometryclass_id).delete()
    for key, value in classes.items():
        geom_equiv_result = value[:]
        if len(geom_equiv_result) == 0:
            geometry_ID = key
            geometry = Geometries.alive_objects.get(pk=geometry_ID)
            if geometry.geometryclass_id <> modeltools.dummyobject.geometryclass.geometryclass_id:
                geometry.geometryclass_id = gt.buildsingolargeometryclass(geometry.geometry)
                geometry.save()
        else:
            geom_equiv_IDs = [g.id_geometry1 for g in geom_equiv_result] + [g.id_geometry2 for g in geom_equiv_result]
            avaragegeometry = gt.buildavaragegeometry(geom_equiv_result)
            
            #update 
            if avaragegeometry:
                #save avaragegeometry in geometryclasses in avggeomcls
                avggeomcls = Geometryclasses(geometry = avaragegeometry.buildcmloutput())
                avggeomcls.calculate_md5()
                avggeomcls.save()
                #update Geometries
                qs = Geometries.alive_objects.filter(pk__in=geom_equiv_IDs)
                for g in qs:
                    g.geometryclass_id = avggeomcls.pk
                    g.save()
                    pass
            else:
                #error
                pass
            
    #f = open("outputclasses.txt", "w")
    #f.write(outputresult)
    #f.close()
    
    #print outputresult
    t = loader.get_template("utility/bugfixutility.html")
    c = Context({"result": outputresult})
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
#@transaction.commit_on_success
def utility_functional_classes(request, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    md5 = qchitool_tools.returnmd5(request.session.session_key)
    session_path = PATH_FOR_MOLECULE_FILES + str(md5) + "/"
    if not os.path.exists(session_path):
        os.makedirs(session_path)
    TheoryLevels.alive_objects.all().update(xcclass = modeltools.dummyobject.xcclasses)
    Xcclasses.objects.exclude(pk = modeltools.dummyobject.xcclasses.pk).delete()
    ft.build_functional_classes()
    t = loader.get_template("utility/bugfixutility.html")
    c = Context({})
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
#@transaction.commit_on_success
def utility_delete_ionizationenergies(request, PATH_FOR_MOLECULE_FILES):
    IonisationEnergies.alive_objects.exclude(pk=0).delete()
    t = loader.get_template("utility/bugfixutility.html")
    c = Context({})
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
#@transaction.commit_on_success
def utility_calculate_ionizationenergies(request, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    md5 = qchitool_tools.returnmd5(request.session.session_key)
    session_path = PATH_FOR_MOLECULE_FILES + str(md5) + "/"
    if not os.path.exists(session_path):
        os.makedirs(session_path)
    iet.build_ionizationenergies(None)
    t = loader.get_template("utility/bugfixutility.html")
    c = Context({})
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

@login_required(login_url='/accounts/login/')
#@transaction.commit_on_success
def utility_rebuild_elementlist_for_molecules(request, PATH_FOR_MOLECULE_FILES):
    post_data = request.POST.copy()
    md5 = qchitool_tools.returnmd5(request.session.session_key)

    index = 0
    calculations = Calculations.objects.all()
    ElementSpeciesBasisSet.objects.all().delete()
    ElementSpecies.objects.all().delete()    
    for c in calculations:
        #if c.pk == 1163:
        #    pass
        pt.p("Next Calculus " + str(c.pk))
        output_file, input_file, database_file = c.returntemporaryfiles()
        tasks_in_DB = c.tasks_set.all()
        pt.p("start calculation objects")
        calculation_objects = check_import_object.return_parsed_object_for_session(None, output_file.name, input_file.name, database_file.name, None, is_dbfile_humanreadable = True)
        pt.p("end calculation objects")

        #for task_index in range(len(calculation_objects["calculation_object"]["Tasks"])):
        #result = findrelatedtasks(tasks_in_DB.filter(name = calculation_objects["calculation_object"]["Tasks"][task_index].name), calculation_objects["calculation_object"]["Tasks"], calculation_objects["calculation_object"]["ElectronicStates"], calculation_objects["calculation_object"]["Geometries"])
        result = findrelatedtasks(tasks_in_DB, calculation_objects["calculation_object"]["Tasks"], calculation_objects["calculation_object"]["ElectronicStates"], calculation_objects["calculation_object"]["Geometries"])
        for (t_db, t_calc) in result.items():
            electronicstates = ElectronicStates.objects.filter(task = t_db)
            if len(electronicstates)==1:
                electronicstate = electronicstates[0]
                create_elementspeciesbasisset(task = Tasks.objects.get(pk = t_db), elementspecies = create_elementspecies(molecularspecie = electronicstate.species, elements =  calculation_objects["calculation_object"]["Elements"][t_calc]), elements =  calculation_objects["calculation_object"]["Elements"][t_calc])
            else:
                #search corresponding electronicstates
                #TO DO
                pass
        output_file.close()
        input_file.close()
        database_file.close()
        index += 1
    t = loader.get_template("utility/bugfixutility.html")
    c = Context({})
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))


def build_tasks_forms(calculation_object, selected_tasks):
    forms = {}
    if "Calculations" in calculation_object:

        calculation_found, first_calculation_found, calculation_found_tasks_list = modeltools.count_calculation_by_md5(calculation_object["Calculations"].output_md5)
        if calculation_found:
            #get Calculation from DB
            forms["Calculations"] = modelforms.CalculationsForm(instance = Calculations.objects.get(pk = first_calculation_found), 
                                                     prefix = "Calculations")
            
        else:
            #complete calculation_object
            calculation_object["Calculations"].qual_index = 0
            forms["Calculations"] = modelforms.CalculationsForm(instance = calculation_object["Calculations"], 
                                                     prefix = "Calculations")

        if  "ChemistryCodes" in calculation_object:
            all_chemistry_codes = ChemistryCodes.alive_objects.filter(name = calculation_object["ChemistryCodes"].name, 
                                                                version = calculation_object["ChemistryCodes"].version).order_by("-qual_index")
            if len(all_chemistry_codes) > 0: 
                #chemistri code found in DB
                calculation_object["ChemistryCodes"] = all_chemistry_codes[0]
            else:
                #complete new chemistrycode_object
                calculation_object["ChemistryCodes"].qual_index = 0
                #return HttpResponse("Complete ChemistriCode Info")
            forms["ChemistryCodes"] = modelforms.ChemistryCodesForm(instance = calculation_object["ChemistryCodes"],
                                                        prefix = "ChemistryCodes")
        
        #if  ("Tasks" in calculation_object) and ("TheoryLevels" in calculation_object):
        forms.update(build_form.BuildTasksForm(selected_tasks))
        for task_index in range(len(selected_tasks)):
            task = selected_tasks[task_index]
            
            build_form.BuildMolecularSpeciesandRelatedForm(task[1].molecular_specie, 
                                                              task[1].elements, 
                                                              task[1].electronicstate, 
                                                              task[1].dipolemoment,
                                                              task_index,
                                                              forms)
            build_form.BuildGeometriesForm(task[1].geometry, 
                                           calculation_object["optimization_convergence"][task[0]],
                                           task_index, 
                                           forms)
        
            #if task.vibrationanalisisarmonic:
            if task[1].rotationalconstants:
                build_form.BuildVibrationalAnalysisHarmonic(task[1].vibrationanalisisarmonic, 
                                                             task[1].rotationalconstants, 
                                                             calculation_object["TabulatedVibrations"][task[0]],
                                                             task_index,
                                                             forms)
            else:
                build_form.BuildVibrationalAnalysisHarmonic(task[1].vibrationanalisisarmonic, 
                                                             None, 
                                                             calculation_object["TabulatedVibrations"][task[0]],
                                                             task_index,
                                                             forms)
                
            
        tasks_forms_obj = commonsobjects.Tasks_Forms(forms["Calculations"], 
                                                     forms["ChemistryCodes"], 
                                                     forms["Tasks"], 
                                                     forms["TheoryLevels"], 
                                                     forms["Elements"], 
                                                     forms["BasisSets"], 
                                                     forms["MolecularSpecies"], 
                                                     forms["ElectronicStates"], 
                                                     forms["DipoleMoments"], 
                                                     forms["Geometries"], 
                                                     forms["RotationalConstants"], 
                                                     forms["VibrationalAnalysesHarmonic"], 
                                                     forms["TabulatedVibrations"])
        
    return forms, tasks_forms_obj

def associate_task_operation(selected_tasks, request):
    results = []
    tasks = request.session["tasks"]
    request.session["selected_tasks"] = []
    for task_sel in selected_tasks:
        task_sel_tmp = qchitool_tools.val(task_sel)
        current_selected_task = tasks[task_sel_tmp]
        request.session["selected_tasks"].append([task_sel_tmp, current_selected_task])
        if current_selected_task.task.name.upper() == "OPTIMIZE":
            #geometry optimization
            results.append([request.session["calculation_object"]["optimization_convergence"][task_sel_tmp], current_selected_task])
        elif current_selected_task.task.name.upper() == "ENERGY":
            #energy calculation
            results.append([request.session["calculation_object"]["energy_analisys"][task_sel_tmp], current_selected_task])
            pass
        elif current_selected_task.task.name.upper() == "FREQUENCIES" or current_selected_task.task.name.upper() == "FREQ":
            #frequency calculation
            results.append([request.session["calculation_object"]["frequency_analisys"][task_sel_tmp], current_selected_task])
    return results

def findrelatedtasks(databasetasks, calculatedtasks, calculatedelectronicstates, calculatedgeometries):
    #databasetasks -> queryset of Tasks extracted from database
    #calculatedtasks -> list of Tasks extracted from calculus
    #calculatedelectronicstates -> list of electronic states extracted from calculus
    #calculatedgeometries -> list of geometries  extracted from calculus
    result = {}
    #group databasetasks by task name
    databasetasks_group = {} #dictionary {}
    for task in databasetasks:
        if task.name in databasetasks_group:
            databasetasks_group[task.name].append(task)
        else:
            databasetasks_group[task.name] = [task]
    #group calculatedtasks by task name    
    calculatedtasks_group = {} #dictionary {database_task.name: [calculatedtasks_index, calculatedtasks_index]
    index = 0
    for task in calculatedtasks:
        if task.name in calculatedtasks_group:
            calculatedtasks_group[task.name].append(index)
        else:
            calculatedtasks_group[task.name] = [index]
        index += 1
        
    for (key,value) in databasetasks_group.items():
        if len(value) == len(calculatedtasks_group[key]):
            if len(value) == 1:
                result[value[0].pk] = calculatedtasks_group[key][0]
            else:
                #more than one task type for this calculation
                #use the geometry to find the correct task
                found2 = False
                for task2_index in calculatedtasks_group[key]:
                    #search the geometry of the first electronicstate 
                    if value[0].electronicstates_set.all()[0].geom.geometry_md5 == calculatedgeometries[task2_index].geometry_md5 and float(value[0].electronicstates_set.all()[0].total_energy) == calculatedelectronicstates[task2_index].total_energy:
                        result[value[0].pk] = task2_index
                        found2 = True
                        break
                if not found2:
                    #Error: Task Not Found
                    pass
        elif len(value) > len(calculatedtasks_group[key]):
            #ERROR more calculus in database than in calculus
            pass
        elif len(value) < len(calculatedtasks_group[key]):
            if len(value) == 1:
                #search the corresponding geometry
                found2 = False
                for task2_index in calculatedtasks_group[key]:
                    #search the geometry of the first electronicstate 
                    electronicstates1 = value[0].electronicstates_set.all()
                    electronicstates2 = calculatedelectronicstates[task2_index]
                    if len(electronicstates1) == 1:
                        if electronicstates1[0].geom.geometry_md5 == calculatedgeometries[task2_index].geometry_md5 and float(electronicstates1[0].total_energy) == electronicstates2.total_energy:
                            result[value[0].pk] = task2_index
                            found2 = True
                            break
                    else:
                        pass
                if not found2:
                    #Error: Task Not Found
                    pass
            else:
                pass
                #
    
    #return dictionary of {id_task: task_index of calculatedtasks}
    return result
    
def show_molecule_information(request, id, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    try:
        moleculeid = int(id)
    except ValueError:
        raise Http404()
    molecule = MolecularSpecies.alive_objects.filter(pk = moleculeid)
    if len(molecule)>0:
        molecules = MolecularSpecies.alive_objects.filter(name = molecule[0].name, formula = molecule[0].formula).order_by("charge")
        for mol in molecules:
            mol.electronicstategroup = mol.electronicstates_set.all().order_by("task__thlevel__name", "task__thlevel__xcclass__name", "description", "multiplicity", "total_energy") #"basissetshort" ,
            for es in mol.electronicstategroup:
                esid = es.pk
                ion_en = IonisationEnergies.alive_objects.filter(start_state = es)
                if ion_en.count()>0:
                    es.startionisationenergies =  ion_en
                
    moleculefilename = "m" + str(id) + str(esid) + ".cml"
    if not os.path.isfile(URL_FOR_MOLECULE_FILES + moleculefilename):
        electronicstate = ElectronicStates.alive_objects.filter(pk = esid)[0]
        electronicstate.geom.returncmlstructurefile(PATH_FOR_MOLECULE_FILES + moleculefilename)
        
    c = Context({"molecule":molecule[0], 
                 "molecules": molecules, 
                 "roundto": 6,
                 "appletsize": 300,
                 "moleculefilename" : moleculefilename,
                 "URL_FOR_MOLECULE_FILES" : URL_FOR_MOLECULE_FILES,
                 "pageorigin": "/"})
    c.update(commons.nwchem_units)
    t = loader.get_template("explore/base_energy.html")
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

def show_elemetbasisset_information(request, mid, esid, tid, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    try:
        moleculeid = int(mid)
        electronicstateid = int(esid)
        taskid = int(tid)
    except ValueError:
        raise Http404()
    task = Tasks.alive_objects.filter(pk = taskid)[0]
    molecule = MolecularSpecies.alive_objects.filter(pk = moleculeid)[0]
    electronicstate = ElectronicStates.alive_objects.filter(pk = esid)[0]
    electronicstate.basissets = task.elementspeciesbasisset_set.all()
    molecule.electronicstategroup = [electronicstate]
    molecules = [molecule]


    moleculefilename = "m" + str(mid) + str(esid) + ".cml"
    if not os.path.isfile(URL_FOR_MOLECULE_FILES + moleculefilename):
        electronicstate = ElectronicStates.alive_objects.filter(pk = esid)[0]
        electronicstate.geom.returncmlstructurefile(PATH_FOR_MOLECULE_FILES + moleculefilename)

    c = Context({"moleculeid": moleculeid, 
                 "molecule":molecules[0],
                 "molecules": molecules, 
                 "roundto": 6,
                 "appletsize": 300,
                 "moleculefilename" : moleculefilename,
                 "URL_FOR_MOLECULE_FILES" : URL_FOR_MOLECULE_FILES,
                 "pageorigin": "/database/molecule/" + str(moleculeid) + "/"})
    c.update(commons.nwchem_units)
    t = loader.get_template("explore/base_basissetselement.html")
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))
            
        
def show_vibrationalanalysesanarmonic_information(request, mid, esid):
    try:
        moleculeid = int(mid)
        electronicstateid = int(esid)
    except ValueError:
        raise Http404()
    vibrationalanalysesanarmonics = VibrationalAnalysesAnarmonic.alive_objects.filter(state = electronicstateid)

    c = Context({"vibrationalanalysesanarmonics": vibrationalanalysesanarmonics, "pageorigin": "/database/molecule/" + str(moleculeid) + "/"})
    t = loader.get_template("explore/vibrationalanalysesanarmonic.html")
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

def show_vibrationalanalysesharmonic_information(request, mid, esid, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    try:
        moleculeid = int(mid)
        electronicstateid = int(esid)
    except ValueError:
        raise Http404()
    
    
    molecule = MolecularSpecies.alive_objects.filter(pk = moleculeid)[0]
    electronicstate = ElectronicStates.alive_objects.filter(pk = esid)[0]
    electronicstate.vibrationalanalysesharmonics = VibrationalAnalysesHarmonic.alive_objects.filter(state = electronicstateid)
    for va in electronicstate.vibrationalanalysesharmonics:
        va.vibrationsgroup = va.tabulatedvibrations_set.all()
    molecule.electronicstategroup = [electronicstate]
    molecules = [molecule]
    
    
    #Load from DB the calculation output file
    
    moleculefilename = "m" + str(mid) + str(esid) + ".out"
    if not os.path.isfile(URL_FOR_MOLECULE_FILES + moleculefilename):
        electronicstate.task.calc.saveoutputfilewithoutinputdeck(PATH_FOR_MOLECULE_FILES + moleculefilename)
    
    
    c = Context({"moleculeid": moleculeid, 
                 "molecule":molecules[0],
                 "molecules": molecules, 
                 "pageorigin": "/database/molecule/" + str(moleculeid) + "/",
                 "roundto": 5, 
                 "appletsize": 300,
                 "moleculefilename" : moleculefilename,
                 "URL_FOR_MOLECULE_FILES" : URL_FOR_MOLECULE_FILES})
    c.update(commons.nwchem_units)
    t = loader.get_template("explore/base_vibrationanalisis.html")
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

def show_vibrationalanalysesharmoniceigenvectors_information(request, mid, esid, vaid, vatid, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    try:
        moleculeid = int(mid)
        electronicstateid = int(esid)
        vibrationalalisisid = int(vaid)
        vibrationid = int(vatid)
    except ValueError:
        raise Http404()
    
    molecule = MolecularSpecies.alive_objects.filter(pk = moleculeid)[0]
    electronicstate = ElectronicStates.alive_objects.filter(pk = esid)[0]
    vibrationanalisisharmonic = VibrationalAnalysesHarmonic.alive_objects.filter(state = electronicstateid)[0]
    vibration = TabulatedVibrations.objects.filter(pk = vibrationid)[0]
    vibrationanalisisharmonic.vibrationsgroup = [vibration]
    for v in vibrationanalisisharmonic.vibrationsgroup:
        v.eigenvectorsgroup = parser_tools.chunks(eval(v.eigenvectors), 3)
    electronicstate.vibrationalanalysesharmonics = [vibrationanalisisharmonic]
    molecule.electronicstategroup = [electronicstate]
    molecules = [molecule]
    

    moleculefilename = "m" + str(mid) + str(esid) + ".out"
    if not os.path.isfile(URL_FOR_MOLECULE_FILES + moleculefilename):
        electronicstate = ElectronicStates.alive_objects.filter(pk = electronicstateid)[0]
        electronicstate.task.calc.saveoutputfilewithoutinputdeck(PATH_FOR_MOLECULE_FILES + moleculefilename)
    
    c = Context({"moleculeid": moleculeid, 
                 "molecule":molecules[0],
                 "molecules": molecules, 
                 "roundto": 5, 
                 "appletsize": 300,
                 "moleculefilename" : moleculefilename,
                 "URL_FOR_MOLECULE_FILES" : URL_FOR_MOLECULE_FILES,
                 "pageorigin": "/database/vibrationalanalysesharmonic/" + str(moleculeid) + "/" + str(electronicstateid) + "/"})
    c.update(commons.nwchem_units)
    t = loader.get_template("explore/base_vector.html")
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))

def show_ionisationenergies_information(request, mid, esid, ieid, URL_FOR_MOLECULE_FILES, PATH_FOR_MOLECULE_FILES):
    try:
        moleculeid = int(mid)
        electronicstateid = int(esid)
        ionisationenergyid = int(ieid)
    except ValueError:
        raise Http404()
    molecule = MolecularSpecies.alive_objects.filter(pk = moleculeid)[0]
    electronicstate = ElectronicStates.alive_objects.filter(pk = esid)[0]
    electronicstate.ionisationenergy = [IonisationEnergies.alive_objects.get(pk=ionisationenergyid)]
    molecule.electronicstategroup = [electronicstate]
    molecules = [molecule]
    

    moleculefilename = "m" + str(mid) + str(esid) + ".out"
    if not os.path.isfile(URL_FOR_MOLECULE_FILES + moleculefilename):
        electronicstate = ElectronicStates.alive_objects.filter(pk = electronicstateid)[0]
        electronicstate.task.calc.saveoutputfilewithoutinputdeck(PATH_FOR_MOLECULE_FILES + moleculefilename)

    
    
    c = Context({"moleculeid": moleculeid, 
                 "molecule":molecules[0],
                 "molecules": molecules, 
                 "roundto": 5, 
                 "appletsize": 300,
                 "moleculefilename" : moleculefilename,
                 "URL_FOR_MOLECULE_FILES" : URL_FOR_MOLECULE_FILES,
                 "pageorigin": "/database/molecule/" + str(moleculeid) + "/"})
    c.update(commons.nwchem_units)
    t = loader.get_template("explore/base_energy.html")
    c.update(csrf(request)) 
    return HttpResponse(t.render(c))
