'''
Created on 16/dic/2011

@author: asaba
'''

import os

from django.db import transaction


from ....models import Calculations, ChemistryCodes, Bibliography, ElectronicStates, ElectronicStatesBibliographies, VibrationalAnalysesHarmonic, VibrationalAnalysesHarmonicBibliographies
from ....parser.nwchem import check_import_object
from ....parser import commonsobjects
from ....modelforms import modelforms, build_form
from .... import performancetest as pt

@transaction.commit_on_success
def re_run_parsing(session_path):
    for calculation in Calculations.alive_objects.all():
        inp_filename = session_path + "tmp_inp.nw"
        inp = calculation.input
        out_filename = session_path + "tmp_out.out"
        out = calculation.output
        oth_filename = session_path + "tmp_oth.out"
        oth = calculation.other_output
        if not os.path.exists(session_path):
            os.makedirs(session_path)
        inp_file = open(inp_filename, "w+")
        inp_file.write(inp)
        inp_file.close()
        out_file = open(out_filename, "w+")
        out_file.write(out)
        out_file.close()
        oth_file = open(oth_filename, "w+")
        oth_file.write(oth)
        oth_file.close()
        session = check_import_object.return_parsed_object_for_session(None, out_filename, inp_filename, oth_filename, 
                                                session_path, 
                                                is_dbfile_humanreadable = True)
        results = []
        forms = {}
        task_sel_tmp = 0
        if "Calculations" in session["calculation_object"]:

                    #complete calculation_object
                    session["calculation_object"]["Calculations"].qual_index = 0
                    forms["Calculations"] = modelforms.CalculationsForm(instance = session["calculation_object"]["Calculations"], 
                                                             prefix = "Calculations")
        
                    if  "ChemistryCodes" in session["calculation_object"]:
                        all_chemistry_codes = ChemistryCodes.alive_objects.filter(name = session["calculation_object"]["ChemistryCodes"].name, 
                                                                            version = session["calculation_object"]["ChemistryCodes"].version).order_by("-qual_index")
                        if len(all_chemistry_codes) > 0: 
                            #chemistri code found in DB
                            session["calculation_object"]["ChemistryCodes"] = all_chemistry_codes[0]
                        else:
                            #complete new chemistrycode_object
                            session["calculation_object"]["ChemistryCodes"].qual_index = 0
                            #return HttpResponse("Complete ChemistriCode Info")
                        forms["ChemistryCodes"] = modelforms.ChemistryCodesForm(instance = session["calculation_object"]["ChemistryCodes"],
                                                                    prefix = "ChemistryCodes")
                    
        forms.update(build_form.BuildTasksForm(session["tasks"]))
        index = 0
        for current_selected_task in session["tasks"]:
            if current_selected_task.task.name.upper() == "OPTIMIZE":
                #geometry optimization
                results.append([session["calculation_object"]["optimization_convergence"][index], current_selected_task])
            elif current_selected_task.task.name.upper() == "ENERGY":
                #energy calculation
                results.append([session["calculation_object"]["energy_analisys"][index], current_selected_task])
                pass
            elif current_selected_task.task.name.upper() == "FREQUENCIES" or current_selected_task.task.name.upper() == "FREQ":
                #frequency calculation
                results.append([session["calculation_object"]["frequency_analisys"][index], current_selected_task])
            
            
            build_form.BuildMolecularSpeciesForm(current_selected_task.task.molecular_specie, 
                                                              current_selected_task.task.elements, 
                                                              current_selected_task.task.electronicstate, 
                                                              current_selected_task.task.dipolemoment,
                                                              index,
                                                              forms)
            build_form.BuildGeometriesForm(current_selected_task.task.geometry, 
                                           session["calculation_object"]["optimization_convergence"][index],
                                           index, 
                                           forms)
        
            if current_selected_task.task.rotationalconstants:
                build_form.BuildVibrationalAnalysisHarmonic(current_selected_task.task.vibrationanalisisarmonic, 
                                                             current_selected_task.task.rotationalconstants, 
                                                             session["calculation_object"]["TabulatedVibrations"][index],
                                                             index,
                                                             forms)
            else:
                build_form.BuildVibrationalAnalysisHarmonic(current_selected_task.task.vibrationanalisisarmonic, 
                                                             None, 
                                                             session["calculation_object"]["TabulatedVibrations"][index],
                                                             index,
                                                             forms)
                
            index += 1
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
        
            
@transaction.commit_on_success
def addBibliographyToElectronicStates(bib_id):
    b = Bibliography.alive_objects.get(bib_id = bib_id)
    for s in ElectronicStates.alive_objects.all():
        if len(ElectronicStatesBibliographies.alive_objects.filter(electronicstates = s, bibliography = b)) == 0:
            eb = ElectronicStatesBibliographies(electronicstates = s, bibliography = b)
            eb.save()
            
@transaction.commit_on_success
def addBibliographyToVibrationalAnalysesHarmonic(bib_id):
    b = Bibliography.alive_objects.get(bib_id = bib_id)
    for v in VibrationalAnalysesHarmonic.alive_objects.all():
        if len(VibrationalAnalysesHarmonicBibliographies.alive_objects.filter(vibrationalanalysesarmonic = v, bibliography = b)) == 0:
            vb =VibrationalAnalysesHarmonicBibliographies(vibrationalanalysesarmonic = v, bibliography = b)
            vb.save()