'''
Created on 02/mag/2011

@author: asaba
'''

import re
import pybel
import shlex, subprocess
import zlib
import tempfile

from ..tools import parser_tools
from ... import models
import commons
from ...tools import qchitool_tools, modeltools, geometry_tools
from ..impnwc import parsing
from ..commonsobjects import MoleculeInfo, TaskCode, myElement
from ... import performancetest as pt
#from ..tools import qchitool_tools

     
     
def parse_basis_input_directive(basis_rows, basis_sets):
    #parsing of basis set section in input file or "echo of input deck"
    re_library = re.compile(r".*library.*")
    for row in basis_rows:
        if re_library.search(row):
            basissetname = row.split("library",1)[1].strip().split()[0].strip()
            if not basissetname.upper() in commons.nwchem_task_basis_sets:
                basissetname = "ERROR"
            basis_sets.append(models.BasisSets(name = basissetname, 
                                   description = commons.nwchem_task_basis_sets[basissetname.upper()]))  

def check_input(input_lines, theorylevel, tasks):
    #get Input Deck rows or input file rows and extract theory levels and tasks lists
    tasks_list = []
    #basis_rows = []
    #started_basis_section = False
    #for InputDeckRow in input_lines:
    #    if re.match(".*BASIS[\s]*.*", InputDeckRow.upper()): 
    #        # BASIS [<string name default "ao basis">] \
    #        #[(spherical || cartesian) default cartesian] \
    #        #[(segment || nosegment) default segment] \
    #        #[(print || noprint) default print]
    #        started_basis_section = True
    #        continue
    #    if started_basis_section:
    #        if re.match("[\s]*end[\s]*", InputDeckRow.lower()):
    #            break
    #        else:
    #            basis_rows.append(InputDeckRow)
    #parse_basis_input_directive(basis_rows, basis_sets)
        

            #basis[-1].fullprint()
    re_task = re.compile(r"[\s]*TASK.*")
    for InputDeckRow in input_lines:
        if re_task.search(InputDeckRow.upper()):
            tasks_list.append(InputDeckRow.split(None,1)[1])
            theory = InputDeckRow.split(None,1)[1].split()[0]
            if len(InputDeckRow.split(None,1)[1].split())>1:
                operation = InputDeckRow.split(None,1)[1].split()[1]
            else:
                operation = "ENERGY" #default operation value
                
            if theory.upper() in commons.nwchem_task_theory:
                if not (operation.upper() in commons.nwchem_task_operation):
                    operation = "ERROR"
            else:
                theory = "ERROR"
                if not (operation.upper() in commons.nwchem_task_operation):
                    operation = "ERROR"
            theorylevel.append(models.TheoryLevels(name=theory, 
                                                   description = commons.nwchem_task_theory[theory.upper()]))
            tasks.append(models.Tasks(name = operation, 
                                      description = commons.nwchem_task_operation[operation.upper()]))
    ##all tasks use the same basis set
    #while len(tasks)>len(basis_sets):
    #    basis_sets.append(basis_sets[-1])
    
def create_human_readable_database_file(restart_name, session_path):
    #create db output file
    #1. build the new input file on temporary directory
    #require a working version of nwchem on server
    if session_path:
        session_path_tmp = session_path
    else:
        session_path_tmp = tempfile.gettempdir() + "/"
    previous_calculus_name = restart_name[:-1] #remove dot
    new_input_file_name = session_path +   "_input_" +  restart_name + "nw"
    new_output_database_file_name = session_path +   "_database_" + restart_name + "out"
    new_input_file = open(new_input_file_name, "w")
    new_output_file = open(new_output_database_file_name, "w")
    new_input_file.writelines(["permanent_dir " + session_path +  "  \n",
                               "scratch_dir " + session_path +  "  \n",
                               "restart " + previous_calculus_name + "\n",
                               "task rtdbprint"])
    new_input_file.close()
    #2. execute nwchem with new input file and save the output in a new database output file
    args = shlex.split(commons.nwchem_bin_path + " " + new_input_file_name)
    subprocess.call(args, stdout = new_output_file)
    new_output_file.close()        
    return new_output_database_file_name


                              

def get_geometry_from_DB(restart_name, geometry_name, basis_set_name, session_path, new_output_database_file_name):
    pt.p("start")
    if session_path:
        session_path_tmp = session_path
    else:
        session_path_tmp = tempfile.gettempdir() + "/"
    
    #3. parse the database output file
    db_parsing_result = parsing(new_output_database_file_name, True, "NwChemDatabase") 
    #4 search geometry and basis set information
    section_GeometryInDatabase = parser_tools.returnsection("GeometriesInTheDatabase", db_parsing_result, "First")
    section_BasisSetInDatabase = parser_tools.returnsection("BasisSetsInTheDatabase", db_parsing_result, "First")
    #use the default geometry in database CHECK
    if section_GeometryInDatabase:
        geometry_name = section_GeometryInDatabase.section_object.values["database_geometry_name"]
    #use the default basis set name in database
    if section_BasisSetInDatabase:
        basis_set_name = section_BasisSetInDatabase.section_object.values["database_basis_set_name"]
    
    section = parser_tools.returnsection("RTDBDatabase", db_parsing_result, "FIRST")
    sections_databaseentry = parser_tools.returnsection("DatabaseEntry", section, "LIST")
    coordinate_dataentry = None
    tags_dataentry = None
    mass_dataentry = None
    charges_dataentry = None
    charge_value = 0
    unit_conversion_value = None
    multiplicity_value = None
    basis_set_tags_number = None
    basis_set_standard_name_list = None
    basis_set_tags_name_list = None
    star_basis_set_type = None
    for section in sections_databaseentry:
        if section.section_object.values["entry_name"]:
            dataentry_name = section.section_object.values["entry_name"]
            if len(dataentry_name)>1:
                #search coordinates
                if dataentry_name[0] == "geometry" and dataentry_name[1] == geometry_name and dataentry_name[2] == "coords":
                    coordinate_dataentry = section
                    continue
                #search tags (atoms names - usually element symbol)
                if dataentry_name[0] == "geometry" and dataentry_name[1] == geometry_name and dataentry_name[2] == "tags":
                    tags_dataentry = section
                    continue
                #search mass
                if dataentry_name[0] == "geometry" and dataentry_name[1] == geometry_name and dataentry_name[2] == "masses":
                    mass_dataentry = section
                    continue
                #search charge (atomic number)
                if dataentry_name[0] == "geometry" and dataentry_name[1] == geometry_name and dataentry_name[2] == "charges":
                    charges_dataentry = section
                    continue
                if dataentry_name[0] == "geometry" and dataentry_name[1] == geometry_name and dataentry_name[2] == "angstrom_to_au":
                    unit_conversion_value = float(section.section_object.values["values"][0])
                    continue
                if dataentry_name[0] == "dft" and dataentry_name[1] == "mult":
                    multiplicity_value = int(section.section_object.values["values"][0])
                    continue
                if not (basis_set_name is None):
                    if dataentry_name[0] == "basis" and dataentry_name[1] == basis_set_name and dataentry_name[2] == "bs_nr_tags":
                        basis_set_tags_number = int(section.section_object.values["values"][0])
                        continue
                    if dataentry_name[0] == "basis" and dataentry_name[1] == basis_set_name and dataentry_name[2] == "bs_stdname":
                        basis_set_standard_name_list = section
                        continue
                    if dataentry_name[0] == "basis" and dataentry_name[1] == basis_set_name and dataentry_name[2] == "bs_tags":
                        basis_set_tags_name_list = section
                        continue
                    if dataentry_name[0] == "basis" and dataentry_name[1] == basis_set_name and dataentry_name[2] == "star bas type":
                        star_basis_set_type = section
                        continue

            elif len(dataentry_name) == 1 and dataentry_name[0] == "charge":
                charge_value = int(section.section_object.values["values"][0])
                    
    coord = zip(*[iter(coordinate_dataentry.section_object.values["values"])]*3)
    
        #This function extract structural information from nwchem output file
    #Moreover it makes inchi and inchikey 
    #Return a MoleculeInfo object
        #build a dictionary for basis set used (element tag: basis set standard name)
    basis_set_dictionary = {}
    if not (basis_set_name is None):
        if basis_set_tags_number == 0:
            #it was used a basis set for each element of molecule
            for tag in tags_dataentry.section_object.values["values"]:
                basis_set_dictionary[tag] = star_basis_set_type.section_object.values["values"][0]
        else:
            #different basis set for molecule tag (usually a basis set for each element tag name - ex. each H same basis set)
            index = 0
            for tag in basis_set_tags_name_list.section_object.values["values"]:
                basis_set_dictionary[tag] = basis_set_standard_name_list.section_object.values["values"][index]
                index +=1
    molecule = pybel.ob.OBMol()      
    #pybelmol = pybel.Molecule(molecule)
    atom_basis_set = {}
    geometry_list = []
    for atom_index in range(len(tags_dataentry.section_object.values["values"])):
        atom = molecule.NewAtom()
        atom.SetVector(float(coord[atom_index][0])/unit_conversion_value, float(coord[atom_index][1])/unit_conversion_value, float(coord[atom_index][2])/unit_conversion_value)
        atom.SetAtomicNum(int(charges_dataentry.section_object.values["values"][atom_index]))
        #build the geometry list with [index, Tag, AtomicMass, X, Y, Z]
        geometry_list.append([atom.GetIdx(), tags_dataentry.section_object.values["values"][atom_index], int(charges_dataentry.section_object.values["values"][atom_index]), float(coord[atom_index][0])/unit_conversion_value, float(coord[atom_index][1])/unit_conversion_value, float(coord[atom_index][2])/unit_conversion_value])
        if not (basis_set_name is None):
            if "*" in basis_set_dictionary: #same basis set for each atom
                atom_basis_set[atom.GetIdx()] = basis_set_dictionary["*"]
            else:
                atom_basis_set[atom.GetIdx()] = basis_set_dictionary[tags_dataentry.section_object.values["values"][atom_index]]
        #CHEKC
        #pybelmol.atoms[atom_index].exactmass = mass_dataentry.section_object.values["values"][atom_index]
    molecule.ConnectTheDots()
    #get first molecule information found
    xyz = pybel.Molecule(molecule).write("xyz")
    cml = pybel.Molecule(molecule).write("cml")
    molecule.PerceiveBondOrders()
    sdf = pybel.Molecule(molecule).write("sdf")
    re_name = re.compile(r".*[<]name[>].*[<][/]name[>]\n")
    re_molecule = re.compile(r".*[<]molecule.*")
    re_molecule2 = re.compile(".*[<]molecule.*[>]")
    for row in cml:
        #remove molecule name
        row = re_name.sub("", row)
        #remove molecule id
        if re_molecule.match(row):
            row = re_molecule2.sub("<molecule>", row)
    xyz_file_name = parser_tools.writefile(xyz, restart_name, session_path_tmp, "xyz")
    cml_file_name = parser_tools.writefile(cml, restart_name, session_path_tmp, "cml")
    sdf_file_name = parser_tools.writefile(sdf, restart_name, session_path_tmp, "sdf")
    
    #molecule = pybel.readfile("xyz", str(xyz_file_name)).next()
    uncharged_molecule = pybel.Molecule(molecule)
    pybel.Molecule(molecule).atoms[0].OBAtom.SetFormalCharge(charge_value)
    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("inchi") 
    conv.SetOptions("K", conv.OUTOPTIONS)
    #pybelmol.GetAtom(1).SetFormalCharge(charge)
    
    inchikey = conv.WriteString(molecule)
    inchi = pybel.Molecule(molecule).write("inchi").split("=")[1]

    molecule_info = MoleculeInfo(inchi, inchikey, uncharged_molecule, 
                                 multiplicity_value, charge_value,
                                 geometry_list,
                                 xyz_file_name, xyz,
                                 sdf_file_name, sdf,
                                 cml_file_name, cml,
                                 atom_basis_set)
    pt.p("end")
    return molecule_info

def check_geometry_optimization_convergence(parsing_result):
    #check geometry optimization convergence
    #return the step information
    optimization_convergence_step = {}
    re_DftGradientModule = re.compile(r"DftGradientModule")
    re_OptimizationConverged = re.compile(r"OptimizationConverged")
    re_StepInfo = re.compile(r"StepInfo")
    re_Geometry = re.compile(r"Geometry")
    re_NwchemScfModule = re.compile(r"NwchemScfModule")
    
    for supersection in parsing_result:
        for subsection in supersection.list_of_subsections:
            #search section DftGradientModule_supersection in supersection Step
            if re_DftGradientModule.search(subsection.section_class):
                for sub_subsection in subsection.list_of_subsections:
                    if re_OptimizationConverged.search(sub_subsection.section_class):
                        for sub_subsection in subsection.list_of_subsections:
                            if re_StepInfo.search(sub_subsection.section_class):
                                optimization_convergence_step = supersection.section_object.values
                                optimization_convergence_step.update(subsection.section_object.values)
                                optimization_convergence_step.update(sub_subsection.section_object.values)
                            elif re_Geometry.search(sub_subsection.section_class):
                                geometry_optimization_convergence_step = sub_subsection
            elif re_NwchemScfModule.search(subsection.section_class):
                for sub_subsection in subsection.list_of_subsections:
                    if re_OptimizationConverged.search(sub_subsection.section_class):
                        for sub_subsection in subsection.list_of_subsections:
                            if re_StepInfo.search(sub_subsection.section_class):
                                optimization_convergence_step = supersection.section_object.values
                                optimization_convergence_step.update(subsection.section_object.values)
                                optimization_convergence_step.update(sub_subsection.section_object.values)
                            elif re_Geometry.search(sub_subsection.section_class):
                                geometry_optimization_convergence_step = sub_subsection
    if optimization_convergence_step:
        for supersection in parsing_result:
            step_section = parser_tools.returnsection("Step", supersection, "LAST")
            if step_section:                
            #if re.search("Step", supersection.section_class):
                if step_section.section_object.values["Step"] == optimization_convergence_step["Step"]: #StepInfo values
                    optimization_convergence_step["section"] = step_section #Step sections
                    optimization_convergence_step["geometry"] = geometry_optimization_convergence_step
                    break
                            
    return optimization_convergence_step

def check_property_analisys(parsing_result):
    #search the dipole information
    property_analisys_section = {}
    NwchemPropertyModule_section = parser_tools.returnsection("NwchemPropertyModule", parsing_result, "LAST")
    if NwchemPropertyModule_section:
        property_analisys_section["section"] = NwchemPropertyModule_section
        #serach basis set information
        SummaryOf_section = parser_tools.returnsection("SummaryOf", NwchemPropertyModule_section, "FIRST")
        if SummaryOf_section:
            property_analisys_section["basis_sets"] = SummaryOf_section.section_object.values["data"]
    return  property_analisys_section


def check_frequency_analisys_armonic(parsing_result):
    #search the frequency analisys information
    #return the section information
    frequency_analisys_section = {}
    VibrationalAnalysisViaTheFXMethod_sections = parser_tools.returnsection("VibrationalAnalysisViaTheFXMethod", parsing_result, "LIST")
    if VibrationalAnalysisViaTheFXMethod_sections:
        for sec in VibrationalAnalysisViaTheFXMethod_sections:
            if sec.section_object.values["project_out"]:
                frequency_analisys_section["data"] = sec.section_object.values
                CenterOfMass_section = parser_tools.returnsection("CenterOfMass", VibrationalAnalysisViaTheFXMethod_sections, "FIRST")
                if CenterOfMass_section:
                    frequency_analisys_section["CenterOfMass"] = CenterOfMass_section.section_object.values
                MomentOfInertia_section = parser_tools.returnsection("MomentOfInertia", VibrationalAnalysisViaTheFXMethod_sections, "FIRST")
                if MomentOfInertia_section:
                    frequency_analisys_section["MomentOfInertia"] = MomentOfInertia_section.section_object.values
                Rot_Const_And_Energy_section = parser_tools.returnsection("Rot_Const_And_Energy", VibrationalAnalysisViaTheFXMethod_sections, "FIRST")
                if Rot_Const_And_Energy_section:
                    frequency_analisys_section["Rot_Const_And_Energy"] = Rot_Const_And_Energy_section.section_object.values
                NormalEigenvalue_DipoleMoment_section = parser_tools.returnsection("NormalEigenvalue_DipoleMoment", VibrationalAnalysisViaTheFXMethod_sections, "LAST")
                if NormalEigenvalue_DipoleMoment_section:
                    frequency_analisys_section["DipoleMoment"] = NormalEigenvalue_DipoleMoment_section.section_object.values
                NormalEigenvalue_InfraRedIntensities_section = parser_tools.returnsection("NormalEigenvalue_InfraRedIntensities", VibrationalAnalysisViaTheFXMethod_sections, "LAST")
                if NormalEigenvalue_InfraRedIntensities_section:
                    frequency_analisys_section["InfraRedIntensities"] = NormalEigenvalue_InfraRedIntensities_section.section_object.values
                frequency_analisys_section["section"] = sec
                break
                
    return frequency_analisys_section

def check_energy_analisys(parsing_result):
    #search the energy information
    #return the section information
    energy_analisys_section = {}
    NwchemDftModule_section = parser_tools.returnsection("NwchemDftModule", parsing_result, "LAST")
    if NwchemDftModule_section:
        energy_analisys_section["section"] = NwchemDftModule_section
        #serach basis set information
        SummaryOf_section = parser_tools.returnsection("SummaryOf", NwchemDftModule_section, "FIRST")
        if SummaryOf_section:
            energy_analisys_section["basis_sets"] = SummaryOf_section.section_object.values["data"]
    return  energy_analisys_section

def return_similar_electronicstates(electronicstate):
    #This function search in the database molecules similar to 'molecule'
    #and return the list of this molecules
    no_ionisation_inchi = re.sub(r"/q.*", "", electronicstate.species.inchi) 
    electronicstates = models.ElectronicStates.alive_objects.filter(species__inchi__istartwith = no_ionisation_inchi)
    if len(electronicstates) > 0:
        return electronicstates
    else:
        return None
    
def build_calculation_object(original_output_file, out_file, inp_file, db_file, parsing_result, session_path, is_dbfile_humanreadable):
    
    pt.p("start")
    if session_path:
        session_path_tmp = session_path
    else:
        session_path_tmp = tempfile.gettempdir() + "/"
    TableDataFound = {}
    
    section_NwchemManifest = parser_tools.returnsection("NwchemManifest", parsing_result, "First")
    section_JobInformation = parser_tools.returnsection("JobInformation", parsing_result, "First")
    
    section_InputDeck = parser_tools.returnsection("InputDeck", parsing_result, "First")
    #section_InputModule = parser_tools.returnsection("NwchemInputModule", parsing_result, "List")
    
    outputfilecontent = None
    inputfilecontent = None
    dbfilecontent = None
    
    if out_file:
        outputfile = open(out_file, "r")
        outputfilecontent = outputfile.read()
        outputfile.close()
    if inp_file:
        inputfile = open(inp_file, "r")
        inputfilecontent = inputfile.read()
        inputfile = open(inp_file, "r")
        inputfilelines = inputfile.readlines()
        inputfile.close()

        
    if db_file:
        if is_dbfile_humanreadable:
            dbfilecontent = db_file
            human_readable_database_file_name = db_file
        else:
            human_readable_database_file_name = create_human_readable_database_file(section_JobInformation.section_object.values["prefix"], session_path_tmp)
            dbfile = open(human_readable_database_file_name, "r")
            dbfilecontent = dbfile.read()
            dbfile.close()
        


    #build calculation
    dummy_ChemistryCode = models.ChemistryCodes.alive_objects.get(name = "DUMMY")
    
    if section_NwchemManifest:
        TableDataFound["ChemistryCodes"] = models.ChemistryCodes(name = "Nwchem",
                                       version = section_NwchemManifest.section_object.values["version"])
    else:
        #Error on out_file
        TableDataFound["ChemistryCodes"] = dummy_ChemistryCode
        
    TableDataFound["Calculations"] = models.Calculations(code = dummy_ChemistryCode)
    TableDataFound["Calculations"].output = outputfilecontent
    TableDataFound["Calculations"].input = inputfilecontent
    TableDataFound["Calculations"].other_output = dbfilecontent #.decode('raw_unicode_escape').encode("utf-8")
    if original_output_file:
        TableDataFound["Calculations"].comments = original_output_file

    TableDataFound["Calculations"].calculate_md5()
    
    #get tasks and theorylevel information
    TableDataFound["Tasks"] = []
    TableDataFound["TheoryLevels"] = []
    TableDataFound["MolecularSpecies"] = []
    TableDataFound["ElectronicStates"] = []
    TableDataFound["DipoleMoments"] = []
    TableDataFound["Elements"] = [] #list of myElement
    TableDataFound["RotationalConstants"] = []
    TableDataFound["VibrationalAnalysesHarmonic"] = []
    TableDataFound["TabulatedVibrations"] = []
    TableDataFound["Geometries"] = []
    TableDataFound["molecule_info"] = []
    TableDataFound["optimization_convergence"] = []
    TableDataFound["energy_analisys"] = []
    TableDataFound["frequency_analisys"] = []
    TableDataFound["property_analisys"] = []
    TableDataFound["SwapOrbitals"] = []
    
    #Search Tasks, Basis and Theory information
    
    #if section_InputDeck:
    #    #read information by Input Deck section
    #    check_input(section_InputDeck.section_object.values["lines"], TableDataFound["TheoryLevels"], TableDataFound["Tasks"])
    #elif inp_file:
    #    #read information by input file
    #    check_input(inputfilelines, TableDataFound["TheoryLevels"], TableDataFound["Tasks"])
    #else:
    
    #infer information by sections present
    input_output_sections = parser_tools.return_inputoutput_sectionschema(parsing_result)
    #last_basis_set = None
    last_molecule_info = None
    last_charge = None
    last_spinmultiplicity = None
    for iosection in input_output_sections:
        iosection.extracttasks()
        t = iosection.theorylevel.upper()
        o = iosection.operation.upper()
        if iosection.xc_name:
            xc_name = iosection.xc_name.upper()
        else:
            xc_name = None
        xc_description = iosection.xc_description
        if not (t in commons.nwchem_task_theory):
            t = "ERROR"
        if not(o in commons.nwchem_task_operation):
            o = "ERROR"

        TableDataFound["TheoryLevels"].append(models.TheoryLevels(name=t, description = commons.nwchem_task_theory[t], xc_name = xc_name, xc_description = xc_description))
        TableDataFound["Tasks"].append(models.Tasks(name = o, description = commons.nwchem_task_operation[o]))
        
        TableDataFound["optimization_convergence"].append(check_geometry_optimization_convergence(iosection.out_section))
        TableDataFound["energy_analisys"].append(check_energy_analisys(iosection.out_section))
        TableDataFound["frequency_analisys"].append(check_frequency_analisys_armonic(iosection.out_section))
        TableDataFound["property_analisys"].append(check_property_analisys(iosection.out_section))
        
        section_GeneralInformation = parser_tools.returnsection("GeneralInformation", iosection.out_section, "First")
        #get charge and spin multiplicity
        if section_GeneralInformation:
            last_charge = qchitool_tools.val(section_GeneralInformation.section_object.values["Charge"])
            last_spinmultiplicity = qchitool_tools.val(section_GeneralInformation.section_object.values["Spin multiplicity"])
        #get orbital swap information
        section_SwapOrbitals = parser_tools.returnsection("SwappingOrbitals", iosection.out_section, "First")
        if section_SwapOrbitals:
            TableDataFound["SwapOrbitals"].append(section_SwapOrbitals.section_object.values)
        else:
            TableDataFound["SwapOrbitals"].append(None)
        
        #get geometry and basis sets
        if TableDataFound["optimization_convergence"][-1]:
            section_Geometry = TableDataFound["optimization_convergence"][-1]["geometry"]
        else:
            section_Geometry = parser_tools.returnsection("Geometry", iosection.inp_section, "Last")
        section_SummaryOf = parser_tools.returnsection("SummaryOf", iosection.inp_section, "Last")
        if section_Geometry is None:
            if section_SummaryOf is None:
                if last_molecule_info is None:
                    #get info from DB
                    section_GeometryInDatabase = parser_tools.returnsection("GeometriesInTheDatabase", parsing_result, "First")
                    section_BasisSetInDatabase = parser_tools.returnsection("BasisSetsInTheDatabase", parsing_result, "First")
                    
                    if section_BasisSetInDatabase is None or section_GeometryInDatabase is None:
                        pass
                        #Error Basis set information not found in database
                    else:
                        #get basis set information from database
                        iosection.molecule_info = get_geometry_from_DB(section_JobInformation.section_object.values["prefix"],
                                             section_GeometryInDatabase.section_object.values["database_geometry_name"],
                                             section_BasisSetInDatabase.section_object.values["database_basis_set_name"],
                                             session_path_tmp,
                                             human_readable_database_file_name) #basis set name unspecified
                else:
                    iosection.molecule_info = last_molecule_info
                last_molecule_info = iosection.molecule_info
            else:
                if last_molecule_info is None:
                    #get info from DB
                    section_GeometryInDatabase = parser_tools.returnsection("GeometriesInTheDatabase", parsing_result, "First")
                    section_BasisSetInDatabase = parser_tools.returnsection("BasisSetsInTheDatabase", parsing_result, "First")
                    
                    if section_BasisSetInDatabase is None or section_GeometryInDatabase is None:
                        pass
                        #Error Basis set information not found in database
                    else:
                        #get basis set information from database
                        iosection.molecule_info = get_geometry_from_DB(section_JobInformation.section_object.values["prefix"],
                                             section_GeometryInDatabase.section_object.values["database_geometry_name"],
                                             section_BasisSetInDatabase.section_object.values["database_basis_set_name"],
                                             session_path_tmp,
                                             human_readable_database_file_name) #basis set name unspecified
                else:
                    iosection.molecule_info = last_molecule_info
                    
                tmp_molecule_info = parser_tools.extract_molecule_info(out_file, 
                                                               last_charge, 
                                                               session_path_tmp,
                                                               section_SummaryOf.section_object.values["data"],
                                                               iosection.molecule_info.geometry_list,
                                                               last_spinmultiplicity,
                                                               1) #section_Geometry.section_object.values["unit_conversion_value"]          
                iosection.molecule_info.atom_basis_set = tmp_molecule_info.atom_basis_set
                last_molecule_info = iosection.molecule_info
            
        else: #section Geometry defined
            if section_SummaryOf is None:
                if last_molecule_info is None:
                    #get info from DB
                    section_GeometryInDatabase = parser_tools.returnsection("GeometriesInTheDatabase", parsing_result, "First")
                    section_BasisSetInDatabase = parser_tools.returnsection("BasisSetsInTheDatabase", parsing_result, "First")
                    
                    if section_BasisSetInDatabase is None or section_GeometryInDatabase is None:
                        pass
                        #Error Basis set information not found in database
                    else:
                        #get basis set information from database
                        iosection.molecule_info = get_geometry_from_DB(section_JobInformation.section_object.values["prefix"],
                                             section_GeometryInDatabase.section_object.values["database_geometry_name"],
                                             section_BasisSetInDatabase.section_object.values["database_basis_set_name"],
                                             session_path_tmp,
                                             human_readable_database_file_name) #basis set name unspecified
                else:
                    iosection.molecule_info = last_molecule_info
                tmp_molecule_info = parser_tools.extract_molecule_info(out_file, 
                                                               last_charge, 
                                                               session_path_tmp,
                                                               None,
                                                               section_Geometry.section_object.values["data"],
                                                               last_spinmultiplicity,
                                                               1) #section_Geometry.section_object.values["unit_conversion_value"]          
                tmp_molecule_info.atom_basis_set = iosection.molecule_info.atom_basis_set
                iosection.molecule_info = tmp_molecule_info
                last_molecule_info = iosection.molecule_info
            else: #section_SummaryOf defined
                iosection.molecule_info = parser_tools.extract_molecule_info(out_file, last_charge, session_path_tmp, section_SummaryOf.section_object.values["data"], section_Geometry.section_object.values["data"], last_spinmultiplicity, 1) #section_Geometry.section_object.values["unit_conversion_value"]
                last_molecule_info = iosection.molecule_info

        TableDataFound["MolecularSpecies"].append(models.MolecularSpecies(charge = iosection.molecule_info.charge))
        es = models.ElectronicStates()
        if TableDataFound["optimization_convergence"][-1]:
            es.multiplicity = iosection.molecule_info.spinmultiplicity
            es.total_energy = TableDataFound["optimization_convergence"][-1]["section"].section_object.values["Energy"]
        elif TableDataFound["energy_analisys"][-1]:
            es.multiplicity = iosection.molecule_info.spinmultiplicity
            es.total_energy = TableDataFound["energy_analisys"][-1]["section"].section_object.values["Total DFT energy"]
        else:
            es.multiplicity = iosection.molecule_info.spinmultiplicity
        #Add conventional description name by multiplicity 
        try: 
            es.description = commons.nwchem_electronicstates_description[es.multiplicity]
        except:
            es.description = "X"

        #Add indication about swapping
        if TableDataFound["SwapOrbitals"][-1]:
            es.description += "SWAP " + TableDataFound["SwapOrbitals"][-1]["spin_orbital"] + " " + TableDataFound["SwapOrbitals"][-1]["first_orbital"] + " " + TableDataFound["SwapOrbitals"][-1]["second_orbital"]
        TableDataFound["ElectronicStates"].append(es)


        section_NwchemGeometryOptimization = parser_tools.returnsection("NwchemGeometryOptimization", iosection.out_section, "First")

        section_DipoleMoment = parser_tools.returnsection("MultipoleAnalysisOfTheDensity", iosection.out_section, "Last")
        #Rotational Costants information
        section_Rot_Const_And_Energy = parser_tools. returnsection("Rot_Const_And_Energy", iosection.out_section, "Last")
        #vibrational analisis information
        section_NormalEigenvalue_DipoleMoment = parser_tools. returnsection("NormalEigenvalue_DipoleMoment", iosection.out_section, "Last")
        section_NormalEigenvalue_InfraRedIntensities = parser_tools. returnsection("NormalEigenvalue_InfraRedIntensities", iosection.out_section, "Last")
        section_NormalModeEigenvectors = parser_tools. returnsection("NormalModeEigenvectors", iosection.out_section, "Last")
        
        #check type of calculus (start or restart)
    
        if section_DipoleMoment:
            TableDataFound["DipoleMoments"].append(models.DipoleMoments(mu_x = section_DipoleMoment.section_object.values["levels"]["x"], 
                                          mu_y = section_DipoleMoment.section_object.values["levels"]["y"], 
                                          mu_z = section_DipoleMoment.section_object.values["levels"]["z"]))
        else:
            TableDataFound["DipoleMoments"].append(None)
        TableDataFound["Geometries"].append(models.Geometries(geometry=iosection.molecule_info.geometries["cml"]["content"]))
        TableDataFound["Geometries"][-1].calculate_md5()   
        TableDataFound["Geometries"][-1].geometryclass_id = geometry_tools.returngeometryclassID(TableDataFound["Geometries"][-1].geometry)
          
        #Check for is_minimun in other geometry in geometry class TODO
        TableDataFound["ElectronicStates"][-1].is_minimum = geometry_tools.electronicstateisminimumbygeometry(None, TableDataFound["Geometries"][-1].geometryclass_id, TableDataFound["MolecularSpecies"][-1].charge, TableDataFound["TheoryLevels"][-1])
              
        elements = [] #tuple of Element and Basis
        for atom in iosection.molecule_info.moleculeobject.atoms:
            elements.append(myElement(models.Elements(atomic_number = atom.atomicnum,
                                     atomic_mass = round(atom.exactmass),
                                     standard_atomic_weight = atom.exactmass), 
                             models.BasisSets(name = iosection.molecule_info.atom_basis_set[atom.idx],
                                              description = commons.nwchem_task_basis_sets[iosection.molecule_info.atom_basis_set[atom.idx].upper()])))
        TableDataFound["Elements"].append(elements)
        
        TableDataFound["MolecularSpecies"][-1].formula = iosection.molecule_info.moleculeobject.formula 
        TableDataFound["MolecularSpecies"][-1].inchi = iosection.molecule_info.inchi
        TableDataFound["MolecularSpecies"][-1].inchikey = iosection.molecule_info.inchikey.strip()
        TableDataFound["MolecularSpecies"][-1].aromatic_cycles = 0
        TableDataFound["MolecularSpecies"][-1].qual_index = 0
        
    
        if section_Rot_Const_And_Energy:
            TableDataFound["RotationalConstants"].append(models.RotationalConstants(a = section_Rot_Const_And_Energy.section_object.values["A"],
                                                                               b = section_Rot_Const_And_Energy.section_object.values["B"],
                                                                               c = section_Rot_Const_And_Energy.section_object.values["C"]))
        else:
            TableDataFound["RotationalConstants"].append(None)
        
        if section_NormalEigenvalue_DipoleMoment and section_NormalEigenvalue_InfraRedIntensities:
            #get frequency information
            
            TableDataFound["VibrationalAnalysesHarmonic"].append(models.VibrationalAnalysesHarmonic())
            TableDataFound["TabulatedVibrations"].append([])
            for row_index in range(len(section_NormalEigenvalue_DipoleMoment.section_object.values["data"])):
                dm =  section_NormalEigenvalue_DipoleMoment.section_object.values["data"][row_index]
                ir = section_NormalEigenvalue_InfraRedIntensities.section_object.values["data"][row_index]
                if dm[1]==0:
                    continue
                TableDataFound["TabulatedVibrations"][-1].append(models.TabulatedVibrations(frequency = dm[1],                                                                    
                                                                                        ir_intensity = ir[5],
                                                                                        diff_mu_x = dm[3],
                                                                                        diff_mu_y = dm[4],
                                                                                        diff_mu_z = dm[5],
                                                                                        eigenvectors = str(section_NormalModeEigenvectors.section_object.values["data"][row_index])))
        else:
            TableDataFound["VibrationalAnalysesHarmonic"].append(None)
            TableDataFound["TabulatedVibrations"].append(None)
        TableDataFound["molecule_info"].append(iosection.molecule_info)
    pt.p("end")
    return TableDataFound

def return_parsed_object_for_session(outputfile, uploaded_output_file, uploaded_input_file, uploaded_database_file, session_path, is_dbfile_humanreadable):
    pt.p("start")
    result = {}
    result["calculation_object"] = build_calculation_object(outputfile,
                                                            uploaded_output_file, 
                                                                      uploaded_input_file, 
                                                                      uploaded_database_file, 
                                                                      parsing(uploaded_output_file, True, "NwChemOutput"), 
                                                                      session_path,
                                                                      is_dbfile_humanreadable)
    
    #last_molecular_species = None
    result["calculation_found"] = False
    result["calculation_found_tasks_list"] = None
    tasks_description = []
    relative_file_path = []
    pt.p("start search calculation")
    if "Calculations" in result["calculation_object"]:
        result["calculation_found"], result["first_calculation_found"], result["calculation_found_tasks_list"] = modeltools.count_calculation_by_md5(result["calculation_object"]["Calculations"].output_md5)
    re_media = re.compile(r".*media")
    pt.p("end search calculation")
    if  "Tasks" in result["calculation_object"]:
        if len(result["calculation_object"]["Tasks"]) > 0:
            for i in range(len(result["calculation_object"]["Tasks"])):
                    
                relative_file_path = re_media.sub("", result["calculation_object"]["molecule_info"][i].geometries["sdf"]["filename"])
                #Delete molecule object because it'snt serializable in session
                result["calculation_object"]["molecule_info"][i].moleculeobject = None
                #Search molecular species in DB
                last_ms = modeltools.return_molecule_specie_by_inchikey(result["calculation_object"]["MolecularSpecies"][i].inchikey)
                if last_ms:
                    result["calculation_object"]["MolecularSpecies"][i] = last_ms
                    
                tasks_description.append(TaskCode(result["calculation_object"]["TheoryLevels"][i], 
                                                                 result["calculation_object"]["Tasks"][i],
                                                                 result["calculation_object"]["molecule_info"][i],
                                                                 result["calculation_object"]["ElectronicStates"][i],
                                                                 result["calculation_object"]["RotationalConstants"][i],
                                                                 result["calculation_object"]["VibrationalAnalysesHarmonic"][i],
                                                                 result["calculation_object"]["DipoleMoments"][i],
                                                                 result["calculation_object"]["MolecularSpecies"][i],
                                                                 result["calculation_object"]["Geometries"][i],
                                                                 relative_file_path,
                                                                 result["calculation_object"]["Elements"][i]))
            result["tasks"] = tasks_description
    pt.p("end")
    return result

