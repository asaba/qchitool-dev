'''
Created on 26/ago/2011

@author: asaba
'''
import re
import string
import pybel
import tempfile

from ...tools import qchitool_tools
from ..commonsobjects import Basis, EnergyForDerivativeDipole, MoleculeInfo, InputOuputSection
from ... import performancetest as pt

#def return_file_row(filename, line_number):
#    #This function return the line 'line_number' in the 'filename' file
#    f = qchitool_tools.openfile(filename)
#    Current_line = None
#    Current_Line_Number = 0
#    while not Current_line is None:
#        Current_line = f.readline
#        if Current_Line_Number == line_number:
#            break
#        Current_Line_Number += 1
#    return Current_line

def return_file_rows(filename, initial_line, final_line):
    #This function return the list of rows between initial_line and final_line in "filename" file
    file_rows = []
    f = qchitool_tools.openfile(filename)
    Current_Line_Number = 0
    if final_line is None:
        Current_line = None
        while not Current_line is None:
            Current_line = f.readline
            if Current_Line_Number == initial_line:
                break
            Current_Line_Number += 1
        file_rows.append(Current_line)
    else:
        while 1: #TO CHECK
            line=f.readline()
            if not line: break
            if Current_Line_Number >= initial_line and Current_Line_Number <= final_line:
                file_rows.append(line)
            Current_Line_Number += 1
    f.close()
    return file_rows


def openbabel_count_hetero_rings(mol):
    count = 0
    for ring in mol.GetSSSR():
        if (ring.IsAromatic() and
               any(mol.GetAtom(atom_id).GetAtomicNum() != 6 for atom_id in ring._path)):
            count += 1
        
    return count


def chunks(l, n):
    #split the list l in chunks of n elements
    return [l[i:i+n] for i in range(0, len(l), n)]


def separate_by_schema(splited_line, schema, description_field, value):
    #This function split the Current line by the schema and save it on value
    #the description_field (0 base) can be composed by more than one token
    index = 0
    if not description_field is None:
        
        if len(splited_line)> len(schema):
            description_field_tokens = len(splited_line) - len(schema) + 1
            description = ""
            for i in splited_line[description_field:description_field + description_field_tokens]:
                description = description + str(i) + " "
            description = description[:-1]
            splited_line = splited_line[: description_field] + [description] + splited_line[description_field + description_field_tokens:]
        else:
            if len(splited_line)< len(schema):
                #the description field is not present
                schema = schema[:description_field]+schema[description_field+1:]
            
    for i in schema:
        if i == "s": #sting
            splited_line[index] = splited_line[index].strip()
        if i == "n": #number
            splited_line[index] = qchitool_tools.val(splited_line[index])
        index += 1
    value.append(tuple(splited_line))

def get_parenteses_value(row, parenteses):
    if len(parenteses) <> 2:
        print "get_parenteses_value ERROR"
        return None
    else:
        if re.search(parenteses[0], row):
            return re.split(parenteses[1], re.split(parenteses[0], row, 1)[-1], 1)[0].strip()
        else:
            print "get_parenteses_value ERROR"
            return None

def get_singlevalue_sep(row, tuple_separator, value_separator, only_one_value, mydictionary):
        #This function split the Current line by the tuple_separator and (after) by values_separator 
    #It add each couple of "variable : value" in mydictionary
    #Ex: something = 23.4; somethingelse = 11 (tuple_separator = ";", value_separator = "="
        Current_line = row
        if Current_line.find(value_separator)!=-1:
            if tuple_separator == "" or tuple_separator is None:
                if only_one_value:
                    tmp = Current_line.split(value_separator, 1)
                    if len(tmp)>1:
                        mydictionary[tmp[0].strip()]=tmp[1].strip()
                else:
                    tmp_1 = Current_line.strip().split(value_separator)
                    tmp = []
                    for s in tmp_1:
                        tmp += s.split()
                    for index in range(0, len(tmp), 2):
                        mydictionary[tmp[index]] = tmp[index+1]
            else:
                tmp1 = Current_line.split(tuple_separator)
                for t in tmp1:
                    tmp = t.split(value_separator, 1)
                    if len(tmp)>1:
                        mydictionary[tmp[0].strip()]=tmp[1].strip()

def get_listvalue_sep(filename, initial_line, end_line, tuple_separator, value_separator, only_one_value, section_rows = None):

    values = {}
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    for row in section_rows[:end_line-initial_line+1]:
        get_singlevalue_sep(row, tuple_separator, value_separator, only_one_value, values)
    return values    
        
def get_listvalue_no_sep(filename, initial_line, end_line, section_rows = None):
        #two field only. The last is a number or a number followed by "local" or "non-local"
        #Ex: No. of variables    17199 
        #Ex: Convergence  1.0E-04 local
        values = {}
        if section_rows is None:
            section_rows = return_file_rows(filename, initial_line, end_line)
        for row in section_rows[:end_line-initial_line+1]:
            tmp = row.split()
            if len(tmp)>0:
                if tmp[-1].strip() == "local" or tmp[-1].strip() == "non-local":
                    values[string.join(tmp[:-2] + ['(', tmp[-1], ')']).strip()]=qchitool_tools.val(tmp[-2].strip())
                else:
                    if string.join(tmp[:-1]).strip() <> "" and qchitool_tools.val(tmp[-1].strip()) <> None:
                        values[string.join(tmp[:-1]).strip()]=qchitool_tools.val(tmp[-1].strip()) 
        return values  
    
def get_value_from_to_no_sep(filename, initial_line, end_line, separator, section_rows = None):
    #Return values that have the schema "variable from x to y" where the couple of value are 
    #separated by separator (ex. "to") only if separator is present
    values = []
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    for row in section_rows[:end_line-initial_line+1]:
        tmp = row.split()
        for i in tmp:
            if i.strip()==separator:
                variable = string.join(tmp[:-3]).strip()
                from_to_value = string.join(tmp[-3:]).split(separator)
                values.append([variable, qchitool_tools.val(from_to_value[0]), qchitool_tools.val(from_to_value[-1])]) #[variable_name, from_value, to_value]
                break;
    return values      

def get_multirow_table(filename, initial_line, end_line, schema, skip_first_row, section_rows = None):
    #return value stored in multiline with vertical schema 'schema' 
    values = []
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    tmp_line = []
    index = 0
    skip_first_row_tmp = skip_first_row
    for row in section_rows[:end_line-initial_line+1]:
        if len(row.strip()) > 0:
            if skip_first_row_tmp:
                skip_first_row_tmp = False
                continue
            tmp = row.split()
            if len(tmp_line)<(index+1):
                tmp_line.append([])
            for value in tmp[1:]: #skip row name
                try:
                    if schema[index]=="s":
                        tmp_line[index].append(value.strip()) 
                    else: #n
                        tmp_line[index].append(qchitool_tools.val(value))
                except:
                    continue
            index+=1
        else:
            index = 0
            skip_first_row_tmp = skip_first_row
            
    for i in range(len(tmp_line[0])):
        tmp_value = []
        for j in range(len(tmp_line)):
            tmp_value.append(tmp_line[j][i])
        values.append(tuple(tmp_value))
    return values

def get_continous_list_no_sep(filename, initial_line, end_line, schema, section_rows = None):

        tmp_values = []
        values = []
        if section_rows is None:
            section_rows = return_file_rows(filename, initial_line, end_line)
        for row in section_rows[: end_line - initial_line + 1]:
            tmp = row.strip().split()
            if len(tmp)>0:
                for v in tmp:
                    if not(v is None):
                        tmp_values.append(v)
                        
        for value_index in range(0, len(tmp_values), len(schema)):
            tmp_output = []
            schema_index = 0
            for j in schema:
                if j=="s":
                    value_found = tmp_values[value_index + schema_index].strip()
                else: # "n"
                    value_found = qchitool_tools.val(tmp_values[value_index + schema_index])
                tmp_output.append(value_found)
                schema_index +=1
            if len(tmp_output) > 1:
                    tmp_output = [l for l in tmp_output if l is not None]
                    values.append(tuple(tmp_output))
            else:
                    values.append(tmp_output[0])
        values = [l for l in values if l is not None]
        return values
    
    
def get_tabulated_data(filename, initial_line, end_line, schema, first_row_value, description_field, value_separator, section_rows = None):
    #Tabulated data are a list of row where the first start with first_row_value, each row follows the schema separated by value_separator (space if value_separator="" or None)
    #the description_field-th field in the schema can be composed by more than one token.
    values = []
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    if first_row_value is None:
        first_row_found = True
    else:
        first_row_found = False
    for row in section_rows[:end_line-initial_line + 1]:
        if value_separator is None or value_separator == "":
            tmp = row.split()
        else:
            tmp = row.split(value_separator)
        if len(tmp)>0:
            if not first_row_found:
                if qchitool_tools.val(tmp[0]) != first_row_value: #find the first row of table
                    continue
                else:
                    first_row_found = True
            if len(tmp)>=len(schema):
                separate_by_schema(tmp, schema, description_field, values)
    return values   
    

    
#def get_basisvalue(filename, initial_line, end_line, schema, description_field, section_rows = None):
#    #Basis values are a list of basis_class where for each basis_class there is a element and a set of gaussians for each orbital 
#    values = []
#    element_pattern = re.compile(r".+ \(.+\)")
#    if section_rows is None:
#        section_rows = return_file_rows(filename, initial_line, end_line)
#    for row in section_rows[:end_line-initial_line+1]:
#        if  element_pattern.search(row): #start of orbitals information for one element
#            b = Basis(row.split()[1].strip()) #element name
#            values.append(b)
#        else:
#            if not len(row.strip)==0: #separator line
#                tmp= row.split()
#                separate_by_schema(tmp, schema, description_field, b.values["Orbitals"])
#    return values


def get_z_matrix(filename, initial_line, end_line, section_rows = None):
    
    values = []
    extra_column = ""
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    prog = re.compile(r".*Type.*")
    for row in section_rows[:end_line-initial_line+1]:
        #Check the number of columns
        if prog.match(row):
            row_len = len(row.strip().split())
            if row_len > 8:
                extra_column = "n" * (row_len - 8)
            continue
        row = row.strip().split()
        if len(row)>4: #separator line
            if row[1].strip() == "Stretch":
                separate_by_schema(row, "nssnnn" + extra_column, 2, values)
            else: 
                if row[1].strip() == "Bend":
                    separate_by_schema(row, "nssnnnn" + extra_column, 2, values)
                else:
                    if row[1].strip() == "Torsion":
                        separate_by_schema(row, "nssnnnnn" + extra_column, 2, values) 
    return values

def get_hessian_derivative_value(filename, initial_line, end_line, section_rows = None):
    #X vector of derivative dipole (au) [debye/angstrom]
    #d_dipole_x/<atom=   1,x> =     0.1818     [    0.8733]
    #d_dipole_x/<atom=   1,y> =    -0.5943     [   -2.8545]
    #d_dipole_x/<atom=   1,z> =     0.0000     [    0.0000]
    
    values = []
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    prog_vector = re.compile(r".*vector of derivative dipole.*")
    prog_dipole = re.compile(r"d_dipole_.*")
    for row in section_rows[:end_line-initial_line+1]:
        if prog_vector.match(row.strip()):
            unit = row.split("[")[1].split("]")[0]
        elif prog_dipole.match(row.strip()):
            splited_row = row.split("=")
            vector = splited_row[0].split("_")[2].split("/")[0]
            atom = qchitool_tools.val(splited_row[1].split(",")[0])
            direction = splited_row[1].split(",")[1].split(">")[0]
            value = qchitool_tools.val(splited_row[2].split("[")[1].split("]")[0])
            values.append(EnergyForDerivativeDipole(vector, atom, direction, value, unit))
    return values

def get_triangular_matrix(filename, initial_line, end_line, section_rows = None):
    matrix = {}
    series = False
    if section_rows is None:
        section_rows = return_file_rows(filename, initial_line, end_line)
    prog = re.compile(r"[-]+")
    for row in section_rows[:end_line-initial_line+1]:
        if row.strip() == "":
            #stop series
            series = False
            continue
        index = row.split()[0]
        if index == "":
            #stop series
            series = False
            continue
        if prog.match(index):
            #start new series
            series = True
            continue
        if series:
            if qchitool_tools.val(index) in matrix.keys():
                matrix[qchitool_tools.val(index)] += row.split()[1:]
            else:
                matrix[qchitool_tools.val(index)] = row.split()[1:]
    return matrix


#def return_inputoutput_sectionschema_OLD(parsing_result):
#    inputsection = None
#    outputsection = None
#    inp_out_sections = []
#    for supersection in parsing_result:
#        if re.search("NwchemInputModule", supersection.section_class):
#            if (not inputsection is None) and (not outputsection is None):
#                inp_out_sections.append(InputOuputSection(inputsection, outputsection))
#                outputsection = None
#            inputsection = supersection
#        else:
#            outputsection = supersection
#    return inp_out_sections

def return_inputoutput_sectionschema(parsing_result):
    inputmoduleindex = []
    inp_out_sections = []
    for index in range(len(parsing_result)):
        supersection = parsing_result[index]
        if re.search("NwchemInputModule", supersection.section_class):
            inputmoduleindex.append(index)
    i = 0
    for index in range(len(inputmoduleindex[:-1])):
        inp_out_sections.append(InputOuputSection(parsing_result[inputmoduleindex[index]], parsing_result[inputmoduleindex[index]+1 : inputmoduleindex[index+1]]))
        i += 1
    return inp_out_sections

def returnsection(section_name, parsing_result, number):
    numberupper = number.upper()
    if numberupper == "ALL": #in sections_result we will find the three level separated 
        sections_result = [[],[],[]]
    elif numberupper == "LIST": #a raw list of sections found in order
        sections_result = [] 
    else: #FIRST or LAST. In the second case sections_result will be override. It is the standard case when the number it'snt different from ALL. LIST or FIRST
        sections_result = None
    exitfunction = False
    if "list_of_subsections" in dir(parsing_result):
        list_sections = parsing_result.list_of_subsections
    else:
        list_sections = parsing_result
    for supersection in list_sections:
        if exitfunction: break
        if re.search(section_name, supersection.section_class):
            if numberupper == "ALL":
                sections_result[0].append(supersection)
            elif numberupper == "LIST":
                sections_result.append(supersection)
            else:
                sections_result = supersection
                if numberupper == "FIRST":
                    return sections_result
                    exitfunction = True
        for subsection in supersection.list_of_subsections:
            if exitfunction: break
            if re.search(section_name, subsection.section_class):
                if numberupper == "ALL":
                    sections_result[1].append(subsection)
                elif numberupper == "LIST":
                    sections_result.append(subsection)
                else:
                    sections_result = subsection
                    if numberupper == "FIRST":
                        return sections_result
                        exitfunction = True
            for subsubsection in subsection.list_of_subsections:
                if exitfunction: break
                if re.search(section_name, subsubsection.section_class):
                    if numberupper == "ALL":
                        sections_result[2].append(subsubsection)
                    elif numberupper == "LIST":
                        sections_result.append(subsubsection)
                    else:
                        sections_result = subsubsection
                        if numberupper == "FIRST":
                            return sections_result
                            exitfunction = True
    if numberupper == "ALL" and sections_result == [[],[],[]]:
        return None
    elif numberupper == "LIST" and sections_result == []:
        return None
    else:
        return sections_result
    
def writefile(lines, restart_name, session_path, fileformat):
    if session_path:
        session_path_tmp = session_path
    else:
        session_path_tmp = tempfile.gettempdir() + "/"
        
    new_file_name = session_path_tmp +   "_" + fileformat + "_" +  restart_name +  fileformat
    new_file = open(new_file_name, "w")
    new_file.writelines(lines)
    new_file.close()
    return new_file_name

def extract_molecule_info(out_file, charge, session_path, 
                          basis_set_list, geometry_list, spin_multiplicity, unit_conversion_value):
    #This function extract structural information in nwchem output file
    #Moreover it makes inchi and inchikey 
    #basis_set_list = [[element tag , basis set standard name], ...]
    #geometry_list = [[index, atom symbol, atom atomic number, X coordinate, Y coordinate, Z coordinate], ...] 
    if session_path:
        session_path_tmp = session_path
    else:
        session_path_tmp = tempfile.gettempdir() + "/"
        
    basis_set_dictionary = {} # {element tag : basis set standard name}
    if basis_set_list:
        index = 0
        for basis_set_row in basis_set_list:
            basis_set_dictionary[basis_set_row[0]] = basis_set_row[1]
            index += 1
    molecule = pybel.ob.OBMol()      
    atom_basis_set = {}

    for atom_by_geom in geometry_list:
        atom = molecule.NewAtom()
        atom.SetVector(float(atom_by_geom[3])/unit_conversion_value, float(atom_by_geom[4])/unit_conversion_value, float(atom_by_geom[5])/unit_conversion_value)
        atom.SetAtomicNum(int(atom_by_geom[2]))
        if basis_set_list:
            if "*" in basis_set_dictionary: #same basis set for each atom
                atom_basis_set[atom.GetIdx()] = basis_set_dictionary["*"]
            else:
                atom_basis_set[atom.GetIdx()] = basis_set_dictionary[atom_by_geom[1]]
        #CHEKC
        #pybelmol.atoms[atom_by_geom].exactmass = mass_dataentry.section_object.values["values"][atom_by_geom]
    molecule.ConnectTheDots()
    xyz = pybel.Molecule(molecule).write("xyz")
    cml = pybel.Molecule(molecule).write("cml")
    molecule.PerceiveBondOrders()
    sdf = pybel.Molecule(molecule).write("sdf")
    uncharged_molecule = pybel.Molecule(molecule)
    
    #Remove the name information in CML file
    for row in cml:
        #remove molecule name
        row = re.sub(".*[<]name[>].*[<][/]name[>]\n", "", row)
        #remove molecule id
        if re.match(".*[<]molecule.*", row):
            row = re.sub(".*[<]molecule.+[>]", "<molecule>", row)
    xyz_file_name = writefile(xyz, "out_file.", session_path_tmp, "xyz")
    cml_file_name = writefile(cml, "out_file.", session_path_tmp, "cml")
    sdf_file_name = writefile(sdf, "out_file.", session_path_tmp, "sdf")
    
    #get first molecule information found
    conv = pybel.ob.OBConversion()
    conv.SetInAndOutFormats("nwo", "inchi") 
    conv.SetOptions("K", conv.OUTOPTIONS)
    mol = pybel.ob.OBMol()
    conv.ReadFile(mol, str(out_file)) 
    
    AtomNum = mol.NumAtoms()
    
    if AtomNum:
        atom = mol.GetAtom(1)
        atom.SetFormalCharge(charge)
        inchikey = conv.WriteString(mol) # we need a conversion because the stantard output isn't inchikey
        #aromatic_cycles = openbabel_count_hetero_rings(mol)
                                    
        first_mol = pybel.readfile("nwo", str(out_file)).next()
        
        first_mol.atoms[0].OBAtom.SetFormalCharge(charge)
        inchi = first_mol.write("inchi").split("=")[1]
        #sdf = re.split("[\r\n]+", first_mol.write("sdf"), 2)[2]
        #save geometry without charge information
    else:
        first_mol = pybel.Molecule(molecule)
        #first_mol.spin = None
        inchikey = None
        inchi = None
        
    return MoleculeInfo(inchi, inchikey, uncharged_molecule, spin_multiplicity, charge, geometry_list, xyz_file_name, xyz, sdf_file_name, sdf, cml_file_name, cml, atom_basis_set)
    