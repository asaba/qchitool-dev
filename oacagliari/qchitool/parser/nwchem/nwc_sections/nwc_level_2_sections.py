'''
Created on 07/mar/2011

@author: andrea saba
@email: asaba@oa-cagliari.inaf.it

classes define subsections  
'''

from ... import commonsobjects
from ...tools import parser_tools
from ....tools import qchitool_tools
from .... import performancetest as pt
from nwc_sec_commons import GenericSection
import re
import qchitool
#from oacagliari.qchitool.tools import qchitool_tools

class JobInformation(GenericSection):
    #In this class there are information about the executed job
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "job_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows).copy() #shallow copy
        
class GeneralInformation(GenericSection):
    #In this class there are general information about the calculus
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "general_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, ";", ":", False, self.section_rows).copy() #shallow copy
        
class XcInformation(GenericSection):
    #In this class there are general information about the calculus
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "xc_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.method = None
        re_method = re.compile(r"Method")
        for row in self.section_rows:
            if re_method.search(row):
                self.method = row.split("Method")[0].strip()
                break
        self.values_list = parser_tools.get_listvalue_no_sep(self.section_filename, self.section_initial_line+3, self.section_final_information_line, self.section_rows)
        self.values = {"method": self.method,
                       "values_list": self.values_list}

        
class GridInformation(GenericSection):
    #In this class there are Grid information about the calculus
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "grid_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        #self.values = tf.(self.section_filename, self.section_initial_line, self.section_final_line, ";", ":").copy() #shallow copy
        sections = GenericSection.split_section(self, [r"Grid used.*", 
                                                           r"Tag[\s]*B.*", 
                                                           r"Grid pruning.*"])
        # sections is a list of (section_name, striped_file_row, [initial_line, final_line])
        self.values_list = parser_tools.get_listvalue_sep(self.section_filename,  sections[0].section_initial_line, sections[0].section_final_line, None, ":", True, self.section_rows[sections[0].section_initial_line - self.section_initial_line : sections[0].section_final_line - self.section_initial_line])
        self.atom_value = parser_tools.get_tabulated_data(self.section_filename, sections[1].section_initial_line, sections[1].section_final_line, "snnnn" ,None, 0, None, self.section_rows[sections[1].section_initial_line - self.section_initial_line +2 : sections[1].section_final_line - self.section_initial_line])
        self.values_list.update(parser_tools.get_listvalue_sep(self.section_filename,  sections[2].section_initial_line, sections[2].section_final_line,  None, ":", True, self.section_rows[sections[2].section_initial_line - self.section_initial_line : sections[2].section_final_line - self.section_initial_line]))
        self.values = {"atom_value": self.atom_value,
                       "values_list": self.values_list}

        
class ConvergenceInformation(GenericSection):
    #In this class there are Convergence information about the calculus
    #To Check for table info
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "convergence_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        sections = GenericSection.split_section(self, [r"Convergence aids.*", 
                                                          r"Damping.*"])
        self.values_list = parser_tools.get_listvalue_sep(self.section_filename,  sections[0].section_initial_line, sections[0].section_final_line, None, ":", True, self.section_rows[sections[0].section_initial_line - self.section_initial_line : sections[0].section_final_line - self.section_initial_line])
        self.de = parser_tools.get_tabulated_data(self.section_filename, sections[1].section_initial_line, sections[1].section_final_line, "ssss" , None, 0, None, self.section_rows[sections[1].section_initial_line - self.section_initial_line +2 : sections[1].section_final_line - self.section_initial_line])
        self.values = {"de": self.de,
                       "values_list": self.values_list}
        
class ScreeningToleranceInformation(GenericSection):
    #In this class there are screening tolerance information about the calculus
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "screening_tolerance_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, ":", True, self.section_rows).copy() #shallow copy
        
class SuperpositionOfAtomicDensityGuess(GenericSection):
    #In this class there are the sum of atomic energies and the renormalizing density
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "superposition_of_atomic_density_guess_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values_list = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, ":", True, self.section_rows)
        self.values_from_to = parser_tools.get_value_from_to_no_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, "to", self.section_rows)
        self.values = {"values_from_to": self.values_from_to,
                       "values_list": self.values_list}


class NonVariationalInitialEnergy(GenericSection):
    #In this class there are information about the non-variational initial energy
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "non_variational_initial_energy_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows).copy() #shallow copy
        
    
class MemoryInformation(GenericSection):
    #In this class there are information about the memory used by nwchem
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "memory_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows).copy()  #shallow copy
        

class DirectoryInformation(GenericSection):
    #In this class there are information about the directory used by nwchem
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "directory_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows).copy()  #shallow copy
        
        
class FiniteDifferenceHessianDelta(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "FiniteDifferenceHessianDelta_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        sections_re = [r"finite difference hessian delta.*"]
        sections = GenericSection.split_section(self, sections_re)
        self.raw_value = []
        for section in sections:
                if section.section_name == sections_re[0]:
                    self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                                 section.section_initial_line, 
                                                                 section.section_initial_line+1, 
                                                                 None, "=", True, 
                                                                 self.section_rows[section.section_initial_line - self.section_initial_line: section.section_initial_line + 1 - self.section_initial_line])
                    numberoflines = sections[0].section_final_line - sections[0].section_initial_line 
                    for i in range(1, numberoflines):
                        if len(self.section_rows[-i].strip()) > 0:
                            hessian_matrix_dimension = qchitool_tools.val(self.section_rows[-i].split()[0])
                            break
                    self.hessian_matrix = parser_tools.get_multirow_table(self.section_filename, 
                                                                          section.section_initial_line+1, 
                                                                          section.section_final_line, 
                                                                          "n"*hessian_matrix_dimension, 
                                                                          True, 
                                                                          self.section_rows[section.section_initial_line + 1 - self.section_initial_line: section.section_final_line - self.section_initial_line])
                    self.values["hessian_matrix"] = self.hessian_matrix
                
        pt.p("end")

class FiniteDifferenceDerivativeDipole(GenericSection):
    #get information from derivative dipole moment separated by atom's components 
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "FiniteDifferenceDerivativeDipole_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        
        sections_re = [r"finite difference derivative dipole.*", 
                       r"triangle hessian written.*", 
                       r"Vibrational analysis.*"]
        sections = GenericSection.split_section(self, sections_re)
        self.raw_value = []
        for section in sections:
                    if section.section_name == sections_re[0]: #r"finite difference derivative dipole.*", 
                        self.values.update(parser_tools.get_listvalue_sep(self.section_filename, 
                                                                          section.section_initial_line, 
                                                                          section.section_initial_line+1, 
                                                                          None, "=", True, 
                                                                          self.section_rows[section.section_initial_line - self.section_initial_line: section.section_initial_line + 1 - self.section_initial_line]))
                        self.vector_of_derivative_dipole = parser_tools.get_hessian_derivative_value(self.section_filename, 
                                                                                                    section.section_initial_line+1, 
                                                                                                    section.section_final_information_line, 
                                                                                                    self.section_rows[section.section_initial_line + 1 - self.section_initial_line: section.section_final_line - self.section_initial_line])
                        self.values["vector_of_derivative_dipole"] = self.vector_of_derivative_dipole
                    elif section.section_name == sections_re[1] or section.section_name == sections_re[2]:
                        for i in self.section_rows[section.section_initial_line - self.section_initial_line: section.section_final_information_line - self.section_initial_line]: self.raw_value.append(i)
        pt.p("end")


class PreviousTaskInformation(GenericSection):
    #In this class there are information about the restarted job
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "previous_task_information_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows).copy()  #shallow copy
        
class GeometriesInTheDatabase(GenericSection):
    #In this class there are information about geometries found in the database of the restarted job
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "geometries_in_the_database_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "nsnssnsn", 1, 1, None, self.section_rows)
        #try to find the restart geometry information name
        database_geometry_name = "UNDEFINED" 
        re_geom_name = re.compile(r"The geometry named.*")
        for row in self.section_rows:
            if re_geom_name.search(row):
                database_geometry_name = parser_tools.get_parenteses_value(row, [r'"', r'"'])
                break
        #'The geometry named "geometry" is the default for restart'
        self.values = {"data": self.data, "database_geometry_name": database_geometry_name}
    
class BasisSetsInTheDatabase(GenericSection):
    #In this class there are information about basis found in the database of the restarted job
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "basis_sets_in_the_database_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "nsnssnsn", 1, 1, None, self.section_rows)
        database_basis_set_name = "UNDEFINED" 
        re_basissetname = re.compile(r"The basis set named.*")
        for row in self.section_rows:
            if re_basissetname.search(row):
                database_basis_set_name = parser_tools.get_parenteses_value(row, [r'"', r'"'])
                break
        self.values = {"data": self.data, "database_basis_set_name": database_basis_set_name}
    
class MultipoleAnalysisOfTheDensity(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "multipole_analysis_of_the_density_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "n"*8, 0 , None, None, self.section_rows)
        self.level = {}
        for row in self.data:
            if row[0] == 0:
                self.charge = row[4]
                continue
            elif row[0] >= 1:
                level_vector = "x"*row[1] + "y"*row[2] + "z"*row[3]
                self.level[level_vector] = row[4] #qchitool_tools.val(row[4])
        self.values = {"data": self.data,
                       "charge": self.charge,
                       "levels": self.level}

class MultipoleAnalysisOfTheDensityWRTOrigini(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "multipole_analysis_of_the_density_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "n"*7, 0 , None, None, self.section_rows)
        self.level = {}
        for row in self.data:
            if row[0] == 0:
                self.charge = row[4]
                continue
            elif row[0] >= 1:
                level_vector = "x"*row[1] + "y"*row[2] + "z"*row[3]
                self.level[level_vector] = row[4] #qchitool_tools.val(row[4])
        self.values = {"data": self.data,
                       "charge": self.charge,
                       "levels": self.level}

class SummaryOf(GenericSection):
    #In this class there are basis sets used for each element
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        self.data = [] # [Tag(element), Description, Shells, Functions, Types] (Shell, Functions and Types can be ['on', 'all', 'atoms'] 
        GenericSection.print_name(self, "summary_of_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        
        if section_rows is None:
            section_rows = parser_tools.return_file_rows(filename, initial_line, end_line)
        #" Summary of "ao basis" -> "" (cartesian)"
        #get "from" and "to" information
        self.summary_from = re.findall("(?<=of \").*(?=\" ->)", GenericSection.go_to_first_line(self, self.section_filename)[1])[0]
        self.summary_to = re.findall("(?<=-> \").*(?=\" \()", GenericSection.go_to_first_line(self, self.section_filename)[1])[0]

        prog = re.compile(r"[-]+")
        for row in self.section_rows[1 : self.section_final_information_line -self.section_initial_line + 1]: #skip first row
            tmp_line = row.split()
            if len(tmp_line) == 0:
                continue
            if tmp_line[0] == "Tag":
                continue
            elif prog.match(tmp_line[0]):
                continue
            if tmp_line[0] == "Caching":
                break
            parser_tools.separate_by_schema(tmp_line, "sssss", 0, self.data)
        self.values = {"summary_from": self.summary_from,
                       "summary_to": self.summary_to,
                       "data": self.data}

        
class CenterOfMass(GenericSection, commonsobjects.CenterOfMass):
    #In this class there is the molecular center of mass coordinates 
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "center_of_mass_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        tmp = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", False, self.section_rows)
        commonsobjects.CenterOfMass.__init__(self, qchitool_tools.val(tmp['x']), qchitool_tools.val(tmp['y']), qchitool_tools.val(tmp['z']))

        
class MomentOfInertia(GenericSection, commonsobjects.MomentOfInertia):
    #In this class there is the molecular moment of inertia 
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "moment_of_inertia_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        Result = GenericSection.go_to_first_line(self, self.section_filename)
        
        Current_line = Result[1]
        #get energy units information
        self.units_of_measurement = parser_tools.get_parenteses_value(Current_line, [r"\(", r"\)"])
        if section_rows is None:
            tmp = parser_tools.get_continous_list_no_sep(self.section_filename, 3 , self.section_final_information_line - self.section_initial_line, "n", self.section_rows)
        else:
            tmp = parser_tools.get_continous_list_no_sep(self.section_filename, 0, self.section_final_information_line - self.section_initial_line, "n", self.section_rows[2:])
        mtx = parser_tools.chunks(tmp, 3) 
        commonsobjects.MomentOfInertia.__init__(self, mtx)

        
class ExpectationValueOfS2(GenericSection):
    #In this class there is the molecular value of S2 (expected and exact)
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        self.values = ({})
        GenericSection.print_name(self, "expectation_value_of_S2_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, "(", "=", False, self.section_rows)
        for k in values.keys():
            self.values[k] = qchitool_tools.val(values[k].replace(")", "")) #remove ')'

class SymmetryAnalysisOfMolecularOrbitals(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "symmetry_analysis_of_molecular_orbitals")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.spin_orbital = GenericSection.go_to_first_line(self, None)[1].split()[-1]
        if self.spin_orbital <> "alpha" and self.spin_orbital <> "beta":
            self.spin_orbital = "both"
        sections = GenericSection.split_section(self, [r"Numbering of irreducible representations:", 
                                                           r"Orbital symmetries:", 
                                                           r"Time (after|prior).*"])
        # sections is a list of (section_name, striped_file_row, [initial_line, final_line])
        self.irreducible_values = parser_tools.get_continous_list_no_sep(self.section_filename, 1, sections[0].section_final_line, "ns", self.section_rows[sections[0].section_initial_line - self.section_initial_line +2 : sections[0].section_final_line - self.section_initial_line])
        self.symmetry_values = parser_tools.get_continous_list_no_sep(self.section_filename, 1, sections[1].section_final_line, "ns", self.section_rows[sections[1].section_initial_line - self.section_initial_line +2 : sections[1].section_final_line - self.section_initial_line])
        self.values = {"spin_orbital": self.spin_orbital, 
                       "irreducible_values": self.irreducible_values, 
                       "symmetry_values": self.symmetry_values}
        #add other data
        
class SwappingOrbitals(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "swapping_orbitals")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        first_row = GenericSection.go_to_first_line(self, None)[1].split()
        self.spin_orbital = first_row[1]
        self.first_orbital = first_row[-2]
        self.second_orbital = first_row[-1]
        if self.spin_orbital <> "alpha" and self.spin_orbital <> "beta":
            self.spin_orbital = "both"
        self.values = {"spin_orbital": self.spin_orbital, 
               "first_orbital": self.first_orbital, 
               "second_orbital": self.second_orbital}
        
class SymmetryInformation(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "symmetry_information")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        sections = GenericSection.split_section(self, [r"Group name.*",
                                                          r"Symmetry unique atoms", 
                                                           r"internuclear distances", 
                                                           r"internuclear angles"]) 
        self.values = parser_tools.get_listvalue_no_sep(self.section_filename, 
                                              sections[0].section_initial_line, 
                                              sections[0].section_final_line, 
                                              self.section_rows[sections[0].section_initial_line - self.section_initial_line-1: sections[0].section_final_line - self.section_initial_line])
        symmetryuniqueatoms = parser_tools.get_continous_list_no_sep(self.section_filename, 
                                                           1, 
                                                           sections[1].section_final_line, 
                                                           "n", 
                                                           self.section_rows[sections[1].section_initial_line - self.section_initial_line: sections[1].section_final_line - self.section_initial_line])
        self.values["symmetryuniqueatoms"] = symmetryuniqueatoms
        if len(sections) >2:
            self.internucleardistances = parser_tools.get_tabulated_data(self.section_filename, 
                                                                sections[2].section_initial_line, 
                                                                sections[2].section_final_line, 
                                                                "ssnn", None, 0, "|", 
                                                                self.section_rows[sections[2].section_initial_line - self.section_initial_line + 3: sections[2].section_final_line - self.section_initial_line])
            self.internuclearangles = parser_tools.get_tabulated_data(self.section_filename, 
                                                             sections[3].section_initial_line, 
                                                             sections[3].section_final_line, 
                                                             "sssn", None, 0, "|", 
                                                             self.section_rows[sections[3].section_initial_line - self.section_initial_line + 3: sections[3].section_final_line - self.section_initial_line])
            self.values["internucleardistances"] = self.internucleardistances
            self.values["internuclearangles"] = self.internuclearangles
            
        
class DftEnergyGradients(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "dft_energy_gradients_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "nsnnnnnn", 1, 1, None, self.section_rows)
        self.values = {"data": self.data}
        
class Geometry(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "geometry_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        Result = GenericSection.go_to_first_line(self, self.section_filename)
        Current_Line_Number = Result[0]
        Current_line = Result[1]
        #get "from" and "to" information
        tmp = Current_line.split('"')
        self.geometry_from=tmp[1]
        self.geometry_to=tmp[3]
        #get coordinates units information
        self.units_of_measurement = None
        self.unit_conversion_value = None
        prog = re.compile(r".*Output coordinates in .*")
        for i in range(self.section_initial_line, self.section_initial_line+4):
            tmp = GenericSection.return_row(self, i) #skip to Output information line
            if prog.match(tmp):
                #print parser_tools.get_parenteses_value(tmp, [r"\(", r"\)"])
                self.units_of_measurement = tmp.split()[3].strip()
                self.unit_conversion_value = qchitool_tools.val(parser_tools.get_parenteses_value(tmp, [r"\(", r"\)"]).split()[2].strip())

        self.values = {"data": parser_tools.get_tabulated_data(self.section_filename, 
                                                      Current_Line_Number, 
                                                      self.section_final_information_line, 
                                                      "nsnnnn", 1, 1, None, self.section_rows),# [index, Tag, AtomicMass, X, Y, Z]
                       "geometry_from": self.geometry_from,
                       "geometry_to": self.geometry_to,
                       "units_of_measurement": self.units_of_measurement, 
                       "unit_conversion_value": self.unit_conversion_value}
                        

class EffectiveNuclearRepulsionEnergy(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "effective_nuclear_repulsion_energy_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        Result = GenericSection.go_to_first_line(self, self.section_filename)
        
        Current_line = Result[1]
        #get energy value
        self.values = {"effective_nuclear_repulsion_energy_class" : qchitool_tools.val(Current_line.split()[-1])}
        #get energy units information
        self.units_of_measurement = parser_tools.get_parenteses_value(Current_line, [r"\(", r"\)"])
        self.values["units_of_measurement"] = self.units_of_measurement

        
class NuclearDipoleMoment(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "nuclear_dipole_moment_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        Result = GenericSection.go_to_first_line(self, self.section_filename)
        
        Current_line = Result[1]
        #get energy units information
        self.units_of_measurement = parser_tools.get_parenteses_value(Current_line, [r"\(", r"\)"])
        
        self.data = parser_tools.get_continous_list_no_sep(self.section_filename, 0, self.section_final_line - self.section_initial_line, "n", self.section_rows)
        self.values = {"data": self.data,
                       "units_of_measurement": self.units_of_measurement}

        
class AtomicMass(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "atomic_mass_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        
        self.values_list= parser_tools.get_listvalue_no_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, self.section_rows)
        self.values = {"values_list":self.values_list}

        
class FinalUhfResult(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "final_uhf_result_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values= parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows)
        
#class energy_minimization_class(tf.generic_section):
#    values =

class ZMatrix(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "z_matrix_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        #get information about units used

        #if filename is None:
        #    tmp = self.section_rows[3]
        #else:
        #    tmp = parser_tools.return_file_row(self.section_filename, self.section_initial_line + 3) #skip to Output information line
        #" Units are Angstrom for bonds and degrees for angles"
        self.units_of_measurement = None
        self.units_for_angles = None
        re_unit = re.compile(r"Unit are .*")
        for row in self.section_rows:
            
            if re_unit.match(row):
                self.units_of_measurement = re.findall("(?<=Unit are).*(?=for bonds)",  row)[0].strip()
                self.units_for_angles = re.findall("(?<=and).*(?=for angles)",  row)[0].strip()
                break

        sections = GenericSection.split_section(self, [r".*Type.*", 
                                                     r".*internuclear distances.*", 
                                                     r".*internuclear angles.*"])
        # sections is a list of (section_name, striped_file_row, [initial_line, final_line])
        #print sections[0].section_rows

        self.value_Z = parser_tools.get_z_matrix(self.section_filename, 
                                       sections[0].section_initial_line, 
                                       sections[0].section_final_line, 
                                       self.section_rows[sections[0].section_initial_line - self.section_initial_line-1: sections[0].section_final_line - self.section_initial_line])
        #don't save internuclear distances and internuclear angles
        self.values = {"units_of_measurement": self.units_of_measurement,
                       "units_for_angles": self.units_for_angles,
                       "value_Z": self.value_Z}  
        
class BasisFrom(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        self.basis = [] #list of basis object
        GenericSection.print_name(self, "basis_from_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        Result = GenericSection.go_to_first_line(self, self.section_filename)
        Current_line = Result[1]
        #get "from" and "to" information
        #"                      Basis "ao basis" -> "ao basis" (cartesian)"
        tmp = Current_line.split('"')
        self.basis_from=re.findall("(?<=Basis \").*(?=\" ->)", Current_line)[0] #  tmp[1]
        self.basis_to=re.findall("(?<=-> \").*(?=\" \()", Current_line)[0] #tmp[3]
        elements = GenericSection.split_section(self, [r".* \(.*\)"])
        for element in elements[1:]:
            new_basis = commonsobjects.Basis(GenericSection.return_row(self, element.section_initial_line).split("(")[0].strip())
            new_basis.values["Orbitals"] = parser_tools.get_tabulated_data(self.section_filename, 0, element.section_final_line - self.section_initial_line, "nsnn", 1, 1, None, self.section_rows[element.section_initial_line - self.section_initial_line : element.section_final_line - self.section_initial_line])
            self.basis.append(new_basis)
        self.values = {"basis": self.basis,
               "basis_from": self.basis_from,
               "basis_to": self.basis_to
               }    
    
class VectorForMolecularOrbitAnalysis(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "vector_for_molecular_orbit_analysis_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.base_function_ns = []
        if self.section_rows is None:
            Current_Line_Number = 0
            f=qchitool_tools.openfile(filename)
        
            while 1:                            #go to initial line
                Current_line = f.readline()
                if Current_Line_Number == self.section_initial_line:
                    break
                Current_Line_Number += 1
                
            tmp_line = Current_line.split("=")
            self.vector = qchitool_tools.val(tmp_line[0].split()[1])
            self.Occupancy = qchitool_tools.val(tmp_line[1].split()[0])
            self.Energy = qchitool_tools.val(tmp_line[2].split()[0])
            if len(tmp_line)>3:
                self.Symmetry = tmp_line[3].strip()
            else:
                self.Symmetry = "Not Used"
            
            Current_line = f.readline()
            Current_Line_Number += 1
            
            tmp_line = Current_line.split("=")
            tmp_mo = tmp_line[1].split(",")
            self.Molucular_Orbital_center = commonsobjects.Point3d(qchitool_tools.val(tmp_mo[0]), qchitool_tools.val(tmp_mo[1]), qchitool_tools.val(tmp_mo[2]))
            self.r2 = qchitool_tools.val(tmp_line[2])
            self.add_vector(Current_Line_Number+3, self.section_final_line,self.section_rows)
        else:
            #"Vector   49  Occ=1.000000D+00  E=-8.174234D-01  Symmetry=ag
            #  MO Center= -2.8D-16,  1.5D-13,  2.2D-20, r^2= 2.2D+01""
            #tmp_line = self.section_rows[0].split("=")
            self.vector = qchitool_tools.val(re.findall("(?<=Vector).*(?=Occ)", self.section_rows[0])[0].strip())    #tf.val(tmp_line[0].split()[1])
            self.Occupancy = qchitool_tools.val(re.findall("(?<=Occ\=).*(?=E\=)", self.section_rows[0])[0].strip()) #  tf.val(tmp_line[1].split()[0])
            if re.search(r"Symmetry", self.section_rows[0]):
                self.Energy = qchitool_tools.val(re.findall("(?<=E\=).*(?=Symmetry)", self.section_rows[0])[0].strip()) # tf.val(tmp_line[2].split()[0])
                self.Symmetry = re.findall("(?<=Symmetry).*", self.section_rows[0])[0].strip() #tmp_line[3].strip()
            else:
                self.Energy = qchitool_tools.val(re.findall("(?<=E\=).*", self.section_rows[0])[0].strip()) # tf.val(tmp_line[2].split()[0])
                self.Symmetry = "Not Used"

            tmp_line = self.section_rows[1].split("=")
            tmp_mo = re.findall("(?<=MO Center\=).*(?=, r\^2\=)", self.section_rows[1])[0].strip().split(",") #tmp_line[1].split(",")
            self.Molucular_Orbital_center = commonsobjects.Point3d(qchitool_tools.val(tmp_mo[0]), qchitool_tools.val(tmp_mo[1]), qchitool_tools.val(tmp_mo[2]))
            self.r2 = qchitool_tools.val(re.findall("(?<=, r\^2\=).*", self.section_rows[1])[0].strip()) #tf.val(tmp_line[2])
            self.add_vector(self.section_initial_line+3, self.section_final_line, self.section_rows[4:])
        self.values = {"vector": self.vector, 
               "Occupancy": self.Occupancy, 
               "Energy": self.Energy, 
               "Symmetry": self.Symmetry, 
               "Molucular_Orbital_center": self.Molucular_Orbital_center.values, 
               "r2": self.r2, 
               "base_function": self.base_function_ns}           
    
    def add_vector(self, initial_line, final_line, vector_rows = None):
        
        schema = "nnss"
        tuple_for_row = 2
        tmp_values = []
        if self.section_rows is None:
            
            Current_Line_Number = 0
            f=qchitool_tools.openfile(self.filename)
            while 1:                            #go to initial line
                f.readline()
                if Current_Line_Number == initial_line:
                    break
                Current_Line_Number += 1
                
            f.readline()
            f.readline()
            f.readline()
            Current_line = f.readline()
            Current_Line_Number += 4
        
            while Current_Line_Number >= final_line:
                tuples = parser_tools.chunks(Current_line, len(Current_line)/tuple_for_row)
                
                for tpl in tuples:
                    tmp = tpl.strip()
                    if len(tmp)>0:
                        tmp_values.append(tmp)
                Current_line = f.readline()
                Current_Line_Number += 1
            f.close()
        else:
            for row in vector_rows:
                if len(row) > 40:
                    tuples = parser_tools.chunks(row, len(row)/tuple_for_row)
                else:
                    tuples = [row]
                for tpl in tuples:
                    tmp = tpl.strip()
                    if len(tmp)>0:
                        tmp_values.append(tmp)
                        
        for value in tmp_values:
            parser_tools.separate_by_schema(value.split(), schema, 2, self.base_function_ns)


class FinalMolecularOrbitalAnalysis(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        self.data = []
        GenericSection.print_name(self, "final_molecular_orbital_analysis_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        if filename is None:
            self.spin_orbital = self.section_rows[0].split()[2]
        else:
            self.spin_orbital = parser_tools.return_file_rows(self.section_filename, self.section_initial_line, None).split()[2]
        if self.spin_orbital <> "Alpha" and self.spin_orbital <> "Beta":
            self.spin_orbital = "Both"
        sections = GenericSection.split_section(self, [r"Vector .*"])
        # values is a list of (section_name, striped_file_row, [initial_line, final_line])
        for v in sections:
            new_vector = VectorForMolecularOrbitAnalysis(self.section_filename, v.section_initial_line, v.section_final_line, v.section_final_line, self.section_rows[v.section_initial_line - self.section_initial_line: v.section_final_line - self.section_initial_line])
            self.data.append(new_vector)
        self.values = {"data": self.data,
                       "spin orbital": self.spin_orbital}
        
class OrbitalOverlaps(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "orbital_overlaps_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_multirow_table(self.section_filename, self.section_initial_line, self.section_final_information_line, "nnn", False, self.section_rows[2:])
        self.values = {"data": self.data}
        
class SummaryOfAllocatedGlobalArray(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "summary_of_allocated_global_array_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        #never used
        
class GaStatisticsForProcess(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "ga_statistics_for_process_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        #uninteresting

        
class TotalElectronDensity(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "total_electron_density_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "nnnn", 1, None, None, self.section_rows)
        self.values = {"data": self.data}

class TotalSpinDensity(GenericSection):

    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "total_spin_density_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_tabulated_data(self.section_filename, self.section_initial_line, self.section_final_information_line, "nnnn", 1, None, None, self.section_rows)
        self.values = {"data": self.data}
        
class EnergyMinimization(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        #TO CHECK
        GenericSection.print_name(self, "energy_minimization_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.data = parser_tools.get_continous_list_no_sep(self.section_filename, 5 , self.section_final_information_line- self.section_initial_line, "n", self.section_rows)
        self.values = {"data": self.data}
        
class NwchemNuclearHessianAndFrequencyAnalysis(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "nwchem_nuclear_hessian_and_frequency_analysis_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        
class NwchemFiniteDifferenceHessian(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "nwchem_finite_difference_hessian_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, self.section_initial_line, self.section_final_information_line, None, "=", True, self.section_rows)
        
class MassWeightedNuclearHessian(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "mass_weighted_nuclear_hessian_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_triangular_matrix(self.section_filename, self.section_initial_line, self.section_final_information_line, self.section_rows[4:])

class StepInfo(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        GenericSection.print_name(self, "step_info_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        #sections = tf.GenericSection.split_section(self, [r"[\s]*Optimization converged[\s]*.*"])
        
        #if len(sections)>0:
            #self.values["optimization converged"] = True
        #else:
            #self.values["optimization converged"] = False
        for row in self.section_rows:
            row = row.split()
            if row[0]=="@":
                self.values["Step"] = qchitool_tools.val(row[1]) 
                self.values["Energy"] = qchitool_tools.val(row[2]) 
                self.values["Delta E"] = qchitool_tools.val(row[3]) 
                self.values["Gmax"] = qchitool_tools.val(row[4]) 
                self.values["Grms"] = qchitool_tools.val(row[5]) 
                self.values["Xrms"] = qchitool_tools.val(row[6]) 
                if len(row)>8: #bug related to 1000000 seconds
                    self.values["Xmax"] = qchitool_tools.val(row[7]) 
                    self.values["Walltime"] = qchitool_tools.val(row[8]) 
                else:
                    self.values["Xmax"] = -1
                    self.values["Walltime"] = -1
                break
        
class OptimizationConverged(GenericSection):
        def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "optimization_converged_class")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            self.values = {"optimization converged": True}
            
class DatabaseEntry(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "database_entry_class")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            value = []
            if self.section_rows:
                
                entry_name, entry_date = re.split(r"[a-z]+[[][0-9]+[]]", self.section_rows[0], 2)
                entry_type = re.findall(r"[a-z]*[[][0-9]*[]]", self.section_rows[0])[0]
                entry_name = entry_name.strip().split(":")
                entry_date = entry_date.strip()
                entry_type = entry_type.split("[")
                entry_type[1] = qchitool_tools.val(entry_type[1].replace("]", ""))
                if entry_type[0] == "char": #string data type
                    for row in self.section_rows[1:]:
                        value.append(row.strip())
                elif entry_type[0] == "logical": #boolena data type
                    for row in self.section_rows[1:]:
                        for token in row.split():
                            if token == "f": #False value
                                value.append(False)
                            elif token == "t": #True value
                                value.append(True)
                else: #numeric data type
                    for row in self.section_rows[1:]:
                        for token in row.split():
                            value.append(qchitool_tools.val(token))
            else:
                entry_name =  entry_date = None
                entry_type = [None, None]
            self.values = {"entry_name" : entry_name,
                           "entry_type" : entry_type[0],
                           "entry_lenght" : entry_type[1],
                           "entry_date" : entry_date,
                           "values": value}
            
class AtomInformation(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "atom_information_class")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            self.values = {}
            #TODO
            
class NormalModeEigenvectors(GenericSection):
    #For each vibrational model(table column) for each atom it is associate a 3 dimensional vector (group of three value)
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "normal_mode_eigenvectors_class")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            value_dictionary = {}
            splited_len = 0
            for row in self.section_rows[3:]:
                splited_row = row.split()
                if len(splited_row)>0:
                    if splited_row[0] == "Frequency" or splited_row[0] == "P.Frequency":
                        if "Frequency" in value_dictionary:
                            value_dictionary["Frequency"] = value_dictionary["Frequency"] + splited_row[1:]
                        else:
                            value_dictionary["Frequency"] = splited_row[1:]
                        splited_len = len(splited_row)
                    elif splited_row[0] == "Normal":
                        break
                    else:
                        if len(splited_row) == splited_len:
                            if splited_row[0] in value_dictionary:
                                value_dictionary[splited_row[0]] = value_dictionary[splited_row[0]] + splited_row[1:]
                            else:
                                value_dictionary[splited_row[0]] = splited_row[1:]
                                
            data = {}
            for index in range(len(value_dictionary["Frequency"])):
                data[index] = []
                for index2 in range(len(value_dictionary)-1):
                    data[index].append(value_dictionary[str(index2 + 1)][index])
            self.values = {"data": data}
            #TODO
    
class NormalEigenvalue_DipoleMoment(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "normal_eigenvalue__dipolemomentclass")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            data = parser_tools.get_tabulated_data(self.section_filename, 
                                                           self.section_initial_line, 
                                                           self.section_final_information_line, 
                                                           "nnsnnn", 1, None, None, 
                                                           self.section_rows)
            self.values = {"data": data}
            
class NormalEigenvalue_InfraRedIntensities(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "normal_eigenvalue_InfraRedIntesities_class")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            data = parser_tools.get_tabulated_data(self.section_filename, 
                                                           self.section_initial_line, 
                                                           self.section_final_information_line, 
                                                           "nnsnnn", 1, None, None, 
                                                           self.section_rows)
            self.values = {"data": data}


class Rot_Const_And_Energy(GenericSection):
        def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
            GenericSection.print_name(self, "rotational_constants_and_energy_class")
            GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
            self.values = {}
            re_A_C = re.compile(r"[A-C]\=.*\(.*\)")
            re_totalentropy = re.compile(r"Total Entropy")
            re_star = re.compile(r".*\-")
            re_constant_volume = re.compile(r"Cv \(constant volume heat capacity\)")
            for row in self.section_rows:
                if re_A_C.search(row):
                    self.values[row.split("=")[0].strip()] = row.split("(")[0].split("=")[1].split()[0]
            sections = GenericSection.split_section(self, [r".*Temperature.*",
                                                           r".*Total Entropy.*",
                                                           r".*constant volume heat capacity.*"])
            
            self.values.update(parser_tools.get_listvalue_sep(self.section_filename, 
                                                              sections[0].section_initial_line, 
                                                              sections[0].section_final_line, 
                                                              None, "=", True, 
                                                              self.section_rows[sections[0].section_initial_line - self.section_initial_line-1: sections[0].section_final_line - self.section_initial_line]))
            
            for row in self.section_rows[sections[1].section_initial_line - self.section_initial_line-1: sections[1].section_final_line - self.section_initial_line]:
                if re_totalentropy.search(row):
                    self.values[row.split("=")[0].strip()] = qchitool_tools.val(row.split("=")[1].split()[0])
                elif re_star.match(row):
                    self.values["Total Entropy " + row.split("=")[0].strip()] = qchitool_tools.val(row.split("=")[1].split()[0])
            if len(sections)>2:
                for row in self.section_rows[sections[2].section_initial_line - self.section_initial_line-1: sections[2].section_final_line - self.section_initial_line]:
                    if re_constant_volume.search(row):
                        self.values[row.split("=")[0].strip()] = qchitool_tools.val(row.split("=")[1].split()[0])
                    elif re_star.match(row):
                        self.values["Cv " + row.split("=")[0].strip()] = qchitool_tools.val(row.split("=")[1].split()[0])
