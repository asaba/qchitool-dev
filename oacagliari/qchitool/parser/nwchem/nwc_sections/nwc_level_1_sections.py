'''
Created on 07/mar/2011

@author: andrea saba
@email: asaba@oa-cagliari.inaf.it

classes define theory_level_sections
'''
from ...tools import parser_tools
from ....tools import qchitool_tools
from .... import performancetest as pt
from nwc_sec_commons import GenericSection  
import re

class NwchemDftModule(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        self.converged = False
        self.values = {}
        GenericSection.print_name(self, "nwchem_dft_module_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        sections = GenericSection.split_section(self, [r"The DFT is already converged"])
        # [(section_name, striped_file_row, [initial_line, final_line]]
        if len(sections)>0:
            self.converged = True
            self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                         sections[0].section_initial_line+2, 
                                                         sections[0].section_final_line, 
                                                         None, "=", True, 
                                                         self.section_rows[sections[0].section_initial_line - self.section_initial_line +2 : sections[0].section_final_line - self.section_initial_line])
        else:
            values_def = [[r"Total DFT energy", "=", re.compile(r"Total DFT energy")],
                          [r"One electron energy", "=", re.compile(r"One electron energy")],
                          [r"Coulomb energy", "=", re.compile(r"Coulomb energy")], 
                          [r"Exchange\-Corr\. energy", "=", re.compile(r"Exchange\-Corr\. energy")], 
                          [r"Exchange\-Corr\. energy", "=", re.compile(r"Exchange\-Corr\. energy")], 
                          [r"Nuclear repulsion energy", "=", re.compile(r"Nuclear repulsion energy")], 
                          [r"Numeric\. integr\. density", "=", re.compile(r"Numeric\. integr\. density")], 
                          [r"Total iterative time", "=", re.compile(r"Total iterative time")] ]
            for values_tmp in values_def:
                tmp_value = GenericSection.find_value(self, values_tmp[0], values_tmp[1], values_tmp[2], True)
                if not tmp_value is None:
                    if qchitool_tools.val(tmp_value[1]) is None:
                        self.values[tmp_value[0]] = tmp_value[1]
                    else:
                        self.values[tmp_value[0]] = qchitool_tools.val(tmp_value[1])
        self.values["converged"]= self.converged
        pt.p("end")

class NwchemPropertyModule(GenericSection):
    #Documentation http://www.nwchem-sw.org/index.php/Properties
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        self.values = {}
        GenericSection.print_name(self, "nwchem_property_module_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        sections = GenericSection.split_section(self, [r".*property selection.*",
                                                          r".*Atomic weights.*",
                                                          r".*Center of charge.*"])
        # [(section_name, striped_file_row, [initial_line, final_line]]
        if len(sections)>0:
            self.properties_calculated = parser_tools.get_listvalue_sep(self.section_filename, 
                                                              sections[0].section_initial_line+1, 
                                                              sections[0].section_final_line, 
                                                              None, "=", False, 
                                                              self.section_rows[sections[0].section_initial_line - self.section_initial_line +1 : sections[0].section_final_line - self.section_initial_line])
            self.atomic_weights_unit = parser_tools.get_parenteses_value(self.section_rows[sections[1].section_initial_line - self.section_initial_line], [r"\(", r"\)"])
            self.atomic_weights = parser_tools.get_tabulated_data(self.section_filename, 
                                                         sections[1].section_initial_line+1,
                                                         sections[1].section_final_line, 
                                                         "nsn", 1, 1, None, 
                                                         self.section_rows[sections[1].section_initial_line - self.section_initial_line +1 : sections[1].section_final_line - self.section_initial_line])
            self.center_of_charge_unit =  parser_tools.get_parenteses_value(self.section_rows[sections[2].section_initial_line - self.section_initial_line], [r"\(", r"\)"])
            tmp = parser_tools.get_listvalue_sep(self.section_filename, 
                                       sections[2].section_initial_line+1,
                                       sections[2].section_final_line, 
                                       None, "=", False, 
                                       self.section_rows[sections[2].section_initial_line - self.section_initial_line +1 : sections[2].section_final_line - self.section_initial_line])
            if "X" in tmp:
                self.center_of_charge = [qchitool_tools.val(tmp['X']), qchitool_tools.val(tmp['Y']), qchitool_tools.val(tmp['Z'])]
            else:
                self.center_of_charge = [qchitool_tools.val(tmp['x']), qchitool_tools.val(tmp['y']), qchitool_tools.val(tmp['z'])]
        
            self.values = {"properties_calculated": self.properties_calculated,
                           "atomic_weights_unit": self.atomic_weights_unit,
                           "atomic_weights": self.atomic_weights,
                           "center_of_charge_unit": self.center_of_charge_unit,
                           "center_of_charge": self.center_of_charge}
            if self.properties_calculated["nodip"] == 0:
                #dipole moment
                property_sections = GenericSection.split_section(self, [r".*Dipole moment.*"])
                if len(property_sections) > 0:
                    unit = self.section_rows[sections[0].section_initial_line - self.section_initial_line].split()[-1]
                    x = qchitool_tools.val(self.section_rows[sections[0].section_initial_line +1 - self.section_initial_line].split()[1])
                    y = qchitool_tools.val(self.section_rows[sections[0].section_initial_line +2 - self.section_initial_line].split()[1])
                    z = qchitool_tools.val(self.section_rows[sections[0].section_initial_line +3 - self.section_initial_line].split()[1])
                    self.values.update({"dipolemoment": {"data": [x, y, z],
                                   "units_of_measurement": unit}})
            if self.properties_calculated["noqdp"] == 0:
                #quadrupole moment
                pass
            if self.properties_calculated["nootp"] == 0:
                #octupole moment
                pass
            if self.properties_calculated["nopop"] == 0:
                #Mulliken population analysis and bond order analysis 
                pass
            if self.properties_calculated["nospin"] == 0:
                pass
            if self.properties_calculated["noloc"] == 0:
                pass
            if self.properties_calculated["nodpl"] == 0:
                pass
            if self.properties_calculated["nodhp"] == 0:
                pass
            if self.properties_calculated["nod2hp"] == 0:
                pass
            if self.properties_calculated["nofpl"] == 0:
                pass
            if self.properties_calculated["nofhp"] == 0:
                pass
            if self.properties_calculated["nof2hp"] == 0:
                pass
            if self.properties_calculated["nosos"] == 0:
                pass
            if self.properties_calculated["noelp"] == 0:
                pass
            if self.properties_calculated["noelf"] == 0:
                #Electric field at nuclei 
                pass
            if self.properties_calculated["noelfg"] == 0:
                #Electric field gradient at nuclei 
                pass
            if self.properties_calculated["noden"] == 0:
                pass
            if self.properties_calculated["nogiao"] == 0:
                #NMR shielding (GIAO method)
                pass
            if self.properties_calculated["noiglo"] == 0:
                pass
            if self.properties_calculated["noston"] == 0:
                pass
        pt.p("end")
            
class NwchemCphfModule(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        self.values = {}
        GenericSection.print_name(self, "nwchem_CPHF_module_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        sections = GenericSection.split_section(self, [r".*scftype.*", 
                                                       r".*Iterative solution of linear equations.*",
                                                       r".*iter[\s]*nsub[\s]*residual[\s]*time.*",
                                                       r".*HESSIAN:.*"])
        if len(sections) > 1:
            
            data1 = parser_tools.get_listvalue_sep(self.section_filename, 
                                       sections[0].section_initial_line,
                                       sections[0].section_final_line, 
                                       None, "=", True, 
                                       self.section_rows[sections[0].section_initial_line - self.section_initial_line : sections[0].section_final_line - self.section_initial_line])
            
            data2 = parser_tools.get_listvalue_no_sep(self.section_filename, 
                                       sections[1].section_initial_line + 1,
                                       sections[1].section_final_line, 
                                       self.section_rows[sections[1].section_initial_line + 1 - self.section_initial_line : sections[1].section_final_line - self.section_initial_line])
        else:
            data1 = None
            data2 = None
        self.values["CPHF_values"] = data1
        self.values["CPHF_values2"] = data2
        pt.p("end")
        
        
class DftGradientModule(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "dft_gradient_module_supersection_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                     self.section_initial_line, 
                                                     self.section_final_information_line, 
                                                     None, "=", True, self.section_rows).copy()  #shallow copy
        pt.p("end")                                  

class VibrationalAnalysisViaTheFXMethod(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "vibrational_analysis_via_the_fx_method_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        trans_and_rot = re.compile(r".*with translations and rotations projected out.*")
        projected = re.compile(r".*Projected Nuclear Hessian.*")
        self.values["project_out"] = False
        for row in self.section_rows:
            if trans_and_rot.search(row):
                self.values["project_out"] = True
            elif projected.search(row):
                parser_tools.get_singlevalue_sep(row, None, ":", True, self.values)
        pt.p("end")

class NwchemScfModule(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "nwchem_scf_module_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        self.converged = False
        sections = GenericSection.split_section(self, [r"The SCF is already converged"])
        if len(sections)>0:
            self.converged = True
            self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                         sections[0].section_initial_line+2, 
                                                         sections[0].section_final_line, 
                                                         None, "=", True, 
                                                         self.section_rows[sections[0].section_initial_line - self.section_initial_line +2 : sections[0].section_final_line - self.section_initial_line])
        else:
            values_def = [[r"Total SCF energy", "=", re.compile(r"Total SCF energy")]]
            for values_tmp in values_def:
                tmp_value = GenericSection.find_value(self, values_tmp[0], values_tmp[1], values_tmp[2], True)
                if not tmp_value is None:
                    if qchitool_tools.val(tmp_value[1]) is None:
                        self.values[tmp_value[0]] = tmp_value[1]
                    else:
                        self.values[tmp_value[0]] = qchitool_tools.val(tmp_value[1])
        self.values["converged"]= self.converged
        pt.p("end")