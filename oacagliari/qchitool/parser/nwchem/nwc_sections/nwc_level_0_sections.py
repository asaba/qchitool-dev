'''
Created on 10/mar/2011

@author: andrea saba
@email: asaba@oa-cagliari.inaf.it

classes define super_sections
'''

from ...tools import parser_tools
from ....tools import qchitool_tools
from .... import performancetest as pt
from nwc_sec_commons import GenericSection  
import re

class Arguments(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "arguments_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                     self.section_initial_line, 
                                                     self.section_final_information_line, 
                                                     None, "=", True, self.section_rows).copy()  #shallow copy
        pt.p("end")

class Step(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "Step")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {'Step': qchitool_tools.val(GenericSection.go_to_first_line(self, None)[1].split()[-1])}
        pt.p("end")

class InputDeck(GenericSection):
    #TO DO
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        self.lines = []
        GenericSection.print_name(self, "input_deck_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        Current_Line_Number = 0
        if self.section_rows is None:
            
            f=qchitool_tools.openfile(self.section_filename)
        
            while 1:                            #go to initial line
                Current_line = f.readline()
                if Current_Line_Number == self.section_initial_line:
                    break
                Current_Line_Number += 1
                
            while Current_Line_Number < self.section_final_line:
                self.lines.append(Current_line.strip())
            f.close()
        else:
            for row in self.section_rows:
                self.lines.append(row.strip())
        self.values = {"lines": self.lines}
        pt.p("end")
        
class NwchemGeometryOptimization(GenericSection):

    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "nwchem_geometry_optimization_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                     self.section_initial_line, 
                                                     self.section_final_information_line, 
                                                     None, "=", True, 
                                                     self.section_rows).copy()  #shallow copy
        pt.p("end")

class NwchemManifest(GenericSection):

    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        self.values = dict({})
        GenericSection.print_name(self, "nwchem_manifest_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values['version'] = GenericSection.go_to_first_line(self, None)[1].split()[-1]
        #this value is present in the DB schema
        pt.p("end")

class NwchemInputModule(GenericSection):

    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "nwchem_input_module_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = parser_tools.get_listvalue_sep(self.section_filename, 
                                                     self.section_initial_line, 
                                                     self.section_final_information_line, 
                                                     None, ":", True, 
                                                     self.section_rows).copy()  #shallow copy
        #add citation section
        pt.p("end")


class HessianStep(GenericSection):
    
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "hessian_step_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        #" atom:  24 xyz: 1(-) wall time:   13348.5      date:  Fri Oct  8 13:01:59 2004"
        first_row = GenericSection.return_row(self, self.section_initial_line)
        self.atom = qchitool_tools.val(re.findall("(?<=atom:).*(?=xyz:)", first_row)[0])
        self.xyz = re.findall("(?<=xyz:).*(?=wall time:)", first_row)[0].strip()
        self.wall_time = qchitool_tools.val(re.findall("(?<=wall time:).*(?=date:)", first_row)[0])
        self.date = re.findall("(?<=date:).*", first_row)[0].strip()
        
        #result = tf.GenericSection.return_row(self, self.section_initial_line).split(":", 4)
        #self.atom = tf.val(result[1].split()[0].strip())
        #self.xyz = result[2].split()[0].strip()
        #self.wall_time = tf.val(result[3].split()[0].strip())
        #self.date = result[-1].strip()
        self.values = {"atom": self.atom,
                       "xyz": self.xyz,
                       "wall_time": self.wall_time,
                       "date": self.date}
        pt.p("end")


        
        
class RTDBDatabase(GenericSection):
    def __init__(self, filename, initial_line, end_line, end_information_line, section_rows = None):
        pt.p("start")
        GenericSection.print_name(self, "RTDB_database_class")
        GenericSection.__init__(self, filename, initial_line, end_line, end_information_line, section_rows)
        self.values = {}
        pt.p("end")