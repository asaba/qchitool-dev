'''
Created on 26/ago/2011

@author: asaba
'''
import re
from ....tools import qchitool_tools
from ...tools import parser_tools
from .... import performancetest as pt



class SectionFound:
   
    def __init__(self, section_name, striped_file_row, section_initial_line, section_final_line, section_final_information_line, section_class):
        self.section_name = section_name
        self.striped_file_row = striped_file_row
        self.section_initial_line = section_initial_line
        self.section_final_line = section_final_line
        self.section_final_information_line = section_final_information_line
        self.section_class = section_class
        self.section_object = None
        self.list_of_subsections = []
        
    def print_me(self):
        print "section_name = " + self.section_name
        print "striped_file_row = " + self.striped_file_row
        print "section_initial_line = " + str(self.section_initial_line)
        print "section_final_information_line = " + str(self.section_final_information_line)
        print "section_final_line = " + str(self.section_final_line)
        print "section_class = " + str(self.section_class)
        for sec in self.list_of_subsections:
            sec.print_me()
            
    def __unicode__(self):
        return "SectionFound Object:" + "section_name = " + self.section_name
    
class GenericSection:
    #This class is a generic section object
    def __init__(self, passed_filename, passed_initial_line, passed_final_line, passed_final_information_line, section_rows = None):
        self.section_final_line = passed_final_line
        self.section_initial_line = passed_initial_line
        if passed_final_information_line == 0:
            self.section_final_information_line = passed_final_line
        else:
            self.section_final_information_line = passed_final_information_line
        if section_rows is None:
            self.section_filename = passed_filename
            self.section_rows = None
        else:
            self.section_filename = None
            self.section_rows = section_rows #list of file row
        
    
    def go_to_first_line(self, file_pointer = None):
        #skip all lines before the initial line on the section
        #return [the initial line number, the initial line]
        if self.section_rows is None:
            
            Current_Line_Number = 0

            while 1: #go to initial line
                Current_line = file_pointer.readline()
                if Current_Line_Number == self.section_initial_line:
                    break
                Current_Line_Number += 1
        
            while Current_line == "": #skip empty lines
                Current_line = file_pointer.readline()
                Current_Line_Number += 1
            
            
        else:
            index = 0
            for row in self.section_rows:
                if row.strip() == "":
                    index += 1
                else:
                    break
            Current_line = self.section_rows[index]
            
        return [self.section_initial_line + index, Current_line]
                    
    def print_section(self):
        #print all lines between the first and the last line of the section
        if self.section_rows is None:
            
            f=qchitool_tools.openfile(self.filename)
            Result = self.go_to_first_line(f)
            Current_Line_Number = Result[0]
        
            print Current_Line_Number & " " & Result[1]
                
            while Current_Line_Number < self.section_final_line: #skip empty lines
                print Current_Line_Number + 1 + " " + f.readline()
                Current_Line_Number += 1
        else:
            Result = self.go_to_first_line()
            Current_Line_Number = Result[0]
            for row in self.section_rows:
                print Current_Line_Number + " " + row
                Current_Line_Number += 1
                
    def print_name(self, name):
        pass
        #print "Creating " + name + " object..."
                
    def __unicode__(self):
        #Return value associated to this section
        if self.values:
            return str(self.values)
        else:
            return "No values"
                
    def split_section(self, list_of_subsection_re):
        #this function split the section by the list_of_subsection_re and return the initial and final row for each subsection found
        sections_found = [] # [list of sectionsplited]
        f = None
        if self.section_rows is None:
            f=qchitool_tools.openfile(self.filename)
        Result = self.go_to_first_line(f)
        Current_Line_Number = Result[0]
        Current_Line = Result[1]
        re_dict = {}
        for subsection_re in list_of_subsection_re:
            re_dict[subsection_re] = re.compile(subsection_re)
        if self.section_rows is None:
            self.section_rows = parser_tools.return_file_rows(self.section_filename, self.section_initial_line, self.section_final_line)
        for row in self.section_rows:
            Current_Line = row.strip()
            for (current_re, compiled_re) in re_dict.items():
                if compiled_re.match(Current_Line):
                    if len(sections_found) > 0:
                        sections_found[-1].section_final_line = Current_Line_Number-1
                        sections_found[-1].section_final_information_line = Current_Line_Number-1
                    new_section_found  = SectionFound(current_re, Current_Line, Current_Line_Number, 0,0, None)
                    sections_found.append(new_section_found)
                    break
            Current_Line_Number += 1
        if len(sections_found)>0:
            sections_found[-1].section_final_line = Current_Line_Number-1
            if sections_found[0].section_initial_line - 1 >  self.section_initial_line:
                self.section_final_line = sections_found[0].section_initial_line - 1
                self.section_final_information_line = sections_found[0].section_initial_line - 1
            else:
                self.section_final_line = sections_found[-1].section_final_line
                self.section_final_information_line = sections_found[-1].section_final_line
        return sections_found

    def return_row(self, line_number):
        #this function return the line-number-th line
        if line_number < self.section_initial_line or line_number > self.section_final_line:
            #line requested is out of range
            return None
        else:
            if self.section_rows is None:
                return parser_tools.return_file_rows(self.section_filename, line_number, None)
            else:
                return self.section_rows[line_number - self.section_initial_line]
    
    def find_value(self, re_valuename, separator, re_compiled, search_in_all_section):
        element_pattern = re_compiled
        if search_in_all_section:
            search_limit = self.section_final_line - self.section_initial_line + 1
        else:
            search_limit = self.section_final_information_line - self.section_initial_line + 1
        if self.section_rows is None:
            self.section_rows = parser_tools.return_file_rows(self.section_filename, self.section_initial_line, self.section_final_line)
        for row in self.section_rows[:search_limit]:
            if  element_pattern.search(row):
                tmp = row.strip().split(separator)
                found = True
                break
            else:
                found = False
        if found:
            return [tmp[0].strip(), tmp[-1].strip()]
        else:
            return None
        
        
