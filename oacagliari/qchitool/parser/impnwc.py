'''
Created on 09/feb/2011
Module for parsing of NWchem output files
@author: andrea saba
@email: asaba@oa-cagliari.inaf.it
'''

import re


from .nwchem.nwc_sections import nwc_sec_commons, nwc_mapping
from .nwchem import commons#, check_import_object
#from .commonsobjects import TaskCode
from ..tools import qchitool_tools
from .. import performancetest as pt
import tempfile



def create_objects(current_section, file_rows, filename, load_file_in_memory):
    for section in current_section.list_of_subsections:
        create_objects(section, file_rows, filename, load_file_in_memory)
        if section.section_class is None:
            print "Object undefined."
        else:
            if load_file_in_memory:
                section.section_object = section.section_class(None, section.section_initial_line, section.section_final_line, section.section_final_information_line, file_rows[section.section_initial_line : section.section_final_line])
            else:
                section.section_object = section.section_class(filename, section.section_initial_line, section.section_final_line, section.section_final_information_line, None)
            #print "Object created."
            section.section_class = str(section.section_class)
            
def set_final_information_line(sections_found):

    #This function set the end of information in each section by the start of the first subsection

    for section in sections_found:
        if len(section.list_of_subsections)>0:
            section.section_final_information_line = section.list_of_subsections[0].section_initial_line # -1
            set_final_information_line(section.list_of_subsections)
        else:
            if section.section_final_information_line != 0:
                section.section_final_information_line = section.section_final_line

    
def clean_output_file_rows(file_content, cutter_extremes):
    #make something
    pt.p("start")
    file_content = file_content.replace("\r\n", "\n")
    #temporaryinputfile1 = tempfile.NamedTemporaryFile(mode="w")
    #temporaryinputfile1.writelines(file_content)
    #temporaryinputfile1.seek(0)
    #result = ""
    for extremes in cutter_extremes:
        #for s in re.split(extremes[0] + "|" + extremes[1], file_content)[::2]:
        #    result = result + s
        file_content = "".join(s for s in re.split(extremes[0] + "|" + extremes[1], file_content)[::2])
        #temporaryinputfile2 = tempfile.NamedTemporaryFile(mode="w")
        #temporaryinputfile2.writelines(result)
        #temporaryinputfile2.seek(0)
    
    #return splited file_content
    #file_content = re.split("\n", result)
    file_content = re.split("\n", file_content)
    pt.p("end")
    return file_content
        

def parsing(filename, load_file_in_memory, filetype):
    #find sections in "filename" with starting and ending line
    pt.p("start")
    Current_row_index = 0
    sections = nwc_mapping.MappingSections(filetype)
    sections_found = []  #list of section_class
    last_found = ""
    file_rows = []
    next_row = False
    sections_schema = "" #it is the schema of sections found (S = super_section, T = theory_level_section, U = sub_sections)

    # [(section_name, striped_file_row, [initial_line, final_line], section_class, [subsections,...]), ]
    f = qchitool_tools.openfile(filename)

    if load_file_in_memory:
        file_rows = clean_output_file_rows(f.read(), commons.cutter_extremes)
        for row in file_rows:
            if Current_row_index == 3888:
                pass
            tmp_row = row.strip()
            if len(tmp_row)==0:
                Current_row_index += 1
                continue
            for super_section_tuple in sections.levels[0]: #in the level 0 sections we can't found a super_section named as the previous super_section
                if super_section_tuple.re_ex_compiled.match(tmp_row):  
                    new_super_section = nwc_sec_commons.SectionFound(super_section_tuple.name, tmp_row, Current_row_index, 0, 0, super_section_tuple.relatedclass)
                    if len(sections_schema) > 0:
                        if sections_found[-1].section_name <> new_super_section.section_name:
                            #close last super_section
                            sections_found[-1].section_final_line = Current_row_index-1 
                            if sections_schema[-1] == "T":
                                sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1  #close last theory_section or last sub_section
                            elif sections_schema[-1] == "U":
                                sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1  #close last theory_section or last sub_section
                                if len(sections_found[-1].list_of_subsections[-1].list_of_subsections)>0:  #close the last, eventually, sub_section
                                    sections_found[-1].list_of_subsections[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                            sections_found.append(new_super_section)
                            sections_schema += "S"
                            next_row = True
                            break
                    else:
                        sections_found.append(new_super_section)
                        sections_schema += "S"
                        next_row = True
                        break


            if next_row:
                next_row = False
                Current_row_index += 1
                continue

            for theory_tuple in sections.levels[1]:
                if theory_tuple.re_ex_compiled.match(tmp_row):
                    new_theory_level_section = nwc_sec_commons.SectionFound(theory_tuple.name, tmp_row, Current_row_index,0,0, theory_tuple.relatedclass)
                    if sections_schema[-1] == "T":
                        sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                    elif sections_schema[-1] == "U":
                        sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                        if len(sections_found[-1].list_of_subsections[-1].list_of_subsections)>0:#close the last, eventually, sub_section
                            sections_found[-1].list_of_subsections[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                    sections_found[-1].list_of_subsections.append(new_theory_level_section)
                    sections_schema += "T"
                    next_row = True
                    break

            if next_row:
                next_row = False
                Current_row_index += 1
                continue

            for section_tuple in sections.levels[2]:
                if section_tuple.re_ex_compiled.match(tmp_row):
                    new_sub_section = nwc_sec_commons.SectionFound(section_tuple.name, tmp_row, Current_row_index, 0,0, section_tuple.relatedclass)
                    if sections_schema[-1] == "S":
                        sections_found[-1].list_of_subsections.append(new_sub_section)
                    elif sections_schema[-1] == "T":
                        sections_found[-1].list_of_subsections[-1].list_of_subsections.append(new_sub_section)
                    elif sections_schema[-1] == "U":
                        #close last sub_section
                        if len(sections_found[-1].list_of_subsections)>0:
                            if len(sections_found[-1].list_of_subsections[-1].list_of_subsections)>0:
                                sections_found[-1].list_of_subsections[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                                sections_found[-1].list_of_subsections[-1].list_of_subsections.append(new_sub_section)
                            else:
                                sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                                sections_found[-1].list_of_subsections.append(new_sub_section)
                    sections_schema += "U"
                    break
            Current_row_index += 1
    else:
        while 1: #find sections in "filename"
            line=f.readline()
            if not line: break
            section_title = line.strip()

            #ADD super_sections
            for super_section_tuple in sections.levels[0]:
                if super_section_tuple.re_ex_compiled.match(section_title):
                    if len(sections_found) > 0:
                        sections_found[-1].section_final_line = Current_row_index # - 1
                        if len(sections_found[-1].list_of_subsections)>0: #list of subsections
                            sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                    new_super_section = nwc_sec_commons.SectionFound(super_section_tuple.name, section_title, Current_row_index, 0,0, super_section_tuple.relatedclass)
                    sections_found.append(new_super_section)
                    last_found = "super_section"
                    sections_schema += "S"
                    break

            for theory_tuple in sections.levels[1]:
                if theory_tuple.re_ex_compiled.match(section_title):
                    if len(sections_found[-1].list_of_subsections) > 0:
                        sections_found[-1].list_of_subsections[-1].section_final_line = Current_row_index # - 1
                        if len(sections_found[-1].list_of_subsections)>0:
                            sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index #- 1
                    new_theory_level_section = nwc_sec_commons.SectionFound(theory_tuple.name, section_title, Current_row_index,0,0, theory_tuple.relatedclass)
                    sections_found[-1].list_of_subsections.append(new_theory_level_section)
                    last_found = "theory_section"
                    sections_schema += "T"
                    break

            for section_tuple in sections.levels[2]:
                if section_tuple.re_ex_compiled.match(section_title):
                    if len(sections_found[-1].list_of_subsections)>0:
                        sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index # - 1
                    new_sub_section = nwc_sec_commons.SectionFound(section_tuple.name, section_title, Current_row_index, 0,0, section_tuple.relatedclass)
                    sections_schema += "U"
                    if last_found == "theory_section":
                        #add this section to the last one theory section
                        sections_found[-1].list_of_subsections[-1].list_of_subsections.append(new_sub_section)
                    else:
                        #'add this section to the last one super section
                        sections_found[-1].list_of_subsections.append(new_sub_section)
                    break

            Current_row_index += 1

        f.close()
        
    sections_found[-1].section_final_line = Current_row_index # - 1

    if len(sections_found[-1].list_of_subsections)>0:
        if sections_found[-1].list_of_subsections[-1].section_final_line == 0:
            sections_found[-1].list_of_subsections[-1].section_final_line= Current_row_index #- 1
        if len(sections_found[-1].list_of_subsections[-1].list_of_subsections)>0:
            if sections_found[-1].list_of_subsections[-1].list_of_subsections[-1].section_final_line == 0:
                sections_found[-1].list_of_subsections[-1].list_of_subsections[-1].section_final_line = Current_row_index # - 1

    set_final_information_line(sections_found)

    #print "Sections founded, creating objects..."
    #print sections_schema
    for section in sections_found:
        #section.print_me()
        create_objects(section, file_rows, filename, load_file_in_memory)
        if not section.section_class is None:
            if load_file_in_memory:
                section.section_object = section.section_class(None, section.section_initial_line, section.section_final_line, section.section_final_information_line, file_rows[section.section_initial_line : section.section_final_line])
            else:
                section.section_object = section.section_class(filename, section.section_initial_line, section.section_final_line, section.section_final_information_line, None)
            #print "Object created."
            section.section_class = str(section.section_class)
        else:
            print "Object undefined."
    pt.p("end")
    return sections_found


