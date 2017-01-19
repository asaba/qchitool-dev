'''
Created on 26/ago/2011

@author: asaba
'''

import re
import sys
import tools.parser_tools
from ..models import Elements, BasisSets
from .. import performancetest as pt
from nwchem.commons import nwchem_task_basis_sets

class Basis:
    #set of Gaussians that define the model of orbital for one element
    def __init__(self, element):
        self.values = {"Element": element,
                       "Orbitals" : []}
    def __unicode__(self):
        return str(self.values)

        
class GeneralValues():
    #generic class for list of value (foo = 99 or foo : something ...)
    def __init__(self, values):
        self.values = dict({})
        self.values.update(values)

class Tensor():
    #this class is a generic matrix
    def __init__(self, rows):
        self.matrix = []
        for i in rows:
            self.matrix.append(i)
        self.values = {"matrix" : self.matrix}

    
class Point3d():
    #This class is a generic point in the 3 dimensional Cartesian space
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.values = {"x":self.x, "y":self.y, "z":self.z}

        
class CenterOfMass(Point3d):
    #The center of mass is a 3 dimensional point
    def __init__(self, x, y, z):
        Point3d.__init__(self, x, y, z)

class MomentOfInertia(Tensor):
    #the moment of inertia is a tensor
    def __init__(self, rows):
        Tensor.__init__(self, rows)
        self.values.update({"dimension" : len(self.matrix)}) 

class EnergyForDerivativeDipole():
    def __init__(self, vector, atom, direction, energy, unit):
        self.vector = vector
        self.atom = atom
        self.direction = direction
        self.energy = energy
        self.unit = unit
        self.values ={"vector": self.vector, 
              "atom": self.atom, 
              "direction": self.direction,
              "energy": self.energy,
              "unit": self.unit }
    
class Atom:
    def __init__(self, atomicmass, atomicnum, coords, exactmass, formalcharge, 
                 heavyvalence, heterovalence, hyb, idx, implicitvalence, isotope,
                 partialcharge, spin, type, valence, vector):
        self.atomicmass = atomicmass
        self.atomicnum = atomicnum
        self.coords = coords
        self.exactmass = exactmass
        self.formalcharge = formalcharge
        self.heavyvalence = heavyvalence
        self.heterovalence = heterovalence
        self.hyb = hyb
        self.idx = idx
        self.implicitvalence = implicitvalence
        self.isotope = isotope
        self.partialcharge = partialcharge
        self.spin = spin
        self.type = type
        self.valence = valence
        self.vector = vector


class Molecule:
    def __init__(self, inchi, inchikey, formula, energy, 
                 exactmass, charge, aromatic_cycles, geometry, atoms, xyz):
        self.inchi = inchi
        self.inchikey = inchikey
        self.formula = formula
        self.energy = energy
        self.exactmass = exactmass
        self.charge = charge
        self.aromatic_cycles = aromatic_cycles
        self.geometry = geometry
        self.atoms = atoms
        self.xyz = xyz
        
class TaskCode: 
    def __init__(self, theorylevel, task, molecule_info, 
                 electronic_state, rotational_constants, 
                 vibration_analisis_armonic, dipolemoment,
                 molecular_specie, geometry, relative_file_path,
                 elements):
        self.task = task
        #calculation info
        self.theorylevel = theorylevel
        #input info
        self.moleculeinfo = molecule_info # basis sets, geometry, charge, spinmultiplicity, inchi , inchikey,...
        #ouput info
        self.electronicstate = electronic_state # energy ...
        self.dipolemoment = dipolemoment
        self.rotationalconstants = rotational_constants
        self.vibrationanalisisarmonic = vibration_analisis_armonic #tabulated vibrations...
        self.molecular_specie = molecular_specie
        self.geometry = geometry
        self.relative_file_path = relative_file_path
        self.elements = elements
    def __unicode__(self):
        returnstring = "Theory level: " + str(self.theorylevel.name) + "; "
        returnstring += "Operation: " + str(self.task.name)
        #returnstring += "Basis Set: " + str(self.basissets)
        return returnstring
    
class MoleculeInfo:
    def __init__(self, inchi, inchikey, 
                 uncharged_molecule, spinmultiplicity, charge,
                 geometry_list, 
                 xyz_file_name, xyz, 
                 sdf_file_name, sdf, 
                 cml_file_name, cml, 
                 atom_basis_set):
        self.inchi = inchi
        self.inchikey = inchikey
        self.moleculeobject = uncharged_molecule #OBMolecule
        self.spinmultiplicity = spinmultiplicity
        self.charge = charge
        self.geometries = {"xyz": {"filename" : xyz_file_name, 
                                   "content" : xyz},
                           "sdf": {"filename" : sdf_file_name, 
                                   "content" : sdf},
                           "cml": {"filename" : cml_file_name, 
                                   "content" : cml}}
        self.atom_basis_set = atom_basis_set #dictionary {OBAtom.index : basis set name}
        self.geometry_list = geometry_list
        #if self.inchikey:
        #    print "MoleculeInfo inchikey: " + self.inchikey
        #else:
        #    print "MoleculeInfo inchikey: None"
   
    def __unicode__(self):
        returnstring = "Inchi: " + str(self.inchi) + "; "
        returnstring += "InchiKey: " + str(self.inchikey) + "; "
        returnstring += "Spin Multiplicity: " + str(self.spinmultiplicity) + "; "
        returnstring += "Charge: " + str(self.charge) + "; "
        returnstring += "xyz geometry: " + str(self.geometries["xyz"]["filename"]) + "; "
        returnstring += "sdf geometry: " + str(self.geometries["sdf"]["filename"]) + "; "
        returnstring += "cml geometry: " + str(self.geometries["cml"]["filename"]) + "; "
        returnstring += "InchiKey: " + str(self.atom_basis_set) + ";"
        returnstring += "Basis Sets: " + str(self.atom_basis_set)
        
        return returnstring

            

class InputOuputSection():
    def __init__(self, inputsection, outputsection):
        self.inp_section = inputsection
        self.out_section = outputsection
        self.task = self.extracttasks()
        self.operation = None
        self.theorylevel = None
        self.xc_name = None
        self.xc_description = None
        self.molecule_info = MoleculeInfo(None, None, None, 
                                          None, None, None,
                                          None, None, None,
                                          None, None, None,
                                          None)
        
        
    def extracttasks(self):
        #change this part to add new task operation

        if re.search("NwchemNuclearHessianAndFrequencyAnalysis", self.out_section[0].section_class):   #NwchemNuclearHessianAndFrequencyAnalysis
            self.operation = "freq"
        elif re.search("NwchemGeometryOptimization", self.out_section[0].section_class):
            self.operation = "optimize"
        elif re.search("NwchemDftModule", self.out_section[0].section_class):
            self.operation = "energy"            
        elif re.search("NwchemPropertyModule", self.out_section[0].section_class):
            self.operation = "property"    
        else:
            self.operation = "ERROR"
            
        #Change this part for add new theory levels
        dft_section = tools.parser_tools.returnsection("NwchemDftModule", self.out_section, "First")
        if dft_section:
            self.theorylevel = "DFT"
            xc_section = tools.parser_tools.returnsection("XcInformation", dft_section, "First")
            if xc_section:
                self.xc_name = xc_section.section_object.values["method"]
                self.xc_description = str(xc_section.section_object.values["values_list"])
        #elif not (tools.parser_tools.returnsection("NwchemCphfModule", self.out_section, "First") is None):
        #    theorylevels = "SCF"
        else:
            self.theorylevel = "ERROR"
        
class Tasks_Forms():
    def __init__(self, calculation, chemistrycode, 
                 tasks, theorylevels, elements, 
                 basissets, molecularspecies, electonicstates, 
                 dipolemoments, geometries, rotationalconstants, 
                 vibrationalanalysesarmonic, tabulatedvibrations ):
        self.calculation = calculation
        self.chemistrycode = chemistrycode
        self.tasks = []
        index = 0
        for task in tasks:
            task.theorylevel = theorylevels[index]
            task.molecularspecie = molecularspecies[index]
            task.electonicstate = electonicstates[index]
            task.dipolemoment = dipolemoments[index]
            task.geometry = geometries[index]
            task.rotationalconstant = rotationalconstants[index]
            task.vibrationalanalysesarmonic = vibrationalanalysesarmonic[index]
            task.tabulatedvibrations = []
            if tabulatedvibrations:
                for tab in tabulatedvibrations[index]:
                    task.tabulatedvibrations.append(tab)
            task.elements = []
            
            for i in range(len(elements[index])):
                task.elements.append({})
                task.elements[-1]["element"] = elements[index][i]
                task.elements[-1]["basisset"] = basissets[index][i]
            self.tasks.append(task)
            index += 1
            
class myElement():
    def __init__(self, element, basisset):
        self.element = element
        self.basisset = basisset
        
    def update_elements_by_database(self, update_elements = True, update_basisset = True, create_basisset = True):
        if update_elements:
            try:
                if not self.element.pk:
                    self.element = Elements.objects.filter(atomic_number = self.element.atomic_number, atomic_mass = round(self.element.atomic_mass))[0]
            except:
                pass
        if update_basisset:
            try:
                if not self.basisset.pk:
                    self.basisset = BasisSets.alive_objects.filter(name = self.basisset.name).order_by("-qual_index")[0]
            except:
                #create it
                if create_basisset:
                    if str(self.basisset.name).upper() in nwchem_task_basis_sets:
                        self.basisset = BasisSets(name = str(self.basisset.name).upper(), description = nwchem_task_basis_sets[str(self.basisset.name).upper()])
                        self.basisset.now()
                    self.basisset.save()