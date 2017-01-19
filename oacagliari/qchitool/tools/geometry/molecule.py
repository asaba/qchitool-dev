'''
Created on 12/mar/2012

@author: asaba
'''

import numpy as np
import pybel
import tempfile
from atom import myAtom
from ... import performancetest as pt

simmetry = np.array([[1,1,1],
                     [1,1,-1],
                     [1,-1,1],
                     [1,-1,-1],
                     [-1, 1, 1],
                     [-1, 1, -1],
                     [-1, -1, 1],
                     [-1, -1, -1]])


class myMolecule:
    def __init__(self, pybel_molecule):
        #extract from molecule the element symbol, element mass, element coordinate
        coords = []
        self.atomicnum = np.array([])
        self.masses = np.array([])
        self.atomnumber = 0
        if pybel_molecule:
            for atom in pybel_molecule.atoms:
                coords.append(atom.coords)
                self.atomicnum = np.concatenate((self.atomicnum, [atom.atomicnum]))
                self.masses = np.concatenate((self.masses, [atom.exactmass]))
                self.atomnumber += 1
        self.coords = np.array(coords)
        
    def copy(self):
        new_molecule = myMolecule(None)
        new_molecule.coords = self.coords.copy()
        new_molecule.atomicnum = self.atomicnum.copy()
        new_molecule.masses = self.masses.copy()
        new_molecule.atomnumber = self.atomnumber
        return new_molecule
        
    def translate(self, translation_vector):
        for coord in self.coords:
            coord -= translation_vector
                
    def rotate(self, rotation_matrix):
        for index in range(self.atomnumber):
            self.coords[index] = np.dot(rotation_matrix, self.coords[index])
        return rotation_matrix
          
    def buildunixyzoutput(self):
        result =[]
        result.append("\n")
        result.append(str(self.atomnumber)+ "\n")
        i = 0
        for atomcoords in self.coords:
            result.append("%s  %s %s %s \n"%(str(self.atomicnum[i]), str(atomcoords[0]), str(atomcoords[1]), str(atomcoords[2])))
            i += 1
        return result
          
    def exportunixyz(self, suffix):
        f = open("geom_" + suffix + ".unixyz", "w")
        f.writelines(self.buildunixyzoutput())
        f.close()
        
    def buildcmloutput(self):
        f1 = tempfile.NamedTemporaryFile(mode="w")
        f1.writelines(self.buildunixyzoutput())
        f1.seek(0)
        mol = pybel.readfile("unixyz", f1.name).next()
        f1.close()
        result = mol.write(format="cml")
        return result
        
    def eigenvalues(self):
        a, tmpeigenvalues, b = self.getmomentofinertia()
        return tmpeigenvalues
    
    def eigenvectors(self):
        a, b, tmpeigenvectors = self.getmomentofinertia()
        return tmpeigenvectors
    
    def rotatebyaxesandatom(self, symmetryaxeslist, atom_index):
        atom_vector = self.coords[atom_index].copy()
        axis_vector = np.array([0,0,0])
        for axis in range(3):
            if not (axis in symmetryaxeslist):
                atom_vector[axis] = 0
                axis_vector[axis] = 1
        atom_vector = atom_vector/np.linalg.norm(atom_vector)
        ortogonal_vector = np.cross(atom_vector, axis_vector)
        rotation_matrix = np.array([atom_vector, ortogonal_vector, axis_vector])
        return self.rotate(rotation_matrix)
        
    def rotatebytwoatoms(self, atom_index1, atom_index2):
        atom_vector1 = self.coords[atom_index1].copy()
        atom_vector2 = self.coords[atom_index2].copy()
        vector1 = atom_vector1/np.linalg.norm(atom_vector1)
        scalar1 = np.vdot(atom_vector2, vector1)
        vector2 = atom_vector2 - scalar1 * vector1
        vector2 = vector2/np.linalg.norm(vector2)
        vector3 = np.cross(vector1, vector2)
        rotation_matrix = np.array([vector1, vector2, vector3])
        return self.rotate(rotation_matrix)
        
        
    def sortbydistance(self):
        distances = np.zeros(self.atomnumber, np.float64)
        for index in range(self.atomnumber):
            distances[index] = myAtom(self, index).distance()
        
        perm = np.argsort(distances)  # sort in ascending order
        self.coords = self.coords[perm]
        self.atomicnum = self.atomicnum[:, perm]
        self.masses = self.masses[:, perm]
        distances = distances[:, perm]
        return  distances
    
    def sortbycoords(self):
        perm = np.lexsort(tuple([self.coords[:,x] for x in range(3)]))# sort in descending order
        self.coords = self.coords[perm]
        self.atomicnum = self.atomicnum[:, perm]
        self.masses = self.masses[:, perm]
            
    def reflect(self, vector):
        self.coords = self.coords * vector
        
    def compatibledistance(self, molecule2, roundingerror):
        if self.atomnumber == molecule2.atomnumber:
            result = True
            skipindex = []
            for index1 in range(self.atomnumber):
                for index2 in range(molecule2.atomnumber):
                    if index2 in skipindex:
                        continue
                    myAtom1 = myAtom(self, index1)
                    myAtom2 = myAtom(molecule2, index2)
                    if myAtom1.sameIsotope(myAtom2, roundingerror) and roundingerror.isabsoluteequal(myAtom1.distance(), myAtom2.distance()):
                        skipindex.append(index2)
                        break
                if len(skipindex) <> index1 +1:
                    result = False
                    break
            return result
        else:
            return False
        
    def isequal(self, molecule, roundingerror, printcoordinatematrix=False, calculate_rotational_vectors = False):
        #perm = self.coords[:, 0].argsort() #sort in ascending order
        result = False
        #self.sortbycoords()
        simmetryvector = np.array([0,0,0])
        rotationalerror = []
        correspondigatoms = {} #dictionary with corresponding atom index in self and molecule
        i = 0
        for s in simmetry:
            i+=1
            tmpmolecule = molecule.copy()
            tmpmolecule.reflect(s)
            #tmpmolecule.sortbycoords()
            result = False
            totalerror = 0.0
            maxerror = 0.0
            skipatomindex = []
            correspondigatoms = {}
            for index in range(self.atomnumber):
                atom1 = myAtom(self, index)
                result, error, atomfoundindex = tmpmolecule.atom_in_molecule(atom1, roundingerror, skipatomindex)

                if result:
                    if calculate_rotational_vectors:
                        rotationalerror.append(atom1.rotationalerror(myAtom(tmpmolecule, atomfoundindex)))                    
                    skipatomindex.append(atomfoundindex)
                    correspondigatoms[index] =  atomfoundindex
                    maxerror = max(maxerror, error)
                    totalerror += error
                else:
                    break
            if result:
                if i>1:
                    pass
                simmetryvector = s
                break
        v = np.array([0.0,0.0,0.0])
        n = 0.0
        for vect in rotationalerror:
            #add check for 0
            v = v + vect[0]/vect[1]
            n = n + vect[1]
        rerror = v/n
        return result, totalerror/self.atomnumber, simmetryvector, rotationalerror, rerror, maxerror, correspondigatoms
    
    def isinline(self, roundingerror, I= None, values= None, vectors= None):
        #count zero eigenvalues
        if I is None:
            I, values, vectors = self.getmomentofinertia()
        result = []
        for index in range(len(values)):
            equal, error = roundingerror.isabsoluteequal(values[index], 0)
            if equal:
                result.append(index)
        return result 
    
    def coordinaxes(self, coord_index, axeslist, roundingerror):
        coord = self.coords[coord_index]
        result = False
        for axis in range(3):
            if not (axis in axeslist):
                result = result and roundingerror.isabsoluteequal(coord[axis], 0)[0]
        return result

    
    def getmomentofinertia(self):
        """
        Calculate and return the moment of inertia tensor for the current 
        geometry. If the coordinates are not at the center of mass,
        they are temporarily shifted there for the purposes of this calculation.
        """
  
        I = np.zeros((3,3), np.float64)
        centerOfMass = self.center_of_mass()
        for atom_index in range(self.atomnumber):
            mass = self.masses[atom_index]
            coord = self.coords[atom_index] - centerOfMass
            I[0,0] += mass * (coord[1] * coord[1] + coord[2] * coord[2])
            I[1,1] += mass * (coord[0] * coord[0] + coord[2] * coord[2])
            I[2,2] += mass * (coord[0] * coord[0] + coord[1] * coord[1])
            I[0,1] -= mass * coord[0] * coord[1]
            I[0,2] -= mass * coord[0] * coord[2]
            I[1,2] -= mass * coord[1] * coord[2]
        I[1,0] = I[0,1]
        I[2,0] = I[0,2]
        I[2,1] = I[1,2]
        values1, vectors1 = np.linalg.eig(I)
        return I, values1, vectors1
    
    def center_of_mass(self):
        cm = np.zeros(3, np.float64)
        m = 0.0
        for atom_index in range(self.atomnumber):
            m += self.masses[atom_index]
            cm += self.coords[atom_index] * self.masses[atom_index]
        cm /= m
        return cm
    
    def returnfirstcloseratom(self, roundingerror):
        #The atoms must be sorted by distance 
        #return the first closer atom to center of mass, not coinciding with the center of mass, and the second one linear independent
        result1 = -1
        result2 = -1
        n1 = None
        for index in range(self.atomnumber):
            if not roundingerror.coord3dequal(self.coords[index], np.array([0,0,0]))[0]:
                result1 = index
                myAtomM1 = myAtom(self, index)
                break
        if result1 > -1:
            for index in range(self.atomnumber):
                if index == result1:
                    continue
                if not roundingerror.coord3dequal(self.coords[index], np.array([0,0,0]))[0]:
                    myAtomM2 = myAtom(self, index)
                    n1 = myAtomM1.linearindipendent(myAtomM2, roundingerror)[0]
                    if not n1:
                        #atom 1 and atom 2 are linear dependent
                        continue
                    result2 = index
                    break
        return result1, result2, n1
    
    def atom_in_molecule(self, atom2, roundingerror, skipindex):
        result = False
        error = 0.0
        i = -1
        for index in range(self.atomnumber):
            if index in skipindex:
                continue
            atom1 = myAtom(self, index)
            result, error = atom1.similar(atom2, roundingerror)
            if result:
                i = index
                break
        return result, error, i   
    
    def averagemolecule(self, molecule, correspondingatoms):
        avgmol = self.copy()
        for i in range(self.atomnumber):
            for k in range(3):
                if i <> correspondingatoms[i]:
                    pass
                avgmol.coords[i][k] = np.average([self.coords[i][k], molecule.coords[correspondingatoms[i]][k]])
        return  avgmol