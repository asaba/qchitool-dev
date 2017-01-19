'''
Created on 12/mar/2012

@author: asaba
'''
import numpy as np
from ... import performancetest as pt

class myAtom:
    def __init__(self, molecule, index):
        self.mass = molecule.masses[index]
        self.atomicnum = molecule.atomicnum[index]
        self.coords = molecule.coords[index]
        
    def similar(self, myatom2, roundingerror):
        if self.sameIsotope(myatom2, roundingerror):
            return roundingerror.coord3dequal(self.coords, myatom2.coords)
        else:
            return False, 0.0
        
    def sameIsotope(self, myatom2, roundingerror):
        return self.atomicnum == myatom2.atomicnum and roundingerror.isabsoluteequal(self.mass, myatom2.mass)[0]
            
    def linearindipendent(self, atom2, roundingerror):
        crossvector = np.cross(self.coords, atom2.coords)
        error = np.linalg.norm(np.inner(self.coords,atom2.coords))
        n = np.linalg.norm(crossvector)
        if roundingerror.isrelativeequal(n, 0)[0]:
            return False, 0.0
        else:
            return n, error
        
    def rotationalerror(self, atom2):
        atom1vector = self.coords.copy()
        atom2vector = atom2.coords.copy()
        delta = atom1vector - atom2vector
        norm = np.linalg.norm((atom1vector + atom2vector)/2)
        vector2 = ((atom1vector + atom2vector)/2)/pow(norm, 2)
        crossvector = np.cross(delta, vector2)
        return [crossvector, norm]
    
    def distance(self):
        #distance from center of mass
        return np.linalg.norm(self.coords)
