'''
Created on 28/feb/2012

@author: asaba
'''


import pybel
import tempfile 
import os
#setting for import matplotlib

#import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import numpy as np



from ..models import Geometries, Geometryclasses, MolecularSpecies
from ..tools import modeltools
from ..tools.qchitool_tools import RoundingError
from ..import performancetest as pt

from geometry.atom import myAtom
from geometry.molecule import myMolecule
#from scipy.cluster.hierarchy import maxdists
from oacagliari.qchitool.models import ElectronicStates


def comparevectors(vector1, vector2, roundingerror, relative=True):
    #This function compare two sorted eigenvalues
    result = True
    
    error = 0.0
    maxvalue = 0.0
    for index in range(len(vector1)):
        if relative:
            tmpresult, tmperror = roundingerror.isrelativeequal(vector1[index], vector2[index])
        else:
            tmpresult, tmperror = roundingerror.isabsoluteequal(vector1[index], vector2[index])
        result = result and tmpresult
        error += tmperror
    return result, error


def test(a, b, debug=False, findmaxerrorvalues=False):
    print a + " "  + b
    
    mol = pybel.readfile("cml", a).next()
    myMol1 = myMolecule(mol)
    mol = pybel.readfile("cml", b).next()
    myMol2 = myMolecule(mol)
    if findmaxerrorvalues:
        relativeerror_b = 0.001 #* 20
        absoluterror_b = 0.0001 #* 2.5
        for i in range(1, 5000):
            relativeerror = relativeerror_b * (0.5 * i) 
            print "Relative Error = %s"%(str(relativeerror))
            for j in range(1,5):
                absoluterror = absoluterror_b * (0.5 * j)
                print "Absolute Error = %s"%(str(absoluterror))
                rounderr= RoundingError(absoluterror = absoluterror, relativeerror = relativeerror)
                if similargeometry(myMol1, myMol2, rounderr, debug)["isgeometrysimilar"]:
                    print "Relative Error = %s"%(str(relativeerror))
                    print "Absolute Error = %s"%(str(absoluterror))
                    break
    else:
        rounderr= RoundingError()
        similargeometry(myMol1, myMol2, rounderr, debug)
        
def rotate(filecml, filexyzrotated, rotation_matrix, traspose=False):
    if not rotation_matrix:
        rotation_matrix = np.array([[ -1.25016671e-07,   1.00000000e+00,   7.38595472e-14],
                                    [ -1.00000000e+00,  -1.25016671e-07,  -5.71157281e-08],
                                    [ -5.71157281e-08,  -8.09825149e-14,   1.00000000e+00]] )
        
    mol = pybel.readfile("cml", filecml).next()
    myMol1 = myMolecule(mol)
    if traspose:
        myMol1.rotate(rotation_matrix.T)
    else:
        myMol1.rotate(rotation_matrix)
    myMol1.exportxyz(filexyzrotated)

def returnOrientedmyMolbyString(mol1_string):
        myMol1 = buildmyMolbystring(mol1_string, None)
        return returnOrientedmyMol(myMol1)
    
def returnOrientedmyMol(myMol1):
        I1, eigenvaluesM1, eigenvectorsM1 = myMol1.getmomentofinertia()
        #sort eigenvalues in molecule 1 (Ascending order) 
        perm = np.argsort(eigenvaluesM1)   
        eigenvaluesM1 = eigenvaluesM1[perm]
        eigenvectorsM1 = eigenvectorsM1[:, perm]
        rotational_matrix_M1 = myMol1.rotate(eigenvectorsM1.T)
        return myMol1
        
def similargeometrybystrings(mol1_string, mol2_string, id_g1= None, id_g2 = None, rounderr=None, tmpdirectory=""):
    myMol1 = buildmyMolbystring(mol1_string, rounderr)
    myMol2 = buildmyMolbystring(mol2_string, rounderr)
    if not rounderr:
        rounderr= RoundingError()
    return similargeometry(myMol1, myMol2, rounderr, debug = True)
        

def similargeometry(myMol1, myMol2, rounderr, debug=False):
    #test getMomentOfInertia
    #add check for 1 and 2 atoms
    result = ""
    isgeometrysimilar = tmpisgeometrysimilar = False
    similarityerror = tmpsimilarityerror = 0.0
    rotational_matrix_M1 = tmprotational_matrix_M1 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    rotational_matrix_M2 = tmprotational_matrix_M2 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    simmetryvector = tmpsimmetryvector = np.array([1, 1, 1])
    rotationalerror = tmprotationalerror = np.array([])
    verror = tmpverror = np.array([])
    tmperror = 0.0
    maxerror = tmpmaxerror = 0.0
    correspondigatoms = tmpcorrespondigatoms = {}
    if myMol1.atomnumber == myMol2.atomnumber:
        #the number of atom in molecule 1 are equal to number of atom in molecule 2
        I1, eigenvaluesM1, eigenvectorsM1 = myMol1.getmomentofinertia()
        I2, eigenvaluesM2, eigenvectorsM2 = myMol2.getmomentofinertia()
        if debug:
            result += printoutput("Eigenvalues", eigenvaluesM1, eigenvaluesM2)
        
        #sort eigenvalues in molecule 1 (Ascending order) 
        perm = np.argsort(eigenvaluesM1)  
        eigenvaluesM1 = eigenvaluesM1[perm]
        eigenvectorsM1 = eigenvectorsM1[:, perm]
    
        #sort eigenvalues in molecule 2 (Ascending order) 
        perm = np.argsort(eigenvaluesM2)  
        eigenvaluesM2 = eigenvaluesM2[perm]
        eigenvectorsM2 = eigenvectorsM2[:, perm]
    
        if debug:
            #printoutput("Center of Mass" , myMol1.center_of_mass(), myMol2.center_of_mass())
            #printoutput("Eigenvalues", eigenvaluesM1, eigenvaluesM2)
            #printoutput("Eigenvectors", eigenvectorsM1, eigenvectorsM2)
            pass
        if comparevectors(eigenvaluesM1, eigenvaluesM2, rounderr)[0]:
            #eigenvectors of molecule 1 are equal to eigenvectors of molecule 2
            if myMol1.compatibledistance(myMol2, rounderr):
                error = 0.0
                #count equal eigenvalues and save their indices in listequaleigenvalues1
                countequaleigenvalues1, listequaleigenvalues1 =  rounderr.countmaxequal(eigenvaluesM1, relative=True)
                countequaleigenvalues2, listequaleigenvalues2 =  rounderr.countmaxequal(eigenvaluesM2, relative=True)
                #aggiungere anmche il controllo della seconda
                if countequaleigenvalues1 > 1:
                    #more then one eigenvalues equal in molecule 1
                    if countequaleigenvalues1 == 3:
                        #exactly 3 eigenvalues equal in molecule 1
                        if debug:
                            result += "three equal Eigenvalues\n"
                        #sort molecule by distance from center of mass
                        distance1 = myMol1.sortbydistance()
                        distance2 = myMol2.sortbydistance()
                        compareresult, error = comparevectors(distance1, distance2, rounderr)
                        if compareresult:
                            
                            tmpMol1 = myMol1.copy()
                            tmpMol2 = myMol2.copy()
                            molecules_equal = False
                            indexAtom1_M1, indexAtom2_M1, n1 = tmpMol1.returnfirstcloseratom(rounderr)
                            #atom 1 in molecule 1
                            #rotate molecule 1 by the base atom1, atom2, atom1^atom2 (atom1, atom2 in molecule 1)
                            rotational_matrix_M1 = tmpMol1.rotatebytwoatoms(indexAtom1_M1, indexAtom2_M1)
                            #myMol1.exportxyz("mol1rotated")
                            
                            for indexAtom1_M2 in range(tmpMol2.atomnumber-1):
                                #atom 1 in molecule 2
                                if not rounderr.isabsoluteequal(distance1[indexAtom1_M1], distance2[indexAtom1_M2])[0]:
                                    continue
                                myAtomM2_1 = myAtom(tmpMol2, indexAtom1_M2)
                                #confronto su tutte
                                for indexAtom2_M2 in range(indexAtom1_M2 + 1, tmpMol2.atomnumber):
                                    if not rounderr.isabsoluteequal(distance1[indexAtom2_M1], distance2[indexAtom2_M2])[0]:
                                        continue
                                    #atom 2 in molecule 2
                                    myAtomM2_2 = myAtom(tmpMol2, indexAtom2_M2)
                                    n2, coefficienterror = myAtomM2_1.linearindipendent(myAtomM2_2, rounderr)
                                    if not n2:
                                        #atom 1 and atom 2 in molecule 2 are linear dependent
                                        continue
                                    else:
                                        #if not rounderr.isabsolute(coefficienterror1, coefficienterror2, absolutecoeffient=1E-05)[0]:
                                        if not rounderr.isrelativeequal(n1, n2, coefficienterror)[0]:
                                            #The angle between vectors atom1 and atom2 in molecule 1 are different from the angle between atom1 and atom2 in molecule 2
                                            continue
                                    #rotate molecule 2 by the base atom1, atom2, atom1^atom2 (atom1, atom2, in molecule 2)
                                    #tmpMol2.exportxyz("mol2NotRotated")
                                    rotational_matrix_M2 = tmpMol2.rotatebytwoatoms(indexAtom1_M2, indexAtom2_M2)
                                    #tmpMol2.exportxyz("mol2rotated")
                                    #calculate difference between molecule 1 and molecule 2
                                    equal, tmperror, tmpsimmetryvector, tmprotationalerror, tmpverror, tmpmaxerror, tmpcorrespondigatoms = tmpMol1.isequal(tmpMol2, rounderr, True, True)
                                    if equal:
                                        molecules_equal = True
                                        break
                                    else:
                                        if tmperror > 0:
                                            if debug:
                                                result += "Error= %s Mol1: Atom %s, Atom %s - Mol2: Atom %s, Atom %s\n"%(str(error), indexAtom1_M1, indexAtom2_M1, indexAtom1_M2, indexAtom2_M2) 
                                    #raw_input("Press Enter to continue...")
                                if molecules_equal:
                                    break
                                else:
                                    #molecule 2 in initial position
                                    tmpMol2 = myMol2.copy()
    
                            if molecules_equal:
                                if debug:
                                    result += printoutputsimilarmolecule(tmperror, tmpmaxerror, tmpsimmetryvector, tmprotationalerror, tmpverror, distance1[-1], rotational_matrix_M1, rotational_matrix_M2)
                            else:
                                if debug:
                                    result +=  "Molecules are different Error= %s\n"%(str(error))
                            isgeometrysimilar = molecules_equal
                            similarityerror = tmperror
                            simmetryvector = tmpsimmetryvector
                            rotationalerror = tmprotationalerror
                            verror = tmpverror
                            maxerror = tmpmaxerror
                            correspondigatoms = tmpcorrespondigatoms
                        else:
                            if debug:
                                result +=  "Molecules are different Error= %s\n"%(str(error))  
                            similarityerror = error
                    else:
                        #two eigenvalue equal in molecule 1
                        if debug:
                            result +=  "two equal Eigenvalues\n"
                        #rotate molecule 1 and molecule 2 by their eigenvectors matrix
                        rotational_matrix_M1 = myMol1.rotate(eigenvectorsM1.T)
                        rotational_matrix_M2 = myMol2.rotate(eigenvectorsM2.T)
                        
                        eigenvaluesM1 = myMol1.eigenvalues()
                        eigenvaluesM2 = myMol2.eigenvalues()
                        
                        perm = np.argsort(eigenvaluesM1)  # sort in descending order
                        eigenvaluesM1 = eigenvaluesM1[perm]
                        eigenvectorsM1 = eigenvectorsM1[:, perm]
                        
                        perm = np.argsort(eigenvaluesM2)  # sort in descending order
                        eigenvaluesM2 = eigenvaluesM2[perm]
                        eigenvectorsM2 = eigenvectorsM2[:, perm]
                        
                        countequaleigenvalues1, listequaleigenvalues1 =  rounderr.countmaxequal(eigenvaluesM1, relative=True)
                        countequaleigenvalues2, listequaleigenvalues2 =  rounderr.countmaxequal(eigenvaluesM2, relative=True)
                        #printoutput("Eigenvalues", eigenvaluesM1, eigenvaluesM2)
                        #in line atoms
    
                        if len(myMol1.isinline(rounderr)) > 0: #atoms are in line
                            equal, tmperror, tmpsimmetryvector, tmprotationalerror, tmpverror, tmpmaxerror, tmpcorrespondigatoms = myMol1.isequal(myMol2, rounderr, calculate_rotational_vectors = True)
                            if equal:
                                if debug:
                                    result += printoutputsimilarmolecule(tmperror, tmpmaxerror, tmpsimmetryvector, tmprotationalerror, tmpverror, distance1[-1], rotational_matrix_M1, rotational_matrix_M2)
                            else:
                                if debug:
                                    result +=  "Molecules are different Error= %s\n"%(str(tmperror))
                            isgeometrysimilar = equal
                            similarityerror = tmperror
                            simmetryvector = tmpsimmetryvector
                            rotationalerror = tmprotationalerror
                            verror = tmpverror
                            maxerror = tmpmaxerror
                            correspondigatoms = tmpcorrespondigatoms
                        else:
                            distance1 = myMol1.sortbydistance()
                            distance2 = myMol2.sortbydistance()
                            molecules_equal = False
                            index1 = myMol1.returnfirstcloseratom(rounderr)[0]
                            if  myMol1.coordinaxes(index1, listequaleigenvalues1, rounderr): #atom in symmetry axes
                                raw_input("ERROR - Press Enter to continue...")
                            rotational_matrix_M1 = myMol1.rotatebyaxesandatom(listequaleigenvalues1, index1)
                            #myMol1.exportxyz("mol1rotated")
                            for index2 in range(0, myMol2.atomnumber):
                                if  myMol2.coordinaxes(index2, listequaleigenvalues2, rounderr): #atom in symmetry axes
                                    continue
                                if rounderr.isabsoluteequal(distance1[index1], distance2[index2])[0]: #
                                    rotational_matrix_M2 = myMol2.rotatebyaxesandatom(listequaleigenvalues2, index2)
                                    #calculate difference between molecule 1 and molecule 2
                                    #myMol2.exportxyz("mol2rotated")
                                    equal, tmperror, tmpsimmetryvector, tmprotationalerror, tmpverror, tmpmaxerror, tmpcorrespondigatoms = myMol1.isequal(myMol2, rounderr, True, calculate_rotational_vectors = True)
                                    #printoutput("Eigenvalues", myMol1.eigenvalues(), myMol2.eigenvalues())
                                    #printoutput("Eigenvectors", myMol1.eigenvectors(), myMol2.eigenvectors())
                                    if equal:
                                        molecules_equal = True
                                        break
                                    else:
                                        if tmperror > 0:
                                            if debug:
                                                result +=  "Error=\n" + str(tmperror)
                                        else:
                                            if debug:
                                                result +=  str(index1) + str(index2) + "\n" 
                            if molecules_equal:
                                if debug:
                                    result += printoutputsimilarmolecule(tmperror, tmpmaxerror, tmpsimmetryvector, tmprotationalerror, tmpverror, distance1[-1], rotational_matrix_M1, rotational_matrix_M2)
                            else:
                                if debug:
                                    result +=  "Molecules are different Error= %s\n"%(str(tmperror))
                            isgeometrysimilar = molecules_equal
                            similarityerror = tmperror
                            simmetryvector = tmpsimmetryvector
                            rotationalerror = tmprotationalerror
                            verror = tmpverror
                            maxerror = tmpmaxerror
                            correspondigatoms = tmpcorrespondigatoms
                            
                else:
                    #three eigenvalues different in molecule 1
                    if debug:
                        result +=  "Three different Eigenvalues\n"
                    #rotate molecule 1 and molecule 2 by their eigenvectors matrix
                    rotational_matrix_M1 = myMol1.rotate(eigenvectorsM1.T)
                    rotational_matrix_M2 = myMol2.rotate(eigenvectorsM2.T)
                    #calculate difference between molecule 1 and molecule 2
                    equal, tmperror, tmpsimmetryvector, tmprotationalerror, tmpverror, tmpmaxerror, tmpcorrespondigatoms = myMol1.isequal(myMol2, rounderr, calculate_rotational_vectors = True)
                    if equal:
                        if debug:
                            distance1 = myMol1.copy().sortbydistance()
                            result += printoutputsimilarmolecule(error, tmpmaxerror, tmpsimmetryvector, tmprotationalerror, tmpverror, distance1[-1], rotational_matrix_M1, rotational_matrix_M2)
                        isgeometrysimilar = True
                        similarityerror = tmperror
                        simmetryvector = tmpsimmetryvector
                        rotationalerror = tmprotationalerror
                        verror = tmpverror
                        maxerror = tmpmaxerror
                        correspondigatoms = tmpcorrespondigatoms
                    else:
                        if debug:
                            result += "Molecules are different Error= %s\n"%(str(tmperror))  
            else:
                if debug:
                    result +=  "Molecules are different (Distances not compatible)\n"
        else:
            if debug:
                result +=  "Molecules are different (Different eigenvalues)\n"
    else:
        if debug:
            result +=  "Molecules are different (Different number of atoms)\n"
    
    return {"isgeometrysimilar" : isgeometrysimilar, 
            "similarityerror" : similarityerror, 
            "simmetryvector" : simmetryvector, 
            "rotationalerror" : rotationalerror, 
            "verror" : verror, 
            "maxerror" : maxerror, 
            "result" : result, 
            "myMol1" : myMol1, 
            "myMol2" : myMol2, 
            "rotational_matrix_M1" : rotational_matrix_M1, 
            "rotational_matrix_M2" : rotational_matrix_M2,
            "correspondigatoms" : correspondigatoms} 
        
    
def printoutputsimilarmolecule(error, maxerror, simmetryvector, rotationalerror, verror, maxdistance, rotational_matrix_M1, rotational_matrix_M2):
    result = "Molecules are similar\n"
    result += "Average Error = %s; Max Error = %s; Simmetry vector = %s \n"%(str(error), str(maxerror), str(simmetryvector))
    result +=  "Rotational Error: %s\n"%(rotationalerror)
    result +=  "Rotational Error vector: %s\n"%(verror)
    result +=  "Max distance from center of mass = %s\n"%(str(maxdistance))
    result +=  "Molecules Rotational Matrix: \n"
    result +=  "Molecule 1: %s \n"%(rotational_matrix_M1)
    result +=  "Molecule 2: %s \n"%(rotational_matrix_M2)
    return result
    
class geometryequivalenceresult():
    def __init__(self, equivalent_result_dict):
        self.equivalent = equivalent_result_dict["isgeometrysimilar"]
        self.id_geometry1 = 0
        self.molecule1 = equivalent_result_dict["myMol1"]
        self.rotationalmatrix1 = equivalent_result_dict["rotational_matrix_M1"]
        self.id_geometry2 = 0
        self.molecule2 = equivalent_result_dict["myMol2"]
        self.rotationalmatrix2 = equivalent_result_dict["rotational_matrix_M2"]
        self.error = equivalent_result_dict["similarityerror"]
        self.simmetryvector = equivalent_result_dict["simmetryvector"]
        self.rotationalerror = equivalent_result_dict["rotationalerror"]
        self.vectorerror = equivalent_result_dict["verror"]
        self.maxerror = equivalent_result_dict["maxerror"]
        self.stringresult = equivalent_result_dict["result"]
        self.correspondigatoms = equivalent_result_dict["correspondigatoms"]
        
    
def buildequivalenceclasses2(tmpdirectory, limit = 0, geometryidlist =[]):
    if limit > 0:
        geoms = Geometries.alive_objects.all()[:limit]
    else:
        if len(geometryidlist) > 0:
            geoms = Geometries.alive_objects.filter(pk__in = geometryidlist)
        else:
            geoms = Geometries.alive_objects.all()
    classes1 = {}
    ctotal = {}
    classes2 = []
    maxerror = 0.0
    result = ""
    for index1 in range(geoms.count()):
        #print str(geoms.count() - index1)
        classes1[geoms[index1].pk]=[]
        if geoms[index1].pk in classes2:
            #it is in one class
            continue
        rerr = RoundingError(absoluterror = 0.001, relativeerror = 0.1)
        #atomsnumber1 = geoms[index1].states.countatoms()
        #print geoms[index1].geometry
        g1 = str(geoms[index1].geometry)
        formula = geoms[index1].electronicstates_set.all()[0].species.formula.replace("+", "").replace("-", "")
        result += "Geometry ID = %s, Formula: %s\n"%(str(geoms[index1].pk), str(formula))
        try:
            file = geoms[index1].electronicstates_set.all()[0].task.calc.comments
            result += "Calculus File: %s\n"%(file)
        except:
            result += "Calculus File: ERROR ON DB\n"
        #search MolecularSpecies with same stoichiometric formula
        species = MolecularSpecies.alive_objects.filter(formula__istartswith=formula.replace("+", "").replace("-", ""))
        geoms2 = []
        for s in species:
            for s2 in s.electronicstates_set.all():
                #select Geometries with same stoichiometric formula
                geoms2.append(s2.geom)
        for index2 in range(len(geoms2)):
            if geoms[index1].pk == geoms2[index2].pk:
                continue
            formula = geoms2[index2].electronicstates_set.all()[0].species.formula.replace("+", "").replace("-", "")
            
            result +=  "Checking...\n"
            result +=  "Geometry ID = %s, Formula: %s\n"%(str(geoms2[index2].pk), str(formula))
            try:
                file = geoms2[index2].electronicstates_set.all()[0].task.calc.comments
                result +=  "Calculus File: %s\n"%(file)
            except:
                result += "Calculus File: ERROR ON DB\n"
            tmpresult = geometryequivalenceresult(similargeometrybystrings(g1, geoms2[index2].geometry, rerr, tmpdirectory))
            result += tmpresult.stringresult
            if tmpresult.equivalent:
                tmpresult.id_geometry1 = geoms[index1].pk
                tmpresult.id_geometry2 = geoms2[index2].pk
                if not geoms[index1].pk in classes1:
                    classes1[geoms[index1].pk]=[]
                classes1[geoms[index1].pk].append(tmpresult)
                classes2.append(geoms2[index2].pk)
                result +=  " similar to %s (error= %s)\n"%(str(geoms2[index2].pk), str(tmpresult.error))
                maxerror = max(maxerror, tmpresult.error)
        
    ctotal = {}
    collisions = 0
    for c in  classes1.keys():
        if len(classes1[c]) in ctotal:
            ctotal[len(classes1[c])] += 1
        else:
            ctotal[len(classes1[c])] = 1
        for c2 in classes1.keys():
            if c2<>c:
                for v in classes1[c]:
                    if v in [k.id_geometry2 for k in classes1[c2]]:
                        collisions +=1
    for c in ctotal.keys():
        result +=  "Class Lenght %d (classes %d)\n"%(c, ctotal[c])
    result +=  "Collisions: %d\n"%(collisions)
    result +=  "MaxError in similar geometries: %.15f\n"%(maxerror)
    return result, classes1

def printoutput(string, mol1text, mol2text):
    result =  string + "\n"
    result +=  "Molecule 1 \n"
    result +=  str(mol1text) + "\n"
    result +=  "Molecule 2 \n" 
    result +=  str(mol2text) + "\n"
    return result

def buildavaragegeometry(geometry_simmetryresult_list):
    avgmol = avgmolresult = None 
    for g in geometry_simmetryresult_list:
        if not avgmol:
            avgmol = g.molecule1.copy()
            avgmol.rotate(g.rotationalmatrix1)
        g.molecule2.reflect(g.simmetryvector)    
        g.molecule2.rotate(g.rotationalmatrix2)
        avgmolresult = avgmol.averagemolecule(g.molecule2, g.correspondigatoms)
    return avgmolresult

def buildavaragegeometry_byID(geometry_ID_list):
    geometry_myMol_list = []
    for g_id in geometry_ID_list:
        geometry_myMol_list.append(buildmyMolbystring(Geometryclasses.objects.get(pk = g_id).geometry))
    return buildavaragegeometry(geometry_myMol_list)

def returngeometryclassID(geometry):
    #This function search the geometry class of geometry
    #if not this function create a new one geometry class based on geometry
    
    geometry_classes = Geometryclasses.objects.all()
    geometryequivalenceresultlist = []
    result = None
    for gc in geometry_classes:
        if gc.pk <> modeltools.dummyobject.geometryclass.pk:
            gre = geometryequivalenceresult(similargeometrybystrings(gc.geometry, geometry))
            if gre.equivalent:
                gre.id_geometry1 = gc.pk
                geometryequivalenceresultlist.append(gre)
        
    if len(geometryequivalenceresultlist) == 0:
        #geometry class not found: create it
        result = buildsingolargeometryclass(geometry)
    elif len(geometryequivalenceresultlist) == 1:
        #geometry class
        result = geometryequivalenceresultlist[0].id_geometry1
    else:
        #union of classes
        listofgeometryID = []
        for gre in geometryequivalenceresultlist:
            listofgeometryID += [g.id_geometry1 for g in gre] + [g.id_geometry2 for g in gre]
        #TO DO
        
    return result

def electronicstatebygeometryclass(electronicstate, geometryclass, charge, thlevel):
    #this function seach in the geometry class of electronistate.geom
    #different electronic states with same charge, task.teorylevel, basisset
    if electronicstate:
        return ElectronicStates.alive_objects.filter(geom__geometryclass = electronicstate.geom.geometryclass, species__charge = electronicstate.species.charge, task__thlevel = electronicstate.task.thlevel)
    else:
        return ElectronicStates.alive_objects.filter(geom__geometryclass = geometryclass, species__charge = charge, task__thlevel = thlevel)

            
def electronicstateisminimumbygeometry(electronicstate, geometryclass, charge, thlevel):
    #select electronic states by electronicstates geometryclass
    #if also one electronicstate is minimum return True
    result = False
    for es in electronicstatebygeometryclass(electronicstate, geometryclass, charge, thlevel):
        if es.is_minimum:
            result = True
            break
    return result


def buildmyMolbystring(geometrystring, rounderr = None):
        f1 = tempfile.NamedTemporaryFile(mode="w")
        f1.write(geometrystring)
        f1.seek(0)
        if not rounderr:
            rounderr= RoundingError()
        mol = pybel.readfile("cml", f1.name).next()
        f1.close()
        return myMolecule(mol)

def buildsingolargeometryclass(geometry):
        avaragegeometry = returnOrientedmyMolbyString(geometry)
        #save avaragegeometry in geometryclasses in avggeomcls
        avggeomcls = Geometryclasses(geometry = avaragegeometry.buildcmloutput())
        avggeomcls.calculate_md5()
        avggeomcls.save()
        return avggeomcls.pk
