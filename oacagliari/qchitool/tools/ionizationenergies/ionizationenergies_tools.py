'''
Created on 07/mag/2012

@author: asaba
'''
from ...models import TheoryLevels, Xcclasses, Geometryclasses, Geometries, ElectronicStates, ElectronicTransitions, Tasks, IonisationEnergies, IonisationTypes
from operator import itemgetter
from ..qchitool_tools import RoundingError
from ... import performancetest as pt

def build_ionizationenergies(electronicstate):
    
    #Group by Inchi without charge
    es_by_inchi = {} #{inchi_without_charge: [Electronicstate, Electronicstate]}
    ionisationenergies_list = []
    if electronicstate:
        es_by_inchi = group_by_inchi(ElectronicStates.alive_objects.filter(pk=electronicstate.pk))
    else:
        es_by_inchi = group_by_inchi(ElectronicStates.alive_objects.all())
    #Group by Basis sets
    for esg_group_by_inchi in es_by_inchi:
        by_basisset = group_by_basisset(es_by_inchi[esg_group_by_inchi]) #{Primarykey Electronicstate : [Primarykey Electronicstate, Primarykey Electronicstate]}
        #Group by Theory Level
        for index in by_basisset:
            by_thlevel = group_by_thlevel(by_basisset[index])
            #Group by Geometry
            for th_pk in by_thlevel:
                #by_geometry = group_by_geometryclass(by_thlevel[th_pk])
                #for geom_pk in by_geometry:
                    #current_geometry = by_geometry[geom_pk]
                        by_energy = group_by_energy(by_thlevel[th_pk])
                        for energy_es in by_energy:
                            current_es_1 = by_energy[energy_es][0] #get olny the first one electronicstate for each energy group
                            if current_es_1.is_minimum:
                                current_es_1_min = True
                            else:
                                current_es_1_min = False
                            for energy_es2 in by_energy:
                                current_es_2 = by_energy[energy_es2][0] #get olny the first one electronicstate for each energy group
                                if current_es_2.is_minimum:
                                    current_es_2_min = True
                                else:
                                    current_es_2_min = False
                                if current_es_1 != current_es_2: #is not the same electronic state
                                    iontype = "TH"
                                    if current_es_1_min and current_es_2_min:
                                        #states are minimum (Adiabatic)
                                        iontype += "AD"
                                        if current_es_1.species.charge < current_es_2.species.charge:
                                            #ionization energy
                                            iontype += "ION"
                                            ionisationenergies_list.append(IonisationEnergies(start_state = current_es_1, ion_state = current_es_2, iontype = IonisationTypes.objects.get(description=iontype), energy = current_es_1.total_energy - current_es_2.total_energy))
                                        elif current_es_1.species.charge > current_es_2.species.charge:
                                            #electron affinity
                                            iontype += "AFF"
                                            ionisationenergies_list.append(IonisationEnergies(start_state = current_es_1, ion_state = current_es_2, iontype = IonisationTypes.objects.get(description=iontype), energy = current_es_1.total_energy - current_es_2.total_energy))
        
                                    elif current_es_1_min:
                                        #same geometry class
                                        if current_es_1.geom.geometryclass.pk == current_es_2.geom.geometryclass.pk:
                                            #one state only is minimum (vertical)
                                            iontype += "VE"
                                            if current_es_1.species.charge < current_es_2.species.charge:
                                                #ionization energy
                                                iontype += "ION"
                                                ionisationenergies_list.append(IonisationEnergies(start_state = current_es_1, ion_state = current_es_2, iontype = IonisationTypes.objects.get(description=iontype), energy = current_es_1.total_energy - current_es_2.total_energy))
                                            elif current_es_1.species.charge > current_es_2.species.charge:
                                                #electron affinity
                                                iontype += "AFF"
                                                ionisationenergies_list.append(IonisationEnergies(start_state = current_es_1, ion_state = current_es_2, iontype = IonisationTypes.objects.get(description=iontype), energy =  current_es_1.total_energy - current_es_2.total_energy))
            #electronicstates_group = {} #{primarykey of Xcclass: [Primarykey Electronicstate, Primarykey Electronicstate]}

    for newionisationenergy in ionisationenergies_list:
        if len(IonisationEnergies.alive_objects.filter(start_state = newionisationenergy.start_state, ion_state = newionisationenergy.ion_state, energy = newionisationenergy.energy)) <= 0:
            #if newionisationenergy not exist create and save it
            newionisationenergy.now()
            newionisationenergy.save()

def group_by_inchi(ess):
    #this function groups the electronicstates by inchi
    es_by_inchi = {} #{inchi_without_charge: [Electronicstate, Electronicstate]}
    tmp_dic= {}
    for es in ess:
        current_es_inchi = es.species.inchi_without_charge()
        if current_es_inchi in tmp_dic:
            tmp_dic[current_es_inchi].append(es.pk)
        else:
            tmp_dic[current_es_inchi] = [es.pk]
    for es_inchi in tmp_dic:
        es_by_inchi[es_inchi] = ElectronicStates.alive_objects.filter(pk__in=tmp_dic[es_inchi])
    return es_by_inchi

def group_by_basisset(ess):
    #this function groups the electronicstates primary key by basissets
    by_basisset = {} #{Primarykey Electronicstate : [Primarykey Electronicstate, Primarykey Electronicstate]}
    indexes = range(len(ess))
    for i in indexes:
        current_es_i = ess[i]
        by_basisset[i] = [current_es_i.pk]
        indexes.remove(i)
        for j in indexes:
            current_es_j = ess[j]
            if current_es_i.equal_basis_set(current_es_j.return_basissets_for_atoms()):
                by_basisset[i].append(current_es_j.pk)
                indexes.remove(j)
    return by_basisset

def group_by_thlevel(ess):
    #this function groups the electronicstates by theory levels
    by_thlevel = {}
    for es_pk in ess:
        es_thlevel = ElectronicStates.alive_objects.get(pk=es_pk).task.thlevel.pk
        if es_thlevel in by_thlevel:
            by_thlevel[es_thlevel].append(es_pk)
        else:
            by_thlevel[es_thlevel] = [es_pk]
    return by_thlevel

def group_by_energy(ess):
    #this function groups the electronicstates by energy
    #return a dictionary of {energy1:[electronicstate1, electronicstate2, ...], energy2:[electronicstate3, electronicstate4, ...],...}
    by_energy = {}
    rerror = RoundingError(absoluterror = 1E-2)
    found = None
    for es_pk in ess:
        es =  ElectronicStates.alive_objects.get(pk=es_pk)
        found = None
        for energy in by_energy:
            if len(by_energy[energy]) > 0:
                result, error = rerror.isabsoluteequal(float(energy), float(es.total_energy))
                if result:
                    if not found:
                        by_energy[energy].append(es)
                        found = energy
                    else:
                        by_energy[found] += by_energy[energy]
                        by_energy[energy] = []
        if not found:
            by_energy[es.total_energy] = [es]
    by_energy = dict((k, v) for k, v in by_energy.iteritems() if len(v) > 0)
    return by_energy

#def group_by_geometryclass(ess):
#    #this function groups the electroncstates by geometry classes
#    by_geometryclass = {}
#    for es_pk in ess:
#        es_thlevel = ElectronicStates.alive_objects.get(pk=es_pk).geom.geometryclass.pk
#        if es_thlevel in by_geometryclass:
#            by_geometryclass[es_thlevel].append(es_pk)
#        else:
#            by_geometryclass[es_thlevel] = [es_pk]
#    return by_geometryclass

