'''
Created on 07/giu/2012

@author: asaba
'''

from .. import models
from .. import performancetest as pt

        
def authorsexpire(a):
    for ba in models.BibliographyAuthors.objects.filter(authors = a):
        bibliographyauthorsexpire(ba)
    a.expire()


def editorsexpire(e):
    for be in models.BibliographyEditors.objects.filter(editors = e):
        bibliographyeditorsexpire(be)
    e.expire()
    
def geometriesexpire(g):
    for es in models.ElectronicStates.objects.filter(geom = g):
        electronicstatesexpire(es)
    g.expire()

def molecularspeciesexpire(ms):
    for es in models.ElectronicStates.objects.filter(species = ms):
        electronicstatesexpire(es)
#    for es in models.ElementSpecies.objects.filter(species = model_id):
#        elementspeciesexpire(es.pk)
    ms.expire()

def publishersexpire(p):
    for b in models.Bibliography.objects.filter(publisher = p):
        bibliographyexpire(b)
    p.expire()

def reftypeexpire(rt):
    for b in models.Bibliography.objects.filter(type = rt):
        bibliographyexpire(b)
    rt.expire()

def bibliographyexpire(b):
    for ba in models.BibliographyAuthors.objects.filter(bibliographies = b):
        bibliographyauthorsexpire(ba)
    for be in models.BibliographyEditors.objects.filter(bibliographies = b):
        bibliographyeditorsexpire(be)
    for tlb in models.TheoryLevelsBibliographies.objects.filter(bibliographies = b):
        theorylevelsbibliographiesexpire(tlb)
    for bsb in models.BasisSetsBibliographies.objects.filter(bibliographies = b):
        basissetsbibliographiesexpire(bsb)
    for ccb in models.ChemistryCodesBibliographies.objects.filter(bibliographies = b):
        chemistrycodesbibliographiesexpire(ccb)
    for vahb in models.VibrationalAnalysesHarmonicBibliographies.objects.filter(bibliography = b):
        vibrationalanalysesharmonicbibliographiesexpire(vahb)
    for vaab in models.VibrationalAnalysesAnarmonicBibliographies.objects.filter(bibliographies = b):
        vibrationalanalysesanarmonicbibliographiesexpire(vaab)
    for etb in models.ElectronicTransitionsBibliographies.objects.filter(bibliographies = b):
        electronictransitionsbibliographiesexpire(etb)
    for ieb in models.IonisationEnergiesBibliographies.objects.filter(bibliography = b):
        ionisationenergiesbibliographiesexpire(ieb)
    for esb in models.ElectronicStatesBibliographies.objects.filter(bibliography = b):
        electronicstatesbibliographiesexpire(esb)
    for rcb in models.RotationalConstantsBibliographies.objects.filter(bibliographies = b):
        rotationalconstantsbibliographiesexpire(rcb)
    for dmb in models.DipoleMomentsBibliographies.objects.filter(bibliographies = b):
        dipolemomentsbibliographiesexpire(dmb)
    for pb in models.PolarisabilitiesBibliographies.objects.filter(bibliographies = b):
        polarisabilitiesbibliographiesexpire(pb)
    for vdwb in models.VanDerWallsBibliographies.objects.filter(bibliographies = b):
        vanderwallsbibliographiesexpire(vdwb)
    b.expire()

def theorylevelsexpire(tl):
    for t in models.Tasks.objects.filter(thlevel = tl):
        tasksexpire(t)
    tl.expire()
        
def bibliographyauthorsexpire(model_id):
    models.BibliographyAuthors.objects.get(pk=model_id).expire()
        
def chemistrycodesexpire(cc):
    for ccb in models.ChemistryCodesBibliographies.objects.filter(chemistrycodes = cc):
        chemistrycodesbibliographiesexpire(ccb)
    cc.expire()



def calculationsexpire(c):
    for t in models.Tasks.objects.filter(calc = c):
        tasksexpire(t)
    c.expire()


def basissetsexpire(bs):
    for bsb in models.BasisSetsBibliographies.objects.filter(basissets = bs):
        basissetsbibliographiesexpire(bsb)
#    for esbs in models.ElementSpeciesBasisSet.objects.filter(basisset = model_id):
#        elementspeciesbasissetexpire(esbs.pk)
    bs.expire()


def tasksexpire(t):
    for es in  models.ElectronicStates.objects.filter(task = t):
        electronicstatesexpire(es)
    for vbaa in models.VibrationalAnalysesAnarmonic.objects.filter(task = t):
        vibrationalanalysesanarmonicexpire(vbaa)
    for vbah in models.VibrationalAnalysesHarmonic.objects.filter(task = t):
        vibrationalanalysesharmonicexpire(vbah)
    for dm in models.DipoleMoments.objects.filter(task = t):
        dipolemomentsexpire(dm)
    for et in models.ElectronicTransitions.objects.filter(task = t):
        electronictransitionsexpire(et)
    for p in models.Polarisabilities.objects.filter(task = t):
        polarisabilitiesexpire(p)
    for vdw in models.VanDerWaals.objects.filter(task = t):
        vanderwaalsexpire(vdw)
#    for esbs in models.ElementSpeciesBasisSet.objects.filter(task = model_id):
#        elementspeciesbasissetexpire(esbs.pk)
    t.expire()
            

def bibliographyeditorsexpire(be):
    be.expire()


def electronicstatesexpire(es):
    for vaa in models.VibrationalAnalysesAnarmonic.objects.filter(state = es):
        vibrationalanalysesanarmonicexpire(vaa)
    for vah in models.VibrationalAnalysesHarmonic.objects.filter(state = es):
        vibrationalanalysesharmonicexpire(vah)
    for dm in models.DipoleMoments.objects.filter(state = es):
        dipolemomentsexpire(dm)
    for et in models.ElectronicTransitions.objects.filter(low_state = es): 
        electronictransitionsexpire(et)
    for et in models.ElectronicTransitions.objects.filter(up_state = es):  
        electronictransitionsexpire(et)
    for ie in models.IonisationEnergies.objects.filter(start_state = es):   
        ionisationenergiesexpire(ie)
    for ie in models.IonisationEnergies.objects.filter(ion_state = es):  
        ionisationenergiesexpire(ie)
    for p in models.Polarisabilities.objects.filter(state = es):
        polarisabilitiesexpire(p)
    for rc in models.RotationalConstants.objects.filter(state = es):
        rotationalconstantsexpire(rc)
    for vdw in models.VanDerWaals.objects.filter(state_id_1 = es): 
        vanderwaalsexpire(vdw)
    for vdw in models.VanDerWaals.objects.filter(state_id_2 = es): 
        vanderwaalsexpire(vdw)
    for esb in models.ElectronicStatesBibliographies.objects.filter(electronicstates = es):
        electronicstatesbibliographiesexpire(esb)
    es.expire()

def vibrationalanalysesanarmonicexpire(vaa):
    for vaab in models.VibrationalAnalysesAnarmonicBibliographies.objects.filter(vibrationalanalysisanarmonic = vaa):
        vibrationalanalysesanarmonicbibliographiesexpire(vaab)
    vaa.expire()

def vibrationalanalysesharmonicexpire(vah):
    for vahb in models.VibrationalAnalysesHarmonicBibliographies.objects.filter(vibrationalanalysesarmonic = vah):
        vibrationalanalysesharmonicbibliographiesexpire(vahb)
    vah.expire()
            
def dipolemomentsexpire(dm):
    for dmb in models.DipoleMomentsBibliographies.objects.filter(dipolemoments = dm):
        dipolemomentsbibliographiesexpire(dmb)
    dm.expire()


def electronictransitionsexpire(et):
    for etb in models.ElectronicTransitionsBibliographies.objects.filter(electronictransitions = et):
        electronictransitionsbibliographiesexpire(etb)
    et.expire()
            

def ionisationenergiesexpire(ie):
    for ieb in models.IonisationEnergiesBibliographies.objects.filter(ionisationenergies = ie):
        ionisationenergiesbibliographiesexpire(ieb)
    ie.expire()

            
def polarisabilitiesexpire(p):
    for pb in models.PolarisabilitiesBibliographies.objects.filter(polarisabilities = p):
        polarisabilitiesbibliographiesexpire(pb)
    p.expire()
            
def rotationalconstantsexpire(rc):
    for rcb in models.RotationalConstantsBibliographies.objects.filter(rotationalcostants = rc):
        rotationalconstantsbibliographiesexpire(rcb)
    rc.expire()

def vanderwaalsexpire(vdw):
    for vdwb in models.VanDerWallsBibliographies.objects.filter(vanderwalls = vdw):
        vanderwallsbibliographiesexpire(vdwb)
    vdw.expire()

def theorylevelsbibliographiesexpire(tlb):
    tlb.expire()
        
def basissetsbibliographiesexpire(bsb):
    bsb.expire()

        
def chemistrycodesbibliographiesexpire(ccb):
    ccb.expire()
        
def vibrationalanalysesharmonicbibliographiesexpire(vahb):
    vahb.expire()
        
def vibrationalanalysesanarmonicbibliographiesexpire(vaab):
    vaab.expire()
        
def electronictransitionsbibliographiesexpire(etb):
    etb.expire()
        
def ionisationenergiesbibliographiesexpire(ieb):
    ieb.expire()

def electronicstatesbibliographiesexpire(esb):
    esb.expire()
        
def rotationalconstantsbibliographiesexpire(rcb):
    rcb.expire()
                
def dipolemomentsbibliographiesexpire(dmb):
    dmb.expire()
        
def polarisabilitiesbibliographiesexpire(pb):
    pb.expire()
        
def vanderwallsbibliographiesexpire(vdwb):
    vdwb.expire()
        