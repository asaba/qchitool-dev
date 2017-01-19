'''
Created on 03/mag/2012

@author: asaba
'''
from ...models import TheoryLevels, Xcclasses
from ... import performancetest as pt

def build_functional_classes():
    #build a dictionary for xc field
    xc_classes = {}
    for tl in TheoryLevels.alive_objects.all():
        if tl.xc_name in xc_classes:
            xc_classes[tl.xc_name].append(tl.pk)
        else:
            xc_classes[tl.xc_name] = [tl.pk]
        
    for xc_class in xc_classes:
        xc_name = TheoryLevels.alive_objects.get(pk = xc_classes[xc_class][0]).xc_name
        if not xc_name:
            xc_name = "NONE"
        xc_class_record = Xcclasses(name = xc_name)
        xc_class_record.save()
        for tl in TheoryLevels.alive_objects.filter(pk__in = xc_classes[xc_class]):
            tl.xcclass = xc_class_record
            tl.save()
            