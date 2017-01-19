from models import *
import performancetest as pt
from tools.expire import molecularspeciesexpire, bibliographyexpire, editorsexpire, reftypeexpire, publishersexpire, authorsexpire, dipolemomentsexpire, electronicstatesexpire, chemistrycodesexpire, calculationsexpire, tasksexpire, basissetsexpire, theorylevelsexpire, geometriesexpire, calculationsexpire, ionisationenergiesexpire

from django.contrib import admin 
#from django.contrib.admin.helpers import AdminForm

#Molecules

#class IsotopesInLine(admin.TabularInline):
#    model = Isotopes
#    extra = 5

class ElementsAdmin(admin.ModelAdmin):
    list_display = ('atomic_number', "atomic_mass", 'symbol', 'name', "element_group", "standard_atomic_weight")
    #exclude = ("isotope_of",)

class IonisationTypesAdmin(admin.ModelAdmin):
    list_display = ("description",)
    
class IonizationenergiesAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for ie in queryset:
            if not ie.expired():
                ionisationenergiesexpire(ie)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("start_state", "ion_state", "iontype", "energy", "expired",)
    list_filter = ("iontype",)
    filter_horizontal = ("bibliographies",)
    actions = (expire, )
    expire.short_description = "Expire"
#class MolecularIsotopologuesInLine(admin.TabularInline):
#    model = MolecularIsotopologues
#    estras = 3
    
class MolecularSpeciesAdmin(admin.ModelAdmin):
    list_display = ("inchikey", "name", "formula", "aromatic_cycles", "charge", "time_stamp", "expired")
    list_filter = ("charge", "formula",)
    search_fields = ("inchikey", "name",)
    isotopologue_of = models.ForeignKey("self", blank=True, null=True)
    def update_aromatic_cycle_to_plus_1(self, request, queryset):
        for ms in queryset:
            ms.aromatic_cycles += 1
            ms.save()
        if len(queryset) == 1:
            message_bit = "1 Molecular Specie was"
        else:
            message_bit = "%s Molecular Specie" % len(queryset)
        self.message_user(request, "%s successfully changed." % message_bit)
    def update_aromatic_cycle_to_minus_1(self, request, queryset):
        for ms in queryset:
            ms.aromatic_cycles -= 1
            ms.save
        if len(queryset) == 1:
            message_bit = "1 Molecular Specie was"
        else:
            message_bit = "%s Molecular Specie" % len(queryset)
        self.message_user(request, "%s successfully changed." % message_bit)
    def expire(self, request, queryset):
        changed = 0
        for ms in queryset:
            if not ms.expired():
                molecularspeciesexpire(ms)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    actions = (update_aromatic_cycle_to_plus_1, update_aromatic_cycle_to_minus_1, expire, )
    update_aromatic_cycle_to_plus_1.short_description = "Update Aromatic Cycles to +1"
    update_aromatic_cycle_to_minus_1.short_description = "Update Aromatic Cycles to -1"
    expire.short_description = "Expire"
    
class GeometriesAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for g in queryset:
            if not g.expired():
                geometriesexpire(g)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("geometry_md5", "sym_group", "sym_elements", "time_stamp", "expired")
    actions = (expire, )
    expire.short_description = "Expire"


class GeometryclassesAdmin(admin.ModelAdmin):
    list_display = ("geometry_md5", "sym_group", "sym_elements")

class XcclassesAdmin(admin.ModelAdmin):
    list_display = ("name", "description")

class TheoryLevelsAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for tl in queryset:
            if not tl.expired():
                theorylevelsexpire(tl)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ('name', 'description', 'xc_name', 'xc_description', 'time_stamp', "expired")
    filter_horizontal = ("bibliographies",)
    actions = (expire, )
    expire.short_description = "Expire"
    
class BasisSetsAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for bs in queryset:
            if not bs.expired():
                basissetsexpire(bs)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("name", "description", "time_stamp", "qual_index", "expired")
    filter_horizontal = ("bibliographies",)
    actions = (expire, )
    expire.short_description = "Expire"
    
class TasksAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for t in queryset:
            if not t.expired():
                tasksexpire(t)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    #fields = (("thlevel", "calc", "name", "description", "time_stamp", "qual_index", "comments"))
    list_display = ("thlevel", "name", "description", "time_stamp", "qual_index", "expired")
    actions = (expire, )
    expire.short_description = "Expire"
    
class ChemistryCodesAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for cc in queryset:
            if not cc.expired():
                chemistrycodesexpire(cc)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("name", "version", "description", "time_stamp", "qual_index", "expired")
    filter_horizontal = ("bibliographies",)
    actions = (expire, )
    expire.short_description = "Expire"
    
class CalculationsAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for c in queryset:
            if not c.expired():
                calculationsexpire(c)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("time_stamp", "qual_index", "comments", "expired")
    actions = (expire, )
    expire.short_description = "Expire"

class ElectronicStatesAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for es in queryset:
            if not es.expired():
                electronicstatesexpire(es)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("total_energy", "is_minimum", "rel_energy", "symmetry", "multiplicity", "description", "time_stamp", "expired")
    list_filter = ("is_minimum", "multiplicity", )
    search_fields = ("description",)
    filter_horizontal = ("bibliographies",)
    electronic_state_energy = models.ForeignKey("self", null=True, blank=True)
    actions = (expire, )
    expire.short_description = "Expire"

class DipoleMomentsAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for dm in queryset:
            if not dm.expired():
                dipolemomentsexpire(dm)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("mu_x", "mu_y", "mu_z", "time_stamp", "qual_index", "expired")
    filter_horizontal = ("bibliographies",)
    actions = (expire, )
    expire.short_description = "Expire"
#Bibliography

class AuthorsAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for a in queryset:
            if not a.expired():
                authorsexpire(a)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ('name', 'address', 'email', 'time_stamp', "expired")
    actions = (expire, )
    expire.short_description = "Expire"
    
class PublishersAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for p in queryset:
            if not p.expired():
                publishersexpire(p)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ('name', 'address', 'time_stamp', "expired")
    actions = (expire, )
    expire.short_description = "Expire"
    
class PublicationSeriesAdmin(admin.ModelAdmin):
    list_display = ('name', 'shortname')
    
class ReftypeAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for rt in queryset:
            if not rt.expired():
                reftypeexpire(rt)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ('description', 'time_stamp', "expired")
    actions = (expire, )
    expire.short_description = "Expire"

class EditorsAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for e in queryset:
            if not e.expired():
                editorsexpire(e)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("name", "address", "expired")
    actions = (expire, )
    expire.short_description = "Expire"
    
class BibliographyAdmin(admin.ModelAdmin):
    def expire(self, request, queryset):
        changed = 0
        for b in queryset:
            if not b.expired():
                bibliographyexpire(b)
                changed += 1
        self.message_user(request, "%s records on %s successfully changed." % changed, queryset.count())
    list_display = ("title", "date", "expired")
    filter_horizontal = ("authors","editors",)
    actions = (expire, )
    expire.short_description = "Expire"
    
#class VibrationalAnalysesAnarmonicAdmin(admin.ModelAdmin):
#    list_display = ("polymode", "time_stamp")
#Polarisability



admin.site.register(Elements, ElementsAdmin)
admin.site.register(TheoryLevels, TheoryLevelsAdmin)
admin.site.register(Authors, AuthorsAdmin)
admin.site.register(Publishers, PublishersAdmin)
admin.site.register(PublicationSeries, PublicationSeriesAdmin)
admin.site.register(Reftype, ReftypeAdmin)
admin.site.register(BasisSets, BasisSetsAdmin)
admin.site.register(Tasks, TasksAdmin)
admin.site.register(ChemistryCodes, ChemistryCodesAdmin)
admin.site.register(Editors, EditorsAdmin)
admin.site.register(Bibliography, BibliographyAdmin)
admin.site.register(Calculations, CalculationsAdmin)
admin.site.register(IonisationTypes, IonisationTypesAdmin)
admin.site.register(IonisationEnergies, IonizationenergiesAdmin)
admin.site.register(MolecularSpecies, MolecularSpeciesAdmin)
admin.site.register(ElectronicStates, ElectronicStatesAdmin)
admin.site.register(Geometries, GeometriesAdmin)
admin.site.register(Geometryclasses, GeometryclassesAdmin)
admin.site.register(Xcclasses, XcclassesAdmin)
admin.site.register(DipoleMoments, DipoleMomentsAdmin)
#admin.site.register(VibrationalAnalysesAnarmonic, VibrationalAnalysesAnarmonicAdmin)