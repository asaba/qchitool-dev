from django.conf.urls.defaults import patterns, include, url
from django.conf.urls.defaults import *
from django.views.generic import RedirectView

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
from django.conf import settings
admin.autodiscover()


urlpatterns = patterns('qchitool.views',
    (r'^info/$', 'global_index', {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
    (r'^nwc_import/$', 'nwc_import_index'),
    (r'^nwc_import/file_info/$', 'nwc_import_file_info', {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
    (r'^nwc_import/task_info/$', 'nwc_import_file_task_info'),
    (r'^nwc_import/atoms_info/$', 'nwc_import_file_atoms_info'),
    (r'^nwc_import/to_save/$', 'nwc_post_selected_task_to_save'),
    (r'^nwc_import/save/$', 'nwc_calculus_information_to_save'),
    (r'^nwc_import/save2/$', 'nwc_electronic_states_similar'),
    (r'^octopus_import/$', 'octopus_import_index'),
    (r'^octopus_import/file_info/$', 'octopus_import_file_info', {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),

)

urlpatterns += patterns('',
    (r'^admin/', include(admin.site.urls)),
    (r'^accounts/login/$', 'django.contrib.auth.views.login'),

)

#urlpatterns += patterns("", 
#                        (r'^jmol/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '/home/asaba/workspace/pah_import_output/src/oacagliari/qchitool/templates/JmolFolder/'}),
#                        (r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.PATH_FOR_MOLECULE_FILES}),
#                        (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.PATH_STATIC_FILES}),
#                        (r'^molecule/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT}),
                        
#                        )


urlpatterns += patterns("qchitool.views",
                        (r"^utility/$", "utility_index"),
                        (r"^utility/geom/$", "utility_geom_classes", {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
                        (r"^utility/functional/$", "utility_functional_classes", {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
                        (r"^utility/deleteionizationenergy/$", "utility_delete_ionizationenergies", {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
                        (r"^utility/ionizationenergy/$", "utility_calculate_ionizationenergies", {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
                        (r"^utility/molecularelements/$", "utility_rebuild_elementlist_for_molecules", {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),)

urlpatterns += patterns("qchitool.views",
                        (r"^$", "explore_index",{"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),
                        (r"^database/molecule/(\d+)/$", "show_molecule_information", {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}), #pass molecule ID
                        (r"^database/elementbasissets/(\d+)/(\d+)/(\d+)/$", "show_elemetbasisset_information", {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}), #pass molecule ID, Electronic state ID, Task ID
                        (r"^database/vibrationalanalysesanarmonic/(\d+)/(\d+)/$", "show_vibrationalanalysesanarmonic_information", {"PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}), #pass molecule ID, Electronic State ID
                        (r"^database/vibrationalanalysesharmonic/(\d+)/(\d+)/$", "show_vibrationalanalysesharmonic_information", {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}), #pass molecule ID, Electronic State ID
                        (r"^database/eigenvectors/(\d+)/(\d+)/(\d+)/(\d+)/$", "show_vibrationalanalysesharmoniceigenvectors_information", {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}),#pass molecule ID, Electronic State ID, vibrational analisis ID
                        (r"^database/ionisationenergies/(\d+)/(\d+)/(\d+)/$", "show_ionisationenergies_information", {"URL_FOR_MOLECULE_FILES": settings.URL_FOR_MOLECULE_FILES, "PATH_FOR_MOLECULE_FILES": settings.PATH_FOR_MOLECULE_FILES}) #pass molecule ID, electronicstateid, ionisation energy ID
                        )
