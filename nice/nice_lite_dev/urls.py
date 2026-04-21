from django.urls import path

from . import views
from . import views_classic
from . import views_stream
from . import api

app_name = 'nice_lite'

urlpatterns = [

    # ------------------------------------------------------------------
    # Global
    # ------------------------------------------------------------------
    path("",                                   views.index,          name="index"),
    path("login",                              views.user_login,     name="login"),
    path("logout",                             views.user_logout,    name="logout"),
    path("newproject/<str:mode>",              views.new_project,    name="new_project"),
    path("createproject",                      views.create_project, name="create_project"),
    path("filebrowser/<str:type>",             views.file_browser,   name="file_browser"),
    path("filebrowser/<str:type>/<path:path>", views.file_browser,   name="file_browser"),
    path("image <path:src>",                   api.image,            name="image"),
    path("api",                                api.index,            name="api"),
    path("api_classic",                        api.index_classic,    name="api_classic"),

    # ------------------------------------------------------------------
    # Classic mode — workspace and job management
    # ------------------------------------------------------------------
    path("classic",                                              views_classic.classic,                        name="classic"),
    path("workspace",                                            views_classic.workspace,                      name="workspace"),
    path("createworkspace/<int:projectid>",                      views_classic.create_workspace,               name="create_workspace"),
    path("deleteworkspace/<int:workspaceid>",                    views_classic.delete_workspace,               name="delete_workspace"),
    path("newjobtype/<int:parentid>",                            views_classic.new_job_type,                   name="new_job_type"),
    path("newjob/<int:parentid>/<str:package>/<str:jobtype>",    views_classic.new_job,                        name="new_job"),
    path("rerunjob/<int:jobid>",                                 views_classic.rerun_job,                      name="rerun_job"),
    path("createjob/<int:parentid>/<str:package>/<str:jobtype>", views_classic.create_job,                     name="create_job"),
    path("updateworkspacename",                                  views_classic.update_workspace_name,          name="update_workspace_name"),
    path("updateworkspacedescription",                           views_classic.update_workspace_description,   name="update_workspace_description"),
    path("markjobcomplete/<int:jobid>",                          views_classic.mark_job_complete,              name="mark_job_complete"),
    path("updatejobdescription/<int:jobid>",                     views_classic.update_job_description,         name="update_job_description"),
    path("viewjob/<int:jobid>",                                  views_classic.view_job,                       name="view_job"),
    path("viewjobmicrographs/<int:jobid>",                       views_classic.view_job_micrographs,           name="view_job_micrographs"),
    path("viewjobmicrographshistogram/<int:jobid>",              views_classic.view_job_micrographs_histogram, name="view_job_micrographs_histogram"),
    path("viewjobmicrographsplot/<int:jobid>",                   views_classic.view_job_micrographs_plot,      name="view_job_micrographs_plot"),
    path("viewjobcls2d/<int:jobid>",                             views_classic.view_job_cls2D,                 name="view_job_cls2D"),
    path("viewjobcls2dhistogram/<int:jobid>",                    views_classic.view_job_cls2D_histogram,       name="view_job_cls2D_histogram"),
    path("viewjobcls2dplot/<int:jobid>",                         views_classic.view_job_cls2D_plot,            name="view_job_cls2D_plot"),
    path("viewjoblogs/<int:jobid>",                              views_classic.view_job_logs,                  name="view_job_logs"),
    path("selectjobmicrographs/<int:jobid>",                     views_classic.select_job_micrographs,         name="select_job_micrographs"),
    path("selectjobcls2d/<int:jobid>",                           views_classic.select_job_cls2D,               name="select_job_cls2D"),
    path("deletejob/<int:jobid>",                                views_classic.delete_job,                     name="delete_job"),

    # ------------------------------------------------------------------
    # Stream mode — dataset lifecycle
    # ------------------------------------------------------------------
    path("stream",                    views_stream.stream_index,              name="stream"),
    path("dataset",                   views_stream.dataset,                   name="dataset"),
    path("deletedataset",             views_stream.delete_dataset,            name="delete_dataset"),
    path("updatedatasetname",         views_stream.update_dataset_name,       name="update_dataset_name"),
    path("updatedatasetdescription",  views_stream.update_dataset_description, name="update_dataset_description"),

    # Stream job lifecycle
    path("newstream",                 views_stream.new_stream,                name="new_stream"),
    path("createstream",              views_stream.create_stream,             name="create_stream"),
    path("termstream",                views_stream.terminate_stream,          name="terminate_stream"),
    path("deletestream",              views_stream.delete_stream,             name="delete_stream"),
    path("updatestreamdescription",   views_stream.update_stream_description, name="update_stream_description"),
    path("updatestream",              views_stream.update_stream_parameters,  name="update_stream_parameters"),

    # Stream sub-process control (terminate / restart)
    path("termstreamprocess",                               views_stream.terminate_stream_process,          name="term_stream_process"),
    path("restartstreamprocess",                            views_stream.restart_stream_process,            name="restart_stream_process"),

    # Stream stage views (panel + zoom)
    path("viewstream/<int:jobid>",          views_stream.view_stream,                        name="view_stream"),
    path("viewlogs/<int:jobid>/<path:log>/<path:error>", views_stream.view_logs,             name="view_logs"),
    path("viewstreammovies",                views_stream.view_stream_movies,                 name="view_stream_movies"),
    path("viewstreampreprocess",            views_stream.view_stream_preprocess,             name="view_stream_preprocess"),
    path("viewstreampreprocesszoom",        views_stream.view_stream_preprocess_zoom,        name="view_stream_preprocess_zoom"),
    path("viewstreamoptics",                views_stream.view_stream_optics,                 name="view_stream_optics"),
    path("viewstreamopticszoom",            views_stream.view_stream_optics_zoom,            name="view_stream_optics_zoom"),
    path("viewstreaminitialpick",           views_stream.view_stream_initial_pick,           name="view_stream_initial_pick"),
    path("viewstreaminitialpickzoom",       views_stream.view_stream_initial_pick_zoom,      name="view_stream_initial_pick_zoom"),
    path("viewstreamgeneratepickrefs",      views_stream.view_stream_generate_pickrefs,      name="view_stream_generate_pickrefs"),
    path("viewstreamgeneratepickrefszoom",  views_stream.view_stream_generate_pickrefs_zoom, name="view_stream_generate_pickrefs_zoom"),
    path("viewstreamreferencepicking",      views_stream.view_stream_reference_picking,      name="view_stream_reference_picking"),
    path("viewstreamreferencepickingzoom",  views_stream.view_stream_reference_picking_zoom, name="view_stream_reference_picking_zoom"),
    path("viewstreamsieveparticles",        views_stream.view_stream_sieve_particles,        name="view_stream_sieve_particles"),
    path("viewstreamsieveparticleszoom",    views_stream.view_stream_sieve_particles_zoom,   name="view_stream_sieve_particles_zoom"),
    path("viewstreamclassification2D",      views_stream.view_stream_classification_2D,      name="view_stream_classification_2D"),
    path("viewstreamclassification2Dzoom",  views_stream.view_stream_classification_2D_zoom, name="view_stream_classification_2D_zoom"),
    path("viewstreamparticlesets",          views_stream.view_stream_particle_sets,          name="view_stream_particle_sets"),

    # Stream user interactions (parameter updates, selections)
    path("updatestreamclassification2Dmskdiam", views_stream.update_classification_2D_mskdiam,     name="update_classification_2D_mskdiam"),
    path("snapshotstreamclassification2D",      views_stream.snapshot_stream_classification_2D,    name="snapshot_stream_classification_2D"),
    path("selectstreamclassification2D",        views_stream.select_stream_classification_2D,      name="select_stream_classification_2D"),
    path("regeneratepickrefs",                  views_stream.regenerate_pickrefs,                  name="regenerate_pickrefs"),
    path("selectpickrefs",                      views_stream.select_pickrefs,                      name="select_pickrefs"),
    path("selectstreamsieveparticles",          views_stream.select_stream_sieve_particles,        name="select_stream_sieve_particles"),
    
    path("linkstreamparticleset/<int:jobid>/<int:setid>/<path:filename>/<path:type>",  views_stream.link_stream_particle_set,             name="link_stream_particle_set"),

]
