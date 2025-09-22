from django.urls import path

from . import views
from . import views_classic
from . import api

app_name = 'nice_lite' # Define an app namespace
urlpatterns = [
    ## Global URLs
    path("",                 views.index,          name="index"),
    path("login",            views.user_login,     name="login"),
    path("logout",           views.user_logout,    name="logout"),
    path("image <path:src>", api.image,            name="image"),
    path("api",              api.index,            name="api"),
    path("newproject/<str:caller>",       views.new_project,    name="new_project"),
    path("createproject",    views.create_project, name="create_project"),

    ## Classic specific URLs
    path("classic",         views_classic.classic,  name="classic"),
    path("workspace",       views_classic.workspace,          name="workspace"),
    path("createworkspace/<int:projectid>", views_classic.create_workspace,   name="create_workspace"),
    path("newjobtype/<int:parentid>",          views_classic.new_job_type,       name="new_job_type"),
    path("newjob/<int:parentid>/<str:package>/<str:jobtype>",          views_classic.new_job,       name="new_job"),
    path("rerunjob/<int:parentid>",          views_classic.rerun_job,       name="rerun_job"),
    path("createjob/<int:parentid>/<str:package>/<str:jobtype>",          views_classic.create_job,       name="create_job"),
    path("updateworkspacename",          views_classic.update_workspace_name,       name="update_workspace_name"),
    path("updateworkspacedescription",   views_classic.update_workspace_description,       name="update_workspace_description"),
    path("markjobcomplete/<int:jobid>",   views_classic.mark_job_complete,       name="mark_job_complete"),
    path("updatejobdescription/<int:jobid>",   views_classic.update_job_description,       name="update_job_description"),

    ## Stream specific URLs
    path("stream",          views.stream,           name="stream"),
    path("createdataset",   views.create_dataset,   name="create_dataset"),
    path("dataset",         views.dataset,          name="dataset"),
    path("newstream",       views.new_stream,       name="new_stream"),
    path("termstream/<int:jobid>",      views.terminate_stream, name="terminate_stream"),
    path("updatedatasetname",          views.update_dataset_name,       name="update_dataset_name"),
    path("updatedatasetdescription",   views.update_dataset_description,       name="update_dataset_description"),
    path("createstream",    views.create_stream,    name="create_stream"),
    path("rerunstream/<int:parentid>",    views.rerun_stream,    name="rerun_stream"),
    path("updatestreamdescription/<int:jobid>",    views.update_stream_description,    name="update_stream_description"),

    path("viewstream/<int:jobid>",      views.view_stream,      name="view_stream"),
    path("viewlogs/<int:jobid>/<path:log>/<path:error>",       views.view_logs,         name="view_logs"),
    path("viewstreammovies/<int:jobid>",         views.view_stream_movies,        name="view_stream_movies"),
    path("viewstreampreprocess/<int:jobid>",         views.view_stream_preprocess,        name="view_stream_preprocess"),
    path("viewstreampreprocesszoom/<int:jobidzoom>", views.view_stream_preprocess,        name="view_stream_preprocess_zoom"),
    path("viewstreamoptics/<int:jobid>",           views.view_stream_optics,            name="view_stream_optics"),
    path("viewstreaminitialpick/<int:jobid>",      views.view_stream_initial_pick,      name="view_stream_initial_pick"),
    path("viewstreamgeneratepickrefs/<int:jobid>", views.view_stream_generate_pickrefs, name="view_stream_generate_pickrefs"),
    path("viewstreamreferencepicking/<int:jobid>",  views.view_stream_reference_picking, name="view_stream_reference_picking"),
    path("viewstreamsieveparticles/<int:jobid>",   views.view_stream_sieve_particles,    name="view_stream_sieve_particles"),
    path("viewstreamclassification2D/<int:jobid>",   views.view_stream_classification_2D,    name="view_stream_classification_2D"),
    path("viewstreamparticlesets/<int:jobid>",   views.view_stream_particle_sets,    name="view_stream_particle_sets"),

    path("termstreampreprocess/<int:jobid>",       views.term_stream_preprocess,        name="term_stream_preprocess"),
    path("termstreamoptics/<int:jobid>",          views.term_stream_optics,             name="term_stream_optics"),
    path("termstreaminitialpick/<int:jobid>",      views.term_stream_initial_pick,        name="term_stream_initial_pick"),
    path("termstreamgeneratepickrefs/<int:jobid>",      views.term_stream_generate_pickrefs,        name="term_stream_generate_pickrefs"),
    path("termstreamreferencepicking/<int:jobid>",      views.term_stream_reference_picking,        name="term_stream_reference_picking"),
    path("termstreamsieveparticles/<int:jobid>",      views.term_stream_sieve_particles,        name="term_stream_sieve_particles"),
    path("termstreamclassification2D/<int:jobid>",      views.term_stream_classification_2D,        name="term_stream_classification_2D"),

    path("updatestreampreprocessctfres/<int:jobid>",      views.update_preprocess_ctfres,          name="update_preprocess_ctfres"),
    path("updatestreampreprocessastig/<int:jobid>",       views.update_preprocess_astig,           name="update_preprocess_astig"),
    path("updatestreampreprocessicescore/<int:jobid>",    views.update_preprocess_icescore,        name="update_preprocess_icescore"),
    path("updatestreamclassification2Dmskdiam/<int:jobid>",    views.update_classification_2D_mskdiam,        name="update_classification_2D_mskdiam"),

    path("selectmoldiamstreaminitialpick/<int:jobid>",      views.select_moldiam_stream_initial_pick,        name="select_moldiam_stream_initial_pick"),
    path("increasemoldiamstreaminitialpick/<int:jobid>",      views.increase_moldiam_stream_initial_pick,        name="increase_moldiam_stream_initial_pick"),
    path("decreasemoldiamstreaminitialpick/<int:jobid>",      views.decrease_moldiam_stream_initial_pick,        name="decrease_moldiam_stream_initial_pick"),
    path("selectrefsstreamgeneratepickrefs/<int:jobid>",      views.select_refs_stream_generate_pickrefs,        name="select_refs_stream_generate_pickrefs"),
    path("snapshotstreamclassification2D/<int:jobid>",      views.snapshot_stream_classification_2D,        name="snapshot_stream_classification_2D"),
    path("selectstreamsieveparticles/<int:jobid>",         views.select_stream_sieve_particles,        name="select_stream_sieve_particles"),
    path("linkstreamparticleset/<int:jobid>/<int:setid>/<path:filename>",         views.link_stream_particle_set,        name="link_stream_particle_set"),
   
]