"""URL routes for the nice_lite app.

This module defines:
- top-level index/auth routes
- API and file-browser routes
- project/workspace routes
- stream lifecycle/stage interaction routes
"""

# django imports
from django.urls import path

# local imports
from .      import api
from .      import views
from .views import stream_views

app_name = "nice_lite"

urlpatterns = [

    # ------------------------------------------------------------------
    # Index
    # ------------------------------------------------------------------
    path("",                                   views.view_index,                                  name="index"),
    path("login",                              views.view_user_login,                             name="login"),
    path("logout",                             views.view_user_logout,                            name="logout"),

    # ------------------------------------------------------------------
    # API
    # ------------------------------------------------------------------
    path("image <path:src>",                   api.image,                                         name="image"),
    path("api",                                api.index,                                         name="api"),
    path("api_classic",                        api.index_classic,                                 name="api_classic"),

    # ------------------------------------------------------------------
    # File browser
    # ------------------------------------------------------------------
    path("filebrowser/<str:type>",             views.view_file_browser,                           name="file_browser"),
    path("filebrowser/<str:type>/<path:path>", views.view_file_browser,                           name="file_browser"),

    # ------------------------------------------------------------------
    # Project
    # ------------------------------------------------------------------
    path("newproject/<str:mode>",              views.view_new_project,                            name="new_project"),
    path("createproject",                      views.view_create_project,                          name="create_project"),

    # ------------------------------------------------------------------
    # Workspace
    # ------------------------------------------------------------------
    path("workspace",                          views.view_workspace,                               name="workspace"),
    path("workspacejobs",                      views.view_workspace_jobs,                          name="workspace_jobs"),
    path("deleteworkspace",                    views.view_delete_workspace,                        name="delete_workspace"),
    path("updateworkspacename",                views.view_update_workspace_name,                   name="update_workspace_name"),
    path("updateworkspacedescription",         views.view_update_workspace_description,            name="update_workspace_description"),

    # Stream job lifecycle
    path("newstream",                          views.view_job_builder,                             name="new_stream"),
    path("createstream",                       stream_views.view_stream_create_stream,             name="create_stream"),
    path("termstream",                         stream_views.view_stream_terminate_stream,          name="terminate_stream"),
    path("deletestream",                       stream_views.view_stream_delete_stream,             name="delete_stream"),
    path("killdeletestream",                   stream_views.view_stream_kill_delete_stream,        name="kill_delete_stream"),
    path("updatestreamdescription",            stream_views.view_stream_update_description,        name="update_stream_description"),
    path("updatestream",                       stream_views.view_stream_update_parameters,         name="update_stream_parameters"),

    # Stream sub-process control (terminate / restart)
    path("termstreamprocess",                  stream_views.view_stream_terminate_stream_process,  name="term_stream_process"),
    path("restartstreamprocess",               stream_views.view_stream_restart_stream_process,    name="restart_stream_process"),

    # Stream stage views (panel + zoom)
    path("viewstream/<int:jobid>",             stream_views.view_stream,                           name="view_stream"),
    path("viewlogs/<int:jobid>/<path:log>/<path:error>", stream_views.view_stream_logs,            name="view_logs"),
    path("viewstreammovies",                   stream_views.view_stream_movies,                    name="view_stream_movies"),
    path("viewstreampreprocess",               stream_views.view_stream_preprocess,                name="view_stream_preprocess"),
    path("viewstreampreprocesszoom",           stream_views.view_stream_preprocess_zoom,           name="view_stream_preprocess_zoom"),
    path("viewstreamoptics",                   stream_views.view_stream_optics,                    name="view_stream_optics"),
    path("viewstreamopticszoom",               stream_views.view_stream_optics_zoom,               name="view_stream_optics_zoom"),
    path("viewstreaminitialpick",              stream_views.view_stream_initial_pick,              name="view_stream_initial_pick"),
    path("viewstreaminitialpickzoom",          stream_views.view_stream_initial_pick_zoom,         name="view_stream_initial_pick_zoom"),
    path("viewstreamgeneratepickrefs",         stream_views.view_stream_generate_pickrefs,         name="view_stream_generate_pickrefs"),
    path("viewstreamgeneratepickrefszoom",     stream_views.view_stream_generate_pickrefs_zoom,    name="view_stream_generate_pickrefs_zoom"),
    path("viewstreamreferencepicking",         stream_views.view_stream_reference_picking,         name="view_stream_reference_picking"),
    path("viewstreamreferencepickingzoom",     stream_views.view_stream_reference_picking_zoom,    name="view_stream_reference_picking_zoom"),
    path("viewstreamsieveparticles",           stream_views.view_stream_sieve_particles,           name="view_stream_sieve_particles"),
    path("viewstreamsieveparticleszoom",       stream_views.view_stream_sieve_particles_zoom,      name="view_stream_sieve_particles_zoom"),
    path("viewstreamclassification2D",         stream_views.view_stream_classification_2D,         name="view_stream_classification_2D"),
    path("viewstreamclassification2Dzoom",     stream_views.view_stream_classification_2D_zoom,    name="view_stream_classification_2D_zoom"),
    path("viewstreamparticlesets",             stream_views.view_stream_particle_sets,             name="view_stream_particle_sets"),

    # Stream user interactions (parameter updates, selections)
    path("updatestreamclassification2Dmskdiam", stream_views.view_stream_update_classification_2D_mskdiam, name="update_classification_2D_mskdiam"),
    path("snapshotstreamclassification2D",      stream_views.view_stream_snapshot_classification_2D,       name="snapshot_stream_classification_2D"),
    path("selectstreamclassification2D",        stream_views.view_stream_select_classification_2D,         name="select_stream_classification_2D"),
    path("selectpickrefs",                      stream_views.view_stream_select_pickrefs,                  name="select_pickrefs"),

    path(
        "linkstreamparticleset/<int:jobid>/<int:setid>/<path:filename>/<path:type>",
        stream_views.view_stream_link_particle_set,
        name="link_stream_particle_set",
    ),

]
