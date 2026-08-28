"""Lifecycle actions for classic SIMPLE and SINGLE batch jobs."""

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect
from django.views.decorators.http import require_POST

from ..data_structures.batchjob import BatchJob
from ..data_structures.project import Project
from ..data_structures.workspace import Workspace
from ..features import batch_job_controls_enabled
from ..helpers import get_job_id, print_error


def _get_accessible_batch_job(request, log_context):
    """Return an owned batch job and model selected by the request."""
    job_id = get_job_id(request)
    if job_id is None:
        print_error(f"{log_context}: missing job id")
        return None, None

    batch_job = BatchJob(id=job_id)
    jobmodel = batch_job.get_jobmodel()
    workspace = getattr(jobmodel, "dset", None)
    owner = (getattr(workspace, "user", "") or "").strip()
    if jobmodel is None or owner != request.user.username:
        print_error(f"{log_context}: access denied for job {job_id}")
        return None, None
    return batch_job, jobmodel


def _controls_available(request):
    """Reject lifecycle writes while the opt-in feature is disabled."""
    if batch_job_controls_enabled():
        return True
    messages.add_message(request, messages.ERROR, "batch job controls are disabled")
    return False


@login_required(login_url="/login")
@require_POST
def view_batch_stop(request):
    """Stop an owned running batch job and return to its workspace."""
    if not _controls_available(request):
        return redirect("nice_lite:workspace")

    batch_job, jobmodel = _get_accessible_batch_job(request, "stop_batch")
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
    elif jobmodel.status != "running":
        messages.add_message(request, messages.ERROR, "batch job is not running")
    elif batch_job.stop():
        messages.add_message(request, messages.INFO, "batch job stopped successfully")
    else:
        messages.add_message(request, messages.ERROR, "failed to stop batch job")
    return redirect("nice_lite:workspace")


@login_required(login_url="/login")
@require_POST
def view_batch_delete(request):
    """Permanently delete an owned terminal batch job and return to its workspace."""
    if not _controls_available(request):
        return redirect("nice_lite:workspace")

    batch_job, jobmodel = _get_accessible_batch_job(request, "delete_batch")
    if batch_job is None:
        messages.add_message(request, messages.ERROR, "invalid batch job selection")
    elif jobmodel.status not in BatchJob.DELETABLE_STATUSES:
        messages.add_message(request, messages.ERROR, "batch job is not complete")
    else:
        project = Project(id=jobmodel.dset.proj_id)
        workspace = Workspace(jobmodel.dset_id)
        if batch_job.delete(project, workspace):
            messages.add_message(request, messages.INFO, "batch job permanently deleted successfully")
        else:
            messages.add_message(request, messages.ERROR, "failed to delete batch job")
    return redirect("nice_lite:workspace")
