"""Shared assertion helpers for legacy NICE Lite tests."""

from django.test import TestCase


def _workspace_field(workspace, field):
  if hasattr(workspace, field):
    return getattr(workspace, field)
  workspacemodel = getattr(workspace, "workspacemodel", None)
  if workspacemodel is not None and hasattr(workspacemodel, field):
    return getattr(workspacemodel, field)
  return None


def assertProject(project, name=None, desc=None, dirc=None, id=None):
  testcase = TestCase()
  if name is not None:
    testcase.assertEqual(project.name, name, "project name is not equal to " + name)
  if desc is not None:
    testcase.assertEqual(project.desc, desc, "project description is not equal to " + desc)
  if dirc is not None:
    testcase.assertEqual(project.dirc, dirc, "project directory is not equal to " + dirc)
  if id is not None:
    testcase.assertEqual(project.id, id, "project id is not equal to " + str(id))

def assertWorkspace(workspace, id=None, user=None, name=None, desc=None):
  testcase = TestCase()
  if id is not None:
    testcase.assertEqual(_workspace_field(workspace, "id"), id, "workspace id is not equal to " + str(id))
  if user is not None:
    testcase.assertEqual(_workspace_field(workspace, "user"), user, "workspace user is not equal to " + user)
  if name is not None:
    testcase.assertEqual(_workspace_field(workspace, "name"), name, "workspace name is not equal to " + name)
  if desc is not None:
    testcase.assertEqual(_workspace_field(workspace, "desc"), desc, "workspace description is not equal to " + desc)

def assertDispatch(dispatch, id=None):
  testcase = TestCase()
  if id is not None:
    testcase.assertEqual(dispatch.id, id, "dispatch id is not equal to " + str(id))