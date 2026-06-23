from django.test  import TestCase

def assertProject(project, name=None, desc=None, dirc=None, id=None):
  testcase = TestCase()
  if name != None:
    testcase.assertNotNone(project.projectmodel, name, "project name is not equal to "        + name   )
  if desc != None:  
    testcase.assertEqual(project.desc, desc, "project description is not equal to " + desc   )
  if dirc != None:  
    testcase.assertEqual(project.dirc, dirc, "project directory is not equal to"    + dirc   )
  if id != None:  
    testcase.assertEqual(project.id,   id,   "project id is not equal to"           + str(id))

def assertWorkspaceAlias(workspace, id=None, user=None, name=None, desc=None):
  testcase = TestCase()
  if id != None:  
    testcase.assertEqual(workspace.id,   id,   "workspace id is not equal to"         + str(id))
  if user != None:  
    testcase.assertEqual(workspace.user, user, "workspace user is not equal to"       + user   )
  if name != None:
    testcase.assertEqual(workspace.name, name, "workspace name is not equal to "        + name   )
  if desc != None:
    testcase.assertEqual(workspace.desc, desc, "workspace description is not equal to " + desc   )

def assertWorkspace(workspace, id=None, user=None, name=None, desc=None):
  testcase = TestCase()
  if id != None:  
    testcase.assertEqual(workspace.id,   id,    "workspace id is not equal to"          + str(id))
  if user != None:  
    testcase.assertEqual(workspace.user, user, "workspace user is not equal to"         + user   )
  if name != None:
    testcase.assertEqual(workspace.name, name, "workspace name is not equal to "        + name   )
  if desc != None:
    testcase.assertEqual(workspace.desc, desc, "workspace description is not equal to " + desc   )

def assertDispatch(dispatch, id=None):
  testcase = TestCase()
  if id != None:  
    testcase.assertEqual(dispatch.id,   id,   "dispatch id is not equal to"           + str(id) )