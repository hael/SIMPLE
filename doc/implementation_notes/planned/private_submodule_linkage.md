# Private helpers called from tester submodules

## Observed failure

After the test-submodule refactor in `a307152e0`, the maintainer reports that
macOS builds successfully. The supplied Linux build output shows the production
executables linking successfully, then `simple_test_exec` failing to link.
There are 30 distinct unresolved module-procedure symbols across 15 parent
modules, all private helpers called by their descendant tester submodules.

GFortran has documented private-parent/submodule linkage defects:
[upstream description](https://gcc.gnu.org/pipermail/gcc-cvs/2026-August/460346.html).
The Linux compiler version was not included in the supplied output. No
compiler-version diagnosis or successful Linux rebuild is claimed.

## Contract and implementation

Keep the helpers private and retain the descendant test submodules. Declare
each affected helper in a private `module` procedure interface in its parent,
and give its existing definition the matching `module` prefix. Separate module
procedures may be implemented in the parent itself; the implementations stay
in their current files. Preserve the argument/result characteristics, purity,
elemental attributes, implementation bodies and existing public API lists.

The patch applies this declaration pattern to every distinct symbol in the
Linux linker output. It introduces no wrappers, extra source files, compiler
options, CMake targets or numerical changes. Existing tester calls exercise
the linkage on the next build; no additional numerical test is needed.

## Review and validation

- Source audit confirmed coverage of all 30 unresolved symbols and verified
  that implementation bodies and public API lists match `HEAD` after removing
  the added declaration prefixes.
- `git diff --check` passed.
- No build or test was run by the agent, following the maintainer's instruction.
- Outstanding: rebuild and link `simple_test_exec` on Linux, then run the
  normal fast gate. The Mac rebuild and Linux test results remain unobserved.
