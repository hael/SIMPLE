# Full Refactor View

This folder is now a clean design pack for the refactor strategy we converged on.

The strategy is:

- keep `params%box` style flat access
- remove `init_strings` by moving defaults onto declarations
- replace the giant `check_*arg` wall with generic parsing over type-based registries
- keep derivation and validation split by semantic phase rather than by storage type

## Files

- [simple_parameters_type_reference.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_type_reference.f90:1)
  Extracted current type/inventory reference from the live
  [`src/main/simple_parameters.f90`](/Users/elmlundho/src/SIMPLE/src/main/simple_parameters.f90:1).
  This is the real surface area we are refactoring.

- [simple_parameters_refactor_module_view.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_refactor_module_view.f90:1)
  Shows how the top-level module and `type(parameters)` would be organized after the refactor.

- [simple_parameters_refactor_bindings_view.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_refactor_bindings_view.f90:1)
  Shows how the giant parsing wall gets replaced by type-based binding procedures.

- [simple_parameters_refactor_phases_view.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_refactor_phases_view.f90:1)
  Shows how the current derivation and validation dump gets split into named semantic phases.

## What to look at first

If you want the shortest path to judging the design:

1. inspect [simple_parameters_type_reference.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_type_reference.f90:1) to re-ground on the real size
2. inspect [simple_parameters_refactor_module_view.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_refactor_module_view.f90:1) for the top-level shape
3. inspect [simple_parameters_refactor_bindings_view.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_refactor_bindings_view.f90:1) to see how `check_*arg` disappears
4. inspect [simple_parameters_refactor_phases_view.f90](/Users/elmlundho/src/SIMPLE/doc/refactoring_notes/parameters_prototype/simple_parameters_refactor_phases_view.f90:1) to see how the semantic code dump gets structured

## Intended edit workflow

With this strategy, a normal new input parameter would mean:

1. add the field once in the flat type declaration
2. add one binding in the matching storage-type binding procedure
3. only touch a semantic phase file if the parameter actually has custom behavior

That directly addresses the three pains you called out:

- no `init_strings`
- no scattered `check_*arg` edits in the constructor
- no monolithic derivation/checking block
