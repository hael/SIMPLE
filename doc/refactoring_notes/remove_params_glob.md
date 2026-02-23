# Refactoring: Remove params_glob Global State
**Date:** February 22, 2026  
**Status:** Complete  
**Scope:** SIMPLE execution contexts  

## Executive Summary
Eliminated module-level global variable `params_glob` from the library by refactoring procedures to explicitly pass `class(parameters)` instances as arguments. This improves code safety, testability, and eliminates hidden dependencies that complicate distributed execution and error diagnostics.

### Global Dependency Issues
- **`params_glob`** was a module-level global in `simple_parameters` module, used everywhere
- Hidden dependency made code difficult to trace and refactor
- Violated encapsulation principles
- Complicated distributed execution where multiple parameter contexts exist simultaneously
- Created implicit couplings between procedures

## Validation Checklist

- [x] All procedure signatures updated consistently
- [x] All call sites passing `params` argument explicitly
- [x] No remaining `params_glob` references

## Benefits

### Safety
✅ **Eliminates hidden global state** — All parameter dependencies explicit and traceable  
✅ **Prevents race conditions** — Multiple concurrent parameter contexts won't interfere  
✅ **Improves testability** — Procedures can be tested with specific parameter configurations  

### Maintainability
✅ **Code clarity** — Dependencies visible in function signatures  
✅ **Refactoring safety** — Grep/IDE tools can identify all affected call sites  
✅ **Future extensibility** — New parameter versions can coexist during migration  

### For Future Development
- All new procedures should follow the explicit parameter-passing pattern
- Do not introduce new global state for parameter management
- Test with multiple parameter contexts to detect hidden dependencies

## Conclusion
This refactoring eliminates a significant source of hidden complexity in the SIMPLE library. By threading parameters explicitly, we gain:
- **Transparency** — All dependencies visible in code
- **Safety** — No implicit global state pollution
- **Debuggability** — Better error traces and diagnostic output
The fix enables safer concurrent execution contexts and provides a foundation for future distributed execution improvements.
