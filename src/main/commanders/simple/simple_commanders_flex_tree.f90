!@descr: flex_pca mode=tree: the FLEX state tree, driving flex_pca and refine3D per node
module simple_commanders_flex_tree
use simple_commanders_api,      only: builder, cmdline, commander_base, parameters, simple_end
use simple_commanders_flex_pca, only: commander_flex_pca, apply_flex_pca_defaults
implicit none

public :: commander_flex_tree
private

type, extends(commander_base) :: commander_flex_tree
  contains
    procedure :: execute => exec_flex_tree
end type commander_flex_tree

contains

    subroutine exec_flex_tree( self, cline )
        use simple_flex_tree_strategy,  only: flex_tree_strategy
        use simple_commanders_refine3D, only: commander_refine3D
        class(commander_flex_tree), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(flex_tree_strategy) :: tree
        type(commander_flex_pca) :: xflex
        type(commander_refine3D) :: xrefine3D
        type(parameters)         :: params
        type(builder)            :: build
        call apply_flex_pca_defaults(cline)
        call tree%initialize(params, build, cline)
        call tree%execute(params, cline, xflex, xrefine3D)
        call tree%cleanup
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_PCA TREE NORMAL STOP ****')
    end subroutine exec_flex_tree

end module simple_commanders_flex_tree
