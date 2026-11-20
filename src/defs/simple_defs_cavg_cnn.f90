!@descr: constants of the program cavg_cnn_predict: its table and the class-average CNN's groups
module simple_defs_cavg_cnn
implicit none
private

public :: CAVG_CNN_PRED_TABLE, CAVG_CNN_GROUP_NAMES, CAVG_CNN_GROUP_EMPTY

character(len=*),  parameter :: CAVG_CNN_PRED_TABLE     = 'cavg_cnn_predictions.txt'
!> the network's groups in the order of its probabilities, then a class with no particles
character(len=19), parameter :: CAVG_CNN_GROUP_NAMES(4) = &
    &[character(len=19) :: 'junk', 'secondary_structure', 'likely_protein', 'empty']
integer,           parameter :: CAVG_CNN_GROUP_EMPTY    = 4

end module simple_defs_cavg_cnn
