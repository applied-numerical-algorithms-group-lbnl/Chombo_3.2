set print pretty on
set print static-members off
set print array
set print asm-demangle
set breakpoint pending on
break exit
break abort
set auto-load safe-path /
break Chombo::ResistivityOp::slowGSRB
break Chombo::ResistivityOp::getFlux

