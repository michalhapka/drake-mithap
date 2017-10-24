module precision
use,intrinsic :: iso_fortran_env, only : int16, int64, real64, real128
implicit none

private int16, int64, real64, real128

integer,parameter :: shortint = int16
integer,parameter :: longint  = int64
integer,parameter :: dble     = real64
integer,parameter :: prec     = 16 !real128

integer,parameter :: log_size  = storage_size(.false.)/8
integer,parameter :: int_size  = storage_size(0)/8
integer,parameter :: real_size = storage_size(0._prec)/8

end module precision
