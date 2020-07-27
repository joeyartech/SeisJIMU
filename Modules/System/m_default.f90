!public entities naming convention
! i : integer
! r : single precision real number (binary32)
! d : double precision real number (binary64)
! c : single precision complex number
! z : double precision complex number
! s : character, string
! m : module
! t : type, class
! o : optional
! is, if : logical

module m_default

!default
parameter :: i_str_len  = 80   !string length
parameter :: i_str_xlen = 256  !string length
parameter :: i_str_xxlen = 512  !string length

!default for CPML
parameter :: r_k_x = 1.
parameter :: r_k_y = 1.
parameter :: r_k_z = 1.
parameter :: r_npower = 2.
parameter :: r_Rcoef = 0.001  
parameter :: r_kappa_max=7. !increase this number if absorbing is not satisfactory at grazing incident angle
                            !(make in/outside PML reflection more separate..)
                            !decrease this number if grid dispersion is not satisfactory

!constant
parameter :: r_pi = 3.1415927
parameter :: d_pi = 3.1415926535897932d0
character(*), parameter :: s_return = achar(13)
    !https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/494946
    !char function: CHAR(10) is linefeed LF. CHAR(13) is carriage return CR. If you are a little paranoid, ACHAR(10) is better - this is a little more robust to the default character kind not being ascii.
    !The NEW_LINE standard intrinsic is even more robust. There's also the C_NEW_LINE and C_CARRIAGE_RETURN constants from the ISO_C_BINDING module for kind C_CHAR characters.

!compilation info
character(*),parameter :: s_commit = commit
character(*),parameter :: s_branch = branch

#ifdef GNU
character(*),parameter :: s_compiler = 'gfortran v' // __VERSION__
character(*),parameter :: s_version = __VERSION__
integer,parameter :: i_endian = __BYTE_ORDER__
#endif

#ifdef INTEL
character(*),parameter :: s_compiler = 'ifort'
integer,parameter :: s_version(2) = [__INTEL_COMPILER, __INTEL_COMPILER_UPDATE]
character(*),parameter :: i_endian = 'unavailable'
#endif

end