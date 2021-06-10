!prefix nomenclature
! i_ : integer,index
! r_ : single precision real number (kind=4)
! d_ : double precision real number (kind=8)
! c_ : single precision complex number (kind=4)
! z_ : double precision complex number (kind=8)
! s_ : character, string, string array
! m_ : module
! t_ : type, class
! o_ : optional
! is_,if_ : logical
! i-l: index
! n,m : number of

module m_global

    !default string length
    integer,parameter :: i_str_slen  = 64
    integer,parameter :: i_str_len  = 128
    integer,parameter :: i_str_xlen = 256
    integer,parameter :: i_str_xxlen = 512
    integer,parameter :: i_str_xxxlen = 1024
    character(*), parameter :: s_return = achar(13)
        !https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/494946
        !char function: CHAR(10) is linefeed LF. CHAR(13) is carriage return CR. If you are a little paranoid, ACHAR(10) is better - this is a little more robust to the default character kind not being ascii.
        !The NEW_LINE standard intrinsic is even more robust. There's also the C_NEW_LINE and C_CARRIAGE_RETURN constants from the ISO_C_BINDING module for kind C_CHAR characters.

    !math constants
    real,parameter :: r_pi = 3.1415927
    double precision,parameter :: d_pi = 3.1415926535897932d0
    !real,parameter :: r_epsilon = epsilon(r_pi)
    real,parameter :: r_epsilon=1e-3
    
    !compilation info
    character(*),parameter :: s_commit = commit
    character(*),parameter :: s_branch = branch

#ifdef GNU
    character(*),parameter :: s_compiler = 'gfortran v' // __VERSION__
    character(*),parameter :: s_version = __VERSION__
    integer,parameter      :: i_endian = __BYTE_ORDER__
#endif

#ifdef INTEL
    character(*),parameter :: s_compiler = 'ifort'
    integer,parameter      :: s_version(2) = [__INTEL_COMPILER, __INTEL_COMPILER_UPDATE]
    character(*),parameter :: i_endian = 'unavailable'
#endif

end