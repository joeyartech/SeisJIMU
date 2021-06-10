module m_format_su
use m_string

!A Fortran version of su format (defined in /SeisUnix/include/segy.h)
!
!Note 1: su header is very similar to segy trace header:
!    BYTES#   1-180: same, all integer types
!    BYTES# 181-240: different. cwp assignment in su whereas unassigned in segy (however they're later defined as e.g. inline, xline numbers in segy rev 1)
!    in addition, su format does not have binary header as in segy format.
!
!Note 2: One potential conflict is that su and segy trace header has unsigned short C type, which is not and will not be provided in Fortran standards.
!    This happens to the two keywords: ns, dt.
!    If the two keywords have values > 32767, then after reading Fortran will provide negative values, by this formula:
!        read_value = true_unsigned_short_value - 32768*2,
!    and we will need integer(kind=4) type to save the true value (integer(kind=2) type does not have more mem space).
!    Fortrunately, ns and dt is usually <= 32767 in practice so we are not concerned with this issue.
!
!Note 3: Most of the keywords are not needed in SeisJIMU..
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                         !!
! 
! !! segy.h - include file for SEGY traces
!  *
!  * declarations for:
!  *	typedef struct {} segy - the trace identification header
!  *	typedef struct {} bhed - binary header
!  *
!  * Note:
!  *	If header words are added, run the makefile in this directory
!  *	to recreate hdr.h.
!  *
!  * Reference:
!  *	K. M. Barry, D. A. Cavers and C. W. Kneale, "Special Report:
!  *		Recommended Standards for Digital Tape Formats",
!  *		Geophysics, vol. 40, no. 2 (April 1975), P. 344-352.
!  *	
!  * $Author: john $
!  * $Source: /usr/local/cwp/src/su/include/RCS/segy.h,v $
!  * $Revision: 1.33 $ ; $Date: 2011/11/11 23:56:14 $
!  !! 

    private

    type,public :: t_header
        sequence
        !BYTE# 1-4
        integer(4) :: tracl=0   !! Trace sequence number within line
                                ! ---numbers continue to increase if the
                                ! same line continues across multiple
                                ! SEG Y files. !!
        !BYTE# 5-8
        integer(4) :: tracr=0   !! Trace sequence number within SEG Y file
                                ! ---each file starts with trace sequence one !!
        !BYTE# 9-12
        integer(4) :: fldr=0    !! Original field record number !!
        !BYTE# 13-16
        integer(4) :: tracf=0   !! Trace number within original field record !!
        !BYTE# 17-20
        integer(4) :: ep=0      !! energy source point number
                                ! ---Used when more than one record occurs
                                ! at the same effective surface location. !!
        !BYTE# 21-24
        integer(4) :: cdp=0     !! Ensemble number (i.e. CDP, CMP, CRP,...) !!
        !BYTE# 25-28
        integer(4) :: cdpt=0    !! trace number within the ensemble
                                ! ---each ensemble starts with trace number one. !!
        !BYTE# 29-30
        integer(2) :: trid=0    !! trace identification code:
                                !  -1 = Other
                                !   0 = Unknown
                                !   1 = Seismic data
                                !   2 = Dead
                                !   3 = Dummy
                                !   4 = Time break
                                !   5 = Uphole
                                !   6 = Sweep
                                !   7 = Timing
                                !   8 = Water break
                                !   9 = Near-field gun signature
                                !  10 = Far-field gun signature
                                !  11 = Seismic pressure sensor
                                !  12 = Multicomponent seismic sensor
                                !          - Vertical component
                                !  13 = Multicomponent seismic sensor
                                !          - Cross-line component 
                                !  14 = Multicomponent seismic sensor
                                !          - in-line component 
                                !  15 = Rotated multicomponent seismic sensor
                                !          - Vertical component
                                !  16 = Rotated multicomponent seismic sensor
                                !          - Transverse component
                                !  17 = Rotated multicomponent seismic sensor
                                !          - Radial component
                                !  18 = Vibrator reaction mass
                                !  19 = Vibrator baseplate
                                !  20 = Vibrator estimated ground force
                                !  21 = Vibrator reference
                                !  22 = Time-velocity pairs
                                !  23 ... N = optional use 
                                !          (maximum N = 32,767)

                                ! Following are CWP id flags:

                                ! 109 = autocorrelation
                                ! 110 = Fourier transformed - no packing
                                !      xr[0],xi[0], ..., xr[N-1],xi[N-1]
                                ! 111 = Fourier transformed - unpacked Nyquist
                                !      xr[0],xi[0],...,xr[N/2],xi[N/2]
                                ! 112 = Fourier transformed - packed Nyquist
                                !      even N:
                                !      xr[0],xr[N/2],xr[1],xi[1], ...,
                                !         xr[N/2 -1],xi[N/2 -1]
                                !         (note the exceptional second entry)
                                !      odd N:
                                !      xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
                                !         xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
                                !         (note the exceptional second & last entries)
                                ! 113 = Complex signal in the time domain
                                !      xr[0],xi[0], ..., xr[N-1],xi[N-1]
                                ! 114 = Fourier transformed - amplitude/phase
                                !      a[0],p[0], ..., a[N-1],p[N-1]
                                ! 115 = Complex time signal - amplitude/phase
                                !      a[0],p[0], ..., a[N-1],p[N-1]
                                ! 116 = real(4) part of complex trace from 0 to Nyquist
                                ! 117 = Imag part of complex trace from 0 to Nyquist
                                ! 118 = Amplitude of complex trace from 0 to Nyquist
                                ! 119 = Phase of complex trace from 0 to Nyquist
                                ! 121 = Wavenumber time domain (k-t)
                                ! 122 = Wavenumber frequency (k-omega)
                                ! 123 = Envelope of the complex time trace
                                ! 124 = Phase of the complex time trace
                                ! 125 = Frequency of the complex time trace
                                ! 130 = Depth-Range (z-x) traces
                                ! 201 = Seismic data packed to BYTE# (by supack1)
                                ! 202 = Seismic data packed to 2 BYTE# (by supack2) !!
        !BYTE# 31-32
        integer(2) :: nvs=0     !! Number of vertically summed traces yielding
                                ! this trace. (1 is one trace, 2 is two summed traces, etc.) !!
        !BYTE# 33-34
        integer(2) :: nhs=0     !! Number of horizontally summed traces yielding
                                ! this trace. (1 is one trace, 2 is two summed traces, etc.) !!
        !BYTE# 35-36
        integer(2) :: duse=0    !! Data use:
                                ! 1 = Production
                                ! 2 = Test !!
        !BYTE# 37-40
        integer(4) :: offset=0  !! Distance from the center of the source point 
                                ! to the center of the receiver group 
                                ! (negative if opposite to direction in which 
                                ! the line was shot). !!
        !BYTE# 41-44
        integer(4) :: gelev=0   !! Receiver group elevation from sea level
                                ! (all elevations above the Vertical datum are 
                                ! positive and below are negative). !!
        !BYTE# 45-48
        integer(4) :: selev=0   !! Surface elevation at source. !!
        !BYTE# 49-52
        integer(4) :: sdepth=0  !! Source depth below surface (a positive number). !!
        !BYTE# 53-56
        integer(4) :: gdel=0    !! Datum elevation at receiver group. !!
        !BYTE# 57-60
        integer(4) :: sdel=0    !! Datum elevation at source. !!
        !BYTE# 61-64
        integer(4) :: swdep=0   !! Water depth at source. !!
        !BYTE# 65-68
        integer(4) :: gwdep=0   !! Water depth at receiver group. !!
        !BYTE# 69-70
        integer(2) :: scalel=0  !! Scalar to be applied to the previous 7 entries
                                ! to give the real(4) value. 
                                ! Scalar = 1, +10, +100, +1000, +10000.
                                ! If positive, scalar is used as a multiplier,
                                ! if negative, scalar is used as a divisor. !!
        !BYTE# 71-72
        integer(2) :: scalco=0  !! Scalar to be applied to the next 4 entries
                                ! to give the real(4) value. 
                                ! Scalar = 1, +10, +100, +1000, +10000.
                                ! If positive, scalar is used as a multiplier,
                                ! if negative, scalar is used as a divisor. !!
        !BYTE# 73-76
        integer(4) :: sx=0      !! Source coordinate - X !!
        !BYTE# 77-80
        integer(4) :: sy=0      !! Source coordinate - Y !!
        !BYTE# 81-84
        integer(4) :: gx=0      !! Group coordinate - X !!
        !BYTE# 85-88
        integer(4) :: gy=0      !! Group coordinate - Y !!
        !BYTE# 89-90
        integer(2) :: counit=0  !! Coordinate units: (for previous 4 entries and
                                ! for the 7 entries before scalel)
                                ! 1 = Length (meters or feet)
                                ! 2 = Seconds of arc
                                ! 3 = Decimal degrees
                                ! 4 = Degrees, minutes, seconds (DMS)

                                ! In case 2, the X values are longitude and 
                                ! the Y values are latitude, a positive value designates
                                ! the number of seconds east of Greenwich
                                !         or north of the equator
         
                                ! In case 4, to encode +-DDDMMSS
                                ! counit = +-DDD*10^4 + MM*10^2 + SS,
                                ! with scalco = 1. To encode +-DDDMMSS.ss
                                ! counit = +-DDD*10^6 + MM*10^4 + SS*10^2 
                                ! with scalco = -100. !!
        !BYTE# 91-92
        integer(2) :: wevel=0   !! Weathering velocity. !!
        !BYTE# 93-94
        integer(2) :: swevel=0  !! Subweathering velocity. !!
        !BYTE# 95-96
        integer(2) :: sut=0     !! Uphole time at source in milliseconds. !!
        !BYTE# 97-98
        integer(2) :: gut=0     !! Uphole time at receiver group in milliseconds. !!
        !BYTE# 99-100
        integer(2) :: sstat=0   !! Source static correction in milliseconds. !!
        !BYTE# 101-102
        integer(2) :: gstat=0   !! Group static correction  in milliseconds.!!
        !BYTE# 103-104
        integer(2) :: tstat=0   !! Total static applied  in milliseconds.
                                ! (Zero if no static has been applied.) !!
        !BYTE# 105-106
        integer(2) :: laga=0    !! Lag time A, time in ms between end of 240-
                                ! byte trace identification header and time
                                ! break, positive if time break occurs after
                                ! end of header, time break is defined as
                                ! the initiation pulse which maybe recorded
                                ! on an auxiliary trace or as otherwise
                                ! specified by the recording system !!
        !BYTE# 107-108
        integer(2) :: lagb=0    !! lag time B, time in ms between the time break
                                ! and the initiation time of the energy source,
                                ! may be positive or negative !!
        !BYTE# 109-110
        integer(2) :: delrt=0   !! delay recording time, time in ms between
                                ! initiation time of energy source and time
                                ! when recording of data samples begins
                                ! (for deep water work if recording does not
                                ! start at zero time) !!
        !BYTE# 111-112
        integer(2) :: muts=0    !! mute time--start !!
        !BYTE# 113-114
        integer(2) :: mute=0    !! mute time--end !!
        !BYTE# 115-116
        integer(2) :: ns=0      !! number of samples in this trace !!
        !BYTE# 117-118
        integer(2) :: dt=0      !! sample interval in micro-seconds
        !BYTE# 119-120
        integer(2) :: gain=0    !! gain type of field instruments code:
                                ! 1 = fixed
                                ! 2 = binary
                                ! 3 = floating point
                                ! 4 ---- N = optional use !!
        !BYTE# 121-122
        integer(2) :: igc=0     !! instrument gain constant !!
        !BYTE# 123-124
        integer(2) :: igi=0     !! instrument early or initial gain !!
        !BYTE# 125-126
        integer(2) :: corr=0    !! correlated:
                                ! 1 = no
                                ! 2 = yes !!
        !BYTE# 127-128
        integer(2) :: sfs=0     !! sweep frequency at start !!
        !BYTE# 129-130
        integer(2) :: sfe=0     !! sweep frequency at end !!
        !BYTE# 131-132
        integer(2) :: slen=0    !! sweep length in ms !!
        !BYTE# 133-134
        integer(2) :: styp=0    !! sweep type code:
                                ! 1 = linear
                                ! 2 = cos-squared
                                ! 3 = other !!
        !BYTE# 135-136
        integer(2) :: stas=0    !! sweep trace length at start in ms !!
        !BYTE# 137-138
        integer(2) :: stae=0    !! sweep trace length at end in ms !!
        !BYTE# 139-140
        integer(2) :: tatyp=0   !! taper type: 1=linear, 2=cos^2, 3=other !!
        !BYTE# 141-142
        integer(2) :: afilf=0   !! alias filter frequency if used !!
        !BYTE# 143-144
        integer(2) :: afils=0   !! alias filter slope !!
        !BYTE# 145-146
        integer(2) :: nofilf=0  !! notch filter frequency if used !!
        !BYTE# 147-148
        integer(2) :: nofils=0  !! notch filter slope !!
        !BYTE# 149-150
        integer(2) :: lcf=0     !! low cut frequency if used !!
        !BYTE# 151-152
        integer(2) :: hcf=0     !! high cut frequncy if used !!
        !BYTE# 153-154
        integer(2) :: lcs=0     !! low cut slope !!
        !BYTE# 155-156
        integer(2) :: hcs=0     !! high cut slope !!
        !BYTE# 157-158
        integer(2) :: year=0    !! year data recorded !!
        !BYTE# 159-160
        integer(2) :: day=0     !! day of year !!
        !BYTE# 161-162
        integer(2) :: hour=0    !! hour of day (24 hour clock) !!
        !BYTE# 163-164
        integer(2) :: minute=0  !! minute of hour !!
        !BYTE# 165-166
        integer(2) :: sec=0     !! second of minute !!
        !BYTE# 167-168
        integer(2) :: timbas=0  !! time basis code:
                                ! 1 = local
                                ! 2 = GMT
                                ! 3 = other !!
        !BYTE# 169-170
        integer(2) :: trwf=0    !! trace weighting factor, defined as 1/2^N
                                ! volts for the least sigificant bit !!
        !BYTE# 171-172
        integer(2) :: grnors=0  !! geophone group number of roll switch
                                ! position one !!
        !BYTE# 173-174
        integer(2) :: grnofr=0  !! geophone group number of trace one within
                                ! original field record !!
        !BYTE# 175-176
        integer(2) :: grnlof=0  !! geophone group number of last trace within
                                ! original field record !!
        !BYTE# 177-178
        integer(2) :: gaps=0    !! gap size (total number of groups dropped) !!
        !BYTE# 179-180
        integer(2) :: otrav=0   !! overtravel taper code:
                                ! 1 = down (or behind)
                                ! 2 = up (or ahead) !!

        !! cwp local assignments !!
        !BYTE# 181-184
        real(4) :: d1=0.        !! sample spacing for non-seismic data !!
        !BYTE# 185-188
        real(4) :: f1=0.        !! first sample location for non-seismic data !!
        !BYTE# 189-192
        real(4) :: d2=0.        !! sample spacing between traces !!
        !BYTE# 193-196
        real(4) :: f2=0.        !! first trace location !!
        !BYTE# 197-200
        real(4) :: ungpow=0.    !! negative of power used for dynamic
                                ! range compression !!
        !BYTE# 201-204
        real(4) :: unscale=0.   !! reciprocal of scaling factor to normalize
                                ! range !!

        !BYTE# 205-206
        integer(2) :: mark=0    !! mark selected traces !!

        !! SLTSU local assignments !! 
        !BYTE# 207-208
        integer(2) :: mutb=0    !! mute time at bottom (start time)
                                ! bottom mute ends at last sample   !!
        !BYTE# 209-212
        real(4) :: dz=0.        !! depth sampling interval in (m or ft)
                                ! if =0.0, input are time samples       !!
        !BYTE# 213-216
        real(4) :: fz=0.        !! depth of first sample in (m or ft)  !!
        !BYTE# 217-218
        integer(2) :: n2=0      !! number of traces per cdp or per shot !!
        !BYTE# 219-220
        integer(2) :: shortpad=0!! alignment padding !!
        !BYTE# 221-224
        integer(4) :: ntr=0     !! number of traces !!

        !! SLTSU local assignments end !! 
        !BYTE# 225-240
        integer(2),dimension(8) :: unass=0  !! unassigned !!
        !end of su header
        
    end type
        
    type,public :: t_format_su
        type(t_header),dimension(:),allocatable :: hdrs
        real(4),dimension(:,:),allocatable :: trs
        integer :: ntr !number of traces
        integer :: ns !assume all traces have same number of samples
        real    :: dt !assume all traces have same time sampling interval
        contains
        procedure :: read => read
        procedure :: write => write
        final     :: finalize
    end type  

    contains
    
    subroutine read(self,file,o_sindex)
        class(t_format_su) :: self
        character(*) :: file
        character(*),optional :: o_sindex

        integer file_size
        integer(2) :: dt
        
        inquire(file=file, size=file_size) !file_size in bytes
        file_size=file_size/4 !file_size in sizeof(float), single precision
        
        open(11,file=file,action='read',access='stream')
        
        read(11,pos=115) ns  !get number of samples from 1st trace
        if(ns<0) ns=ns+32768*2  !ns keyword is unsigned short in C
        self%ns=ns

        read(11,pos=117) dt
        if(dt<0) dt=dt+32768*2  !dt keyword is unsigned short in C
        self%dt=dt/1e-6
        
        close(11)
        
        self%ntr=file_size/(self%ns+60)

        if(self%ntr*(self%ns+60)/=file_size) call error('ntr*(ns+60) /= file_size. Possibly ns is not same for all traces..')

        if(present(o_sindex)) write(*,*) 'Shot# '//o_sindex//' will read '//num2str(self%ntr)//' traces, each trace has '//num2str(self%ns)//' samples.'
        
        allocate(self%hdrs(self%ntr))
        allocate(self%trs(self%ns,self%ntr))
        
        open(11,file=file,action='read',access='stream') !access='direct',recl=4*(ns+60))
        read(11) (self%hdrs(i),self%trs(:,i), i=1,self%ntr)
        close(11)
        
        call hud('Data read success.')
        
    end subroutine
    
    subroutine write(self,file,o_sindex)
        class(t_format_su) :: self
        character(*) :: file
        character(*),optional :: o_sindex

        if(present(o_sindex)) write(*,*) 'Shot# '//o_sindex//' will write '//num2str(self%ntr)//' traces, each trace has '//num2str(self%ns)//' samples.'

        open(12,file=file,action='read',access='stream') !access='direct',recl=4*(ns+60))
        write(12) (self%hdrs(i),self%trs(:,i), i=1,self%ntr)
        close(12)
        
        call hud('Data write success.')

    end subroutine
    
    subroutine finalize(self)
        type(t_format_su) :: self
        deallocate(self%hdrs,self%trs)
    end subroutine

end