module m_suformat
use m_sysio

    private
    public :: t_suformat, init_suheader, read_sudata

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
!Note 3: Most of the keywords are not needed in lego..
!
! !! Copyright (c) Colorado School of Mines, 2011.!!
! !! All rights reserved.                       !!
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


    type t_suhdr
        sequence
        !BYTE# 1-4
        integer(4) tracl   !! Trace sequence number within line
!                             ---numbers continue to increase if the
!                             same line continues across multiple
!                             SEG Y files. !!
        !BYTE# 5-8
        integer(4) tracr   !! Trace sequence number within SEG Y file
!                             ---each file starts with trace sequence one !!
        !BYTE# 9-12
        integer(4) fldr    !! Original field record number !!
        !BYTE# 13-16
        integer(4) tracf   !! Trace number within original field record !!
        !BYTE# 17-20
        integer(4) ep      !! energy source point number
!                             ---Used when more than one record occurs
!                             at the same effective surface location. !!
        !BYTE# 21-24
        integer(4) cdp     !! Ensemble number (i.e. CDP, CMP, CRP,...) !!
        !BYTE# 25-28
        integer(4) cdpt    !! trace number within the ensemble
!                           ---each ensemble starts with trace number one. !!
        !BYTE# 29-30
        integer(2) trid    !! trace identification code:
!                              -1 = Other
!                               0 = Unknown
!                               1 = Seismic data
!                               2 = Dead
!                               3 = Dummy
!                               4 = Time break
!                               5 = Uphole
!                               6 = Sweep
!                               7 = Timing
!                               8 = Water break
!                               9 = Near-field gun signature
!                              10 = Far-field gun signature
!                              11 = Seismic pressure sensor
!                              12 = Multicomponent seismic sensor
!                                      - Vertical component
!                              13 = Multicomponent seismic sensor
!                                      - Cross-line component 
!                              14 = Multicomponent seismic sensor
!                                      - in-line component 
!                              15 = Rotated multicomponent seismic sensor
!                                      - Vertical component
!                              16 = Rotated multicomponent seismic sensor
!                                      - Transverse component
!                              17 = Rotated multicomponent seismic sensor
!                                      - Radial component
!                              18 = Vibrator reaction mass
!                              19 = Vibrator baseplate
!                              20 = Vibrator estimated ground force
!                              21 = Vibrator reference
!                              22 = Time-velocity pairs
!                              23 ... N = optional use 
!                                      (maximum N = 32,767)
!
!                             Following are CWP id flags:
!
!                             109 = autocorrelation
!                             110 = Fourier transformed - no packing
!                                  xr[0],xi[0], ..., xr[N-1],xi[N-1]
!                             111 = Fourier transformed - unpacked Nyquist
!                                  xr[0],xi[0],...,xr[N/2],xi[N/2]
!                             112 = Fourier transformed - packed Nyquist
!                                  even N:
!                                  xr[0],xr[N/2],xr[1],xi[1], ...,
!                                     xr[N/2 -1],xi[N/2 -1]
!                                     (note the exceptional second entry)
!                                  odd N:
!                                  xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
!                                     xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
!                                     (note the exceptional second & last entries)
!                             113 = Complex signal in the time domain
!                                  xr[0],xi[0], ..., xr[N-1],xi[N-1]
!                             114 = Fourier transformed - amplitude/phase
!                                  a[0],p[0], ..., a[N-1],p[N-1]
!                             115 = Complex time signal - amplitude/phase
!                                  a[0],p[0], ..., a[N-1],p[N-1]
!                             116 = real(4) part of complex trace from 0 to Nyquist
!                             117 = Imag part of complex trace from 0 to Nyquist
!                             118 = Amplitude of complex trace from 0 to Nyquist
!                             119 = Phase of complex trace from 0 to Nyquist
!                             121 = Wavenumber time domain (k-t)
!                             122 = Wavenumber frequency (k-omega)
!                             123 = Envelope of the complex time trace
!                             124 = Phase of the complex time trace
!                             125 = Frequency of the complex time trace
!                             130 = Depth-Range (z-x) traces
!                             201 = Seismic data packed to BYTE# (by supack1)
!                             202 = Seismic data packed to 2 BYTE# (by supack2) !!
        !BYTE# 31-32
        integer(2) nvs     !! Number of vertically summed traces yielding
!                             this trace. (1 is one trace, 2 is two summed traces, etc.) !!
        !BYTE# 33-34
        integer(2) nhs     !! Number of horizontally summed traces yielding
!                             this trace. (1 is one trace, 2 is two summed traces, etc.) !!
        !BYTE# 35-36
        integer(2) duse    !! Data use:
!                             1 = Production
!                             2 = Test !!
        !BYTE# 37-40
        integer(4) offset  !! Distance from the center of the source point 
!                             to the center of the receiver group 
!                             (negative if opposite to direction in which 
!                             the line was shot). !!
        !BYTE# 41-44
        integer(4) gelev   !! Receiver group elevation from sea level
!                             (all elevations above the Vertical datum are 
!                             positive and below are negative). !!
        !BYTE# 45-48
        integer(4) selev   !! Surface elevation at source. !!
        !BYTE# 49-52
        integer(4) sdepth  !! Source depth below surface (a positive number). !!
        !BYTE# 53-56
        integer(4) gdel    !! Datum elevation at receiver group. !!
        !BYTE# 57-60
        integer(4) sdel    !! Datum elevation at source. !!
        !BYTE# 61-64
        integer(4) swdep   !! Water depth at source. !!
        !BYTE# 65-68
        integer(4) gwdep   !! Water depth at receiver group. !!
        !BYTE# 69-70
        integer(2) scalel  !! Scalar to be applied to the previous 7 entries
!                             to give the real(4) value. 
!                             Scalar = 1, +10, +100, +1000, +10000.
!                             If positive, scalar is used as a multiplier,
!                             if negative, scalar is used as a divisor. !!
        !BYTE# 71-72
        integer(2) scalco  !! Scalar to be applied to the next 4 entries
!                             to give the real(4) value. 
!                             Scalar = 1, +10, +100, +1000, +10000.
!                             If positive, scalar is used as a multiplier,
!                             if negative, scalar is used as a divisor. !!
        !BYTE# 73-76
        integer(4) sx      !! Source coordinate - X !!
        !BYTE# 77-80
        integer(4) sy      !! Source coordinate - Y !!
        !BYTE# 81-84
        integer(4) gx      !! Group coordinate - X !!
        !BYTE# 85-88
        integer(4) gy      !! Group coordinate - Y !!
        !BYTE# 89-90
        integer(2) counit  !! Coordinate units: (for previous 4 entries and
!                             for the 7 entries before scalel)
!                             1 = Length (meters or feet)
!                             2 = Seconds of arc
!                             3 = Decimal degrees
!                             4 = Degrees, minutes, seconds (DMS)
!
!                             In case 2, the X values are longitude and 
!                             the Y values are latitude, a positive value designates
!                             the number of seconds east of Greenwich
!                                     or north of the equator
!      
!                             In case 4, to encode +-DDDMMSS
!                             counit = +-DDD*10^4 + MM*10^2 + SS,
!                             with scalco = 1. To encode +-DDDMMSS.ss
!                             counit = +-DDD*10^6 + MM*10^4 + SS*10^2 
!                             with scalco = -100. !!
        !BYTE# 91-92
        integer(2) wevel   !! Weathering velocity. !!
        !BYTE# 93-94
        integer(2) swevel  !! Subweathering velocity. !!
        !BYTE# 95-96
        integer(2) sut     !! Uphole time at source in milliseconds. !!
        !BYTE# 97-98
        integer(2) gut     !! Uphole time at receiver group in milliseconds. !!
        !BYTE# 99-100
        integer(2) sstat   !! Source static correction in milliseconds. !!
        !BYTE# 101-102
        integer(2) gstat   !! Group static correction  in milliseconds.!!
        !BYTE# 103-104
        integer(2) tstat   !! Total static applied  in milliseconds.
!                             (Zero if no static has been applied.) !!
        !BYTE# 105-106
        integer(2) laga    !! Lag time A, time in ms between end of 240-
!                             byte trace identification header and time
!                             break, positive if time break occurs after
!                             end of header, time break is defined as
!                             the initiation pulse which maybe recorded
!                             on an auxiliary trace or as otherwise
!                             specified by the recording system !!
        !BYTE# 107-108
        integer(2) lagb    !! lag time B, time in ms between the time break
!                             and the initiation time of the energy source,
!                             may be positive or negative !!
        !BYTE# 109-110
        integer(2) delrt   !! delay recording time, time in ms between
!                             initiation time of energy source and time
!                             when recording of data samples begins
!                             (for deep water work if recording does not
!                             start at zero time) !!
        !BYTE# 111-112
        integer(2) muts    !! mute time--start !!
        !BYTE# 113-114
        integer(2) mute    !! mute time--end !!
        !BYTE# 115-116
        integer(2) ns      !! number of samples in this trace !!
        !BYTE# 117-118
        integer(2) dt      !! sample interval! in micro-seconds
        !BYTE# 119-120
        integer(2) gain    !! gain type of field instruments code:
!                             1 = fixed
!                             2 = binary
!                             3 = floating point
!                             4 ---- N = optional use !!
        !BYTE# 121-122
        integer(2) igc     !! instrument gain constant !!
        !BYTE# 123-124
        integer(2) igi     !! instrument early or initial gain !!
        !BYTE# 125-126
        integer(2) corr    !! correlated:
!                             1 = no
!                             2 = yes !!
        !BYTE# 127-128
        integer(2) sfs     !! sweep frequency at start !!
        !BYTE# 129-130
        integer(2) sfe     !! sweep frequency at end !!
        !BYTE# 131-132
        integer(2) slen    !! sweep length in ms !!
        !BYTE# 133-134
        integer(2) styp    !! sweep type code:
!                             1 = linear
!                             2 = cos-squared
!                             3 = other !!
        !BYTE# 135-136
        integer(2) stas    !! sweep trace length at start in ms !!
        !BYTE# 137-138
        integer(2) stae    !! sweep trace length at end in ms !!
        !BYTE# 139-140
        integer(2) tatyp   !! taper type: 1=linear, 2=cos^2, 3=other !!
        !BYTE# 141-142
        integer(2) afilf   !! alias filter frequency if used !!
        !BYTE# 143-144
        integer(2) afils   !! alias filter slope !!
        !BYTE# 145-146
        integer(2) nofilf  !! notch filter frequency if used !!
        !BYTE# 147-148
        integer(2) nofils  !! notch filter slope !!
        !BYTE# 149-150
        integer(2) lcf     !! low cut frequency if used !!
        !BYTE# 151-152
        integer(2) hcf     !! high cut frequncy if used !!
        !BYTE# 153-154
        integer(2) lcs     !! low cut slope !!
        !BYTE# 155-156
        integer(2) hcs     !! high cut slope !!
        !BYTE# 157-158
        integer(2) year    !! year data recorded !!
        !BYTE# 159-160
        integer(2) day     !! day of year !!
        !BYTE# 161-162
        integer(2) hour    !! hour of day (24 hour clock) !!
        !BYTE# 163-164
        integer(2) minute  !! minute of hour !!
        !BYTE# 165-166
        integer(2) sec     !! second of minute !!
        !BYTE# 167-168
        integer(2) timbas  !! time basis code:
!                             1 = local
!                             2 = GMT
!                             3 = other !!
        !BYTE# 169-170
        integer(2) trwf    !! trace weighting factor, defined as 1/2^N
!                             volts for the least sigificant bit !!
        !BYTE# 171-172
        integer(2) grnors  !! geophone group number of roll switch
!                             position one !!
        !BYTE# 173-174
        integer(2) grnofr  !! geophone group number of trace one within
!                             original field record !!
        !BYTE# 175-176
        integer(2) grnlof  !! geophone group number of last trace within
!                             original field record !!
        !BYTE# 177-178
        integer(2) gaps    !! gap size (total number of groups dropped) !!
        !BYTE# 179-180
        integer(2) otrav   !! overtravel taper code:
!                             1 = down (or behind)
!                             2 = up (or ahead) !!

!       !! cwp local assignments !!
        !BYTE# 181-184
        real(4) d1         !! sample spacing for non-seismic data !!
        !BYTE# 185-188
        real(4) f1         !! first sample location for non-seismic data !!
        !BYTE# 189-192
        real(4) d2         !! sample spacing between traces !!
        !BYTE# 193-196
        real(4) f2         !! first trace location !!
        !BYTE# 197-200
        real(4) ungpow     !! negative of power used for dynamic
!                             range compression !!
        !BYTE# 201-204
        real(4) unscale    !! reciprocal of scaling factor to normalize
!                             range !!
!
        !BYTE# 205-206
        integer(2) mark    !! mark selected traces !!

!       !! SLTSU local assignments !! 
        !BYTE# 207-208
        integer(2) mutb    !! mute time at bottom (start time)
!                             bottom mute ends at last sample   !!
        !BYTE# 209-212
        real(4) dz         !! depth sampling interval in (m or ft)
!                             if =0.0, input are time samples       !!
        !BYTE# 213-216
        real(4) fz         !! depth of first sample in (m or ft)  !!
        !BYTE# 217-218
        integer(2) n2      !! number of traces per cdp or per shot !!
        !BYTE# 219-220
        integer(2) shortpad  !! alignment padding !!
        !BYTE# 221-224
        integer(4) ntr     !! number of traces !!

!       !! SLTSU local assignments end !! 
        !BYTE# 225-240
        integer(2),dimension(8) ::  unass     !! unassigned !!
        !end of su header
        
    end type
    
    
    type t_suformat
        type(t_suhdr) :: hdr
        
        real(4),dimension(:),allocatable :: trace
        
    end type
    
    contains
    
    subroutine init_suheader(sudata,ntr)
        type(t_suformat),dimension(ntr) :: sudata
        sudata(:)%hdr%tracl=0
        sudata(:)%hdr%tracr=0
        sudata(:)%hdr%fldr=0
        sudata(:)%hdr%tracf=0
        sudata(:)%hdr%ep=0
        sudata(:)%hdr%cdp=0
        sudata(:)%hdr%cdpt=0
        sudata(:)%hdr%trid=0
        sudata(:)%hdr%nvs=0
        sudata(:)%hdr%nhs=0
        sudata(:)%hdr%duse=0
        sudata(:)%hdr%offset=0
        sudata(:)%hdr%gelev=0
        sudata(:)%hdr%selev=0
        sudata(:)%hdr%sdepth=0
        sudata(:)%hdr%gdel=0
        sudata(:)%hdr%sdel=0
        sudata(:)%hdr%swdep=0
        sudata(:)%hdr%gwdep=0
        sudata(:)%hdr%scalel=0
        sudata(:)%hdr%scalco=0
        sudata(:)%hdr%sx=0
        sudata(:)%hdr%sy=0
        sudata(:)%hdr%gx=0
        sudata(:)%hdr%gy=0
        sudata(:)%hdr%counit=0
        sudata(:)%hdr%wevel=0
        sudata(:)%hdr%swevel=0
        sudata(:)%hdr%sut=0
        sudata(:)%hdr%gut=0
        sudata(:)%hdr%sstat=0
        sudata(:)%hdr%gstat=0
        sudata(:)%hdr%tstat=0
        sudata(:)%hdr%laga=0
        sudata(:)%hdr%lagb=0
        sudata(:)%hdr%delrt=0
        sudata(:)%hdr%muts=0
        sudata(:)%hdr%mute=0
        sudata(:)%hdr%ns=0
        sudata(:)%hdr%dt=0
        sudata(:)%hdr%gain=0
        sudata(:)%hdr%igc=0
        sudata(:)%hdr%igi=0
        sudata(:)%hdr%corr=0
        sudata(:)%hdr%sfs=0
        sudata(:)%hdr%sfe=0
        sudata(:)%hdr%slen=0
        sudata(:)%hdr%styp=0
        sudata(:)%hdr%stas=0
        sudata(:)%hdr%stae=0
        sudata(:)%hdr%tatyp=0
        sudata(:)%hdr%afilf=0
        sudata(:)%hdr%afils=0
        sudata(:)%hdr%nofilf=0
        sudata(:)%hdr%nofils=0
        sudata(:)%hdr%lcf=0
        sudata(:)%hdr%hcf=0
        sudata(:)%hdr%lcs=0
        sudata(:)%hdr%hcs=0
        sudata(:)%hdr%year=0
        sudata(:)%hdr%day=0
        sudata(:)%hdr%hour=0
        sudata(:)%hdr%minute=0
        sudata(:)%hdr%sec=0
        sudata(:)%hdr%timbas=0
        sudata(:)%hdr%trwf=0
        sudata(:)%hdr%grnors=0
        sudata(:)%hdr%grnofr=0
        sudata(:)%hdr%grnlof=0
        sudata(:)%hdr%gaps=0
        sudata(:)%hdr%otrav=0
        
        sudata(:)%hdr%d1=0.
        sudata(:)%hdr%f1=0.
        sudata(:)%hdr%d2=0.
        sudata(:)%hdr%f2=0.
        sudata(:)%hdr%ungpow=0.
        sudata(:)%hdr%unscale=0.
        sudata(:)%hdr%mark=0
        
        sudata(:)%hdr%mutb=0
        sudata(:)%hdr%dz=0.
        sudata(:)%hdr%fz=0.
        sudata(:)%hdr%n2=0
        sudata(:)%hdr%shortpad=0
        sudata(:)%hdr%ntr=0
        
        sudata(:)%hdr%unass(8)=0.
        
    end subroutine
    
    subroutine read_sudata(survey,cshot,sudata)
        character(4) :: survey
        character(4) :: cshot
        type(t_suformat),dimension(:),allocatable :: sudata

        character(:),allocatable :: tmp,data_file
        logical alive
        
        integer file_size
        integer(2) :: ns
        
        if(survey=='base') then
            tmp=get_setup_char('FILE_DATA')
        else
            tmp=get_setup_char('FILE_DATA_2')
        endif
        
        data_file=tmp//cshot//'.su'
                
        inquire(file=data_file, size=file_size, exist=alive) !file_size in bytes
        file_size=file_size/4 !file_size in sizeof(float)
        
        if(.not.alive) then
            if(mpiworld%is_master) write(*,*) 'ERROR: data file '//data_file//'does NOT exist!'
            stop
        endif
            
        open(11,file=data_file,action='read',access='stream')
        read(11,pos=115) ns  !get number of samples from 1st trace
        !if(ns<0) ns=ns+32768*2  !ns keyword is unsigned short in C
        close(11)
                
        ntr=file_size/(ns+60)
        if(mpiworld%is_master) write(*,*) 'Shot# '//cshot//' will read',ntr,' traces, each trace has',ns,'samples.'
        
        if(ntr*(ns+60)/=file_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: ntr*(ns+60) /= file_size. Possibly ns is not same for all traces..'
            stop
        endif
        
        if(allocated(sudata))deallocate(sudata)
        allocate(sudata(ntr))
        do itr=1,ntr
            allocate(sudata(itr)%trace(ns))
        enddo
        
        open(11,file=data_file,action='read',access='stream') !access='direct',recl=4*(ns+60))
        read(11) (sudata(itr)%hdr,sudata(itr)%trace, itr=1,ntr)
        close(11)
        
        sudata(:)%hdr%ntr=ntr
        
        call hud('Data read success.')
        
    end subroutine
    
    
!     subroutine write_sudata(cshot,nt,ntr,gather)
!         integer cshot,ntr
!         real,dimension(nt,ntr) :: gather
!         
!         type(t_suformat),dimension(ntr),allocatable :: sudata
!         
!         character(4) :: cshot
!         character,dimension(:),allocatable :: data_file
!         
!         data_file=tmp//cshot//'.su'
!         
!         open(12,file=data_file,action='write',access='direct',recl=4*(nt+60))
!         do itr=1,ntr
!             sudata(itr)=gather(itr)
!             sudata%ns=nt
!             write(11,rec=itr) sudata(itr)
!         enddo
!         
!         call hud('Data write success.')
!         
!         close(12)
!         
!     end subroutine

end
