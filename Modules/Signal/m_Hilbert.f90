module m_hilbert
use m_System
use m_math
use singleton

    contains

    subroutine hilbert_transform(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(:,:),allocatable :: dout

        complex(fftkind),dimension(nt,ntr) :: dfft

        dfft = fft(dcmplx(din),dim=[1])  !may require "ulimit -s unlimited"

        !90deg phase shift == multiply by i*sgn(freq)
        !ref: https://github.com/yanhuay/seisDD/blob/master/seisDD/lib/src/m_hilbert_transform.f90
        !created by Yanhua O. Yuan (yanhuay@princeton.edu)
        ! data_fft(1:imid-1    ,:) =-c_i*data_fft(1:imid-1    ,:) ! pos. spectrum (-i)
        ! data_fft(  imid      ,:) = 0.0                          ! d.c. component
        ! data_fft(  imid+1:npt,:) = c_i*data_fft(  imid+1:npt,:) ! neg. spectrum (i)

        if(mod(nt,2)==0) then !if nt is even
            dfft(1          ,:)= 0.0                     ! DC
            dfft(2:nt/2     ,:)=-c_i*dfft(2:nt/2     ,:) ! pos. spectrum (-i)
            dfft(  nt/2+1:nt,:)= c_i*dfft(  nt/2+1:nt,:) ! neg. spectrum (i)
        else !if nt is odd, 1:(nt+1)/2 are DC & positive freq, (nt+1)/2+1:nt are negative freq
            dfft(1              ,:)= 0.0                         ! DC
            dfft(2:(nt+1)/2     ,:)=-c_i*dfft(2:(nt+1)/2     ,:) ! pos. spectrum (-i)
            dfft(  (nt+1)/2+1:nt,:)= c_i*dfft(  (nt+1)/2+1:nt,:) ! neg. spectrum (i)
        endif

        !inverse fourier transform
        dout=real(fft(dfft,dim=[1],inv=.true.),kind=4)

    end subroutine

    ! !this is way problematic, esp when ntr>1
    subroutine hilbert_envelope(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(:,:),allocatable :: dout

        call hilbert_transform(din,dout,nt,ntr)
        
        dout=sqrt(din**2+dout**2)

    end subroutine

    ! subroutine hilbert_envelope(datain,dataout,nt,ntr)
    !     real,dimension(nt,ntr) :: datain
    !     real,dimension(:,:),allocatable :: dataout

    !     type(t_suformat) :: seismo
        
    !     call suformat_write('tmp',datain,nt,ntr)
    !     call execute_command_line('suenv < '//dir_out//'/tmp.su > tmp1.su')
    !     ! call seismo%init(nt,ntr)
    !     call seismo%read('tmp1.su')
    !     dataout=seismo%trs
    ! end subroutine

    subroutine hilbert_phase(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(:,:),allocatable :: dout
        
        call hilbert_transform(din,dout,nt,ntr)
        
        dout=atan2(dout,din)

    end subroutine

end module m_hilbert