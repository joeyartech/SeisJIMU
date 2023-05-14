module m_hilbert
use m_System
use m_math
use singleton

    contains

    subroutine hilbert_transform(datain, dataout, nt, ntr)
        real,dimension(nt,ntr) :: datain
        real,dimension(:,:),allocatable :: dataout
        !! return hilbert transform of real signal datain
        !! and the analytic signal = cmplx(datain, dataout)

        real,dimension(:,:),allocatable :: data
        complex(fftkind),dimension(:,:),allocatable :: data_fft

        npt = 2**(int(log10(real(nt))/0.30104)+1)
        ! if ( npt /= nt) print*,'pad trace from length ', nt, ' to ',npt
        if (npt > 16784) call error('npt in Hilbert transform exceeds 16784 !')
        imid = npt/2

        !input
        allocate(data(npt,ntr))

        data(1:nt,:)=datain

        data_fft = fft(dcmplx(data),dim=[1])  !may require "ulimit -s unlimited"

        ! ! fourier transform 
        ! call cfft(c,npt,1)
        ! ! scaling 
        ! c=c/npt

        !90deg phase shift == multiply by i*sgn(freq)
        !ref: https://github.com/yanhuay/seisDD/blob/master/seisDD/lib/src/m_hilbert_transform.f90
        !created by Yanhua O. Yuan (yanhuay@princeton.edu)
        ! data_fft(1:imid-1    ,:) =-c_i*data_fft(1:imid-1    ,:) ! pos. spectrum (-i)
        ! data_fft(  imid      ,:) = 0.0                          ! d.c. component
        ! data_fft(  imid+1:npt,:) = c_i*data_fft(  imid+1:npt,:) ! neg. spectrum (i)
        data_fft(1       ,:) = 0.0                      ! DC
        data_fft(2:imid-1,:) =-c_i*data_fft(2:imid-1,:) ! pos. spectrum (-i)
        data_fft(imid:npt,:) = c_i*data_fft(imid:npt,:) ! neg. spectrum (i)

        ! inverse fourier transform
        data=real(fft(data_fft,dim=[1],inv=.true.),kind=4)

        !output
        dataout=data(1:nt,:)

        deallocate(data,data_fft)

    end subroutine

    ! !this is way problematic, esp when ntr>1
    ! subroutine hilbert_envelope(datain,dataout,nt,ntr)
    !     real,dimension(nt,ntr) :: datain
    !     real,dimension(:,:),allocatable :: dataout

    !     call hilbert_transform(datain,dataout,nt,ntr)
        
    !     dataout=sqrt(datain**2+dataout**2)

    ! end subroutine

    subroutine hilbert_envelope(datain,dataout,nt,ntr)
        real,dimension(nt,ntr) :: datain
        real,dimension(:,:),allocatable :: dataout

        type(t_suformat) :: seismo
        
        call suformat_write('tmp',datain,nt,ntr)
        call execute_command_line('suenv < '//dir_out//'/tmp.su > tmp1.su')
        ! call seismo%init(nt,ntr)
        call seismo%read('tmp1.su')
        dataout=seismo%trs
    end subroutine

    subroutine hilbert_phase(datain,dataout,nt,ntr)
        real,dimension(nt,ntr) :: datain
        real,dimension(:,:),allocatable :: dataout
        
        call hilbert_transform(datain,dataout,nt,ntr)
        
        dataout=atan2(dataout,datain)

    end subroutine

end module m_hilbert