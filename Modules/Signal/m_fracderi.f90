module m_fracderi
use m_System
use m_math
use singleton

    contains

    subroutine fracderi(datain, dataout, nt, ntr, dt, order)
        real,dimension(nt,ntr) :: datain
        real,dimension(:,:),allocatable :: dataout
        
        real,dimension(:,:),allocatable :: data
        complex(fftkind),dimension(:,:),allocatable :: data_fft

        real,dimension(:),allocatable :: w

        npt = 2**(int(log10(real(nt))/0.30104)+1)
        ! if ( npt /= nt) print*,'pad trace from length ', nt, ' to ',npt
        if (npt > 16784) call error('npt in fracderi exceeds 16784 !')
        imid = npt/2

        !input
        call alloc(data,npt,ntr)

        data(1:nt,:)=datain

        data_fft = fft(dcmplx(data),dim=[1])  !may require "ulimit -s unlimited"

        dw=1./(nt-1)/dt
        allocate(w(npt))
        do i=1,imid-1
            w(i)=(i-1)*dw
        enddo
        do i=imid,npt
            w(i)=(i-1-npt)*dw
        enddo

        do i=1,ntr
        	data_fft(:,i) = data_fft(:,i) * (c_i*w)**order
        enddo
        
        ! inverse fourier transform
        data=real(fft(data_fft,dim=[1],inv=.true.),kind=4)

        !output
        dataout=data(1:nt,:)

        deallocate(data,data_fft,w)

    end subroutine



end