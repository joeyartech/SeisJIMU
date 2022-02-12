module m_matchfilter
use m_System
use m_math
use singleton

    private
    public matchfilter_estimate,matchfilter_apply_to_wavelet,matchfilter_apply_to_data,matchfilter_correlate_filter_residual,matchfilter_adjointsrc
        
    integer :: nt, ntr

    complex(fftkind),dimension(:),allocatable :: filter
    
    contains
    
    subroutine matchfilter_estimate(dsyn,dobs,nt_,ntr_,o_index,oif_stack,o_filter_time)
        real,dimension(nt_,ntr_) :: dsyn,dobs
        integer,optional :: o_index
        logical,optional :: oif_stack
        real,dimension(nt_),optional :: o_filter_time
        
        complex(fftkind),dimension(nt_,ntr_) :: dsynfft,dobsfft   !fftkind=double precision
        complex(fftkind),dimension(nt_)      :: numer, denom

        nt=nt_; ntr=ntr_
        
        dsynfft = fft(dcmplx(dsyn),dim=[1])  !may require "ulimit -s unlimited"
        dobsfft = fft(dcmplx(dobs),dim=[1])
        
        numer=sum(conjg(dsynfft)*dobsfft,2)
        denom=sum(conjg(dsynfft)*dsynfft,2)

        if(allocated(filter)) deallocate(filter)
        allocate(filter(nt))
        
        ! filter = numer/(denom+r_eps)

        filter = numer/(denom+1e-5*maxval(abs(denom)))

        ! if(present(o_index)) then
        !     call sysio_write('matchfilters_amp',real(abs (filter),kind=4),nt) !for purpose of quality control of results
        !     call sysio_write('matchfilters_ph' ,real(atan(filter),kind=4),nt) !for purpose of quality control of results
        ! endif

        if(either(oif_stack,.false.,present(oif_stack))) then
            call mpi_allreduce(MPI_IN_PLACE, filter, nt, MPI_DOUBLE_COMPLEX, MPI_SUM, mpiworld%communicator, mpiworld%ierr)
            filter =filter /mpiworld%nproc
        endif

        if(present(o_filter_time)) then
            o_filter_time=real(fft(dcmplx(filter),dim=[1],inv=.true.),kind=4)
        endif

    end subroutine
    
    subroutine matchfilter_apply_to_wavelet(wavelet)
        real,dimension(nt) :: wavelet
        
        complex(fftkind),dimension(nt) :: wlfft

        wlfft=fft(dcmplx(wavelet),dim=[1])
        
        wlfft=wlfft*filter
        
        wavelet=real(fft(wlfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    subroutine matchfilter_apply_to_data(dsyn)
        real,dimension(nt,ntr) :: dsyn
        
        complex(fftkind),dimension(nt,ntr) :: dsynfft
        
        dsynfft = fft(dcmplx(dsyn),dim=[1])
        
        do i=1,ntr
            dsynfft(:,i)=dsynfft(:,i)*filter(:)
        enddo
        
        dsyn=real(fft(dsynfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    subroutine matchfilter_correlate_filter_residual(dres)
        real,dimension(nt,ntr) :: dres
            
        complex(fftkind),dimension(nt,ntr) :: dresfft
        
        dresfft=fft(dcmplx(dres),dim=[1])
        
        do i=1,ntr
            dresfft(:,i)=conjg(filter)*dresfft(:,i) !crosscorrelation
        enddo
        
        dres= real(fft(dresfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine

    subroutine matchfilter_adjointsrc(dobs,dadj,scaled_w)
        real,dimension(nt) :: dobs,dadj
        real,dimension(nt) :: scaled_w
        
        complex(fftkind),dimension(nt) :: dobsfft   !fftkind=double precision
        complex(fftkind),dimension(nt) :: dadjfft
        complex(fftkind),dimension(nt) :: numer, denom
        
        dobsfft = fft(dcmplx(dobs),dim=[1])
        
        numer=fft(dcmplx(scaled_w),dim=[1])*dobsfft
        denom=conjg(dobsfft)*dobsfft

        dadjfft = numer/(denom+1e-5*maxval(abs(denom)))

        dadj = real(fft(dadjfft,dim=[1],inv=.true.),kind=4)

    end subroutine
    
end
