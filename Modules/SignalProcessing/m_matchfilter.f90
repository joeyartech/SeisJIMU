module m_matchfilter
use m_mpienv
use singleton

    private
    public matchfilter_estimate,matchfilter_apply_to_wavelet,matchfilter_apply_to_data,matchfilter_correlate_filter_residual
    
    real,parameter :: eps=1e-15
    
    complex(fftkind),dimension(:),allocatable :: filter
    
    contains
    
    subroutine matchfilter_estimate(nt,nd,dsyn,dobs,if_stack)
        integer nt,nd
        real,dimension(nt,nd) :: dsyn,dobs
        logical if_stack
        
        complex(fftkind),dimension(nt,nd) :: dsynfft,dobsfft   !fftkind=double precision
        complex(fftkind),dimension(nt)    :: nom, denom
        
        dsynfft = fft(dcmplx(dsyn),dim=[1])  !maybe this array is so large that require "ulimit -s unlimited"
        dobsfft = fft(dcmplx(dobs),dim=[1])
        
        nom  =sum(conjg(dsynfft)*dobsfft,2)
        denom=sum(conjg(dsynfft)*dsynfft,2)
        
        if(allocated(filter)) deallocate(filter)
        allocate(filter(nt))
        
        filter = nom/(denom+eps)
        
        if(if_stack) then
            call mpi_allreduce(MPI_IN_PLACE, filter, nt, MPI_DOUBLE_COMPLEX, MPI_SUM, mpiworld%communicator, mpiworld%ierr)
            filter =filter /mpiworld%nproc
        endif
        
    end subroutine
    
    subroutine matchfilter_apply_to_wavelet(nt,wavelet)
        integer nt
        real,dimension(nt) :: wavelet
        
        complex(fftkind),dimension(nt) :: wlfft
        
        wlfft=fft(dcmplx(wavelet),dim=[1])
        
        wlfft=wlfft*filter
        
        wavelet=real(fft(wlfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    subroutine matchfilter_apply_to_data(nt,nd,dsyn)
        integer nt,nd
        real,dimension(nt,nd) :: dsyn
        
        complex(fftkind),dimension(nt,nd) :: dsynfft
        
        dsynfft = fft(dcmplx(dsyn),dim=[1])
        
        do i=1,nd
            dsynfft(:,i)=dsynfft(:,i)*filter(:)
        enddo
        
        dsyn=real(fft(dsynfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    subroutine matchfilter_correlate_filter_residual(nt,nd,dres)
        integer nt,nd
        real,dimension(nt,nd) :: dres
        
        complex(fftkind),dimension(nt,nd) :: dresfft
        
        dresfft=fft(dcmplx(dres),dim=[1])
        
        do i=1,nd
            dresfft(:,i)=conjg(filter)*dresfft(:,i) !crosscorrelation
        enddo
        
        dres= real(fft(dresfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    
end
