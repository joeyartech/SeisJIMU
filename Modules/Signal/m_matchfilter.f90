module m_matchfilter
use m_System
use singleton

    private
    public matchfilter_estimate,matchfilter_apply_to_wavelet,matchfilter_apply_to_data,matchfilter_correlate_filter_residual
    
    real,parameter :: eps=1e-15
    
    integer :: nt, ntr

    complex(fftkind),dimension(:),allocatable :: filter
    
    contains
    
    subroutine matchfilter_estimate(dsyn,dobs,nt,ntr,o_index,oif_stack)
        real,dimension(nt,ntr) :: dsyn,dobs
        integer,optional :: o_index
        logical,optional :: oif_stack
        
        complex(fftkind),dimension(:,:),allocatable :: dsynfft,dobsfft   !fftkind=double precision
        complex(fftkind),dimension(:),allocatable   :: nom, denom
        
        allocate(dsynfft(nt,ntr))
        allocate(dobsfft(nt,ntr))

        allocate(nom(nt), denom(nt))

        dsynfft = fft(dcmplx(dsyn),dim=[1])  !maybe this array is so large that require "ulimit -s unlimited"
        dobsfft = fft(dcmplx(dobs),dim=[1])
        
        nom  =sum(conjg(dsynfft)*dobsfft,2)
        denom=sum(conjg(dsynfft)*dsynfft,2)
        
        if(allocated(filter)) deallocate(filter)
        allocate(filter(nt))
        
        filter = nom/(denom+eps)
        
        if(present(o_index)) then
            open(12,file=dir_out//'matchfilters_amp',access='direct',recl=4*nt) !for purpose of quality control of results
            write(12,rec=o_index) real(abs(filter),kind=4)
            close(12)
            open(12,file=dir_out//'matchfilters_phase',access='direct',recl=4*nt) !for purpose of quality control of results
            write(12,rec=o_index) real(atan(filter),kind=4)
            close(12)
        endif
        
        if(either(oif_stack,.false.,present(oif_stack))) then
            call mpi_allreduce(MPI_IN_PLACE, filter, nt, MPI_DOUBLE_COMPLEX, MPI_SUM, mpiworld%communicator, mpiworld%ierr)
            filter =filter /mpiworld%nproc
        endif
        
    end subroutine
    
    subroutine matchfilter_apply_to_wavelet(wavelet)
        real,dimension(nt) :: wavelet
        
        complex(fftkind),dimension(:),allocatable :: wlfft

        allocate(wlfft(nt))
        
        wlfft=fft(dcmplx(wavelet),dim=[1])
        
        wlfft=wlfft*filter
        
        wavelet=real(fft(wlfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    subroutine matchfilter_apply_to_data(dsyn)
        real,dimension(nt,ntr) :: dsyn
        
        complex(fftkind),dimension(:,:),allocatable :: dsynfft
        
        allocate(dsynfft(nt,ntr))

        dsynfft = fft(dcmplx(dsyn),dim=[1])
        
        do i=1,ntr
            dsynfft(:,i)=dsynfft(:,i)*filter(:)
        enddo
        
        dsyn=real(fft(dsynfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    subroutine matchfilter_correlate_filter_residual(dres)
        real,dimension(nt,ntr) :: dres
        
        complex(fftkind),dimension(:,:),allocatable :: dresfft

        allocate(dresfft(nt,ntr))
        
        dresfft=fft(dcmplx(dres),dim=[1])
        
        do i=1,ntr
            dresfft(:,i)=conjg(filter)*dresfft(:,i) !crosscorrelation
        enddo
        
        dres= real(fft(dresfft,dim=[1],inv=.true.),kind=4)
        
    end subroutine
    
    
end
