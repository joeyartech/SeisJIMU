!References
!for SOCIG:
!M Giboli, R Baina, L Nicoletis & B Duquet, 2012, "Reverse Time Migration surface offset gathers part 1: a new method to produce ‘classical’ common image gathers", SEG expanded abstracts.
!R Baina et al, 2017, Method of processing seismic data by providing surface offset common image gathers. US Patent 9,632,192 B2
!https://github.com/slimgroup/ImageGather.jl/tree/main
!
subroutine modeling_imaging
use mpi
use m_System
use m_Modeling
use m_weighter
use m_smoother_laplacian_sparse

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u,fld_a
    type(t_correlate) :: a_star_u
    type(t_correlate) :: u_star_u
    real,dimension(:,:),allocatable :: tmp

    call alloc(correlate_image,m%nz,m%nx,m%ny,ppg%nimag)
    call alloc(correlate_energy,m%nz,m%nx,m%ny,ppg%nengy)
    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer


        call hud('----  Solving Au=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)

        if(setup%get_str('JOB')=='forward modeling') cycle

        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)

        call hud('----  Preparing adjoint source  ----')
        call wei%update
        call alloc(shot%dadj,shot%nt,shot%nrcv)

        if(setup%get_str('RTM_ADJSRC',o_default='dobs')=='dobs') then
            shot%dadj=shot%dobs
        else
            shot%dadj=shot%dobs-shot%dsyn
        endif

        call fld_a%ignite(o_wavelet=shot%dadj)
        call shot%write('dadj_',shot%dadj)

        call hud('----  Solving adjoint eqn & xcorrelate  ----')
        !Aᴴa = -Rᴴd
        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%init_correlate(u_star_u,'u_star_u')

        call ppg%adjoint(fld_a,fld_u,a_star_u,u_star_u)

        call hud('----  Assemble  ----')
        call ppg%assemble(a_star_u)
        call ppg%assemble(u_star_u)

        call hud('---------------------------------')
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    if(setup%get_str('JOB')=='forward modeling') then
        call mpiworld%final
        stop
    endif


    if(mpiworld%is_master) call a_star_u%write

    !collect global correlations
    call mpi_allreduce(mpi_in_place,correlate_image, m%n*ppg%nimag, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place,correlate_energy, m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    if(mpiworld%is_master) call sysio_write('correlate_image',correlate_image,m%n*ppg%nimag)
    if(mpiworld%is_master) call sysio_write('correlate_energy',correlate_energy,m%n*ppg%nengy)

    !write correlate
    if(mpiworld%is_master) then
        call a_star_u%write
        call u_star_u%write
    endif

    !allreduce energy, gradient
    call mpi_allreduce(mpi_in_place, correlate_energy, m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place, correlate_image, m%n*ppg%nimag, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    
    !scale by shotlist
!    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)
    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_image)

    if(mpiworld%is_master) call sysio_write('correlate_image', correlate_image, m%n*ppg%nimag)
    if(mpiworld%is_master) call sysio_write('correlate_energy',correlate_energy,m%n*ppg%nengy)

    call mpiworld%barrier
    
    is_first_in=.false.

    call mpiworld%barrier

end subroutine


subroutine modeling_imagegathers
use mpi
use m_System
use m_Modeling
use m_weighter
use m_smoother_laplacian_sparse

    logical,save :: is_first_in=.true.

    type(t_field) :: fld_u,fld_a
    type(t_correlate) :: a_star_u,a2_star_u
    type(t_correlate) :: u_star_u
    
    integer,dimension(:),allocatable :: ix_CIG
    real,dimension(:,:,:),allocatable :: CIGs

    off_shift = 5e3 !https://github.com/slimgroup/ImageGather.jl/blob/main/src/surface_gather.jl#L74
    call hud('off_shift = '//num2str(off_shift))

    call alloc(correlate_image,m%nz,m%nx,m%ny,ppg%nimag)
    call alloc(correlate_energy,m%nz,m%nx,m%ny,ppg%nengy)

    !max offset
    ! amax_offset=setup%get_real('CIG_MAX_OFFSET',o_default=num2str(m%nx*m%dx/2))
    amax_offset=setup%get_real('CIG_MAX_OFFSET',o_default=num2str(m%nx*m%dx))
    nh_CIG=floor(amax_offset/m%dx)

    !CIG positions
    ix_CIG=nint(setup%get_reals('CIG_X',o_default=num2str(shot%src%x))/m%dx-m%ox)+1
    call hud('ix_CIG:'//strcat(nums2strs(ix_CIG)))
    nx_CIG=size(ix_CIG)

    call hud('nh_CIG, nx_CIG = '//num2str(nh_CIG)//', '//num2str(nx_CIG))

    call alloc(CIGs,m%nz,nh_CIG,nx_CIG)

    
    call hud('===== START LOOP OVER SHOTS =====')
    
    do i=1,shls%nshots_per_processor
    
        call shot%init(shls%yield(i))
        call shot%read_from_data
        call shot%set_var_time
        call shot%set_var_space(index(ppg%info,'FDSG')>0)

        call hud('Modeling Shot# '//shot%sindex)
        
        call cb%init(ppg%nbndlayer)
        call cb%project
        
        call ppg%check_discretization
        call ppg%init
        call ppg%init_abslayer


        call hud('----  Solving Au=s  ----')
        call ppg%init_field(fld_u, name='fld_u');    call fld_u%ignite
        call ppg%forward(fld_u)
        call fld_u%acquire; call shot%write('Ru_',shot%dsyn)

        if(setup%get_str('JOB')=='forward modeling') cycle


        call hud('----  Preparing adjoint source  ----')
        call wei%update
        call alloc(shot%dadj,shot%nt,shot%nrcv)


        call hud('----  migration: Solving adjoint eqn & xcorrelate  ----')
        !Aᴴa = -Rᴴd
        if(setup%get_str('RTM_ADJSRC',o_default='dobs')=='dobs') then
            shot%dadj=shot%dobs*wei%weight
        else
            shot%dadj=(shot%dobs-shot%dsyn)*wei%weight
        endif
        
        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)
        call shot%write('dadj_',shot%dadj)
        call fld_a%ignite(o_wavelet=shot%dadj)

        call ppg%init_correlate(a_star_u,'a_star_u')
        call ppg%init_correlate(u_star_u,'u_star_u')

        call ppg%adjoint(fld_a,fld_u,a_star_u,u_star_u)


        call hud('----  attribute migration: Solving adjoint eqn & xcorrelate  ----')
        !Aᴴa = -Rᴴd
        call add_attribute(shot%dadj)

        call ppg%init_field(fld_a,name='fld_a',ois_adjoint=.true.)
        call shot%write('dadj2_',shot%dadj)
        call fld_a%ignite(o_wavelet=shot%dadj)

        call ppg%init_correlate(a2_star_u,'a2_star_u')

        call ppg%adjoint(fld_a,fld_u,a2_star_u)

        call compute_cig() !,ext='surface_offset')

        call hud('----  Assemble  ----')
        ! call ppg%assemble(a_star_u)
        ! call ppg%assemble(u_star_u)
        correlate_image(:,:,:,1) =correlate_image(:,:,:,1) +a_star_u%ipp
        correlate_energy(:,:,:,1)=correlate_energy(:,:,:,1)+u_star_u%epp
        
        call hud('---------------------------------')
        
    enddo
    
    call hud('        END LOOP OVER SHOTS        ')

    if(setup%get_str('JOB')=='forward modeling') then
        call mpiworld%final
        stop
    endif

    !write correlate
    if(mpiworld%is_master) then
        call a_star_u%write
        call a2_star_u%write
        call u_star_u%write
    endif

    !collect global correlations
    call mpi_allreduce(mpi_in_place,correlate_image,  m%n*ppg%nimag, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place,correlate_energy, m%n*ppg%nengy, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)
    call mpi_allreduce(mpi_in_place,CIGs,        m%nz*nx_CIG*nh_CIG, mpi_real, mpi_sum, mpiworld%communicator, mpiworld%ierr)

    !postprocess correlation image
    call laplacian_filter(correlate_image)
    correlate_image(:,:,:,1) = correlate_image(:,:,:,1) / (correlate_energy(:,:,:,1)+r_eps)
    
    call sysio_write('image',correlate_image, size(correlate_image))
    call sysio_write('illum',correlate_energy,size(correlate_energy))
    call sysio_write('CIGs',CIGs,size(CIGs))

    !scale by shotlist
!    call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_gradient)
    ! call shls%scale(m%n*ppg%ngrad,o_from_sampled=correlate_image)
    is_first_in=.false.

    call mpiworld%barrier

    contains

    subroutine add_attribute(dadj)
use m_math
        real,dimension(shot%nt,shot%nrcv) :: dadj
        do ir=1,shot%nrcv    
            ! dadj(:,ir)=dadj(:,ir)*(shot%rcv(ir)%aoffset +off_shift)
            dadj(:,ir)=dadj(:,ir)*(shot%rcv(ir)%x - shot%src%x +off_shift)
        enddo

        ! do ir=0,9
        !     dadj(:,ir)          = dadj(:,ir)          *exp(-(9.-ir))!cos((9-ir)/9.*r_pi/2)
        ! enddo

        ! do ir=0,9
        !     dadj(:,shot%nrcv-ir)= dadj(:,shot%nrcv-ir)*exp(-(9.-ir))!cos((9-ir)/9.*r_pi/2)
        ! enddo

! ! print*,(shot%rcv(:)%aoffset+off_shift)
! print*, (shot%rcv(:)%x - shot%src%x)

        ! do ir=1,shot%nrcv
        !     offset=shot%rcv(ir)%x - shot%src%x
        !     if(offset>0) then
        !         offset=floor(offset/m%dx/3)*3*m%dx
        !     else
        !         offset=ceiling(offset/m%dx/3)*3*m%dx
        !     endif
        !     print*, ir, (shot%rcv(ir)%x - shot%src%x), offset
        !     dadj(:,ir)=dadj(:,ir)*offset
        ! enddo

    end subroutine

    subroutine compute_cig()
    use m_math
        real,dimension(:,:),allocatable :: imag1, imag2, offset_map

        !compute the regularized LS division of a_star_u (RTM image) & a2_star_u (attribute image)
        !regularization consists of the envelope and moving average
        call alloc(imag1,m%nz,m%nx); call alloc(imag2,m%nz,m%nx)
        imag1=proc_imag( a_star_u%ipp)
        imag2=proc_imag(a2_star_u%ipp)

call sysio_write('imag1',imag1,size(imag1))
call sysio_write('imag2',imag2,size(imag2))

        thres1=1e-6*maxval(imag1)
        thres2=1e-6*maxval(imag2)

        offset_map = (imag1*imag2 +1e-5) / (imag1*imag1 +1e-5) -off_shift

call sysio_write('offset_map',offset_map,size(offset_map))

        ! call dealloc(imag1,imag2)

        !transform RTM image into surface-offset CIGs
        do ix=1,nx_CIG
            do iz=1,m%nz

                if(imag1(iz,ix_CIG(ix))>thres1 .and. imag2(iz,ix_CIG(ix))>thres2) then

                    ! if(abs(offset_map(iz,ix)-shot%rcv(ir)%aoffset) <= 2*avg_daoffset) then !found binned offsets
                    ih=nint(abs(offset_map(iz,ix_CIG(ix)))/m%dx)+1
                    if(ih<=nh_CIG) then
print*,iz,ix,offset_map(iz,ix_CIG(ix)),ih
                        CIGs(iz,ih,ix) = CIGs(iz,ih,ix) + a_star_u%ipp(iz,ix_CIG(ix),1)

                    endif

                endif
            
            enddo
        enddo

    end subroutine

    function proc_imag(ipp) result(res)
    use m_hilbert
        real,dimension(m%nz,m%nx) :: ipp
        real,dimension(:,:),allocatable :: env,res

        integer,parameter :: jfx=-2,jlx=2 !window of the moving average in x dir, window length=5
        integer,parameter :: jfz=-2,jlz=2 !window of the moving average in z dir, window length=5
        
        env=ipp
        call hilbert_envelope(ipp,env,m%nz,m%nx)

        !moving average
        scal=1./(jlx-jfx+1)/(jlz-jfz+1)

        res=env
        do ix=1-jfx, m%nx-jlx
        do iz=1-jfz, m%nz-jlz
            res(iz,ix) = sum(env(iz+jfz:iz+jlz, ix+jfx:ix+jlx))*scal
        enddo
        enddo

        call dealloc(env)

    end function

    subroutine laplacian_filter(image)
        real,dimension(m%nz,m%nx) :: image, tmp

        do ix=2,m%nx-1
        do iz=2,m%nz-1
            tmp(iz,ix) = image(iz+1,ix)+image(iz-1,ix)+image(iz,ix+1)+image(iz,ix-1)-4*image(iz,ix)
        enddo
        enddo

        image = tmp/(m%dz*m%dx)

    end subroutine

end subroutine
