!copyright 2013-2015 seiscopeii project, all rights reserved.

subroutine fcost_source_adjoint_ot(pbdir,acqui,inv,ot,flag)

    implicit none
    include 'mpif.h'
    include 'pbdirect.h'
    include 'acqui.h'
    include 'inversion.h'
    include 'optimaltransport.h'
    include 'common.h'


    !pbdirect
    type(pbdirect) :: pbdir
    !acqui
    type(acquisition) :: acqui
    !inversion
    type(fullwaveinv) :: inv
    !optimal transport
    type(optimaltransport) :: ot

    character(4) :: flag

    integer :: irec,it,i, itwind
    real :: dist
    integer :: ntaper
    real,dimension(:),allocatable :: taper
    real :: tmp
    real :: fcost_tmp
    real :: const
    real :: fcost_comm(3),tmp_comm(3)

    character(4) :: str_mype

    double precision :: dot_prod(2), dot_prod_tmp(2), angle

#ifdef time_profiling
    real(kind=8) :: time_fcost,tmp_time,time_comm
#endif

#ifdef time_profiling
    tmp_time=mpi_wtime()
#endif

    fcost_tmp=0.

    const = 0.5* pbdir%dt

  
    !compute L2 residual

    inv%residual_p(it,irec)=(inv%data_p(it,irec)-pbdir%data_p(it,irec))*acqui%weight_data_time(it,irec)
    !remultiply by the weight for the source of the adjoint !
    inv%residual_p(it,irec)=inv%residual_p(it,irec)*acqui%weight_data_time(it,irec)

    !write L2 residual & fcost
    inv%residual_l2=inv%residual_p
    call write_residu(inv,pbdir%nt,acqui%nrec,inv%residual_l2,'residu_l2_')
    call write_fcost(inv,1,[fcost_tmp],'fcost_l2_')

    !zw check mass conservation condition
    if( sum(inv%residual_p) > 1e-2*maxval(inv%residual_p) ) then
        write(*,*)'************************************************************ '
        write(*,*)'*** warning : mass conservation may not satisfacted !!! **** '
        write(*,*)'*** shot #',mype,' has non-negligible sum of data residuals: ', sum(inv%residual_p)/maxval(inv%residual_p)
        write(*,*)'************************************************************ '
    end if

    !compute KR residual
    inv%residual_p(:,:)=inv%residual_p(:,:)*ot%scale_residual
    select case (ot%ichoice_dmm)
    case (7);  call ot_sdmm_fishpack(ot,inv,pbdir%nt,acqui%nrec)    !zw fft solver for poisson eq, recommended
    case (8);  call ot_sdmm_mudpack(ot,inv,pbdir%nt,acqui%nrec)     !zw multigraid solver for poisson eq.
    case (11); call ot_sdmm_mudpack_3d(ot,inv,pbdir%nt,acqui%nrec)  !zw multigrid solver for poisson eq
    end select

    inv%residual_kr=inv%residual_p
    

    !ZW weight after sdmm
    if(ot%ifpost_sdmm_weight==1) then
        inv%residual_p=inv%residual_p*acqui%weight_data_time
        ntaper=nint(ot%taper_rfl/pbdir%dt)
        allocate(taper(ntaper))
        do i=1,ntaper
            taper(i)=(cos(  (i-1)*3.14/ntaper + 3.14  )+1)/2.
        enddo

        do irec=1,acqui%nrec
            dist=sqrt((acqui%s1(1)-acqui%r1(irec))**2+(acqui%s2(1)-acqui%r2(irec))**2+(acqui%s3(1)-acqui%r3(irec))**2)
            i=1  !ntaper
            if (dist>=acqui%offset_cut) then
                itwind= int(acqui%twind_boundary(irec)/pbdir%dt)+1

                do it=1, itwind-ntaper-1, 1 !before reflection & transition
                    inv%residual_p(it,irec)=0.
                end do

                do it=itwind-ntaper , itwind-1, 1 !before reflection,transition, need tapering to maintain dc=0 for ot
                    if (it<1 .or. it>pbdir%nt) then
                        i=i+1
                        cycle
                    end if
                    inv%residual_p(it,irec)=inv%residual_p(it,irec)*taper(i)
                    i=i+1
                end do
            else
                inv%residual_p(:,irec)=0.
            end if
        end do

        deallocate(taper)
    end if
    !end ZW

    call mpi_barrier(mpi_comm_world, infompi)

!     !dot product of l2 residual and unscaled wasserstein residual
!     dot_prod(1)=sum(dble(inv%residual_l2)*dble(inv%residual_kr))
!     call mpi_reduce (dot_prod(1), dot_prod_tmp(1), 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world,  infompi )
!     if(mype==0) print*, 'dot product of l2 and unscaled w residuals',dot_prod_tmp(1)
! 
!     !angle between l2 residual and unscaled wasserstein residual
!     dot_prod(1)=sum(dble( inv%residual_l2 * inv%residual_l2 ))
!     dot_prod(2)=sum(dble( inv%residual_kr * inv%residual_kr ))
!     call mpi_allreduce (dot_prod, dot_prod_tmp, 2, mpi_double_precision,  mpi_sum, mpi_comm_world,  infompi )
! 
!     dot_prod(1)=sum(dble( inv%residual_l2/dsqrt(dot_prod_tmp(1)) * inv%residual_kr/dsqrt(dot_prod_tmp(2)) ))
! 
!     call mpi_reduce (dot_prod(1), angle, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world,  infompi )
! 
!     angle=acos(angle)*180/3.1415927
!     if(mype==0) print*,'angle between l2 and unscaled kr residuals',angle
    
  
  
    !communicate fcost
    fcost_comm(1)=inv%fcost

    !write w residual & fcost   
    call write_residu(inv,pbdir%nt,acqui%nrec,inv%residual_kr,'residu_kr_')
    call write_fcost(inv,1,[fcost_comm(1)],'fcost_kr_')


    !sum on all processors
    tmp_comm=0.
    call mpi_allreduce (  fcost_comm(1),  tmp_comm,   1,   mpi_real,  mpi_sum,  mpi_comm_world,  infompi )

    inv%fcost=tmp_comm(1)


end subroutine fcost_source_adjoint_ot
