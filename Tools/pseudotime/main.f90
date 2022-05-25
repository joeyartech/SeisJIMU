program main

    character(120) :: fin, fvel
    logical is_vel

    integer :: file_size1, file_size2
    logical :: exist

    write(*,*) 'Filename of the input model:'
    read(*,*) fin

    inquire(file=fin,size=file_size1,exist=exist)
    if(file_size1==0) exist=.false.
    if(.not.exist) then
        write(*,*) 'ERROR: file '//trim(adjustl(fin))//' does NOT exist or has 0 length.'
        stop
    endif

    write(*,*) 'Is it velocity?'
    read(*,*) is_vel

    if(is_vel) then
        fvel=fin
    else
        write(*,*) 'Filename of the (background) velocity model:'
        read(*,*) fvel
        inquire(file=fvel,size=file_size2,exist=exist)
        if(file_size2==0) exist=.false.
        if(.not.exist) then
            write(*,*) 'ERROR: file '//trim(adjustl(fvel))//' does NOT exist or has 0 length.'
            stop
        endif
        if(file_size2/=file_size1) then
            write(*,*) 'ERROR: Size of velcity model /= Size of input model.'
            stop
        endif
    endif

    write(*,*) 'Enter job (=1: depth to time; =2: time to depth):'
    read(*,*) job

    if(job==1) call depth2time(fin,is_vel,fvel)

    if(job==2) call time2depth(fin,is_vel,fvel)

end program

subroutine depth2time(fin,is_vel,fvel)
use m_pseudotime

    character(*) :: fin, fvel
    logical is_vel

    integer file_size
    real,dimension(:,:,:),allocatable :: min, xout
    real,dimension(:,:,:),allocatable :: v_z


    write(*,*) '1st dim of input depth-domain model (nz, dz):'
    read(*,*) nz,dz

    inquire(file=fin,size=file_size)
    nx=file_size/4/nz
    write(*,*) 'Shape of the model (nz x nx):'
    write(*,*) nz,nx,nz*nx

    allocate(min(nz,nx,1))
    open(11,file=fin,access='direct',recl=4*nz*nx)
    read(11,rec=1)min
    close(11)

    allocate(v_z(nz,nx,1))
    open(11,file=fvel,access='direct',recl=4*nz*nx)
    read(11,rec=1)v_z
    close(11)

    call pseudotime_init('z->t',vmin=minval(min),vmax=maxval(min),nx_=nx,ny_=1, &
        nz_=nz, dz_=dz, &
        nt_=nt, dt_=dt)

    write(*,*) '1st dim of output time-domain model (nt, dt)'
    write(*,*) nt, dt

    write(*,*) 'Shape of the output depth-domain model (nt x nx):'
    write(*,*) nt,nx,nt*nx


    if(is_vel) then
        call pseudotime_convert('z->t',min,xout)
    else
        call pseudotime_convert('z->t',min,xout,o_v=v_z)
    endif
    
    !write(*,*) 'Enter filename of the output depth-domain model:'
    !read(*,*) filename

    call execute_command_line('rm output_time_model',wait=.true.)
    open(12,file='output_time_model',access='direct',recl=4*size(xout))
    write(12,rec=1) xout
    close(12)

end subroutine

subroutine time2depth(fin,is_vel,fvel)
use m_pseudotime

    character(*) :: fin, fvel
    logical is_vel

    integer file_size
    real,dimension(:,:,:),allocatable :: xin, mout
    real,dimension(:,:,:),allocatable :: v_t


    write(*,*) '1st dim of input time-domain model (nt, dt):'
    read(*,*) nt,dt

    inquire(file=fin,size=file_size)
    nx=file_size/4/nt
    write(*,*) 'Shape of the model (nt x nx):'
    write(*,*) nt,nx,nt*nx

    allocate(xin(nt,nx,1))
    open(11,file=fin,access='direct',recl=4*nt*nx)
    read(11,rec=1)xin
    close(11)

    allocate(v_t(nt,nx,1))
    open(11,file=fvel,access='direct',recl=4*nt*nx)
    read(11,rec=1)v_t
    close(11)

    call pseudotime_init('t->z',vmin=minval(xin),vmax=maxval(xin),nx_=nx,ny_=1, &
        nz_=nz, dz_=dz, &
        nt_=nt, dt_=dt)

    write(*,*) '1st dim of output depth-domain model (nz, dz)'
    write(*,*) nz, dz

    write(*,*) 'Shape of the output depth-domain model (nz x nx):'
    write(*,*) nz,nx,nz*nx


    if(is_vel) then
        call pseudotime_convert('t->z',xin,mout)
    else
        call pseudotime_convert('t->z',xin,mout,o_v=v_t)
    endif
    
    !write(*,*) 'Enter filename of the output depth-domain model:'
    !read(*,*) filename

    call execute_command_line('rm output_time_model',wait=.true.)
    open(12,file='output_depth_model',access='direct',recl=4*size(mout))
    write(12,rec=1) mout
    close(12)

end subroutine
