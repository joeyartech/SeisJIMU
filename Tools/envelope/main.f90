program main
use m_sysio, only: dir_in
use m_suformat
use m_hilbert
    
    character(80) :: file
    type(t_suformat) :: sudata
    real,dimension(:,:),allocatable :: datain, dataout

    !read(*,*) file
    file='dsyn_Shot0001.su'

    dir_in='./'
    call sudata%read(file)


    call alloc(datain,sudata%ns,sudata%ntr) ; datain=sudata%trs
    call alloc(dataout,sudata%ns,sudata%ntr)

    open(12,file='hilbert_transformed',access='direct',recl=4*sudata%ns*sudata%ntr)
    write(12,rec=1) datain
    call hilbert_transform(datain,dataout,sudata%ns,sudata%ntr)
    write(12,rec=2) dataout
    call hilbert_envelope(datain,dataout,sudata%ns,sudata%ntr)
    write(12,rec=3) dataout
    call hilbert_phase(datain,dataout,sudata%ns,sudata%ntr)
    write(12,rec=4) dataout
    close(12)
        

end program
