  character(80) :: fname_in, fname_out  
  real,dimension(:,:),allocatable :: data
  real vmax

  read(*,*) fname_in, fname_out
  read(*,*) n1,n2
  read(*,*) vmax

  allocate(data(n1,n2))
  
  open(1,file=fname_in, access='direct',recl=n1*n2*4) 
  read(1,rec=1) data
  close(1)
  
  where (data>vmax) data=vmax
  
  open(2,file=fname_out,access='direct',recl=n1*n2*4)
  write(2,rec=1) data
  close(2)


  end
