real,allocatable :: vp(:,:)
character(8) :: str

dx=20.; fpeak=7.; dt=0.004

CALL get_command_argument(1, str)
read(str,*) n
print*,'read n=',n

allocate(vp(n,n))

call random_number(vp)

!vpmin/(2.5*fpeak) = lamda >= 5dx
vpmin=2.5*fpeak *5*dx

!dt <= 0.9*0.6*dx/vpmax
vpmax = 0.9*0.6*dx/dt

!scaling
vp=vp*(vpmax-vpmin)+vpmin

open(12,file='model',access='direct',recl=4*n*n)
write(12,rec=1) vp
close(12)


end
