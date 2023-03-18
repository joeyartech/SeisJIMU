program main

integer,parameter :: nz=124, nx=124
real,parameter :: inv_2h = 1./2./20.
real,dimension(nz,nx) :: u,v, uu, E, ph, dph_dz, dph_dx
character(80) :: file_u, file_v

read(*,*) file_u, file_v, i

open(11,file=file_u,access='direct',recl=4*nz*nx)
read(11,rec=i) u
close(11)

open(11,file=file_v,access='direct',recl=4*nz*nx)
read(11,rec=i) v
close(11)

E=sqrt(u*u+v*v)

ph=atan2(v,u)

do ix=2,nx-1
do iz=2,nz-1
    dph_dz(iz,ix) = asin(sin(ph(iz+1,ix) - ph(iz-1,ix)))*inv_2h
    dph_dx(iz,ix) = asin(sin(ph(iz,ix+1) - ph(iz,ix-1)))*inv_2h
enddo
enddo

open(12,file='E',access='direct',recl=4*nz*nx)
write(12,rec=1) E
close(12)

open(12,file='ph',access='direct',recl=4*nz*nx)
write(12,rec=1) ph
close(12)

open(12,file='dph_dz',access='direct',recl=4*nz*nx)
write(12,rec=1) dph_dz
close(12)

open(12,file='dph_dx',access='direct',recl=4*nz*nx)
write(12,rec=1) dph_dx
close(12)

end
