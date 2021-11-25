program main

integer,parameter :: nz=349, nx=349, nt=51
real,dimension(nz,nx,nt) :: snap1,snap2,snap3

character(80) :: file1,file2

read(*,*) file1,file2

open(11,file=file1,access='direct',recl=4*nz*nx*nt)
read(11,rec=1) snap1
close(11)

open(12,file=file2,access='direct',recl=4*nz*nx*nt)
read(12,rec=1) snap2
close(12)

amaxval1=maxval(abs(snap1))
amaxval2=maxval(abs(snap2))

snap3=abs(snap1/amaxval1) - abs(snap2/amaxval2)

open(13,file='snaps',access='direct',recl=4*nz*nx*nt)
write(13,rec=1) snap3
close(13)

end
