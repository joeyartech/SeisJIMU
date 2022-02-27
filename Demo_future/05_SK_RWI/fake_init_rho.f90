character(64) :: fname
character(2) :: ck
integer,parameter :: n=301
real,dimension(n,n) :: rho_in, rho_out

CALL GET_COMMAND_ARGUMENT(1, fname)
CALL GET_COMMAND_ARGUMENT(2, ck); read(ck,*) k

open(11,file=fname,access='direct',recl=4*n*n)
read(11,rec=k) rho_in
close(11)

rho_out(:,1) = sum(rho_in,2)/n

do i=2,n
    rho_out(:,i)=rho_out(:,1)
enddo

open(12,file='init_rho',access='direct',recl=4*n*n)
write(12,rec=1) rho_out
close(12)
    
end
