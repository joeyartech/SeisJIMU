program main

use m_resamp


real,dimension(100,4)  :: in
real,dimension(200,4)  :: out

do it=1,100
    in(it,:)=sin((it-1)*0.2)
enddo

open(11,file='in',access='direct',recl=4*size(in))
write(11,rec=1) in
close(11)

call resamp(ntr=4, o_fin=0.1,din=0.002, nin=100, datain=in,&
                  o_fout=0.02,dout=0.001,nout=200,dataout=out)

open(12,file='out',access='direct',recl=4*size(out))
write(12,rec=1) out
close(12)

end
