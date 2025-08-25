module m_hilbert
use m_System
use m_math
use singleton

    contains

    subroutine hilbert_transform(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(nt,ntr) :: dout

        complex(fftkind),dimension(nt,ntr) :: dfft

        dfft = fft(dcmplx(din),dim=[1])  !may require "ulimit -s unlimited"

        !90deg phase shift == multiply by i*sgn(freq)
        !1st implementation
        !ref: https://github.com/yanhuay/seisDD/blob/master/seisDD/lib/src/m_hilbert_transform.f90
        !created by Yanhua O. Yuan (yanhuay@princeton.edu)
        ! data_fft(1:imid-1    ,:) =-c_i*data_fft(1:imid-1    ,:) ! pos. spectrum (-i)
        ! data_fft(  imid      ,:) = 0.0                          ! d.c. component
        ! data_fft(  imid+1:npt,:) = c_i*data_fft(  imid+1:npt,:) ! neg. spectrum (i)

        ! !2nd implementation
        ! if(mod(nt,2)==0) then !if nt is even
        !     dfft(1          ,:)= 0.0                     ! DC
        !     dfft(2:nt/2     ,:)=-c_i*dfft(2:nt/2     ,:) ! pos. spectrum (-i)
        !     dfft(  nt/2+1:nt,:)= c_i*dfft(  nt/2+1:nt,:) ! neg. spectrum (i)
        ! else !if nt is odd, 1:(nt+1)/2 are DC & positive freq, (nt+1)/2+1:nt are negative freq
        !     dfft(1              ,:)= 0.0                         ! DC
        !     dfft(2:(nt+1)/2     ,:)=-c_i*dfft(2:(nt+1)/2     ,:) ! pos. spectrum (-i)
        !     dfft(  (nt+1)/2+1:nt,:)= c_i*dfft(  (nt+1)/2+1:nt,:) ! neg. spectrum (i)
        ! endif

        !3rd implementation, set DC and neg freq to 0
        if(mod(nt,2)==0) then !if nt is even
            dfft(1          ,:)= 0.0
            dfft(  nt/2+1:nt,:)= 0.0
        else !if nt is odd, 1:(nt+1)/2 are DC & positive freq, (nt+1)/2+1:nt are negative freq
            dfft(1              ,:)= 0.0
            dfft(  (nt+1)/2+1:nt,:)= 0.0
        endif


        !inverse fourier transform
        ! dout=real(fft(dfft,dim=[1],inv=.true.),kind=4)
        dout=2.*aimag(fft(dfft,dim=[1],inv=.true.))

    end subroutine

    subroutine hilbert_transform2(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(nt,ntr) :: dout

        complex(fftkind),dimension(nt,ntr) :: dfft

        dfft = fft(dcmplx(din),dim=[2])  !may require "ulimit -s unlimited"

        !90deg phase shift == multiply by i*sgn(freq)
        !1st implementation
        !ref: https://github.com/yanhuay/seisDD/blob/master/seisDD/lib/src/m_hilbert_transform.f90
        !created by Yanhua O. Yuan (yanhuay@princeton.edu)
        ! data_fft(1:imid-1    ,:) =-c_i*data_fft(1:imid-1    ,:) ! pos. spectrum (-i)
        ! data_fft(  imid      ,:) = 0.0                          ! d.c. component
        ! data_fft(  imid+1:npt,:) = c_i*data_fft(  imid+1:npt,:) ! neg. spectrum (i)

        ! !2nd implementation
        ! if(mod(nt,2)==0) then !if nt is even
        !     dfft(1          ,:)= 0.0                     ! DC
        !     dfft(2:nt/2     ,:)=-c_i*dfft(2:nt/2     ,:) ! pos. spectrum (-i)
        !     dfft(  nt/2+1:nt,:)= c_i*dfft(  nt/2+1:nt,:) ! neg. spectrum (i)
        ! else !if nt is odd, 1:(nt+1)/2 are DC & positive freq, (nt+1)/2+1:nt are negative freq
        !     dfft(1              ,:)= 0.0                         ! DC
        !     dfft(2:(nt+1)/2     ,:)=-c_i*dfft(2:(nt+1)/2     ,:) ! pos. spectrum (-i)
        !     dfft(  (nt+1)/2+1:nt,:)= c_i*dfft(  (nt+1)/2+1:nt,:) ! neg. spectrum (i)
        ! endif

        !3rd implementation, set DC and neg freq to 0
        if(mod(ntr,2)==0) then !if nt is even
            dfft(:,1          )= 0.0
            dfft(:,ntr/2+1:ntr)= 0.0
        else !if nt is odd, 1:(nt+1)/2 are DC & positive freq, (nt+1)/2+1:nt are negative freq
            dfft(:,1              )= 0.0
            dfft(:,  (nrt+1)/2+1:ntr)= 0.0
        endif


        !inverse fourier transform
        ! dout=real(fft(dfft,dim=[1],inv=.true.),kind=4)
        dout=2.*aimag(fft(dfft,dim=[2],inv=.true.))

    end subroutine

    subroutine hilbert_transform_1d(din, dout, nt)
        implicit none
        integer, intent(in) :: nt
        real, dimension(nt), intent(in) :: din
        real, dimension(nt), intent(out) :: dout

        complex, dimension(nt) :: dfft
        integer :: i, mid

        ! 进行 FFT，转换为复数频域表示
        dfft = fft(dcmplx(din), dim=[1])

        ! 设置 Hilbert 滤波器：去掉 DC 和负频率成分
        if (mod(nt, 2) == 0) then
            ! nt 是偶数
            dfft(1)        = 0.0
            dfft(nt/2+1:)  = 0.0
        else
            ! nt 是奇数
            dfft(1)            = 0.0
            dfft((nt+1)/2+1:)  = 0.0
        endif

        ! IFFT 并提取虚部，乘以2 得到 Hilbert 变换结果
        dout = 2.0 * aimag(fft(dcmplx(dfft), dim=[1], inv=.true.))

    end subroutine


    ! !this is way problematic, esp when ntr>1
    subroutine hilbert_envelope(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(:,:),allocatable :: dout

        call hilbert_transform(din,dout,nt,ntr)
        
        dout=sqrt(din**2+dout**2)

    end subroutine

    ! subroutine hilbert_envelope(datain,dataout,nt,ntr)
    !     real,dimension(nt,ntr) :: datain
    !     real,dimension(:,:),allocatable :: dataout

    !     type(t_suformat) :: seismo
        
    !     call suformat_write('tmp',datain,nt,ntr)
    !     call execute_command_line('suenv < '//dir_out//'/tmp.su > tmp1.su')
    !     ! call seismo%init(nt,ntr)
    !     call seismo%read('tmp1.su')
    !     dataout=seismo%trs
    ! end subroutine

    subroutine hilbert_phase(din,dout,nt,ntr)
        real,dimension(nt,ntr) :: din
        real,dimension(:,:),allocatable :: dout
        
        call hilbert_transform(din,dout,nt,ntr)
        
        dout=atan2(dout,din)

    end subroutine

end module m_hilbert