module m_medfilt2d

    ! !2d median filter to denoise
    ! call alloc(tmp,cb%mz,cb%mx)
    ! call medfilt2d(o_F2_star_E0%drp_dt_dsp_dt(:,:,1,1),tmp,cb%mz,cb%mx); o_F2_star_E0%drp_dt_dsp_dt(:,:,1,1)=tmp
    ! call medfilt2d(o_F2_star_E0%nab_rp_nab_sp(:,:,1,1),tmp,cb%mz,cb%mx); o_F2_star_E0%nab_rp_nab_sp(:,:,1,1)=tmp

    contains

    subroutine medfilt2d(a,b,m,n)
        real,dimension(m,n) :: a,b
        real :: x(3,3)
        do j=2,n-1
        do i=2,m-1
            x=a(i-1:i+1,j-1:j+1)
            call bubble_sort(x,3*3)
            b(i,j)=x(2,2)
        enddo;enddo
    end subroutine

    subroutine bubble_sort(x,n)
    !resort x from min to max
        real x(n)

        do j=2,n
            do i=j,n
                if(x(i)<x(j-1))then
                    exchange=x(j-1)
                    x(j-1)=x(i)
                    x(i)=exchange
                endif
            enddo
        enddo
    end subroutine

end