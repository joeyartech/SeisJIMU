module m_arrayop
use m_either
use m_mpienv

    integer,parameter :: max_array_size=2147483647  !=huge(integer(kind=4))

    interface alloc
        module procedure alloc_int1_ubound
        module procedure alloc_int2_ubound
        module procedure alloc_real1_ubound
        module procedure alloc_real2_ubound
        module procedure alloc_real3_ubound
        module procedure alloc_real4_ubound
        module procedure alloc_real5_ubound
        module procedure alloc_real1
        module procedure alloc_real2
        module procedure alloc_real3
    end interface

    interface dealloc
        module procedure dealloc_real1
        module procedure dealloc_real2
        module procedure dealloc_real3
        module procedure dealloc_real4
    end interface

    interface add
        module procedure add_int1
    end interface

    interface rm
        module procedure rm_int1
    end interface

    interface rsign
        module procedure rsign_real2
    end interface

    interface total_size
        module procedure total_size_real2
        module procedure total_size_real3
        module procedure total_size_real4
    end interface
    
    contains
    
    subroutine alloc_int1_ubound(a,n1,old2,oif_protect,o_init)
        integer n1
        integer,dimension(:),allocatable :: a
        integer,dimension(:),optional :: old2
        logical,optional :: oif_protect
        integer,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1),source=either(o_init,0,present(o_init)))
        
    end subroutine
    
    subroutine alloc_int2_ubound(a,n1,n2,old2,oif_protect,o_init)
        integer n1,n2
        integer,dimension(:,:),allocatable :: a
        integer,dimension(:,:),optional :: old2
        logical,optional :: oif_protect
        integer,optional :: o_init

        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2<1.or.n2>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif
        
        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1,n2),source=either(o_init,0,present(o_init)))
        
    end subroutine
    
    subroutine alloc_real1_ubound(a,n1,old2,oif_protect,o_init)
        integer n1
        real,dimension(:),allocatable :: a
        real,dimension(:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
         
        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1),source=either(o_init,0.,present(o_init)))
        
    end subroutine
    
    subroutine alloc_real2_ubound(a,n1,n2,old2,oif_protect,o_init)
        integer n1,n2
        real,dimension(:,:),allocatable :: a
        real,dimension(:,:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2<1.or.n2>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif
        
        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1,n2),source=either(o_init,0.,present(o_init)))
        
    end subroutine
    
    subroutine alloc_real3_ubound(a,n1,n2,n3,old2,oif_protect,o_init)
        integer n1,n2,n3
        real,dimension(:,:,:),allocatable :: a
        real,dimension(:,:,:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2<1.or.n2>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif
        
        if(n3<1.or.n3>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif

        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1,n2,n3),source=either(o_init,0.,present(o_init)))
        
    end subroutine
    
    subroutine alloc_real4_ubound(a,n1,n2,n3,n4,old2,oif_protect,o_init)
        integer n1,n2,n3,n4
        real,dimension(:,:,:,:),allocatable :: a
        real,dimension(:,:,:,:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2<1.or.n2>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif
        
        if(n3<1.or.n3>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n3=',n3
            error stop
        endif
        
        if(n4<1.or.n4>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n4=',n4
            error stop
        endif

        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1,n2,n3,n4),source=either(o_init,0.,present(o_init)))
        
    end subroutine

    subroutine alloc_real5_ubound(a,n1,n2,n3,n4,n5,old2,oif_protect,o_init)
        integer n1,n2,n3,n4,n5
        real,dimension(:,:,:,:,:),allocatable :: a
        real,dimension(:,:,:,:,:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2<1.or.n2>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif
        
        if(n3<1.or.n3>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n3=',n3
            error stop
        endif
        
        if(n4<1.or.n4>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n4=',n4
            error stop
        endif

        if(n5<1.or.n5>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n4=',n4
            error stop
        endif

        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1,n2,n3,n4,n5),source=either(o_init,0.,present(o_init)))
        
    end subroutine

    subroutine alloc_real1(a,n1,old2,oif_protect,o_init)
        integer,dimension(2) :: n1
        real,dimension(:),allocatable :: a
        real,dimension(:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1(1)>n1(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1(1):n1(2)),source=either(o_init,0.,present(o_init)))
        
    end subroutine
    
    subroutine alloc_real2(a,n1,n2,old2,oif_protect,o_init)
        integer,dimension(2) :: n1,n2
        real,dimension(:,:),allocatable :: a
        real,dimension(:,:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1(1)>n1(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2(1)>n2(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif

        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return

            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1(1):n1(2),n2(1):n2(2)),source=either(o_init,0.,present(o_init)))
        
    end subroutine
    
    subroutine alloc_real3(a,n1,n2,n3,old2,oif_protect,o_init)
        integer,dimension(2) :: n1,n2,n3
        real,dimension(:,:,:),allocatable :: a
        real,dimension(:,:,:),optional :: old2
        logical,optional :: oif_protect
        real,optional :: o_init
        
        if(n1(1)>n1(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if(n2(1)>n2(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n2=',n2
            error stop
        endif
        
        if(n3(1)>n3(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n3=',n3
            error stop
        endif
        
        if (allocated(a)) then
            if(either(oif_protect,.false.,present(oif_protect))) return
            
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif

        allocate(a(n1(1):n1(2),n2(1):n2(2),n3(1):n3(2)),source=either(o_init,0.,present(o_init)))
        
    end subroutine

    subroutine dealloc_real1(a,b,c,d,e,f,g)
        real,dimension(:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) deallocate(a)
        if(present(b)) then
            if(allocated(b)) deallocate(b)
        endif
        if(present(c)) then
            if(allocated(c)) deallocate(c)
        endif
        if(present(d)) then
            if(allocated(d)) deallocate(d)
        endif
        if(present(e)) then
            if(allocated(e)) deallocate(e)
        endif
        if(present(f)) then
            if(allocated(f)) deallocate(f)
        endif
        if(present(g)) then
            if(allocated(g)) deallocate(g)
        endif

    end subroutine

    subroutine dealloc_real2(a,b,c,d,e,f,g)
        real,dimension(:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) deallocate(a)
        if(present(b)) then
            if(allocated(b)) deallocate(b)
        endif
        if(present(c)) then
            if(allocated(c)) deallocate(c)
        endif
        if(present(d)) then
            if(allocated(d)) deallocate(d)
        endif
        if(present(e)) then
            if(allocated(e)) deallocate(e)
        endif
        if(present(f)) then
            if(allocated(f)) deallocate(f)
        endif
        if(present(g)) then
            if(allocated(g)) deallocate(g)
        endif

    end subroutine

    subroutine dealloc_real3(a,b,c,d,e,f,g)
        real,dimension(:,:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) deallocate(a)
        if(present(b)) then
            if(allocated(b)) deallocate(b)
        endif
        if(present(c)) then
            if(allocated(c)) deallocate(c)
        endif
        if(present(d)) then
            if(allocated(d)) deallocate(d)
        endif
        if(present(e)) then
            if(allocated(e)) deallocate(e)
        endif
        if(present(f)) then
            if(allocated(f)) deallocate(f)
        endif
        if(present(g)) then
            if(allocated(g)) deallocate(g)
        endif

    end subroutine

    subroutine dealloc_real4(a,b,c,d,e,f,g)
        real,dimension(:,:,:,:),allocatable :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) deallocate(a)
        if(present(b)) then
            if(allocated(b)) deallocate(b)
        endif
        if(present(c)) then
            if(allocated(c)) deallocate(c)
        endif
        if(present(d)) then
            if(allocated(d)) deallocate(d)
        endif
        if(present(e)) then
            if(allocated(e)) deallocate(e)
        endif
        if(present(f)) then
            if(allocated(f)) deallocate(f)
        endif
        if(present(g)) then
            if(allocated(g)) deallocate(g)
        endif
        
    end subroutine
    
    ! subroutine flatten_alias(a,n1,p)
    !     integer,intent(in),dimension(2) :: n1
    !     real,intent(in),dimension(n1(1):n1(2)),target :: a
    !     real,intent(out),dimension(:),pointer :: p
    !     p=>a(1:)
    ! end subroutine
    
    ! subroutine inflate_alias(a,n1,n2,n3,p)
    !     integer,intent(in),dimension(2) :: n1,n2,n3
    !     real,intent(in),dimension(n1(1):n1(2),n2(1):n2(2),n3(1):n3(2)),target :: a
    !     real,intent(out),dimension(:,:,:),pointer :: p
    !     p=>a
    ! end subroutine

    subroutine add_int1(a,item)
        integer,dimension(:),allocatable :: a
        integer item

        integer,dimension(:),allocatable :: b
        
        n=size(a)

        allocate(b(n+1))
        b(1:n)=a(1:n)
        b(n+1)=item
        
        deallocate(a)
        a=b
        deallocate(b)

    end subroutine

    subroutine rm_int1(a,item)
        integer,dimension(:),allocatable :: a
        integer :: item

        integer,dimension(:),allocatable :: b

        n=size(a)
        allocate(b(n))

        j=0
        do i=1,n
            if(a(i)/=item) then
                j=j+1
                b(j)=a(i)
            endif
        enddo

        deallocate(a)
        allocate(a(j))
        a(:)=b(1:j)
        deallocate(b)

    end subroutine

    ! function unify_int1
    ! end function


    pure function rsign_real2(a) result(res)
        real,dimension(:,:),intent(in) :: a
        real,dimension(:,:),allocatable :: res
        
        allocate(res,source=a)

        where (a>0.)
            res=1.
        elsewhere (a<0.)
            res=-1.
        elsewhere
            res=0.
        endwhere
    end function

    pure function total_size_real2(a,b,c,d,e,f,g) result(n)
        real,dimension(:,:),allocatable,intent(in) :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) n=size(a)
        if(present(b)) then
            if(allocated(b)) n=n+size(b)
        endif
        if(present(c)) then
            if(allocated(c)) n=n+size(c)
        endif
        if(present(d)) then
            if(allocated(d)) n=n+size(d)
        endif
        if(present(e)) then
            if(allocated(e)) n=n+size(e)
        endif
        if(present(f)) then
            if(allocated(f)) n=n+size(f)
        endif
        if(present(g)) then
            if(allocated(g)) n=n+size(g)
        endif

    end function

    pure function total_size_real3(a,b,c,d,e,f,g) result(n)
        real,dimension(:,:,:),allocatable,intent(in) :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) n=size(a)
        if(present(b)) then
            if(allocated(b)) n=n+size(b)
        endif
        if(present(c)) then
            if(allocated(c)) n=n+size(c)
        endif
        if(present(d)) then
            if(allocated(d)) n=n+size(d)
        endif
        if(present(e)) then
            if(allocated(e)) n=n+size(e)
        endif
        if(present(f)) then
            if(allocated(f)) n=n+size(f)
        endif
        if(present(g)) then
            if(allocated(g)) n=n+size(g)
        endif

    end function

    pure function total_size_real4(a,b,c,d,e,f,g) result(n)
        real,dimension(:,:,:,:),allocatable,intent(in) :: a,b,c,d,e,f,g
        optional :: b,c,d,e,f,g
        
            if(allocated(a)) n=size(a)
        if(present(b)) then
            if(allocated(b)) n=n+size(b)
        endif
        if(present(c)) then
            if(allocated(c)) n=n+size(c)
        endif
        if(present(d)) then
            if(allocated(d)) n=n+size(d)
        endif
        if(present(e)) then
            if(allocated(e)) n=n+size(e)
        endif
        if(present(f)) then
            if(allocated(f)) n=n+size(f)
        endif
        if(present(g)) then
            if(allocated(g)) n=n+size(g)
        endif

    end function


end
