module m_arrayop
use m_mpienv
use m_setup

    integer,parameter :: max_array_size=2147483647  !=huge(integer(kind=4))

    interface alloc
        module procedure alloc_int1_ubound
        module procedure alloc_int2_ubound
        module procedure alloc_real1_ubound
        module procedure alloc_real2_ubound
        module procedure alloc_real3_ubound
        module procedure alloc_real4_ubound
        module procedure alloc_real1
        module procedure alloc_real2
        module procedure alloc_real3
    end interface

    interface add
        module procedure add_int1
    end interface

    interface rm
        module procedure rm_int1
    end interface
    
    contains
    
    subroutine alloc_int1_ubound(a,n1,old2,o_init)
        integer n1
        integer,dimension(:),allocatable :: a
        integer,dimension(:),optional :: old2
        integer,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if (allocated(a)) then
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0
        endif
        
    end
    
    subroutine alloc_int2_ubound(a,n1,n2,old2,o_init)
        integer n1,n2
        integer,dimension(:,:),allocatable :: a
        integer,dimension(:,:),optional :: old2
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
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1,n2))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0
        endif
        
    end
    
    subroutine alloc_real1_ubound(a,n1,old2,o_init)
        integer n1
        real,dimension(:),allocatable :: a
        real,dimension(:),optional :: old2
        real,optional :: o_init
        
        if(n1<1.or.n1>max_array_size) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
         
        if (allocated(a)) then
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end
    
    subroutine alloc_real2_ubound(a,n1,n2,old2,o_init)
        integer n1,n2
        real,dimension(:,:),allocatable :: a
        real,dimension(:,:),optional :: old2
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
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1,n2))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end
    
    subroutine alloc_real3_ubound(a,n1,n2,n3,old2,o_init)
        integer n1,n2,n3
        real,dimension(:,:,:),allocatable :: a
        real,dimension(:,:,:),optional :: old2
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
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1,n2,n3))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end
    
    subroutine alloc_real4_ubound(a,n1,n2,n3,n4,old2,o_init)
        integer n1,n2,n3,n4
        real,dimension(:,:,:,:),allocatable :: a
        real,dimension(:,:,:,:),optional :: old2
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
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1,n2,n3,n4))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end

    subroutine alloc_real1(a,n1,old2,o_init)
        integer,dimension(2) :: n1
        real,dimension(:),allocatable :: a
        real,dimension(:),optional :: old2
        real,optional :: o_init
        
        if(n1(1)>n1(2)) then
            if(mpiworld%is_master) write(*,*) 'ERROR: invalid required array size! n1=',n1
            error stop
        endif
        
        if (allocated(a)) then
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1(1):n1(2)))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end
    
    subroutine alloc_real2(a,n1,n2,old2,o_init)
        integer,dimension(2) :: n1,n2
        real,dimension(:,:),allocatable :: a
        real,dimension(:,:),optional :: old2
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
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1(1):n1(2),n2(1):n2(2)))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end
    
    subroutine alloc_real3(a,n1,n2,n3,old2,o_init)
        integer,dimension(2) :: n1,n2,n3
        real,dimension(:,:,:),allocatable :: a
        real,dimension(:,:,:),optional :: old2
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
            if(present(old2)) then
                old2=a
            endif
            deallocate(a)
        endif
        allocate(a(n1(1):n1(2),n2(1):n2(2),n3(1):n3(2)))
        
        if(present(o_init)) then
            a=o_init
        else
            a=0.
        endif
        
    end
    
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

    subroutine add_int1(a,n,item)
        integer,dimension(:),allocatable :: a
        integer item

        integer,dimension(:),allocatable :: b
        
        allocate(b(n+1)); b(1:n)=a(1:n)
        b(n+1)=item
        
        deallocate(a)
        a=b
        deallocate(b)

    end subroutine

    subroutine rm_int1(a,n,item)
        integer,dimension(:),allocatable :: a
        integer :: item

        integer,dimension(:),allocatable :: b

        allocate(b(n))

        j=0
        do i=1,n
            if(a(i)/=item) then
                j=j+1
                b(j)=a(i)
            endif
        enddo

        deallocate(a)
        allocate(a(j)); a(:)=b(1:j)

        deallocate(b)

    end subroutine

    ! function unify_int1
    ! end function

end
