module m_querypoint
use m_System
use m_parametrizer

    private

    type,public :: t_querypoint
        real,dimension(:,:,:,:),allocatable :: x  !point in parameter space X (unitless)
        real,dimension(:,:,:,:),allocatable :: g  !gradient, unit in [1/Nm]
        real,dimension(:,:,:,:),allocatable :: pg !preconditioned gradient
        real,dimension(:,:,:,:),allocatable :: d  !descent direction
        real :: f !objective function value, unit in [Nm]
        real :: gdotd  !projection of g onto d

        contains
        procedure :: init
        final :: fin

    end type

    contains

    subroutine init(self)
        class(t_querypoint) :: self

        call alloc(self%x, param%n1,param%n2,param%n3,param%npars)
        ! call alloc(self%g, param%n1,param%n2,param%n3,param%npars)
        ! call alloc(self%pg,param%n1,param%n2,param%n3,param%npars)
        ! call alloc(self%d, param%n1,param%n2,param%n3,param%npars)

        self%f=0.
        ! self%gdotd=0.

        call param%transform('m->x',o_x=self%x)
        ! call param%transform(o_g=self%g)

    end subroutine

    subroutine fin(self)
        type(t_querypoint) :: self

        call dealloc(self%x,self%g,self%pg,self%d)

    end subroutine

end