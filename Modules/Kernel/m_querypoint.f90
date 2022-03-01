module m_querypoint
use m_System
use m_parametrizer

    private

    type,public :: t_querypoint
        character(:),allocatable :: name

        real,dimension(:,:,:,:),allocatable :: x  !point in parameter space X, unit in [1]
        real,dimension(:,:,:,:),allocatable :: g  !gradient, in [Nm]
        real,dimension(:,:,:,:),allocatable :: pg !preconditioned gradient, in [Nm]
        real,dimension(:,:,:,:),allocatable :: d  !descent direction, same unit as x
        real :: f !objective function value, unit in [Nm]
        real :: g_dot_d !projection of g onto d, in [Nm]

        logical :: is_fitting_data

        contains
        procedure :: init
        procedure :: set_sign
        final :: fin

        procedure :: is_registered
        procedure :: register

    end type

    contains

    subroutine init(self,name)
        class(t_querypoint) :: self
        character(*) :: name

        self%name=name

        call alloc(self%x, param%n1,param%n2,param%n3,param%npars)
        ! call alloc(self%g, param%n1,param%n2,param%n3,param%npars)
        ! call alloc(self%pg,param%n1,param%n2,param%n3,param%npars)
        ! call alloc(self%d, param%n1,param%n2,param%n3,param%npars)

        self%f=0.
        ! self%gdotd=0.

        call param%transform('m->x',o_x=self%x)
        ! call param%transform(o_g=self%g)

    end subroutine

    subroutine set_sign(self,o_sign,o_reference)
        class(t_querypoint) :: self
        character,optional :: o_sign
        real,optional :: o_reference
        
        logical :: if_flip
        
        if_flip=.false.
        
        if(present(o_sign)) then
            if(o_sign=='+' .and. self%f<0.) then
                if_flip=.true.
            elseif(o_sign=='-' .and. self%f>0.) then
                if_flip=.true.
            endif
        endif
        
        if(if_flip) then
            self%f  = -self%f
            if(allocated(self%g )) self%g  = -self%g
            if(allocated(self%pg)) self%pg = -self%pg
            if(allocated(self%d )) self%d  = -self%d
        endif
        
    end subroutine

    subroutine fin(self)
        type(t_querypoint) :: self

        call dealloc(self%x,self%g,self%pg,self%d)

    end subroutine


    
    logical function is_registered(self,chp)
        class(t_querypoint) :: self
        type(t_checkpoint) :: chp

        is_registered=chp%check(self%name)
        if(.not.is_registered) return

        call chp%open(self%name)
        call chp%read(self%x,self%g,self%pg)
        call chp%read(self%f)
        call chp%close
        call hud('Read '//self%name//' from '//chp%name//', size='//num2str(total_size(self%x,self%g,self%pg)+4.))

    end function

    subroutine register(self,chp)
        class(t_querypoint) :: self
        type(t_checkpoint) :: chp

        call chp%open(self%name)
        call chp%write(self%x,self%g,self%pg)
        call chp%write(self%f)
        call chp%close

    end subroutine

end
