module m_freesurface
use m_setup
use m_computebox

	!Levandar & Roberttson's stress image for free surface boundary condition
	!free surface is located at [1,ix,1] level

	character(:),allocatable :: FS_method

	contains

	subroutine freesurface_init()

		FS_method=setup%get_str('FS_METHOD',o_default='stress_image')

	end subroutine


	subroutine freesurface_velocity(vz,vx)
		real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: vz,vx

		select case(FS_method)
		case('stress_vel_image') !Roberttson's 2nd method
			!symmetric mirroring: vz[0.5]=vz[1.5], ie. vz(1,ix,iy)=vz(2,ix,iy) -> p(1,ix,iy)=0.
            vz(1,:,:)=vz(2,:,:)

        case('stress_image') !Roberttson's 3rd method
            vz(cb%ifz:1,:,:)=0.
            vx(cb%ifz:0,:,:)=0.

        endselect

	end subroutine

	subroutine freesurface_stress(sz,o_ss)
		real,dimension(cb%ifz:cb%ilz,cb%ifx:cb%ilx,cb%ify:cb%ily) :: sz, o_ss
		optional :: o_ss
		
        !image szz
        sz( 1,:,1)=0.
        sz(0:cb%ifz:-1, :,1)=-sz(2:2+0-cb%ifz, :,1)

        !image ss
        if(present(o_ss)) o_ss(1:cb%ifz:-1, :,1)= -o_ss(2:2+1-cb%ifz, :,1)

    end subroutine

end


!other FS ethods

!zero_stress on velocities
    ! !Δsz=0 : -(λ+2μ)∂_z vz = λ ∂ₓvx
    ! !         -(λ+2μ)[1,ix]*(vz[1.5,ix]-vz[0.5,ix])/dz = λ[1,ix]*(vx[1,ix-0.5]-vx[1,ix+0.5])/dx
    ! !         -(λ+2μ)(1,ix)*(vz(2  ,ix)-vz(1  ,ix))/dz = λ(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
    ! !          (λ+2μ)(1,ix)*(vz(1  ,ix)-vz(2  ,ix))/dz = λ(1,ix)*(vx(1,ix    )-vx(1,ix+1  ))/dx
    ! !Δss=0 : -∂_z vx = ∂ₓvz
    ! !         -(vx[2,ix+0.5]-vx[0,ix+0.5])/2dz = ( (vz[0.5,ix]+vz[1.5,ix])/2 - (vz[0.5,ix+1]+vz[1.5,ix+1])/2 )/dx
    ! !         -(vx(2,ix+1  )-vx(0,ix+1  ))/2dz = ( (vz(1  ,ix)+vz(2  ,ix))/2 - (vz(1  ,ix+1)+vz(2  ,ix+1))/2 )/dx
    ! !          (vx(0,ix+1  )-vx(2,ix+1  ))/ dz = ( (vz(1  ,ix)+vz(2  ,ix))   -  vz(1  ,ix+1)-vz(2  ,ix+1)    )/dx
    ! dz_dx = m%dz/m%dx

    ! do ix=ifx,ilx
    !     f%vz(1,ix,1)= f%vz(2,ix,1) + self%lda(1,ix)*(f%vx(1,ix,1)-f%vx(1,ix+1,1))*dz_dx/self%ldap2mu(1,ix)
    !     !f%vx(0,ix,1)= f%vx(2,ix,1) + (f%vz(1,ix,1)+f%vz(2,ix,1)-f%vz(1,ix+1,1)-f%vz(2,ix+1,1))*dz_dx/self%ldap2mu(1,ix) !toy2del says this condition is not needed
    ! enddo
!
!zero_stress on stresses
    ! !so explicit boundary condition: sz(1,ix)=0
    ! !and antisymmetric mirroring: ss[0.5,ix-0.5]=-ss[1.5,ix-0.5] -> ss(1,ix)=-ss(2,ix)
    ! sz(1,:,1)=0.
    ! ss(1,:,1)=-ss(2,:,1)


!effective medium on velocities !Mittet, Cao & Chen. but not yet working
!     !required for high-ord FD
!     f%vz(cb%ifz:1,:,1)=0.
!     f%vx(cb%ifz:0,:,1)=0.

!     do ix=ifx,ilx
!         dsz_dz_= (f%sz(2,ix,1)                )/m%dz !c1z*(f%sz(2,ix,1)-f%sz(1,ix  ,1)) +c2z*(f%sz(3,ix  ,1)-f%sz(0,ix,1))
!         dsx_dx_= (f%sx(1,ix,1)-f%sx(1,ix-1,1))/m%dx !c1x*(f%sx(1,ix,1)-f%sx(1,ix-1,1)) +c2x*(f%sx(1,ix+1,1)-f%sx(1,ix-2,1))

!         ! dss_dz_= c1z*(f%ss(2,ix,1)-f%ss(1,ix,1)) +c2z*(f%ss(3,ix,1)-f%ss(0,ix,1))
!         dss_dx_= (f%ss(2,ix+1,1)-f%ss(2,ix,1))/m%dx !c1x*(f%ss(2,ix+1,1)-f%ss(2,ix,1)) +c2x*(f%ss(2,ix+2,1)-f%ss(2,ix-1,1))
                        
!         !velocity
!         f%vz(2,ix,1)=f%vz(2,ix,1) + self%dt*   self%buoz(2,ix)*(dsz_dz_           +dss_dx_)

!         !f%vx(i)=f%vx(i) + self%dt*2.*self%buox(i)*(dss_dz_+dsx_dx_)
!         f%vx(1,ix,1)=f%vx(1,ix,1) + self%dt*2.*self%buox(1,ix)*(f%ss(2,ix,1)/m%dz +dsx_dx_)

!     enddo
!effective medium on stresses
!     !required for high-ord FD
!     ! f%ss(cb%ifz:1,:,1)=0.
!     ! f%sz(cb%ifz:0,:,1)=0.
!     f%sz( 1,:,1)=0.
!     nnz=0-cb%ifz
!     f%sz(cb%ifz:0, :,1)=-f%sz(2+nnz:2:-1, :,1)
!     nnz=1-cb%ifz
!     f%ss(cb%ifz:1, :,1)=-f%ss(2+nnz:2:-1, :,1)
!     f%sx(cb%ifz:0,:,1)=0.

!     do ix=ifx,ilx
    
!         dvx_dx_= (f%vx(1,ix+1,1)-f%vx(1,ix,1))/m%dx !c1x*(f%vx(1,ix+1,1)-f%vx(1,ix,1))  +c2x*(f%vx(1,ix+2,1)-f%vx(1,ix-1,1))
!         !dvx_dx_= c1x*(f%vx(1,ix+1,1)-f%vx(1,ix,1))  +c2x*(f%vx(1,ix+2,1)-f%vx(1,ix-1,1))
        
!         !normal stresses
!         f%sz(1,ix,1) = 0.

!         factor=-self%lda(1,ix)**2/self%ldap2mu(1,ix) + self%ldap2mu(1,ix)
!         factor=factor/2.
!         f%sx(1,ix,1) = f%sx(1,ix,1) + time_dir*self%dt * factor*dvx_dx_


!         dvz_dx_= (f%vz(2,ix,1)-f%vz(2,ix-1,1))/m%dx
!         !dvz_dx_= c1x*(f%vz(2,ix,1)-f%vz(2,ix-1,1))  +c2x*(f%vz(2,ix+1,1)-f%vz(2,ix-2,1))
!         dvx_dz_= (f%vx(2,ix,1)-f%vx(1,ix  ,1))/m%dz
!         !dvx_dz_= c1z*(f%vx(2,ix,1)-f%vx(1,ix  ,1))  +c2z*(f%vx(3,ix  ,1)-f%vx(0,ix  ,1))
        
!         !shear stress
!         f%ss(2,ix,1) = f%ss(2,ix,1) + time_dir*self%dt * self%mu(2,ix)*(dvz_dx_+dvx_dz_)
            
!     enddo
