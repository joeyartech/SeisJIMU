!2D median filter from ChatGPT
!using quickselect, more efficient for sorting
!within large windows (size>=10)
module m_median2d

	contains

	function median2d(in,hw) result(out)
		integer :: hw !half window width
		real,dimension(:,:) :: in

		real,dimension(:,:),allocatable :: out

		integer :: i,j,ii,jj
		integer :: i1,i2,j1,j2
		integer :: n,k
		real,allocatable :: work(:)

		nx=size(in,1)
		ny=size(in,2)
		allocate(out(nx,ny))

		allocate(work((2*hw+1)*(2*hw+1)))

		do j=1,ny
	    do i=1,nx

	        i1=max(1,i-hw)
	        i2=min(nx,i+hw)

	        j1=max(1,j-hw)
	        j2=min(ny,j+hw)

	        n=0
	        do jj=j1,j2
            do ii=i1,i2
                n=n+1
                work(n)=in(ii,jj)
            end do
	        end do

	        k=(n+1)/2
	        out(i,j)=quickselect(work,n,k)

	    end do
		end do

		deallocate(work)

	end function

	real function quickselect(a,n,k)
		integer,intent(in) :: n,k
		real,intent(inout) :: a(n)

		integer :: left,right,pivot

		left=1
		right=n

		do

		    if(left==right) then
		        quickselect=a(left)
		        return
		    endif

		    pivot=(left+right)/2

		    call partition(a,left,right,pivot,pivot)

		    if(k==pivot) then
		        quickselect=a(k)
		        return
		    elseif(k<pivot) then
		        right=pivot-1
		    else
		        left=pivot+1
		    endif

		end do

	end function

	subroutine partition(a,left,right,pivot,newpivot)
		integer,intent(in) :: left,right,pivot
		integer,intent(out):: newpivot
		real,intent(inout):: a(:)

		real :: piv,tmp
		integer :: i,store

		piv=a(pivot)

		tmp=a(pivot)
		a(pivot)=a(right)
		a(right)=tmp

		store=left

		do i=left,right-1
		    if(a(i)<piv) then
		        tmp=a(store)
		        a(store)=a(i)
		        a(i)=tmp
		        store=store+1
		    endif
		enddo

		tmp=a(store)
		a(store)=a(right)
		a(right)=tmp

		newpivot=store

	end subroutine

end