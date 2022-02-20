!Quick select algorithm
!https://rosettacode.org/wiki/Quickselect_algorithm#Fortran
!modified to Fortran 99

!find the k'th element in order of an array of n elements, not necessarily in order.
integer function findelement(k,a,n)
    integer k,n      
    integer a(n),hope,pesty
    integer l,r,l2,r2

    l = 1
    r = n             !the bounds of the work area within which the k'th element lurks.
    do while (l<r)    !keep going until it is clamped.
        hope = a(k)   !if array a is sorted, this will be rewarded.
        l2 = l        !but it probably isn't sorted.
        r2 = r        !so prepare a scan.
        do while (l2<=r2)          !keep squeezing until the inner teeth meet.
            do while (a(l2)<hope)  !pass elements less than hope.
                l2 = l2 + 1        !note that at least element a(k) equals hope.
            end do                 !raising the lower jaw.
            do while (hope<a(r2))  !elements higher than hope
            r2 = r2 - 1            !are in the desired place.
            end do                 !and so we speed past them.

            if(l2<r2) then
                pesty = a(l2)      !on grit. a(l2) > hope and a(r2) < hope.
                a(l2) = a(r2)      !so swap the two troublemakers.
                a(r2) = pesty      !to be as if they had been in the desired order all along.
                l2 = l2 + 1        !advance my teeth.
                r2 = r2 - 1        !as if they hadn't paused on this pest.
            endif
            if(l2==r2) then
                l2 = l2 + 1        !advance my teeth.
                r2 = r2 - 1        !as if they hadn't paused on this pest.
            endif

        end do        !and resume the squeeze, hopefully closing in k.

        if (r2 < k ) l = l2 !the end point gives the order position of value hope.
        if (k  < l2) r = r2 !but we want the value of order position k.

    end do            !have my teeth met yet?
        
    findelement = a(k) !yes. a(k) now has the k'th element in order.

end function !remember! array a has likely had some elements moved!