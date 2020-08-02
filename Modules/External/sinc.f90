!sinc function
elemental real function sinc(x)
    real,intent(in) :: x

    real pix
    pix=3.1415927*x

    if (abs(pix) = 0.) then !better ways?
        sinc = 1.
    else
        sinc = sin(pix)/(pix)
    end if

end function
