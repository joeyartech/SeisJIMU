module m_math

    !math constants
    real,parameter :: r_pi = 3.1415927
    double precision,parameter :: d_pi = 3.1415926535897932d0
    !real,parameter :: r_epsilon = epsilon(r_pi)
    real,parameter :: r_epsilon=1e-3
    real,parameter :: r_eps=1e-15
    
    contains

    !sinc function
    elemental real function sinc(x)
        real,intent(in) :: x

        real pix
        pix=r_pi*x

        if (abs(pix) == 0.) then !better ways?
            sinc = 1.
        else
            sinc = sin(pix)/(pix)
        end if

    end function

    !Zero-order Modified Bessel function of the first kind
    !http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tbessi_f90.txt
    elemental real function bessi0(x)
        real, intent(in) :: x

        real ax
        double precision :: y

        double precision,parameter :: p1=1.0d0,p2=3.5156229d0,p3=3.0899424d0,p4=1.2067492d0,p5=0.2659732d0,p6=0.360768d-1,p7=0.45813d-2
        double precision,parameter :: q1=0.39894228d0,q2=0.1328592d-1,q3=0.225319d-2,q4=-0.157565d-2,q5=0.916281d-2,q6=-0.2057706d-1,q7=0.2635537d-1,q8=-0.1647633d-1,q9=0.392377d-2

        if (abs(x)<3.75) then
            y=x*x/14.0625  !=(x/3.75)^2
            bessi0=real( p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))) )
        else
            ax=abs(x)
            y=3.75/ax
            bessi0=real( (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9)))))))) )
        endif

    end function
    
end