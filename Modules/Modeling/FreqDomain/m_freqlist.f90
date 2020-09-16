module m_freqlist
use m_string
use m_mpienv
use m_message
use m_arrayop

    integer :: nfreq=1
    complex,dimension(:),allocatable :: freq

    contains

    subroutine build_freqlist
        type(t_string),dimension(:),allocatable :: s_freq
        logical :: if_complex_freq

        s_freq = partition(get_setup_char('FREQUENCIES','FREQ'))
        if_complex_freq=.false.

        if(ask_setup('COMPLEX_FREQUENCIES','CFREQ')) then
            s_freq = partition(get_setup_char('COMPLEX_FREQUENCIES','CFREQ'),';')
            if_complex_freq=.true.
        endif

        nfreq = size(s_freq)

        call alloc(freq,nfreq)

        if(if_complex_freq) then
            do i=1,nfreq
                read(s_freq(i)%string,*) freq(i)
            enddo
        else
            do i=1,nfreq
                read(s_freq(i)%string,*) tmp
                freq(i)=cmplx(tmp,0.)
            enddo
        endif

    end subroutine
end