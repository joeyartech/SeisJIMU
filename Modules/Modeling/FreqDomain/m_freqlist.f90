module m_freqlist
use m_mpienv
use m_message
use m_arrayop
use m_string

    integer :: nfreq=1
    real,dimension(:),allocatable :: freqs

    contains

    subroutine build_freqlist
        type(t_string),dimension(:),allocatable :: cfreqs

        cfreqs = partition(get_setup_char('FREQUENCIES'))

        nfreq = size(freq_array)

        call alloc(freqs,nfreq)
        do i=1,nfreq
            freqs(i)=str2real(str_array(i))
        end
    end subroutine
end