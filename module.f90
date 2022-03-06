module method_fourie
    implicit none
    contains

    subroutine tetta_series_tau(tetta,x_dl,Fo,Bi,Fom,dots,k)
        implicit none
        real :: x_dl, tetta, tettaN, Fo, dFo, e1 = 1E-6, Bi, mun, Fom, An
        real,dimension(k,2), intent(out) :: dots
        integer :: n, k, i = 1
        dFo = Fom/(k-1)
        Fo = -dFo + 0.005
        do while(Fo < Fom)
            Fo = Fo + dFo
            n = 0
            tettaN = e1 + 1
            tetta = 0
            do while(abs(tettaN) > e1)
                n = n + 1
                mun = binary_search_(n,Bi)
                An = fmun_(mun)
                tettaN = An * exp(-(mun**2)*Fo)*cos(mun*x_dl)
                tetta = tetta + tettaN
            end do
            dots(i,1) = tetta
            dots(i,2) = Fo
            i = i + 1
        end do
    end subroutine

    subroutine tetta_series_x(tetta,x_dl,Fo,Bi,dots,k)
        implicit none
        real :: x_dl, tetta, tettaN, Fo, dx, e1 = 1E-6, Bi, mun, An
        real,dimension(k,2), intent(out) :: dots
        integer :: n = 0,k,i = 1
        dx = 1.0 / (k-1)
        x_dl = -dx
        do while(x_dl < 1.0)
            x_dl = x_dl + dx
            n = 0
            tettaN = e1 + 1
            tetta = 0
            do while(abs(tettaN) > e1)
                n = n + 1
                mun = binary_search_(n,Bi)
                An = fmun_(mun)
                tettaN = An * exp(-(mun**2)*Fo)*cos(mun*x_dl)
                tetta = tetta + tettaN
            end do
            dots(i,1) = tetta
            dots(i,2) = x_dl
            i = i + 1
        end do
    end subroutine

    real function binary_search_(n,Bi)
        implicit none
        integer :: n
        real :: a , b, pi = 4 * atan(1.0), e2, c, Bi
        e2 = 1E-5
        a = pi * (n - 1)
        b = pi * (n - 0.5)
        do while (abs(b - a) > e2)
            c = ( a + b ) / 2
            if( f_(c,Bi) * f_(b,Bi) > 0) then
                b = c
            else
                a = c
            endif
        end do
        binary_search_ = c
    end function

    real function fmun_(mun)
        implicit none
        real :: mun
        fmun_ = 2 * sin(mun) / ( mun + (sin(mun) * cos(mun)))
    end function

    real function f_(x, Bi)
        real:: x,Bi
        f_ = 1/tan(x) - x/Bi
    end function

    subroutine to_dim_t(dots, k, T0, Te, lambda, delta, c, ro)
        real,dimension(k,2),intent(out) :: dots
        real :: T0, Te, lambda, delta, c, ro
        integer :: i = 0,k
        do i = 1,k,1
            dots(i,1) = dots(i,1) * (T0 - Te) + Te
            dots(i,2) = c*ro*dots(i,2)*(delta/2)**2/lambda
        end do
    end subroutine

    subroutine to_dim_x(dots, k, T0, Te, delta)
        real,dimension(k,2),intent(out) :: dots
        real :: T0,Te,delta
        integer :: i = 0,k
        do i=1,k,1
            dots(i,1) = dots(i,1) * (T0 - Te) + Te
            dots(i,2) = dots(i,2) * delta / 2.0
        end do
    end subroutine

end module
