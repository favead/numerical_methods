program first_method_solution
    use method_fourie
    implicit none

    integer,parameter :: k = 100
    integer :: i = 0,n = 0
    real :: Bi, Fom, lambda, delta, alpha, T0, Te, c, ro, tetta = 0.0, Fo = 0.0, x_dl = 0.0
    real,dimension(100,2) :: dots1, dots2

    open(1,file='data.txt')
    read(1,*)lambda, delta, alpha, T0, Te, c, ro
    close(1)

    open(2,file='output1.txt')
    open(3,file='output2.txt')
    open(4,file='output3.txt')

    Bi = alpha * delta / 2.0 / lambda ! 0.09 1.09 68.18
    Fom = 2.0

    call tetta_series_tau(tetta, x_dl, Fo, Bi, Fom, dots1, k)
    call tetta_series_x(tetta, x_dl, Fo, Bi, dots2, k)

    call to_dim_t(dots1, k, T0, Te, lambda, delta, c, ro)
    call to_dim_x(dots2, k, T0, Te, delta)

    do i = 1,k,1
        write(2,*)dots1(i,1),dots1(i,2)
        write(3,*)dots2(i,1),dots2(i,2)
    end do
end program
