subroutine glasso(S, X, W, L, maxIt, warm, thr, n)
    implicit none
    integer, intent(in) :: n, warm, maxIt
    real(8), intent(in) :: thr
    real(8), intent(inout), dimension(1:n, 1:n) :: S, X, W, L
    real(8), dimension(1:n) :: Vd, WXj

    integer :: iter, i, j
    real(8) :: sp, thrLasso, dw, dx, a, b, c, delta, tmp, EPS

    ! f2py intent(in) :: n, thr, maxIt, warm
    ! f2py intent(inout) :: S, X, W, L

    parameter (EPS = 1.1e-16)

    sp = sum(abs(S))

    do i = 1,n
        sp = sp - abs(S(i, i))
    end do

    if (sp .eq.  0.0) then
        ! Special occasion when S is diagonal
        W = 0.0
        X = 0.0
        do i = 1,n
            W(i, i) = W(i, i) + L(i, i)
        end do

        X = 0.0
        do i = 1,n
            X(i, i) = 1.0 / max(W(i,i), EPS)
        end do
        return
    end if

    ! A more general occasion
    sp = thr * sp / (n-1)
    thrLasso = sp / n
    thrLasso = min(thrLasso, 2*EPS)

    if (warm .eq. 0) then
        ! cold strat
        W = S
        X = 0.0
    else
        ! warm start
        do i = 1,n
            X(1:n, i) = -X(1:n, i) / X(i, i)
            X(i, i) = 0
        end do
    end if

    do i = 1,n
        Vd(i) = S(i, i) + L(i, i)
        W(i, i) = Vd(i)
    end do
        
    do iter = 1, maxIt
        dw = 0.0
        do j = 1,n
            WXj(1:n) = 0.0
            do i = 1,n
                if (X(i, j) .ne. 0.0) then
                    WXj = WXj + W(:, i) * X(i ,j)
                end if
            end do

            do
                dx = 0.0
                do i = 1,n
                    if (i .ne. j) then
                        a = S(i, j) - WXj(i) + Vd(i) * X(i, j)
                        b = abs(a) - L(i, j)
                        if (b .gt. 0.0) then
                            c = sign(b, a) / Vd(i)
                        else
                            c = 0.0
                        end if
                        delta = c - X(i,j)
                        if (delta .ne. 0.0) then
                            X(i ,j) = c
                            WXj(1:n) = WXj(1:n) + W(:, i) * delta
                            dx = max(dx, abs(delta))
                        end if
                    end if
                end do
                if (dx .lt. thrLasso) then
                    exit
                end if
            end do

            WXj(j) = Vd(j)
            dw = max(dw, sum(abs(WXj(1:n) - w(:, j))))
            W(:, j) = WXj(1:n)
            W(j, :) = WXj(1:n)
        end do

        if (dw .le. sp) then
            exit
        end if
    end do

    do i = 1,n
        tmp = 1.0 / (Vd(i) - sum(X(:, i) * W(:, i)))
        X(1:n, i) = -tmp * X(1:n, i)
        X(i, i) = tmp
    end do

    do i = 1, n-1
        X(i+1:n, i) = (X(i+1:n,i) + X(i,i+1:n)) / 2
        X(i,i+1:n) = X(i+1:n,i)
    end do

    return

end subroutine glasso