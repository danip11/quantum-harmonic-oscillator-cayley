program Obligatorio3
  implicit none

  ! Parámetros del espacio y tiempo
  integer, parameter :: S = 1000                       ! N° puntos espaciales
  real(8), parameter :: L = 1.d0                       ! Longitud de la caja
  real(8), parameter :: dx = L/S                       ! Paso espacial
  real(8), parameter :: dt = dx*dx                     ! Paso temporal (dt = dx^2)
  real(8), parameter :: pi = acos(-1.d0)               ! Pi mediante intrínseca

  ! Número de onda fijo (para E teórica ~157)
  real(8), parameter :: k0 = 1.d0/(2.d0*sqrt(1.d-3))    ! 1/(2*sqrt(1e-3)) ~ 15.8114

  ! Variables principales
  real(8) :: t, dx2, r, norma
  real(8) :: x(0:S), P(0:S), xp, x2, pp, deltax, deltap
  real(8) :: dk, mom(0:S), P_mom(0:S)
  integer :: j, n, k
  complex(16), allocatable :: phi(:), phi2(:), q(:)
  complex(16) :: dphi, ep, phi_pk
  complex(16), allocatable :: alpha(:), beta(:), gamma(:), a2(:)

  ! Asignar vectores dimensiónales
  allocate(phi(0:S), phi2(0:S), q(0:S))
  allocate(alpha(0:S-1), beta(0:S-1), gamma(1:S-1), a2(1:S-1))

  ! Coeficiente adimensional
  dx2 = dx*dx
  r = dt/dx2

  ! Inicializar posiciones
  do j = 0, S
    x(j) = j*dx
  end do

  ! Inicializar espacio de momentos (para salida)
  dk = 2.d0*pi/(S*dx)
  do k = 0, S
    mom(k) = (k - S/2)*dk
  end do

  ! Condiciones de contorno
  phi(0) = (0.d0,0.d0)
  phi(S) = (0.d0,0.d0)

  ! Onda inicial: plana modulada por gaussiana
  do j = 1, S-1
    phi(j) = exp(cmplx(0.d0, k0*x(j))) * exp(-8.d0*((4.dble(j)-S)**2)/S**2)
  end do

  ! Normalizar estado inicial
  call normalize(phi, dx, S)
  t = 0.d0

  ! Preálculo de coeficientes Crank--Nicolson
  a2 = cmplx(-2.d0, 2.d0/r)
  alpha(S-1) = (0.d0, 0.d0)
  do j = S-1, 1, -1
    gamma(j)    = 1.d0/(a2(j) + alpha(j))
    alpha(j-1) = -gamma(j)
  end do

  ! Abrir ficheros
  open(unit=8,  file='xreim.dat')
  open(unit=9,  file='prob.dat')
  open(unit=10, file='norma.dat')
  open(unit=11, file='x_esperado.dat')
  open(unit=12, file='p_esperado.dat')
  open(unit=13, file='energia.dat')
  open(unit=14, file='incertidumbre.dat')
  open(unit=15, file='probespmom.dat')

  ! Bucle temporal
  do n = 0, 1000
    ! Resolver paso temporal con Cayley
    beta(S-1) = (0.d0, 0.d0)
    do j = S-1, 1, -1
      beta(j-1) = gamma(j) * (cmplx(0.d0,4.d0/r)*phi(j) - beta(j))
    end do

    q(0) = (0.d0,0.d0)
    do j = 1, S-1
      q(j) = alpha(j-1)*q(j-1) + beta(j-1)
    end do
    q(S) = (0.d0,0.d0)

    do j = 1, S-1
      phi2(j) = q(j) - phi(j)
    end do
    phi2(0) = (0.d0,0.d0)
    phi2(S) = (0.d0,0.d0)
    phi = phi2

    ! Normalizar para controlar desviaciones
    call normalize(phi, dx, S)

    ! Cálculos de observables y salida
    call output_observables(phi, x, S, dx, mom, 15, 8, 9, 10, 11, 12, 13, 14, t)

    t = t + dt
  end do

  ! Cerrar ficheros
  close(8); close(9); close(10); close(11)
  close(12); close(13); close(14); close(15)

contains

  subroutine normalize(u, dx, N)
    ! Normaliza vector u(0:N) de manera que sum |u|^2 dx = 1
    complex(16), intent(inout) :: u(0:N)
    real(8), intent(in) :: dx
    integer, intent(in) :: N
    real(8) :: norm2
    integer :: j
    norm2 = 0.d0
    do j = 0, N
      norm2 = norm2 + abs(u(j))**2
    end do
    norm2 = sqrt(norm2*dx)
    do j = 0, N
      u(j) = u(j)/norm2
    end do
  end subroutine normalize

  subroutine output_observables(phi, x, N, dx, mom, unit_mom, &
                                unit_xreim, unit_prob, unit_norm, &
                                unit_xexp, unit_pexp, unit_E, unit_unc, t)
    ! Calcula y escribe norma, <x>, <p>, E, Delta x * Delta p, densidad real y momentum space
    complex(16), intent(in) :: phi(0:N)
    real(8), intent(in) :: x(0:N), dx, mom(0:N), t
    integer, intent(in) :: N, unit_mom, unit_xreim, unit_prob
    integer, intent(in) :: unit_norm, unit_xexp, unit_pexp, unit_E, unit_unc
    real(8) :: P(0:N), norma, xp, x2, pp, deltax, deltap
    complex(16) :: ep, dphi, phi_pk
    integer :: j, k
    ! Momentum-space density
    do k = 0, N
      phi_pk = (0.d0,0.d0)
      do j = 0, N
        phi_pk = phi_pk + exp(-cmplx(0.d0,1.d0)*mom(k)*x(j))*phi(j)
      end do
      phi_pk = phi_pk * dx/sqrt(2.d0*pi)
      write(unit_mom,*) t, mom(k), real(phi_pk*conjg(phi_pk))
    end do
    ! Real-space amplitudes & probability
    xp = 0.d0; x2 = 0.d0; pp = 0.d0; norma = 0.d0
    do j = 0, N
      write(unit_xreim,*) x(j), real(phi(j)), aimag(phi(j))
      P(j) = abs(phi(j))**2
      write(unit_prob,*) x(j), P(j)
      norma = norma + P(j)
      xp = xp + x(j)*P(j)
      x2 = x2 + x(j)**2*P(j)
      if (j>0 .and. j<N) then
        dphi = (phi(j+1)-phi(j-1))/(2.d0*dx)
        pp = pp - aimag(conjg(phi(j))*dphi)*dx
      end if
    end do
    norma = norma*dx
    write(unit_norm,*) t, norma
    xp = xp*dx; x2 = x2*dx
    write(unit_xexp,*) t, xp
    write(unit_pexp,*) t, pp
    ep = (0.d0,0.d0)
    do j = 1, N-1
      ep = ep + conjg(phi(j))*(phi(j+1)-2.d0*phi(j)+phi(j-1))*(-1.d0)/(2.d0*dx2)*dx
    end do
    write(unit_E,*) t, real(ep)
    deltax = sqrt(x2 - xp**2)
    deltap = sqrt(real(ep) - pp**2)
    write(unit_unc,*) t, deltax*deltap
  end subroutine output_observables

end program Obligatorio3
