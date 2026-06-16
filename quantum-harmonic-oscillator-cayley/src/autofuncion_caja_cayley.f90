program Obligatorio3
  implicit none
  ! Declaración de variables
  real(8) :: dx, dt, w, r, x(0:1000), pi, t, P(0:1000), norma, P_mom(0:1000), mom(0:1000)
  real(8) :: xp, pp, deltax, deltap, x2, dk
  integer :: L, S, j, n, k
  complex(16) :: phi(0:1000), alpha(0:999), beta(0:999), gamma(1:999)
  complex(16) :: phi2(0:1000), q(0:1000), dphi, ep, phi_p(0:1000)
  complex(16) :: a2(1:999), phi_pk

  ! Abrimos los ficheros necesarios
  open(8,  file='xreim.dat')
  open(9,  file='prob.dat')
  open(10, file='norma.dat')
  open(11, file='x_esperado.dat')
  open(12, file='p_esperado.dat')
  open(13, file='energia.dat')
  open(14, file='incertidumbre.dat')
  open(15, file='probespmom.dat')

  ! Asignamos los parámetros iniciales
  S = 1000
  L = 1
  dx = 1.d-3
  dt = 1.d-3
  pi = 3.1415926535d0
  t = 0.d0

  ! Parámetro adimensional
  w = dt/(dx**2)
  r = w

  ! Vector de posiciones x(j)
  do j = 0, S
     x(j) = dble(j)*dx
  end do

  ! Espacio de momentos
     dk = 2.d0*pi/(S*dx)
     do k = 0, S
        mom(k) = (dble(k) - dble(S)/2.d0) * dk
     end do

  ! Condiciones de contorno en la onda
  phi(0) = (0.d0, 0.d0)
  phi(S) = (0.d0, 0.d0)

  ! Función de onda inicial (n = 1)
  do j = 0, S
     if (j.eq.0 .or. j.eq.S) then
         phi(j) = (0.d0, 0.d0)
     else
         phi(j) = cmplx(sqrt(2.d0)*sin(1.d0*pi*x(j)), 0.d0)
     endif
  end do

  ! Pre-cálculo de coeficientes para Cayley
  do j = 1, S-1
     a2(j) = cmplx(-2.d0, 2.d0/r)   
  end do
  alpha(S-1) = (0.d0, 0.d0)
  do j = S-1, 1, -1
     gamma(j)    = 1.d0 / (a2(j) + alpha(j))
     alpha(j-1) = -gamma(j)
  end do

  ! Evolución temporal (Cayley)
  do n = 0, 1000
     beta(S-1) = (0.d0, 0.d0)
     do j = S-1, 1, -1
         beta(j-1) = gamma(j) * (cmplx(0.d0, 4.d0/w) * phi(j) - beta(j))
     end do

     q(0) = (0.d0, 0.d0)
     do j = 1, S-1
         q(j) = alpha(j-1)*q(j-1) + beta(j-1)
     end do
     q(S) = (0.d0, 0.d0)

     ! Actualización de la onda phi
     do j = 1, S-1
         phi2(j) = q(j) - phi(j)
     end do
     phi2(0) = (0.d0, 0.d0)
     phi2(S) = (0.d0, 0.d0)
     phi = phi2


     ! Cálculo de probabilidad en espacio de momentos
     do k = 0, S
        phi_pk = (0.d0, 0.d0)
        do j = 0, S
           phi_pk = phi_pk + exp(-cmplx(0.d0,1.d0)*mom(k)*x(j)) * phi(j)
        end do
        phi_pk = phi_pk * dx / sqrt(2.d0*pi)
        P_mom(k) = real(phi_pk * conjg(phi_pk))
        write(15,*) t, mom(k), P_mom(k)
     end do

     ! Escritura de x, Re(phi), Im(phi)
     do j = 0, S
         write(8,*) x(j), real(phi(j)), aimag(phi(j))
     end do

     ! Densidad de probabilidad
     do j = 0, S
         P(j) = real(phi(j))**2 + aimag(phi(j))**2
         write(9,*) x(j), P(j)
     end do

     ! Norma 
     norma = 0.5d0*((P(0)) + (P(S)))
     do j = 1, S-1
         norma = norma + P(j)
     end do
     norma = norma*dx
     write(10,*) t, norma

     ! Valores esperados <x>, <p>
     xp = 0.5d0*(x(0)*P(0) + x(S)*P(S))
     x2 = 0.5d0*(x(0)**2*P(0) + x(S)**2*P(S))
     do j = 1, S-1
         xp = xp + x(j)*P(j)
         x2 = x2 + x(j)**2*P(j)
     end do
     xp = xp*dx
     x2 = x2*dx
     write(11,*) t, xp

     pp = 0.d0
     do j = 1, S-1
         dphi = (phi(j+1)-phi(j-1))/(2.d0*dx)
         pp = pp - aimag(conjg(phi(j))*dphi)*dx
     end do
     write(12,*) t, pp

     ! Energía <H> = ∫ φ* φ'' dx
     ep = (0.d0, 0.d0)
     do j = 1, S-1
         ep = ep + conjg(phi(j))*(phi(j+1)+phi(j-1)-2.d0*phi(j))*(-1.d0)/(dx**2)*dx
     end do
     write(13,*) t, real(ep)

     ! Incertidumbre Δx·Δp
     deltax = sqrt(x2 - xp**2)
     deltap = sqrt(real(ep) - pp**2)
     write(14,*) t, deltax*deltap

     ! Incremento de tiempo e inserción de bloques
     t = t + dt
     write(8,*)
     write(8,*)
     write(9,*)
     write(9,*)
  end do

  ! Cierre de ficheros
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
end program Obligatorio3
