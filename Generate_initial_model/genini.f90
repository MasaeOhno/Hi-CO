
    program GENINI

    implicit none

    integer :: i, j, k, m, n, ih

    real(8), parameter :: rDNA = 1.0d0                      ! radius of DNA (nm)
    real(8), parameter :: bp_bead = 2.0d0*rDNA/0.34d0       ! ~5.88 bp/bead
    real(8), parameter :: ncp = 133.0d0                     ! DNA length actually wrapping histones
    real(8), parameter :: wrap = 1.65d0                     ! turns wrapping histone
    integer, parameter :: bead_core = nint(ncp/bp_bead)     ! 23 beads/(histone core DNA)
    real(8), parameter :: pitch = 2.39                      ! pitch of DNA wrapping histone is 2.39 nm (JPCB, 2009, 113, 2639-2646)
    real(8), parameter :: rHD = 4.18d0                      ! Distance between center of histone and wrapping DNA 
                                                            ! (Luger K, Mader AW, Richmond RK, Sargent DF, Richmond TJ (1997) Crystal structure of
                                                            ! the nucleosome core particle at 2.8 angstrom resolution. Nature 389:251-260)
    real(8), parameter :: rHCP = 3.0d0                      ! Radius of bead for histone core protein
    real(8), parameter :: rHC = rHD - rDNA - rHCP + 1.0d0   ! Distance between center of histone and histone core protein
    real(8), parameter :: PI = acos(-1.0d0)

    integer, parameter :: nHCD = bead_core ! Number of beads for histone core DNA
    integer, parameter :: nHCP = 4         ! Number of beads for histone core protein

    integer :: nnuc                        ! Number of nucleosomes
    integer :: np                          ! Number of particles

    real(8), allocatable :: x(:,:)
    real(8) :: xi(3), xj(3), xk(3), x0(3)
    real(8) :: dt, dz, ct, cz
    real(8) :: dist, theta, phi
    real(8) :: beta = 0.0d0
    real(8) :: dHC(nHCP), tHC(nHCP), pHC(nHCP)

    character(80) :: argc


    !--- Read input

    i = iargc()
    if ( i == 1 ) then
      call getarg(1,argc)
      read(argc,*) nnuc
    else
      write(6,*) "[USAGE] ./genini.exe <N>"
      write(6,*) "where <N> is the number of nucleosomes"
      stop
    endif


    !--- Calculate coordinates

    np = nnuc*(nHCD+nHCP) 
    allocate(x(3,-7:np))

    x = 0.0d0
    dt = 2.0d0*PI/(ncp/wrap) * bp_bead
    dz = pitch/(ncp/wrap) * bp_bead     ! pitch of DNA wrapping histone is 2.39 nm (JPCB, 2009, 113, 2639-2646)
    x(1:3,-3) = (/ rHD*sin(0.0d0*dt), rHD*cos(0.0d0*dt), 0.0d0*dz /)
    x(1:3,-2) = (/ rHD*sin(1.0d0*dt), rHD*cos(1.0d0*dt), 1.0d0*dz /)
    x(1:3,-1) = (/ rHD*sin(2.0d0*dt), rHD*cos(2.0d0*dt), 2.0d0*dz /)
    x(1:3, 0) = (/ rHD*sin(3.0d0*dt), rHD*cos(3.0d0*dt), 3.0d0*dz /)
    xi = x(:,-3)
    xj = x(:,-2)
    call get_bond(xi,xj,dist)
    xk = x(:,-1)
    call get_angle(xi,xj,xk,theta)
    x0 = x(:, 0)
    call get_torsion(xi,xj,xk,x0,phi)

    cz = (bead_core-1)*dz/2.0d0
    ct = (4.0d0*PI-(bead_core-1)*dt)/2.0d0 + (bead_core-1)*dt
    x0(1:3) = (/ rHC*sin(0.00d0*PI + ct), rHC*cos(0.00d0*PI + ct), cz /)
    call get_bond(xk,x0,dHC(1))
    call get_angle(xj,xk,x0,tHC(1))
    call get_torsion(xi,xj,xk,x0,pHC(1))
    x0(1:3) = (/ rHC*sin(0.50d0*PI + ct), rHC*cos(0.50d0*PI + ct), cz /)
    call get_bond(xk,x0,dHC(2))
    call get_angle(xj,xk,x0,tHC(2))
    call get_torsion(xi,xj,xk,x0,pHC(2))
    x0(1:3) = (/ rHC*sin(1.00d0*PI + ct), rHC*cos(1.00d0*PI + ct), cz /)
    call get_bond(xk,x0,dHC(3))
    call get_angle(xj,xk,x0,tHC(3))
    call get_torsion(xi,xj,xk,x0,pHC(3))
    x0(1:3) = (/ rHC*sin(1.50d0*PI + ct), rHC*cos(1.50d0*PI + ct), cz /)
    call get_bond(xk,x0,dHC(4))
    call get_angle(xj,xk,x0,tHC(4))
    call get_torsion(xi,xj,xk,x0,pHC(4))

    do i=-3, 0
      x(:,i-nHCP) = x(:,i)
    enddo
    x(1:3,-nHCP) = (/ 0.0d0, 0.0d0, 1.0d0 /)

    m = 0

    do n=1, nnuc

      ih = 1
      m = m + 1
      i = m - nHCP - 1
      j = m - nHCP - 2
      k = m - nHCP - 3
      xi(:) = x(:,i)
      xj(:) = x(:,j)
      xk(:) = x(:,k)
      x0(:) = 0.0d0
      call build(xi,xj,xk,2.0d0*rDNA,180.0d0,-999.0d0,x0)
      x(:,m) = x0(:)

      ih = 2
      m = m + 1
      i = m - 1
      j = m - nHCP - 2
      k = m - nHCP - 4
      xi(:) = x(:,i)
      xj(:) = x(:,j)
      xk(:) = x(:,k)
      x0(:) = 0.0d0
      call build(xi,xj,xk,dist,theta,beta,x0)
      x(:,m) = x0(:)

      ih = 3
      m = m + 1
      i = m - 1
      j = m - 2
      k = m - nHCP - 3
      xi(:) = x(:,i)
      xj(:) = x(:,j)
      xk(:) = x(:,k)
      x0(:) = 0.0d0
      call build(xi,xj,xk,dist,theta,phi,x0)
      x(:,m) = x0(:)

      do ih=4, nHCD

        m = m + 1
        i = m - 1 
        j = m - 2
        k = m - 3
        xi(:) = x(:,i)
        xj(:) = x(:,j)
        xk(:) = x(:,k)
        x0(:) = 0.0d0
        call build(xi,xj,xk,dist,theta,phi,x0)
        x(:,m) = x0(:)

      enddo

      do ih=1, nHCP

        m = m + 1
        i = m - nHCD + 3 - ih
        j = m - nHCD + 2 - ih
        k = m - nHCD + 1 - ih
        xi(:) = x(:,i)
        xj(:) = x(:,j)
        xk(:) = x(:,k)
        x0(:) = 0.0d0
        call build(xi,xj,xk,dHC(ih),tHC(ih),pHC(ih),x0)
        x(:,m) = x0(:)

      enddo

    enddo


    !--- Shift center of mass
    x0 = 0.0d0
    do i=1, np
      x0(:) = x0(:) + x(:,i)
    enddo
    x0 = x0/np
    do i=1, np
      x(:,i) = x(:,i)-x0(:)
    enddo


    !--- Output PQR file
    open(1,file='out.pqr')
    i = 0
    do n=1, nnuc
      do j=1, nHCD
        i = i + 1
        write(1,'(a6,i5,X,a4,a,a3,X,a,i4,4X,5f11.4)') 'ATOM  ', i, " CA ", " ", "DNA", "A", n, x(:,i), 0.0, rDNA
      enddo
      i = i + 1
      write(1,'(a6,i5,X,a4,a,a3,X,a,i4,4X,5f11.4)') 'ATOM  ', i, " CB ", " ",&
                                    "HIS", "A", n, x(:,i), 0.0, rHCP
      i = i + 1
      write(1,'(a6,i5,X,a4,a,a3,X,a,i4,4X,5f11.4)') 'ATOM  ', i, " CG ", " ",& 
                                    "HIS", "A", n, x(:,i), 0.0, rHCP
      i = i + 1
      write(1,'(a6,i5,X,a4,a,a3,X,a,i4,4X,5f11.4)') 'ATOM  ', i, " CD ", " ",& 
                                    "HIS", "A", n, x(:,i), 0.0, rHCP
      i = i + 1
      write(1,'(a6,i5,X,a4,a,a3,X,a,i4,4X,5f11.4)') 'ATOM  ', i, " CE ", " ",& 
                                    "HIS", "A", n, x(:,i), 0.0, rHCP
    enddo
    write(1,'(a)') "ENDMDL"
    write(1,*)
    close(1)


    stop

    end

!------------------------------------------------------------------------------

    subroutine cross_product(v1,v2,v3) 

    implicit none
   
    real(8), intent(in)  :: v1(3), v2(3)
    real(8), intent(out) :: v3(3)

    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

    end subroutine

!------------------------------------------------------------------------------

    subroutine build(xb,xa,xt,b,a,t,xi)

    implicit none

    real(8), intent(in)  :: xb(3), xa(3), xt(3)
    real(8), intent(in)  :: b, a, t
    real(8), intent(out) :: xi(3)

    real(8), parameter :: pi = acos(-1.0d0)

    real(8) :: R(3,3), v(3), qx(3), qy(3), qz(3), s

    if ( ( a == 180.0d0 ) .and. ( t == -999.0d0 ) ) then

      qz(:) = xb(:) - xa(:)
      qz(:) = qz(:)/sqrt(dot_product(qz,qz))
      v = b*qz
      xi(:) = v(:) + xb(:)

    else

      ! Vector along Z-axis
      v(1:3) = (/ 0.0d0, 0.0d0, b /)
  
      ! Rotate a along Y-axis 
      s = pi*a/180.0d0
      R(1,1:3) = (/  cos(s), 0.0d0, sin(s) /)
      R(2,1:3) = (/   0.0d0, 1.0d0,  0.0d0 /)
      R(3,1:3) = (/ -sin(s), 0.0d0, cos(s) /)
      v = matmul(R,v)

      ! Rotate t along Z-axis 
      s = pi*t/180.0d0
      R(1,1:3) = (/ cos(s), -sin(s), 0.0d0 /)
      R(2,1:3) = (/ sin(s),  cos(s), 0.0d0 /)
      R(3,1:3) = (/  0.0d0,   0.0d0, 1.0d0 /)
      v = matmul(R,v)

      ! Set internal coordinate
      qz(:) = xb(:) - xa(:)
      s     = sqrt(dot_product(qz,qz))
      if ( s == 0.0d0 ) stop "s1 == 0"
      qz(:) = qz(:)/s
      qz(:) = -qz(:)

      qx(:) = xt(:) - xa(:)
      s     = sqrt(dot_product(qx,qx))
      if ( s == 0.0d0 ) stop "s2 == 0"
      qx(:) = qx(:)/s
      qy(1) = qx(2)*qz(3) - qx(3)*qz(2)
      qy(2) = qx(3)*qz(1) - qx(1)*qz(3)
      qy(3) = qx(1)*qz(2) - qx(2)*qz(1)
      s     = sqrt(dot_product(qy,qy))
      if ( s == 0.0d0 ) stop "s3 == 0"
      qy(:) = qy(:)/s

      qx(1) = qy(2)*qz(3) - qy(3)*qz(2)
      qx(2) = qy(3)*qz(1) - qy(1)*qz(3)
      qx(3) = qy(1)*qz(2) - qy(2)*qz(1)
      s     = sqrt(dot_product(qx,qx))
      if ( s == 0.0d0 ) stop "s4 == 0"
      qx(:) = qx(:)/s
      qx(:) = -qx(:)

      ! Find rotation matrix and rotate v
      R(1,:) = qx(:)
      R(2,:) = qy(:)
      R(3,:) = qz(:)
      call inv(R)
      v = matmul(R,v)

      ! Set xi(:) coordinate
      xi(:) = v(:) + xb(:)

    endif

    return

    end subroutine

!------------------------------------------------------------------------------

    subroutine inv(A)

    implicit none

    real(8), intent(inout)  :: A(3,3)

    real(8) :: detA, B(3,3)

    detA =   A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
           - A(1,1)*A(2,3)*A(3,2) - A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3)

    B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)

    B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)

    B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
    B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
    B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    if ( detA == 0.0d0 ) stop 'detA = 0'

    A = B/detA

    return

    end subroutine

!------------------------------------------------------------------------------

    subroutine get_bond(ri,rj,d)

    implicit none

    real(8), intent(in)  :: ri(3), rj(3)
    real(8), intent(out) :: d

    real(8) :: xi(3)

    xi = ri - rj
    d  = sqrt(dot_product(xi,xi))

    return

    end subroutine

!------------------------------------------------------------------------------

    subroutine get_angle(ri,rj,rk,theta)

    implicit none

    real(8), intent(in)  :: ri(3), rj(3), rk(3)
    real(8), intent(out) :: theta

    real(8), parameter :: PI = acos(-1.0d0)
    real(8) :: xi(3), xj(3), di, dj
   
    xi = ri - rj
    xj = rk - rj 
    di = sqrt(dot_product(xi,xi))
    dj = sqrt(dot_product(xj,xj))
    theta = dot_product(xi,xj)/(di*dj)
    theta = acos(theta)*180.0d0/PI

    return

    end subroutine

!------------------------------------------------------------------------------

    subroutine get_torsion(ri,rj,rk,rl,phi)

    implicit none
  
    real(8), intent(in)  :: ri(3), rj(3), rk(3), rl(3)
    real(8), intent(out) :: phi

    real(8), parameter :: PI = acos(-1.0d0)
    real(8) :: rij(3)
    real(8) :: rkj(3), dkj
    real(8) :: rkl(3)
    real(8) :: n1(3), n1d
    real(8) :: n2(3), n2d
    real(8) :: n3(3)
    real(8) :: n1dn2, n3drkj
    real(8) :: cosphi

    rij = ri - rj
    rkj = rk - rj
    dkj = sqrt(dot_product(rkj,rkj))
    rkl = rk - rl

    call cross_product(rij,rkj,n1)
    n1d = sqrt(dot_product(n1,n1))

    call cross_product(rkj,rkl,n2)
    n2d = sqrt(dot_product(n2,n2))

    call cross_product(n1,n2,n3)

    n1dn2 = dot_product(n1,n2)
    n3drkj = dot_product(n3,rkj)

    if ( n1d*n2d /= 0.0d0 ) then

       cosphi = n1dn2/(n1d*n2d)

       if ( cosphi > 1.0d0 ) then
          phi = 0.0d0
       else if ( cosphi < -1.0d0 ) then
          phi = PI
       else
          phi = acos(cosphi)
       endif

       if ( n3drkj < 0.0d0 ) then
          phi = -phi
       endif

       phi = phi*180.0d0/PI

    else

      phi = -9999.0d0

    endif

    return

    end subroutine


