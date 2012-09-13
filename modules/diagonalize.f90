! Modulo autosuficiente para diagonalizar matrices. usa F77, no esta optimizado. 
! 
! puede pedir valores y vectores propios, o solo valores propios. la diferencia
! es como un factor 2.
!
! Aca se puede ver donde quedan almacenados los vectores propios....
!  do j1=1,n
!     print*,norma(matmul(matrix_in,eigenvectors(:,j1))-eigenvalues2(j1)*eigenvectors(:,j1))
!  enddo

! 
!! todo se basa en rutinas de este estilo 
!  http://netlib.cs.utk.edu/eispack/cbabk2.f

module diagslow
  implicit none
  interface matrix_slow_diagonalize
     module procedure matrix_slow_diagonalize_all
     module procedure matrix_slow_diagonalize_part
  end interface
  interface matrix_slow_unitary
     module procedure matrix_slow_unitary_all
     module procedure matrix_slow_unitary_part
  end interface
contains
  subroutine matrix_slow_unitary_part(matrix_in,eigenvalues)
    implicit none
    complex(kind(1d0)), intent(in)            :: matrix_in(:,:)
    real(kind(1d0)), intent(out)              :: eigenvalues(:)
    complex(kind(1d0)), allocatable           :: eigenvectors(:,:)
    integer                                   :: k_size
    k_size=size(matrix_in(:,1))
    allocate(eigenvectors(0:k_size-1,0:k_size-1))
    call matrix_slow_unitary_all(matrix_in,eigenvectors,eigenvalues)
    deallocate(eigenvectors)
  end subroutine matrix_slow_unitary_part
  subroutine matrix_slow_unitary_all(matrix_in,eigenvectors,eigenvalues)
    use math
    implicit none
    complex(kind(1d0)), intent(in)            :: matrix_in(0:,0:)
    complex(kind(1d0)), intent(out)           :: eigenvectors(0:,0:)
    real(kind(1d0)), intent(out)              :: eigenvalues(0:)

    complex(kind(1d0)), allocatable           :: EOb(:,:),h(:,:)
    real(kind(1d0)), allocatable              :: aux(:)
    real(kind(1d0))                           :: randomphase,tiny
    real(kind(1d0)), pointer                  :: evalues(:)
    integer                                   :: k_size,j1

    tiny=10d0**(-8)
    k_size=size(matrix_in(:,0))
    if ( (k_size.ne.size(matrix_in(0,:))).or.(k_size.ne.size(eigenvectors(0,:)))&
         .or.(k_size.ne.size(eigenvectors(:,0))).or.(k_size.ne.size(eigenvalues)) ) then
       print*,"Error en la rutina matrix_slow_unitary, llamada con tamanos incorrectos"
       print*,size(matrix_in(0,:)),size(eigenvectors(0,:)),size(eigenvectors(:,0)),size(eigenvalues)
    end if
    randomphase=2*pi*0.742265331869212d0 
!    randomphase=0d0
    allocate(EOb(0:k_size-1,0:k_size-1))
    allocate(h(0:k_size-1,0:k_size-1),evalues(0:k_size-1),aux(0:k_size-1))
    EOb=matrix_in*exp(I*randomphase)     ! Esto lo hago por si existen eigenfases +- phi, evadir el error
    h=transpose(conjg(EOb))+EOb
    call matrix_slow_diagonalize(h,eigenvectors,evalues)
    do j1=0,k_size-1
       aux=arg(MatMul(EOb,eigenvectors(:,j1))/eigenvectors(:,j1))
       if (standard_deviation_real(aux)**2>tiny) then
          print*,"Algo falla en el procedimiento  matrix_slow_unitary"
          print*,"Se puede probar cambiar randomphase, o aumentar "
          print*,"la tolerancia, aumentando el parametro tiny"
          stop
       end if
       eigenvalues(j1)=mean(aux)-randomphase+pi
    end do
    eigenvalues=modulo(eigenvalues,2*pi)-pi
    deallocate(h,evalues,EOb,aux)
  end subroutine matrix_slow_unitary_all
  subroutine matrix_slow_diagonalize_part(matrix_in,eigenvalues)
    implicit none
    complex(kind(1d0)), intent(in)   :: matrix_in(:,:)
    real(kind(1d0)), intent(out)     :: eigenvalues(:)
    integer                          :: n,nm,matz,ierr
    real (kind(1d0)), allocatable    :: ar(:,:),ai(:,:),w(:),zr(:,:),zi(:,:)
    real (kind(1d0)), allocatable    :: fv1(:),fv2(:),fm1(:,:)
    n=size(matrix_in(1,:))
    if ((size(matrix_in(:,1)).ne.n).or.(size(eigenvalues).ne.n)) then
      print*,"Error en el tamano de los argumentos de matrix_slow_diagonalize"
      print*,n,size(matrix_in(:,1)),size(eigenvalues)
      stop
    end if
    nm=n
    matz=0
    allocate(ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n))
    allocate(fv1(n),fv2(n),fm1(2,n))
    ar=real(matrix_in,kind(1d0))
    ai=aimag(matrix_in)
    call ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
    if (ierr.ne.0) then
       print*,"Algo malo en matrix_slow_diagonalize el codigo error es ",ierr
       stop
    end if
    eigenvalues=w
    deallocate(ar,ai,w,zr,zi)
    deallocate(fv1,fv2,fm1)
  end subroutine matrix_slow_diagonalize_part
  subroutine matrix_hermitian_real(matrix_in,eigenvectors,eigenvalues)
    real(kind(1d0)), intent(in)      :: matrix_in(:,:)
    complex(kind(1d0)), allocatable  :: matrix_in_complex(:,:)
    complex(kind(1d0)), intent(out)  :: eigenvectors(:,:)
    real(kind(1d0)), intent(out)     :: eigenvalues(:)
    allocate(matrix_in_complex(size(matrix_in(1,:)),size(matrix_in(1,:))))
    matrix_in_complex=matrix_in
    call matrix_slow_diagonalize_all(matrix_in_complex,eigenvectors,eigenvalues)
  end subroutine matrix_hermitian_real
  subroutine matrix_slow_diagonalize_all(matrix_in,eigenvectors,eigenvalues)
    implicit none
    complex(kind(1d0)), intent(in)   :: matrix_in(:,:)
    complex(kind(1d0)), intent(out)  :: eigenvectors(:,:)
    real(kind(1d0)), intent(out)     :: eigenvalues(:)
    integer                          :: n,nm,matz,ierr
    real (kind(1d0)), allocatable    :: ar(:,:),ai(:,:),w(:),zr(:,:),zi(:,:)
    real (kind(1d0)), allocatable    :: fv1(:),fv2(:),fm1(:,:)
    n=size(matrix_in(1,:))
    if ((size(matrix_in(:,1)).ne.n).or.(size(eigenvectors(1,:)).ne.n).or.&
         (size(eigenvectors(:,1)).ne.n).or.(size(eigenvalues).ne.n)) then
      print*,"Error en el tamano de los argumentos de matrix_slow_diagonalize"
      print*,n,size(matrix_in(:,1)),size(eigenvectors(1,:)),&
           size(eigenvectors(:,1)),size(eigenvalues)
      stop
    end if
    nm=n
    matz=1
    allocate(ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n))
    allocate(fv1(n),fv2(n),fm1(2,n))
    ar=real(matrix_in,kind(1d0))
    ai=aimag(matrix_in)
    call ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
    if (ierr.ne.0) then
       print*,"Algo malo en matrix_slow_diagonalize el codigo error es ",ierr
       stop
    end if
    eigenvalues=w
    eigenvectors=zr+(0d0,1d0)*zi
    deallocate(ar,ai,w,zr,zi)
    deallocate(fv1,fv2,fm1)
  end subroutine matrix_slow_diagonalize_all
  subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
    integer         :: i,j,n,nm,ierr,matz
    real(kind(1d0)) :: ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),fv1(n),fv2(n),fm1(2,n)
    
    if (n .le. nm) go to 10
    ierr = 10 * n
    go to 50
10  call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
    if (matz .ne. 0) go to 20
    call  tqlrat(n,w,fv2,ierr)
    go to 50
20  do 40 i = 1, n
       do 30 j = 1, n
          zr(j,i) = 0.0d0
30     continue
       zr(i,i) = 1.0d0
40  continue
    call  tql2(nm,n,w,fv1,zr,ierr)
    if (ierr .ne. 0) go to 50
    call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
50  return
   end subroutine ch
   function epslon (x)
     double precision epslon, x
     !     estimate unit roundoff in quantities of size x.
     double precision a,b,c,eps
     a = 4.0d0/3.0d0
10   b = a - 1.0d0
     c = b + b + b
     eps = dabs(c-1.0d0)
     if (eps .eq. 0.0d0) go to 10
     epslon = eps*dabs(x)
     return
   end function epslon
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si
      if (m .eq. 0) go to 200
      do 50 k = 1, n
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
      if (n .eq. 1) go to 200
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
            s = (s / h) / h
            si = (si / h) / h
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
  130    continue
  140 continue
  200 return
      end subroutine htribk
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
      integer i,j,k,l,n,ii,nm,jp1
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
      do 100 i = 1, n
  100 d(i) = ar(i,i)
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
         do 120 k = 1, l
  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0d0
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
            jp1 = j + 1
            if (l .lt. jp1) go to 220
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
         hh = f / (h + h)
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) &
                                + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) &
                                - fi * e(k) - gi * ar(i,k)
  260    continue
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
      return
      end subroutine htridi
      subroutine tqlrat(n,d,e2,ierr)
      integer i,j,l,m,n,ii,l1,mml,ierr
      double precision d(n),e2(n)
      double precision b,c,f,g,h,p,r,s,t
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
  100 e2(i-1) = e2(i)
      f = 0.0d0
      t = 0.0d0
      e2(n) = 0.0d0
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dsqrt(e2(l))
         if (t .gt. h) go to 105
         t = h
         b = epslon(t)
         c = b * b
!     .......... look for small squared sub-diagonal element ..........
  105    do 110 m = l, n
            if (e2(m) .le. c) go to 120
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         s = dsqrt(e2(l))
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * s)
         r = pythag(p,1.0d0)
         d(l) = s / (p + dsign(r,p))
         h = g - d(l)
         do 140 i = l1, n
  140    d(i) = d(i) - h
         f = f + h
!     .......... rational ql transformation ..........
         g = d(m)
         if (g .eq. 0.0d0) g = b
         h = g
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            p = g * h
            r = p + e2(i)
            e2(i+1) = s * r
            s = e2(i) / r
            d(i+1) = h + s * (h + d(i))
            g = d(i) - e2(i) / g
            if (g .eq. 0.0d0) g = b
            h = g * p / r
  200    continue
         e2(l) = s * g
         d(l) = h
!     .......... guard against underflow in convergence test ..........
         if (h .eq. 0.0d0) go to 210
         if (dabs(e2(l)) .le. dabs(c/h)) go to 210
         e2(l) = h * e2(l)
         if (e2(l) .ne. 0.0d0) go to 130
  210    p = d(l) + f
!     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
  250    i = 1
  270    d(i) = p
  290 continue
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end  subroutine tqlrat

      function pythag(a,b)
      double precision pythag, a,b
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end function pythag

      subroutine tql2(nm,n,d,e,z,ierr)
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
!      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
  100 e(i-1) = e(i)
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         do 140 i = l2, n
  140    d(i) = d(i) - h
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
  200    continue
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
  300 continue
      go to 1001
 1000 ierr = l
 1001 return
      end subroutine tql2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
      do 110 i = low, igh
         s = scale(i)
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
  110 continue
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
  140 continue
  200 return
      end subroutine cbabk2
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
      radix = 16.0d0
      b2 = radix * radix
      k = 1
      l = n
      go to 100
   20 scale(m) = j
      if (j .eq. m) go to 50
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
   50 go to (80,130), iexc
   80 if (l .eq. 1) go to 280
      l = l - 1
  100 do 120 jj = 1, l
         j = l + 1 - jj
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
         m = l
         iexc = 1
         go to 20
  120 continue
      go to 140
  130 k = k + 1
  140 do 170 j = k, l
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
         m = k
         iexc = 2
         go to 20
  170 continue
      do 180 i = k, l
  180 scale(i) = 1.0d0
  190 noconv = .false.
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
  270 continue
      if (noconv) go to 190
  280 low = k
      igh = l
      return
      end subroutine cbal
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
      double precision s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end subroutine cdiv
      subroutine general_complex_eigenvalues(rho,eigenvalues)
        complex(kind(1d0)),intent(in)  :: rho(:,:)
        complex(kind(1d0)),intent(out) :: eigenvalues(:)
        integer                        :: dim,ierr
        real(kind(1d0)),allocatable    :: eigReal(:),eigAimag(:),zr(:,:),zi(:,:),fv1(:),fv2(:),fv3(:),rhoR(:,:),rhoI(:,:)
        dim=size(eigenvalues)
        allocate(eigReal(dim),eigAimag(dim),zr(dim,dim),zi(dim,dim),fv1(dim),fv2(dim),fv3(dim))
        allocate(rhoR(dim,dim),rhoI(dim,dim));        rhoR=real(rho); rhoI=aimag(rho)
        if ((size(rho(1,:)).ne.dim).or.(size(rho(:,1)).ne.dim)) then
           print*,"problema en general_complex_eigenvalues",size(rho(1,:)),dim,size(rho(:,1))
           stop
        end if
!      outine cg(nm, n,   ar,       ai,         wr,      wi,matz,zr,zi,fv1,fv2,fv3,ierr)
        call cg(dim,dim,rhoR,rhoI,eigReal,eigAimag,0,zr,zi,fv1,fv2,fv3,ierr)
        eigenvalues=eigReal+(0d0,1d0)*eigAimag
        deallocate(eigReal,eigAimag,zr,zi,fv1,fv2,fv3)
      end subroutine general_complex_eigenvalues
!      subroutine general_complex_eigenvalues_FOR_CONCURRENCE(rho,eigenvalues)
!    oooooooojo, toca mirar las restricciones de esta joda que diagonaliza
      subroutine eigenvalues_FOR_CONCURRENCE(rho,eigenvalues)
        complex(kind(1d0)),intent(in)  :: rho(:,:)
        real(kind(1d0)),intent(out)    :: eigenvalues(:)
        integer                        :: dim,ierr
        real(kind(1d0)),allocatable    :: eigReal(:),eigAimag(:),zr(:,:),zi(:,:),fv1(:),fv2(:),fv3(:),rhoR(:,:),rhoI(:,:)
        dim=size(eigenvalues)
        if ((size(rho(1,:)).ne.dim).or.(size(rho(:,1)).ne.dim)) then
           print*,"problema en general_complex_eigenvalues",size(rho(1,:)),dim,size(rho(:,1))
           stop
        end if
        allocate(eigReal(dim),eigAimag(dim),zr(dim,dim),zi(dim,dim),fv1(dim),fv2(dim),fv3(dim))
        allocate(rhoR(dim,dim),rhoI(dim,dim));        rhoR=real(rho); rhoI=aimag(rho)
        call cg(dim,dim,rhoR,rhoI,eigReal,eigAimag,0,zr,zi,fv1,fv2,fv3,ierr)
        eigenvalues=eigReal
        deallocate(eigReal,eigAimag,zr,zi,fv1,fv2,fv3,rhoR,rhoI)
      end subroutine eigenvalues_FOR_CONCURRENCE
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for comqr
!           and comqr2.  the normal completion code is zero.
!
!        fv1, fv2, and  fv3  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
       integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),fv1(n),fv2(n),fv3(n)
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end subroutine cg
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2
!      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,pythag
      ierr = 0
      if (low .eq. igh) go to 180
      l = low + 1
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
  170 continue
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1)) + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
      lp1 = l + 1
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
  500 continue
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
  600 continue
      if (si .eq. 0.0d0) go to 240
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
      go to 240
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
 1000 ierr = en
 1001 return
      end subroutine comqr
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2
!      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,pythag
      ierr = 0
      do 101 j = 1, n
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
      iend = igh - low - 1
      if (iend) 180, 150, 105
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
            sr = sr / norm
            si = si / norm
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
  130    continue
  140 continue
  150 l = low + 1
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
  170 continue
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))+ dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
      lp1 = l + 1
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
  500 continue
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
  600 continue
      if (si .eq. 0.0d0) go to 240
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
      go to 240
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
  680 norm = 0.0d0
      do 720 i = 1, n
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
  780    continue
  800 continue
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
  840 continue
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
      go to 1001
 1000 ierr = en
 1001 return
      end subroutine comqr2
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale
!      double precision f,g,h,fi,fr,scale,pythag
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
  103    ortr(m) = g
         ar(m,m-1) = scale
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
            fr = fr / h
            fi = fi / h
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
  130    continue
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
            fr = fr / h
            fi = fi / h
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
  160    continue
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
  200 return
      end subroutine corth
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi
      double precision s,tr,ti!,pythag
      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end subroutine csroot
!       function pythag1(a,b)
!       double precision pythag1, a,b
! !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!       double precision p,r,s,t,u
!       p = dmax1(dabs(a),dabs(b))
!       if (p .eq. 0.0d0) go to 20
!       r = (dmin1(dabs(a),dabs(b))/p)**2
!    10 continue
!          t = 4.0d0 + r
!          if (t .eq. 4.0d0) go to 20
!          s = r/t
!          u = 1.0d0 + 2.0d0*s
!          p = u*p
!          r = (s/u)**2 * r
!       go to 10
!    20 pythag1 = p
!       return
!       end function pythag1
end module diagslow
