! FDMNES subroutines
! From B. Ravel
! Calculation of the anomalous atomic f' and f".

      subroutine fprime(numat,ephoton,fpp,fp)

      use declarations
      implicit real(kind=db) (a-h,o-z)
      parameter( nkm = 2)

      real(kind=db), dimension(nkm):: sumfp, sumfpp, xk
      real(kind=db), dimension(103):: atom_weight

      data atom_weight/ 1.00797,   4.0026,   6.9390,   9.0122,  10.8110,
     &         12.0111,  14.0067,  15.9994,  18.9984,  20.1830,
     &         22.9898,  24.3120,  26.9815,  28.0860,  30.9738,
     &         32.0640,  35.4530,  39.9480,  39.1020,  40.0800,
     &         44.9560,  47.9000,  50.9420,  51.9960,  54.9380,
     &         55.8470,  58.9330,  58.7100,  63.5400,  65.3700,
     &         69.7200,  72.5900,  74.9920,  78.9600,  79.9090,
     &         83.8000,  85.4700,  87.6200,  88.9050,  91.2200,
     &         92.9060,  95.9400,  99.0000, 101.0700, 102.9050,
     &        106.4000, 107.8700, 112.4000, 114.8200, 118.6900,
     &        121.7500, 127.6000, 126.9040, 131.3000, 132.9050,
     &        137.3400, 138.9100, 140.1200, 140.9070, 144.2400,
     &        147.0000, 150.3500, 151.9600, 157.2500, 158.9240,
     &        162.5000, 164.9300, 167.2600, 168.9340, 173.0400,
     &        174.9700, 178.4900, 180.9480, 183.8500, 186.2000,
     &        190.2000, 192.2000, 195.0900, 196.9670, 200.5900,
     &        204.3700, 207.1900, 208.9800, 210.0000, 210.0000,
     &        222.0000, 223.0000, 226.0000, 227.0000, 232.0380,
     &        231.0000, 238.0400, 237.0000, 242.0000, 243.0000,
     &  247.0, 247.0, 251.0, 254.0, 253.0, 256.0, 254.0, 257.0/
      data xkev / 12.397639 /

! NW : number of xray values to be interpolated 
! Change Angstroms(KeV) To KeV(Angstroms)
      xk(1) = ephoton * rydb / 1000
      nw = 1
      wt = atom_weight(numat)

      call calc(nw, numat, xk, sumfp, sumfpp, wt)

      fpp = sumfpp(1)
      fp  = sumfp(1)

      end

!***********************************************************************

      subroutine calc(nw, iz, xk, sumfp, sumfpp, wt)

      use declarations
      implicit real(kind=db) (a-h,o-z)
      parameter( nkm = 2, nshm = 24 )

      character(len=132) file_name
      character isym*2, nat*2, natom*2, nshel(24)*8

      real(kind=db), dimension(5):: eg
      real(kind=db), dimension(11):: el, ew, sig, sl
      real(kind=db), dimension(nkm):: cxb, sumfp, sumfp0, sumfpp,  !!ERR
     &                               xjensn, xk
      real(kind=db), dimension(nw,nshm):: fp, fpp

      common/gaus/ cx, bb, sigg(5), rx, icount
      common/edge/ sedge

      data c / 137.0367 /
      data is / 1 /
      data mx / 5 /
      data n_interp / 2 / 
      data au / 2.80022e+7 /
      data c1 / 0.02721 /

c-----------------------------------------------------------------------
c Open "XSECT.DAT" File And Get # Of Orbitals, "NO", "ETERM" And
c Starting Record For Element, "IREC"
c-----------------------------------------------------------------------

      do i = 1,3
        select case(i)
          case(1)
            file_name = 'xsect.dat'
          case(2)
            file_name = 'C:/fdmnes/xsect.dat'
          case(3)
            file_name = '/Share/public/FDMNES/xsect.dat'
        end select    
        open(is, file = file_name, status='old', iostat=istat)
        if( istat == 0 ) exit
      end do 
      if( istat /= 0 ) then
        file_name = 'xsect.dat'
        call write_open_error(file_name,istat,1)
      endif

      if( iz < 1 .or. iz > 103 ) goto 220

      do i = 1,iz-1
        read(1,*)
      end do
      read(1,'(i4,1x,a2,i6,i3,f7.3)') iiz, isym, iirec, no, eterm

      if( iiz /= iz ) goto 220

      irec = iirec
      natom = isym
      nat = natom

      cxb(1:nw) = 0
      fp(1:nw,1:no) = 0._db
      fpp(1:nw,1:no) = 0._db

c-----------------------------------------------------------------------
c NAT     - Atomic Symbol
c NJ      - Orbital Sequence Number
c NSHEL   - Orbital Type. 1S1/2 Etc.
c XW      - Wavelength In Angstroms
c EW      - Energy In KeV
c SIG     - Cross Section In Barns (10^-24cm^2)
c BE      - Binding Energy(KeV)
c IFtype  - 0,1,2 Function Type
c Read MX Energies And X-Sections For Orbital NJ
c-----------------------------------------------------------------------
c Start Orbital Loop  **************************************************
c-----------------------------------------------------------------------

      do 160 j = 1, no

      rewind(is)
      do i = 1,irec-1
        read(is,*)
      end do

      do k = 1, mx
   21   read(unit=is, fmt=280) nat, nj, nshel(j), xw, ew(k), sig(k)
        irec = irec + 1
        if( nat /= natom ) goto 21
        if( nj /= j ) goto 220
      end do

! Read 5 Energies And X-Sections For Orbital J For The GAUSS Integration 
! Points. Binding_Energy BE and Function Type Iftype also.
      do k = 1, 5
        read(unit=is, fmt=280)
     &             nat, nj, nshel(j), xw, eg(k), sigg(k), be, iftype
        irec = irec + 1
        if( nat /= natom .or. nj /= j .or. be <= 0._db ) goto 220
        ew(k + mx) = eg(k)
        sig(k + mx) = sigg(k)
      end do

      nx = mx + 5

! Iftype = 0 So Read X-Section at Energy = 1.001*Binding_Energy
      if( iftype == 0 ) then
        nx = nx + 1
        read(unit=is, fmt=280) nat, nj, nshel(j), xw, ew(nx), sig(nx)
        irec = irec + 1
        if( nat /= natom ) goto 220
        sigedg = sig(nx)
      endif

! Input Completed

      bb = be / c1

      sigg(1:5) = sigg(1:5) / au

! Sort All X-sections
      call sort(nx, 11, ew, sig)

! Sort The Five X-Sections At Integration Points
      call sort(5, 5, eg, sigg)

! Change X-Section And Energy To Logs For Interpolation
      do k = 1, nx
        el(k) = 0._db
        sl(k) = 0._db
        el(k) = log(ew(k))
        if( abs(sig(k)) < eps10 ) cycle
        sl(k) = log(sig(k))
      end do

! Start # Of Desired Energies Loop 

      do 150 k = 1, nw
      mf = 0

! ZX=LOG Of Xray(KeV) Energy
      zx = log(xk(k))
      cx = 0._db

      if( be <= xk(k) ) then

        if (n_interp == 0) then
          call xsect(zx, el, sl, cx, nx)
        else
          do m = 1, nx
            n1 = m
            if( abs(sl(m)) > eps10 ) exit
          end do
          mm = (nx - n1) + 1
          cx = aknint(zx,mm,n_interp,el(n1),sl(n1))
! CX Is The Interpolated X-section In Barns
          cx = exp(cx)
        endif

! CXB Is Sum To Get MU/RHO
! Change CX To Atomic Units
        cxb(k) = cxb(k) + cx
        cx = cx / au

      endif

! RX=Xray Energy In Atomic Units
  
      icount = 6
      rx = xk(k) / c1

      deux_pi2 = 2.0 * ( pi ** 2) 

      if( iftype /= 0 .or. be < xk(k) ) then
        if( iftype == 0 .or. iftype == 1 .or. iftype == 2 )
     &      fp(k,j) = gauss(1,iftype,5) * c / deux_pi2
      else
! SEDGE Is X-section In Atomic Units At Energy=1.001*BE
        sedge = sigedg / au
        cx = 0._db
        fp(k,j) = gauss(1,3,5) * c / deux_pi2
        mf = 3
      endif

      fpp(k,j) = 0._db

      if( abs(cx) > eps10 ) fpp(k,j) = ( c * cx * rx) / (4.0 * pi)

      corr = 0._db

      if( abs(cx) > eps10 ) corr = - (((( cx * rx * 0.5)  
     &       * log( (rx + bb) / (rx - bb))) * c) / deux_pi2 )
      if (mf == 3 ) then
        if( abs(-bb+rx) > 0.0001 ) then
          bbrx = -bb + rx
        else
          bbrx = -0.0001
        endif
        corr = (((((0.5 * sedge) * (bb**2)) * log( bbrx
     &       / ((- bb) - rx))) / rx) * c) / (2. * (pi**2))
      endif
      fp(k,j) = fp(k,j) + corr

  150 continue    ! End # Of Desired Energies (Wavelength) Loop

  160 continue    ! End Orbital Loop

! Calc f', f", ETERM, JENSEN Term, MU/RHO
      do k = 1, nw
        sumfp0(k) = sum( fp(k,1:no) ) 
        sumfpp(k) = sum( fpp(k,1:no) ) 
        xjensn(k) = - ((0.5*iz) * (( (xk(k) / c1) / (137.0367**2))**2))
        sumfp(k) = (sumfp0(k) + eterm) + xjensn(k)
        cxb(k) = (cxb(k) * 0.602472) / wt
      end do

      close(unit=is) 

      return 

  220 write(unit=*, fmt=370) 

      close(unit=is) 

      return 
  280 format(a2,2x,i4,a6,4e15.8,i2)
  370 format('DATA FILE ERROR')
      end

!***********************************************************************

      subroutine xsect(zx, el, sl, cx, nx)

      use declarations
      implicit real(kind=db) (a-h,o-z)
      dimension el(11), sl(11)

! FIND EL(K) CLOSEST TO ZX
      er = 1000000.
      do l = 1, nx
        p = abs(zx - el(l))
        if( p > er ) cycle
        er = p
        ll = l
      end do

      ll = ll - 1

      if( ll == 0 ) ll = 1
      if( ll == 12 ) ll = 11
      if( abs(sl(ll)) < eps10  ) ll = ll + 1

      det = (((el(ll + 2) ** 2) * (el(ll + 1) - el(ll))) + ((el(ll + 1)
     & ** 2) * (el(ll) - el(ll + 2)))) + ((el(ll) ** 2) * (el(ll + 2) - 
     &el(ll + 1)))

      a0 = ((((el(ll) ** 2) * ((sl(ll + 1) * el(ll + 2)) - (sl(ll + 2)
     & * el(ll + 1)))) + ((el(ll + 1) ** 2) * ((sl(ll + 2) * el(ll)) - (
     &sl(ll) * el(ll + 2))))) + ((el(ll + 2) ** 2) * ((sl(ll) * el(ll + 
     &1)) - (sl(ll + 1) * el(ll))))) / det

      a1 = ((((el(ll) ** 2) * (sl(ll + 2) - sl(ll + 1))) + ((el(ll + 1)
     & ** 2) * (sl(ll) - sl(ll + 2)))) + ((el(ll + 2) ** 2) * (sl(ll + 1
     &) - sl(ll)))) / det

      a2 = (((sl(ll) * (el(ll + 2) - el(ll + 1))) + (sl(ll + 1) * (el(ll
     &) - el(ll + 2)))) + (sl(ll + 2) * (el(ll + 1) - el(ll)))) / det

      cx = exp((a0 + (a1 * zx)) + (a2 * (zx**2)))

      return 
      end

!***********************************************************************

      function sigma(iftype, x)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db) sigma 

      common/edge/ sedge
      common/gaus/ cx, bb, ssg(5), rx, icount

      icount = icount - 1

      select case(iftype)
        case(0)
          sigma = (((ssg(icount) * (bb**3)) / (x**2))
     &           / (((rx**2) * (x**2))
     &           - (bb**2))) - (((bb*cx) * (rx**2)) 
     &           / (((rx**2) * (x**2))- (bb**2)))
        case(1)
          sigma = ((0.5 * (bb**3)) * ssg(icount)) / (sqrt(x) * (((rx**2)
     &       * (x**2)) - ((bb**2) * x)))

        case(2)
          denom = ((x**3) * (rx**2)) - ((bb**2) / x)
          sigma = (((((2.0 * bb) * ssg(icount)) * (bb**2)) / (x**4)) 
     &           / denom) - ((((2.0 * bb) * cx) * (rx**2)) / denom)
        case(3)
          sigma = ( (bb**3) * (ssg(icount) - (sedge * (x**2))) )
     &           / ( (x**2) * (((x**2) * (rx**2)) - (bb**2)) )
      end select

      return 
      end

!***********************************************************************

      function gauss(n,iftype,mm)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db) gauss
      character(len=8), dimension(6):: cmt
      integer, dimension(6):: m
      real(kind=db), dimension(6):: a, g, z

      data cmt / 'GAUSS   ','N NOT IN',' (1,6). ',' RESULTS', 
     &           'SET TO Z','RO.     ' /

      if( n > 6 .or. n < 1 ) then
        call labrt(1, cmt, 1)
        gauss = 0._db
        return 
      endif

! YJ Fonctionne car n est toujours egal a 1
      m(1) = mm
! YJ
      nn = n

      goto (110, 90, 70, 50, 30, 10), nn

   10 j = 1
      g(6) = 0._db

   20 call ltbl(m(6), j, a(6), z(6))

   30 k = 1
      g(5) = 0._db

   40 call ltbl(m(5), k, a(5), z(5))

   50 l = 1
      g(4) = 0._db

   60 call ltbl(m(4), l, a(4), z(4))

   70 jj = 1
      g(3) = 0._db

   80 call ltbl(m(3), jj, a(3), z(3))

   90 kk = 1
      g(2) = 0._db

  100 call ltbl(m(2), kk, a(2), z(2))

  110 ll = 1
      g(1) = 0._db

  120 call ltbl(m(1), ll, a(1), z(1))

      g(1) = g(1) + a(1) * sigma(iftype,z(1))
      ll = ll + 1

      if (ll <= m(1)) goto 120
      if (nn == 1) goto 130

      g(2) = g(2) + a(2) * sigma(iftype,z(2)) * g(1)
      kk = kk + 1

      if (kk <= m(2)) goto 100
      if (nn == 2) goto 130

      g(3) = g(3) + a(3) * sigma(iftype,z(3)) * g(2)
      jj = jj + 1

      if (jj <= m(3)) goto 80
      if (nn == 3) goto 130

      g(4) = g(4) + a(4) * sigma(iftype,z(4)) * g(3)
      l = l + 1

      if (l <= m(4)) goto 60
      if (nn == 4) goto 130

      g(5) = g(5) + a(5) * sigma(iftype,z(5)) * g(4)
      k = k + 1

      if (k <= m(5)) goto 40
      if (nn == 5) goto 130

      g(6) = g(6) + a(6) * sigma(iftype,z(6)) * g(5)
      j = j + 1

      if (j <= m(6)) goto 20

  130 gauss = g(nn)

      return 
      end

!***********************************************************************

      subroutine labrt(isw, cmt, inx)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      character(len=8), dimension(6):: cmt
      logical ps, ts

      data np / 10 /
      data ps / .true. /
      data ts / .false. /

      select case(isw)
        case(1)
          if( ps .and. (np > 0) ) write(6,100) cmt(1:5), inx
          np = np - 1
        case(2)
          ps = .false.
        case(3)
          ps = .true.
          np = inx
        case(4)
          ts = .true.
        case(5)
          ts = .false.
      end select

      return 
  100 format('0',9x,5(a8,2x),i9)
      end

!***********************************************************************

      function aknint(xbar,in,im,x,y)

! AITKEN REPEATED INTERPOLATION

!   XBAR = ABSCISSA AT WHICH INTERPOLATION IS DESIRED
!   IABS(IN) = NO. OF VALUES IN TABLE
!              IF IN.GT.0, CHK ORDERING OF X(I).
!              IF IN.LT.0, SKIP PRECEEDING TEST.
!   IM   = DEGREE OF APPROXIMATING POLYNOMIAL
!   X    = VECTOR OF IABS(IN) VALUES OF ABSCISSA
!   Y    = VECTOR OF IABS(IN) VALUES OF ORDINATE
!   T    = TEMPORARY STORAGE VECTOR OF 4*(M+1) LOCATIONS)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db) aknint 

      dimension t(80), x(11), y(11)

!     DATA MES1 /47HAKNINT WARNING ORDER OF INTERPOLATION TOO LARGE/
!     DATA MES3 /35HAKNINT N.LT.2 YBAR RETURNED AS Y(1)/
!     DATA MES4 /34HAKNINT X(I) NOT SEQUENCED PROPERLY/

      dxbar = xbar
      n = iabs(in)
      m = im

      if( m >= n ) goto 120

   10 k = n - 1

      if( n < 2 ) goto 110

      s = x(2) - x(1)

! CHK IF ORDER MONOTONIC
      if( in >= 0 .and. n /= 2 ) then
        do i = 3, n
          z = (x(i) - x(i - 1)) * s
          if( z <= 0._db ) goto 100
        end do
      endif

! INCREASING ORDER
      if( s < 0._db ) goto 50

      do j = 1, n
        if( xbar <= x(j) ) goto 70
      end do

      j = n

      goto 70

! DECREASING ORDER
   50 do j = 1, n
        if( xbar >= x(j) ) goto 70
      end do

      j = n

   70 k = m
      m = m + 1
      j = j - (m / 2)
      j = max0(j,1)
      j = min0(j,n - k)
      mend = j + k

      do i = j, mend
        kk = (i - j) + 1
        t(kk) = y(i)
        t(kk + m) = x(i) - dxbar
      end do

      do i = 1, k
        kk = i + 1
        do jj = kk, m
          t(jj) = ( ( t(i) * t(jj + m) ) - ( t(jj) * t(i + m) ) )  
     &          / ( x((jj + j) -1) - x((i + j) - 1) )
        end do
      end do

      aknint = t(m)

      return 

  100 write(unit=*, fmt=1) 'AKNINT X(I) NOT SEQUENCED PROPERLY'
  110 write(unit=*, fmt=1) 'AKNINT N.LT.2 YBAR RETURNED AS Y(1)'

      aknint = y(1)

      return 

  120 write(unit=*, fmt=1) 
     &'AKNINT WARNING ORDER OF INTERPOLATION TOO LARGE !!'
    1 format(A)

      m = n - 1

      goto 10

      end

!***********************************************************************

      subroutine sort(n, nm, a, b)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db), dimension(nm):: a(nm), b(nm)

      m = n - 1

      do i = 1, m
        i1 = i + 1
        do j = i1, n
          if( a(j) > a(i) ) cycle
          x = a(j)
          y = a(i)
          a(i) = x
          a(j) = y
          x = b(j)
          y = b(i)
          b(i) = x
          b(j) = y
        end do
      end do

      return 
      end

!***********************************************************************

      subroutine ltbl(m, k, aa, z)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      dimension a(68), x(62)

      data x(1) / .06943184420297 /
      data x(2) / .33000947820757 /
      data x(3) / .04691007703067 /
      data x(4) / .23076534494716 /
      data x(5) / .03376524289992 /
      data x(6) / .16939530676687 /
      data x(7) / .38069040695840 /
      data x(8) / .02544604382862 /
      data x(9) / .12923440720030 /
      data x(10) / .29707742431130 /
      data x(11) / .01985507175123 /
      data x(12) / .10166676129319 /
      data x(13) / .23723379504184 /
      data x(14) / .40828267875217 /
      data x(15) / .01591988024619 /
      data x(16) / .08198444633668 /
      data x(17) / .19331428364971 /
      data x(18) / .33787328829809 /
      data x(19) / .01304673574141 /
      data x(20) / .06746831665551 /
      data x(21) / .16029521585049 /
      data x(22) / .28330230293537 /
      data x(23) / .42556283050918 /
      data x(24) / .01088567092697 /
      data x(25) / .05646870011595 /
      data x(26) / .13492399721298 /
      data x(27) / .24045193539659 /
      data x(28) / .36522842202382 /
      data x(29) / .00921968287664 /
      data x(30) / .04794137181476 /
      data x(31) / .11504866290285 /
      data x(32) / .20634102285669 /
      data x(33) / .31608425050091 /
      data x(34) / .43738329574426 /
      data x(35) / .00790847264071 /
      data x(36) / .04120080038851 /
      data x(37) / .09921095463335 /
      data x(38) / .17882533027983 /
      data x(39) / .27575362448178 /
      data x(40) / .38477084202243 /
      data x(41) / .00685809565159 /
      data x(42) / .03578255816821 /
      data x(43) / .08639934246512 /
      data x(44) / .15635354759416 /
      data x(45) / .24237568182092 /
      data x(46) / .34044381553605 /
      data x(47) / .44597252564632 /
      data x(48) / .600374098758e-2 /
      data x(49) / .31363303799647e-1 /
      data x(50) / .75896708294787e-1 /
      data x(51) / .13779113431991 /
      data x(52) / .21451391369574 /
      data x(53) / .30292432646121 /
      data x(54) / .39940295300128 /
      data x(55) / .00529953250417 /
      data x(56) / .02771248846338 /
      data x(57) / .06718439880608 /
      data x(58) / .12229779582250 /
      data x(59) / .19106187779868 /
      data x(60) / .27099161117138 /
      data x(61) / .35919822461038 /
      data x(62) / .45249374508118 /

      data a(1) / .17392742256873 /
      data a(2) / .32607257743127 /
      data a(3) / .11846344252810 /
      data a(4) / .23931433524968 /
      data a(5) / .28444444444444 /
      data a(6) / .85662246189585e-1 /
      data a(7) / .18038078652407 /
      data a(8) / .23395696728635 /
      data a(9) / .06474248308443 /
      data a(10) / .13985269574464 /
      data a(11) / .19091502525256 /
      data a(12) / .20897959183674 /
      data a(13) / .05061426814519 /
      data a(14) / .11119051722669 /
      data a(15) / .15685332293894 /
      data a(16) / .18134189168918 /
      data a(17) / .04063719418079 /
      data a(18) / .09032408034743 /
      data a(19) / .13030534820147 /
      data a(20) / .15617353852000 /
      data a(21) / .16511967750063 /
      data a(23) / .07472567457529 /
      data a(24) / .10954318125799 /
      data a(25) / .13463335965500 /
      data a(26) / .14776211235738 /
      data a(27) / .02783428355809 /
      data a(28) / .06279018473245 /
      data a(29) / .09314510546387 /
      data a(30) / .11659688229599 /
      data a(31) / .13140227225512 /
      data a(32) / .13646254338895 /
      data a(33) / .02358766819326 /
      data a(34) / .05346966299766 /
      data a(35) / .08003916427167 /
      data a(36) / .10158371336153 /
      data a(37) / .11674626826918 /
      data a(38) / .12457352290670 /
      data a(39) / .02024200238266 /
      data a(40) / .04606074991886 /
      data a(41) / .06943675510989 /
      data a(42) / .08907299038097 /
      data a(43) / .10390802376845 /
      data a(44) / .11314159013145 /
      data a(45) / .11627577661544 /
      data a(46) / .01755973016588 /
      data a(47) / .04007904357988 /
      data a(48) / .06075928534395 /
      data a(49) / .07860158357910 /
      data a(50) / .09276919873897 /
      data a(51) / .10259923186065 /
      data a(52) / .10763192673158 /
      data a(53) / .01537662099806 /
      data a(54) / .03518302374405 /
      data a(55) / .05357961023359 /
      data a(56) / .06978533896308 /
      data a(57) / .08313460290850 /
      data a(58) / .09308050000778 /
      data a(59) / .09921574266356 /
      data a(60) / .10128912096278 /
      data a(61) / .01357622970588 /
      data a(62) / .03112676196932 /
      data a(63) / .04757925584125 /
      data a(64) / .06231448562777 /
      data a(65) / .07479799440829 /
      data a(66) / .08457825969750 /
      data a(67) / .09130170752246 /
      data a(68) / .09472530522754 /

      kk = k

      if( m > 16 .or. m < 4 ) kk = 4

      is = 0
      ih = (m + 1) / 2
      z = .5

      if( mod(m,2) == 1 ) is = -1

      ip = kk
      t = 0._db

      if( ip > ih ) then
        ip = m + 1 - ip
        t = -1
      endif

      i4 = m - 4
      ia = (( (i4 * (m + 4)) + is) / 4) + ip
      aa = a(ia)

      if( ip /= ih .or. is >= 0 ) then 
        ia = ia - ( (i4 + is) / 2 )
        z = - t + sign(x(ia),t)
      endif

      return 
      end
