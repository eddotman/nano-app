! FDMNES subroutines

! Post treatment of the spectra calculated by FDMNES
! 1) Sommation with energy shift of varius spectra
! 2) Integration over the energy, that is convolution by
! a lorentzienne for XANES
!              L(x) = (1/(pi*b)) * 1 / ( 1 + ( (x-a)/b )**2 )

! or something similar for Dafs.

      subroutine convolution(bav_open,Bormann,Conv_done,
     &      convolution_out,Delta_edge,E_cut_imp,E_Fermi_man,Ecent,
     &      Elarg,Estart,
     &      Fit_cal,Gamma_hole,Gamma_hole_imp,Gamma_max,Green_plus,ical,
     &      icheck,indice_par,iscratchconv,itape1,kw_conv,length_line,
     &      Length_word,ngamh,ngroup_par,nkw_conv,nomfich,nomfichbav,
     &      npar,nparm,param,Scan_a,typepar,ncal)

      use declarations
      implicit none

      integer:: i, ical, icheck, ie, ie1, ie2, ifich, ifichref, igr, ii,
     &  initl, initlref, ip, ipar, 
     &  ipl, ipr, ipr1, ipr2, is, iscr, iscratchconv, istop,  
     &  istat, itape1, j, j0, je, jfich, jfichref, jpl, js, jseuil, 
     &  k, l, Length_line, Length_word, long,
     &  longf, mfich, n, n_col, n_selec_core, n1, n2, n3, n4, natomsym, 
     &  ncal, ne2, nef, nelor, nemax, nen2, nenerg, nenerge, nes, nfich, 
     &  nfich_tot, ngamh, ngroup_par, ninit, ninit1, ninitlm, nkw_conv,
     &  nnombre, nparm, nphim, npldafs, npldafs_b, npldafs_th, 
     &  npldafs1, ns, nseuil, nt, numat, nxan, nxan1, nw

      integer, parameter:: nassm = 103   
! njp : Points au dela de la gamme pour diminuer les effets de bord de
! la convolution  
      integer, parameter:: njp = 500   

      character(len=1):: rep
      character(len=9):: keyword, mot9
      character(len=9), dimension(nkw_conv) :: kw_conv
      character(len=9), dimension(ngroup_par,nparm) :: typepar
      character(len=50):: com_date, com_time, Revision 
      character(len=Length_word):: nomab
      character(len=132):: chemin, convolution_out, fdmnes_error,  
     &         fichscanout, identmot, mot, mots, nomfich, nomfichbav
      character(len=Length_word), dimension(:), allocatable:: nom_col
      character(len=132), dimension(:), allocatable:: fichin, fichscanin

      complex(kind=db):: cf
      complex(kind=db), dimension(1):: cdum
      complex(kind=db), dimension(:), allocatable :: dampl, dph, dpht,
     &                                 f0, f0_th
      complex(kind=db), dimension(:,:), allocatable :: f0scan, phdtscan
      complex(kind=db), dimension(:,:,:), allocatable :: Ad, Adafs, As    
      complex(kind=db), dimension(:,:,:,:), allocatable :: mu, mus    

      integer, dimension(0) :: idum
      integer, dimension(3) :: hkl_S
      integer, dimension(10) :: num_core
      integer, dimension(ngroup_par) :: npar
      integer, dimension(ngroup_par,nparm) :: indice_par
      integer, dimension(:), allocatable:: i_done, indf, ne, ninitl,
     &                                     nphi
      integer, dimension(:,:), allocatable:: hkl_dafs, nsup, ne_initl

      logical:: Another_one, Arc, bav_open, Bormann, Dafs, Dafs_bio,
     &  Check_conv, chem, Circular, Conv_done, Cor_abs, decferm,
     &  Deuxieme, Double_cor, E_Fermi_man, Energphot, Epsii_ref_man, 
     &  Extrap, Fermip, Fit_cal, Forbidden, fprim, fprime_atom, 
     &  Full_self_abs, Gamma, Gamma_hole_imp, Gamma_var,
     &  Gaussian_default, Green_int, Green_plus, Magn, new_format,
     &  no_extrap,
     &  nxan_lib, Photo_emission, Scan_a, scan_true, Seah, Self_abs,
     &  Shift_auto, Signal_Sph, Tenseur, Tenseur_car, Thomson
      logical, dimension(:), allocatable:: fichdone, run_done, 
     &                                     Skip_run

      real(kind=db):: a, a1, a2, a3, a4, alambda, Asea, b, b1, b2, b3,
     &  b4, bba, bbb, c, Cal_Volume_maille, conv_mbarn_nelec,
     &  ct_epsilon, ct_nelec,
     &  d, d_dead, de_obj, de1, de2, Delta_edge, Deltar, Densite_atom,
     &  E, E_cut_imp, E_obj, e1m, Ecent, Efermi, EFermi_orig, Eintmax,
     &  Emax, Emin, Elarg, Eph, Epsii_moy, Epsii_ref, Esmin, 
     &  Estart, Estart_l, fac, fpp_avantseuil, fpp0,  
     &  Gamm, Gamma_h, Gamma_max, Im_pi, Im_sig, Ip_pi, Ip_sig, p1, p2,
     &  pasdeb, Pdt, S0_2, V0muf, Vibration, Volume_maille

      real(kind=db), dimension(0):: rdum
      real(kind=db), dimension(3):: angxyz, axyz
      real(kind=db), dimension(10):: Gamma_hole 
      real(kind=db), dimension(ngroup_par,nparm) :: param 
      real(kind=db), dimension(:), allocatable:: angle, bb, 
     &      betalor, decal, e1, e2, Efermip, Elor, En_fermi, Energ,
     &      Eph1, Ephoton, Es, Eseuil, fi, fr, lori, lorix, lorr,   
     &      Pds, p1f, p2f, Tens, Yr, Yi
      real(kind=db), dimension(:,:), allocatable:: decal_initl, Ef,
     &      Epsii, mua_r, mua_i, Xa, Xanes, Xs
      real(kind=db), dimension(:,:,:), allocatable:: Icirc, Icirccor,
     &      Icor, Icircdcor, Idcor 

      common/com_out/ com_date, com_time, fdmnes_error, Revision

      Another_one = .false.

 1000 continue  ! point de retour pour convolution des fichiers TDDFT

      Arc = .true.
      Asea = 0.2_db  ! pente de gamma a l'origine
      chem = .false.
      Check_conv = .false.
      Circular = .false.
      convolution_out = ' '
      d_dead = 0._db
      Dafs = .false.
      decferm = .false.
      Deltar = 0._db
      Densite_atom = 0._db
      Deuxieme = .false.
      Double_cor = .false.
      Efermi = E_cut_imp * rydb
      Epsii_ref_man = .false.
      eintmax = 1000000._db
      Fermip = .false.
      fichscanout = ' '
      Forbidden = .false.
      fpp_avantseuil = 0._db
      fprim = .false.
      fprime_atom = .false.
      Full_self_abs = .false.
      Gamma_var = .false.
      Gaussian_default = .false.
      Green_int = .false.
      hkl_S(:) = 0
      jseuil = 1
      Magn = .false.
      nelor = 0
      nfich = 0
      npldafs = 0
      no_extrap = .false.
      nseuil = -1
      do i = 1,10
        num_core(i) = i
      end do 
      nxan_lib = .false.
      Photo_emission = .false.
      S0_2 = 1._db
      scan_true = .false.
      seah = .false.
      self_abs = .false.
      Shift_auto = .true.
      Signal_Sph = .false.
      tenseur = .false.
      tenseur_car = .false.
      Thomson = .false.
      Vibration = 0._db

      if( bav_open .or. icheck > 1 ) Check_conv = .true.

! -- Lecture --------------------------------------------

      Rewind(itape1)

      boucle_ii: do ii = 1,1000

        n = nnombre(itape1,132)
        read(itape1,'(A)',end=1010) mot

        keyword = identmot(mot,9)

        if( keyword == 'calculati' ) then
          nfich = 0
          do i = 1,1000
            n = nnombre(itape1,132)
            read(itape1,'(A)',end=1010) mots
            mot9 = identmot(mots,9)
            do j = 1,nkw_conv
              if( mot9 == kw_conv(j) ) exit boucle_ii
            end do
            if( mot9 == 'run_done' ) exit boucle_ii
            if( n == 0 ) nfich = nfich + 1
          end do
        endif

      end do boucle_ii

 1010 rewind(itape1)

      allocate( run_done(nfich) )
      allocate( skip_run(nfich) )
      allocate( i_done(nfich) )
      run_done(:) = .false.
      Skip_run(:) = .false.
      do ifich = 1,nfich
        i_done(ifich) = ifich
      end do
      mfich = nfich

      do ii = 1,1000
        n = nnombre(itape1,132)
        read(itape1,'(A)',err=1015,end=1015) mot
        keyword = identmot(mot,9)
        k = 0
        if( keyword == 'run_done' ) then
          do i = 1,mfich
            read(itape1,*,err=1015,end=1015) j, is
            if( j == 0 ) then
              run_done(i) = .true.
            else
              k = k + 1
              if( is == 0 ) Skip_run(k) = .true.
            endif
            if( run_done(i) ) then
              nfich = nfich - 1
              do j = i+1,mfich
                i_done(j) = i_done(j) - 1
              end do
            endif
          end do
          exit
        endif
      end do

 1015 rewind(itape1)

      if( nfich == 0 ) then
        nfich = 1
        mfich = 1
        allocate( fichin(nfich) )
        mot = ' '
        mot = nomfich
        l = len_trim(mot)
        mot(l+1:l+4) = '.txt'
        fichin(1) = mot
      else
        allocate( fichin(nfich) )
      endif

      allocate( decal(nfich) )
      allocate( Efermip(nfich) )
      allocate( En_fermi(nfich) )
      allocate( Eseuil(nfich) )
      allocate( fichscanin(nfich) )
      allocate( ne(nfich) )
      allocate( ninitl(nfich) ) 
      allocate( Pds(nfich) ) 

      Pds(1) = 1._db
      decal(:) = 0._db
      Efermip(1) = 0._db
      ninitl(:) = 1
      ninitlm = 1
      n_selec_core = 0

      if( ncal > 1 .and. ical == 1 ) then
        do jfich = 1,mfich
          iscr = 100 + jfich
          open( iscr, status = 'scratch' )
        end do
      endif

      long = 0
      do ii = 1,1000

        n = nnombre(itape1,132)
        read(itape1,'(A)',end=1020) mot

        keyword = identmot(mot,9)

        select case(keyword)

          case('run_done')
            exit

          case('calculati')

            ifich = 0
            do i = 1,mfich
              n = nnombre(itape1,132)
              if( run_done(i) ) then
                read(itape1,*)                
                n = nnombre(itape1,132)
                if( n /= 0 ) read(itape1,*)
                cycle
              endif
              ifich = ifich + 1
              read(itape1,'(A)') fichin(ifich)
              fichin(ifich) = adjustl( fichin(ifich) )
              n = nnombre(itape1,132)
              if( n == 0 ) then
                Pds(ifich) = 1._db
                decal(ifich) = 0._db
              elseif( n == 1 ) then
                read(itape1,*,err=9999) Pds(ifich)
                decal(ifich) = 0._db
              elseif( n == 2 ) then
                read(itape1,*,err=9999) Pds(ifich), decal(ifich)
                Shift_auto = .false.
              else
                Fermip = .true.
                read(itape1,*,err=9999) Pds(ifich), decal(ifich),
     &                                  Efermip(ifich)
                Shift_auto = .false.
              endif
            end do

          case('cal_tddft')

            ifich = 0
            do i = 1,mfich
              n = nnombre(itape1,132)
              if( run_done(i) ) then
                read(itape1,*)                
                n = nnombre(itape1,132)
                if( n /= 0 ) read(itape1,*)
                cycle
              endif
              ifich = ifich + 1
              if( Another_one ) then
                read(itape1,'(A)') fichin(ifich)
                fichin(ifich) = adjustl( fichin(ifich) )
                Deuxieme = .true.
              else
                read(itape1,'(A)') mot
              endif
              Pds(ifich) = 1._db
              decal(ifich) = 0._db
            end do
            Another_one = .not. Another_one

          case('conv_out')

            n = nnombre(itape1,132)
            read(itape1,'(A)') convolution_out
            convolution_out = adjustl( convolution_out )

          case('gaussian')

            n = nnombre(itape1,132)
            if( n == 1 ) then
              read(itape1,*,err=9999) deltar
              vibration = 0._db
            elseif( n == 2 ) then
              read(itape1,*,err=9999) deltar, vibration
            else
              Gaussian_default = .true.
            endif

          case('fprime')

            fprim = .true.

          case('forbidden')

            Forbidden = .true.

          case('scan_file')

            scan_true = .true.
            do ifich = 1,nfich
              n = nnombre(itape1,132)
              read(itape1,'(A)') fichscanin(ifich)
              fichscanin(ifich) = adjustl( fichscanin(ifich) )
            end do

          case('scan')

            Scan_a = .true.

          case('scan_conv')

            n = nnombre(itape1,132)
            read(itape1,'(A)') fichscanout
            fichscanout = adjustl( fichscanout )

          case('directory')

            chem = .true.
            n = nnombre(itape1,132)
            read(itape1,'(A)') chemin
            chemin = chemin

          case('seah')

            Arc = .false.
            seah = .true.
            n = nnombre(itape1,132)
            select case(n)
              case(0)
                continue
              case(1)
                read(itape1,*,err=9999) asea
              case(2)
                read(itape1,*,err=9999) asea, Gamma_max
                Gamma_max = Gamma_max / rydb
              case(3)
                read(itape1,*,err=9999) asea, Gamma_max, Gamma_hole(1)
                Gamma_hole(1) = Gamma_hole(1) / rydb
                Gamma_max = Gamma_max / rydb
                Gamma_hole_imp = .true.
              case default
                read(itape1,*,err=9999) asea, Gamma_max, Gamma_hole(1),
     &                                  Efermi
                Gamma_hole(1) = Gamma_hole(1) / rydb
                Gamma_max = Gamma_max / rydb
                Gamma_hole_imp = .true.
                E_Fermi_man = .true.
            end select

          case('convoluti')

            Arc = .true.
            n = nnombre(itape1,132)
            select case(n)
              case(0)
                continue
              case(1)
                read(itape1,*,err=9999) Ecent
                Ecent = Ecent / rydb
                Elarg = Elarg / rydb
              case(2)
                read(itape1,*,err=9999) Ecent, Elarg
                Ecent = Ecent / rydb
                Elarg = Elarg / rydb
              case(3)
                read(itape1,*,err=9999) Ecent, Elarg, Gamma_max
                Ecent = Ecent / rydb
                Elarg = Elarg / rydb
                Gamma_max = Gamma_max / rydb
              case(4)
                read(itape1,*,err=9999) Ecent, Elarg, Gamma_max,
     &                                  Gamma_hole(1)
                Ecent = Ecent / rydb
                Elarg = Elarg / rydb
                Gamma_hole(1) = Gamma_hole(1) / rydb
                Gamma_max = Gamma_max / rydb
                Gamma_hole_imp = .true.
              case default
                read(itape1,*,err=9999) Ecent, Elarg, Gamma_max,
     &                                  Gamma_hole(1), Efermi
                Ecent = Ecent / rydb
                Elarg = Elarg / rydb
                Gamma_hole(1) = Gamma_hole(1) / rydb
                Gamma_max = Gamma_max / rydb
                Gamma_hole_imp = .true.
                E_Fermi_man = .true.
            end select

          case('table')

            Arc = .false.
            n = nnombre(itape1,132)
            if( n == 1 ) read(itape1,*)
            nelor = 0
            do ie = 1,10000
              n = nnombre(itape1,132)
              if( n == 0 ) exit
              read(itape1,*)
              nelor = nelor + 1 
            end do
            rewind(itape1)
            do i = 1,10000
              read(itape1,'(A)') mots
              mot9 = identmot(mots,9) 
              if( mot9 == 'table' ) exit
            end do 
            n = nnombre(itape1,132)
            if( n == 1 ) then
              read(itape1,*,err=9999) Efermi
              E_Fermi_man = .true.
            endif
            allocate( Elor(nelor) )
            allocate( betalor(nelor) )
            do ie = 1,nelor
              n = nnombre(itape1,132)
              read(itape1,*,err=9999) Elor(ie), betalor(ie) 
            end do
            Elor(:) = Elor(:) / rydb
            betalor(:) = betalor(:) / rydb

          case('eintmax')

            n = nnombre(itape1,132)
            read(itape1,*,err=9999) eintmax

          case('gamma_fix')
            gamma_var = .false.

          case('gamma_var')
            gamma_var = .true.

! Pour faire le decallage avant la convolution.
          case('dec')

            decferm = .true. 

          case('no_extrap')

            no_extrap = .true.

          case('nxan_lib')

            nxan_lib = .true. 

          case('thomson')

            Thomson = .true. 
            npldafs_th = nnombre(itape1,Length_line) / 2
            allocate( fr(npldafs_th) )
            allocate( fi(npldafs_th) )
            allocate( f0_th(npldafs_th) )
            read(itape1,*,err=9999) (fr(ipl),fi(ipl),ipl = 1,npldafs_th)
            f0_th(:) = cmplx( fr(:), fi(:),db) 
            deallocate( fr )
            deallocate( fi )

          case('s0_2')

            read(itape1,*,err=9999) S0_2
 
          case('photo_emi')

            Photo_emission = .true. 

          case('epsii')
            n = nnombre(itape1,132)
            read(itape1,*,err=9999) Epsii_ref
            Epsii_ref_man = .true.

          case('selec_cor')
            n_selec_core = nnombre(itape1,132)
            read(itape1,*,err=9999) num_core(1:n_selec_core)

          case('surface')
            n = nnombre(itape1,132)
            read(itape1,*,err=9999) hkl_S(:)

          case('circular')
            Circular = .true.

          case('dead_laye')
            n = nnombre(itape1,132)
            read(itape1,*,err=9999) d_dead
            d_dead = d_dead / 10000

          case('double_co')
            Double_cor = .true.

          case default

            call write_error
            do ipr = 6,9,3
              write(ipr,110) mot
            end do
            stop

        end select

      end do
 1020 continue

      if( convolution_out == ' ' ) then
        if( nomfich == 'fdmnes_out' ) then
          mot = fichin(1)
          l = len_trim(mot) - 3
          if( nfich > 1 ) then
            if( mot(l-2:l-2) == '_' ) then
              l = l - 2
            elseif( mot(l-3:l-3) == '_' ) then
              l = l - 3
            endif
          endif
        else
          mot = nomfich
          l = len_trim(mot) + 1
        endif
        if( Deuxieme ) then
          mot(l:l+14) = '_tddft_conv.txt'
        else
          mot(l:l+8) = '_conv.txt'
        endif
        convolution_out = mot
      endif

      if( .not. bav_open .and. Check_conv ) then
        l = len_trim(convolution_out)
        nomfichbav = ' ' 
        nomfichbav(1:l-4) = convolution_out(1:l-4) 
        nomfichbav(l-3:l+4) = '_bav.txt' 
        open(3, file = nomfichbav, status='unknown',iostat=istat)
        if( istat /= 0 ) call write_open_error(nomfichbav,istat,1)
        bav_open = .true.
        write(3,'(A/A/A)') Revision, com_date, com_time
      endif

      if( Check_conv .and. .not. Deuxieme ) write(3,115)

      if( Scan_a .and. .not. Scan_true ) then
        do ifich = 1,nfich
          mots = ' '
          mot = fichin(ifich)
          l = len_trim(mot) - 4
          if( nfich == 1 ) then
            ns = 0
          else
            if( mot(l-2:l-2) == '_' ) then
              ns = 3
            elseif( mot(l-1:l-1) == '_' ) then
              ns = 2
            else
              ns = 0
            endif
          endif
          mots(1:l-ns) = mot(1:l-ns)
          mots(l-ns+1:l-ns+5) = '_scan'
          if( ns > 0 ) mots(l-ns+6:l+5) = mot(l-ns+1:l)
          mots(l+6:l+9) = '.txt'
          fichscanin(ifich) = mots
        end do
        Scan_true = .true.
      endif

      Dafs_bio = .false.
      if( Scan_true ) then
        open(2, file = fichscanin(1), status='old',iostat=istat) 
        if( istat /= 0 )
     &    call write_open_error(fichscanin(1),istat,1)
        n = nnombre(2,Length_line)
        read(2,*) n
        if( n == 4 ) Dafs_bio = .true.
        close(2)
      endif

      if( scan_true .and. fichscanout == ' ' .and. .not. Dafs_bio ) then
        mot = convolution_out 
        l = len_trim(mot)
        mot(l-7:l+5) = 'scan_conv.txt'
        fichscanout = mot
      endif
     
      if( chem ) then
        long = len_trim(chemin)
        do ifich = 1,nfich
          mot = fichin(ifich)
          longf = len_trim(mot)
          fichin(ifich) = chemin(1:long) // mot(1:longf)
        end do
        longf = len_trim(convolution_out)
        convolution_out = chemin(1:long) // convolution_out(1:longf)
        if( scan_true ) then
          do ifich = 1,nfich
            mot = fichscanin(ifich)
            longf = len_trim(mot)
            fichscanin(ifich) = chemin(1:long) // mot(1:longf)
          end do
          longf = len_trim(fichscanout)
          fichscanout = chemin(1:long) // fichscanout(1:longf)
        endif
      endif

      if( .not. ( seah .or. Arc ) .and. nelor == 0 ) then
        nelor = 1
        allocate( Elor(nelor) )
        allocate( betalor(nelor) )
        Elor(1) = 0._db
        betalor(1) = 0._db
      endif

      do ifich = 1,nfich
        if( convolution_out /= fichin(ifich) ) cycle
        write(6,120) 
        read(5,*) rep
        if( rep /= 'y' .and. rep /= 'Y' .and. rep /= 'o' .and.
     &      rep /= 'O' ) stop
      end do
      
! -- Dimensionnement des tableaux -------------------------------------

      ninitlm = 1
      do ifich = 1,nfich
        open(2, file = fichin(ifich), status='old', iostat=istat) 
        if( istat /= 0 ) call write_open_error(fichin(ifich),istat,1)
        n = nnombre(2,13200)
        if( n > 8 ) then
          read(2,*) Eseuil(ifich), numat, nseuil, jseuil,
     &                fpp_avantseuil, v0muf, En_fermi(ifich), ninit
          if( n == ninit + 10 ) then
            new_format = .true.
            ninitl(ifich) = ninit
          else
            new_format = .false.
            ninitl(ifich) = n - 8
          endif
          ninitlm = max( ninitlm, ninitl(ifich) )
        endif
        Close(2)
      end do

      allocate( Eph1(nfich) )
      allocate( Epsii(ninitlm,nfich) )
      allocate( decal_initl(ninitlm,nfich) )
      allocate( nsup(ninitlm,nfich) )
      allocate( ne_initl(ninitlm,nfich) )

      Epsii(:,:) = 0._db
      decal_initl(:,:) = 0._db

      do ifich = 1,nfich
        open(2, file = fichin(ifich), status='old', iostat=istat) 
        n = nnombre(2,132)
        v0muf = 10._db 
        Select case(n)
          Case(1)
            read(2,*) Eseuil(ifich)
            numat = 1
          Case(2)
            read(2,*) Eseuil(ifich), numat
          Case(4)
            read(2,*) Eseuil(ifich), numat, nseuil
          Case(5)
            read(2,*) Eseuil(ifich), numat, nseuil, jseuil, v0muf
          Case(6)
            read(2,*) Eseuil(ifich), numat, nseuil, jseuil,
     &                fpp_avantseuil, v0muf
          Case(7)
            read(2,*) Eseuil(ifich), numat, nseuil, jseuil,
     &                fpp_avantseuil, v0muf, En_fermi(ifich)
          Case(8)
            read(2,*) Eseuil(ifich), numat, nseuil, jseuil,
     &                fpp_avantseuil, v0muf, En_fermi(ifich),
     &                Epsii(1,ifich)
          Case default
            if( new_format ) then
              read(2,*) Eseuil(ifich), numat, nseuil, jseuil,
     &          fpp_avantseuil, v0muf, En_fermi(ifich), ninitl(ifich), 
     &          ninit1, Epsii(1:ninitl(ifich),ifich), Densite_atom
            else
              read(2,*) Eseuil(ifich), numat, nseuil, jseuil,
     &                fpp_avantseuil, v0muf, En_fermi(ifich), ninit1,
     &                Epsii(1:ninitl(ifich),ifich)
            endif
            if( ninit1 < 0 ) then
              Green_int = .true.
              ninit1 = abs( ninit1 )
            else
              Green_int = .false.
            endif
        end select
        if( abs(fpp_avantseuil) > 1.e-10_db ) then
          if( fpp_avantseuil < - 1.e-10_db ) then
            Full_self_abs = .true.
            fpp_avantseuil = - fpp_avantseuil
          else
            self_abs = .true.
          endif
        endif
        Cor_abs = Self_abs .or. Full_self_abs
        n = nnombre(2,Length_line)
        if( n > 0 ) then
          npldafs = n / 2
          read(2,*)
          read(2,*)
          if( Cor_abs ) read(2,*) 
        endif
        read(2,'(10x,a13)') nomab
        nomab = adjustl( nomab )
        nt = ( nnombre(2,Length_line) - 1 ) / ninitl(ifich)
        if( n > 0 ) then
          if( Full_self_abs ) then
            nxan = nt - 6 * npldafs
          elseif( Self_abs ) then
            nxan = nt - 4 * npldafs
          else
            nxan = nt - 2 * npldafs
	    endif
          if( nomab(1:8) == 'Sum_dd_r' ) Signal_Sph = .true.
        elseif( fprim .and. nomab(1:1) == 'D' ) then
          Tenseur = .true.
          nxan = 0
          if( nomab(1:4) == 'D_xx' ) Tenseur_car = .true.
          Magn = nomab(1:6) == 'D_xx_r'
          if( Magn ) then
            npldafs = nt / 2
          else
            npldafs = nt
          endif
        else
          nxan = nt
        endif

        if( ifich == 1 ) then
          nxan1 = nxan
          npldafs1 = npldafs
        else
          if( nxan_lib ) nxan = min(nxan1,nxan)
          if( nxan1 /= nxan .or. npldafs1 /= npldafs ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,130) nxan1, npldafs1, ifich, nxan, npldafs
            end do
            stop
          endif
        endif

        do ie = 1,10000
          Read(2,*,err=1030,end=1030) eph
          if( ie == 1 ) eph1(ifich) = eph
          if( eph > eintmax ) exit
        end do
 1030   ne(ifich) = ie - 1

        if( ne(ifich) < 2 ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,140) fichin(ifich)
          end do
          stop
        endif

        Close(2)
      end do

      if( .not. Full_self_abs ) Circular = .false.
      if( .not. Cor_abs ) Double_cor = .false.

      if( .not. ( jseuil == 1 .or. ( jseuil < 4 .and. numat > 20 ) ) )
     &   no_extrap = .true.

! Modification en cas de fit.
      if( Fit_cal ) then

        do igr = 2,ngroup_par
          istop = 0
          do ipar = 1,npar(igr)
            if( typepar(igr,ipar) /= 'shift' .and.
     &          typepar(igr,ipar) /= 'weight' ) cycle
            if( indice_par(igr,ipar) > nfich ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,145) typepar(igr,ipar), indice_par(igr,ipar),
     &                         nfich
              end do
              istop = 1
             endif
          end do
        end do
        if( istop == 1 ) stop

        do igr = 2,ngroup_par
          do ipar = 1,npar(igr)
            select case( typepar(igr,ipar) )
              case('aseah')
                asea = param(igr,ipar)
              case('ecent')
                Ecent = param(igr,ipar)
                Ecent = Ecent / rydb
              case('gaussian')
                deltar = param(igr,ipar)
              case('vibr')
                vibration = param(igr,ipar)
              case('elarg')
                Elarg = param(igr,ipar)
                Elarg = Elarg / rydb
              case('gamma_max')
                Gamma_max = param(igr,ipar)
                Gamma_max = Gamma_max / rydb
              case('gamma_hol')
                Gamma_hole(1) = param(igr,ipar)
                Gamma_hole(1) = Gamma_hole(1) / rydb
                Gamma_hole_imp = .true.
              case('efermi')
                Efermi = param(igr,ipar)
                E_Fermi_man = .true.
              case('shift')
                if( .not. run_done( indice_par(igr,ipar) ) ) then
                  ifich = i_done( indice_par(igr,ipar) )
                  decal( ifich ) = param(igr,ipar)
                  Shift_auto = .false.
                endif
              case('weight')
                if( .not. run_done( indice_par(igr,ipar) ) ) then
                  ifich = i_done( indice_par(igr,ipar) )
                  Pds( ifich ) = param(igr,ipar)
                endif
            end select
          end do
        end do
      endif

      deallocate( i_done )

      Esmin = Eseuil(1)
      do ifich = 2,nfich
        Esmin = min( Esmin, Eseuil(ifich) )
      end do
      do ifich = 1,nfich
        decal(ifich) = decal(ifich) + Eseuil(ifich) - Esmin
      end do

      if( Epsii_ref_man ) then
        Epsii_moy = Epsii_ref
      elseif( Shift_auto ) then
        Epsii_moy = 0._db
        do ifich = 1,nfich
          Epsii_moy = Epsii_moy + sum( Epsii(1:ninitl(ifich),ifich) )
     &                          / ninitl(ifich)
        end do
        Epsii_moy = Epsii_moy / nfich
      endif
      if( Epsii_ref_man .or. Shift_auto ) then
        do ifich = 1,nfich
          decal(ifich) = decal(ifich) + Epsii_moy
     &             - sum( Epsii(1:ninitl(ifich),ifich) ) / ninitl(ifich)
        end do
      endif

      decal_initl(:,:) = 0._db
      do ifich = 1,nfich
        select case( ninitl(ifich) )
          case(1)
            Epsii_moy = Epsii(1,ifich)
          case(2)
             if( ninit1 == 2 ) then
               Epsii_moy = sum( Epsii(1:ninitl(ifich),ifich) )
     &                          / ninitl(ifich)
             else
               Epsii_moy = Epsii(2,ifich)
             endif 
          case(4,6,10)
             if( ninit1 /= ninitl(ifich) ) then
               Epsii_moy = sum( Epsii(ninit1+1:ninitl(ifich),ifich) )
     &                   / ( ninitl(ifich) - ninit1 ) 
             else
               Epsii_moy = sum( Epsii(1:ninitl(ifich),ifich) )
     &                   / ninitl(ifich)
             endif
        end select 
        decal_initl(1:ninitl(ifich),ifich)
     &                  = Epsii(1:ninitl(ifich),ifich) - Epsii_moy
     &                  + decal(ifich)
        if( abs(Delta_edge) > eps10 .and. ( ninitl(ifich) /= ninit1 ) )
     &        decal_initl(1:ninit1,ifich)
     &             = decal_initl(1:ninit1,ifich) + Delta_edge * Rydb
      end do

      deallocate( decal) 

      if( npldafs > 0 ) then
        Dafs = .true.
        allocate( dph(npldafs) )
        allocate( dpht(npldafs) )
        dpht(:) = (0._db,0._db)            
        allocate( fr(npldafs) ); allocate( fi(npldafs) )
        allocate( f0(npldafs) )
        allocate( nphi(npldafs) )
        if( Cor_abs .or. Dafs_bio ) then
          npldafs_b = npldafs
        else
          npldafs_b = 0
        endif
        allocate( hkl_dafs(3,npldafs_b) ) 

        if( Thomson ) then
          n = min( npldafs, npldafs_th )
          f0(1:n) = f0_th(1:n)
          deallocate( f0_th )
        else  
          f0(:) = (0._db, 0._db )
          Pdt = 0._db
          do ifich = 1,nfich
            if( Skip_run(ifich) ) cycle
            open(2, file = fichin(ifich), status='old', iostat=istat) 
            if( istat /= 0) call write_open_error(fichin(ifich),istat,1)
            read(2,*)
            n = nnombre(2,Length_line)
            if( n > 0 ) then
              read(2,*) ( fr(ipl), fi(ipl), ipl = 1,npldafs )
              f0(:) = f0(:) + Pds(ifich) * cmplx( fr(:), fi(:),db) 
              Pdt = Pdt + Pds(ifich) 
            endif
            Close(2)
          end do
          if( abs(Pdt) > 1e-10_db ) f0(:) = f0(:) / Pdt 
        endif

      endif

      if( E_Fermi_man ) En_Fermi(:) = Efermi

      if( Gaussian_default ) deltar = Eseuil(1) / 10000

      deltar = deltar / sqrt( 8 * log(2._db) )
      vibration = vibration / sqrt( 8 * log(2._db) )

      decal_initl(:,:) = decal_initl(:,:) / rydb
      En_fermi(:) = En_fermi(:) / rydb
      Esmin = Esmin / rydb
      Eph1(:) = Eph1(:) / rydb
      Epsii(:,:) = Epsii(:,:) / rydb
      Eseuil(:) = Eseuil(:) / rydb
      v0muf = v0muf / rydb
      deltar = deltar / rydb
      if( Fermip ) Efermip(1:nfich) = Efermip(1:nfich) / rydb

      if( .not. Arc .and. .not. seah ) Elor(:) = Elor(:) - En_fermi(1)

! Elaboration de la grille en energie

      Estart_l = Estart
      do ifich = 1,nfich
        do initl = 1,ninitl(ifich)
          E = Eph1(ifich) + decal_initl(initl,ifich)
          Estart_l = Min( Estart_l, E )
        end do          
      end do 

      pasdeb = 0.5_db / rydb
      do ifich = 1,nfich
        do initl = 1,ninitl(ifich)
          E = Eph1(ifich) + decal_initl(initl,ifich)
          if( E > Estart_l - 1.e-10_db ) then
            nsup(initl,ifich) = nint( ( E - Estart_l ) / pasdeb )
            ne_initl(initl,ifich) = ne(ifich) + nsup(initl,ifich)
          else
            ne_initl(initl,ifich) = ne(ifich)
            nsup(initl,ifich) = 0
          endif
        end do
      end do 

      nemax = maxval( ne_initl )

      n = sum( ninitl(:) ) 
      allocate( Ef(nemax,n) )

      jfich = 0

      do ifich = 1,nfich

        open(2, file = fichin(ifich), status='old', iostat=istat) 
        if( istat /= 0 ) call write_open_error(fichin(ifich),istat,1)

        read(2,*)
        n = nnombre(2,Length_line)
        if( n > 0 ) then
          read(2,*)
          read(2,*)
          if( Cor_abs ) read(2,*) 
        endif
        read(2,*)  

        jfich = jfich + 1

        do ie = nsup(1,ifich)+1,ne_initl(1,ifich)
          Read(2,*) Ef(ie,jfich)
        end do
        j0 = jfich
        Ef(:,j0) = Ef(:,j0) / rydb

        do initl = 2,ninitl(ifich)
          jfich = jfich + 1
          n1 = nsup(initl,ifich) + 1 
          n2 = ne_initl(initl,ifich) 
          n3 = nsup(1,ifich) + 1 
          n4 = ne_initl(1,ifich) 
          Ef(n1:n2,jfich) = Ef(n3:n4,j0)
        end do

        jfich = j0 - 1
        do initl = 1,ninitl(ifich)
          jfich = jfich + 1
          do ie = nsup(initl,ifich),1,-1
            Ef(ie,jfich) = Ef(ie+1,jfich) - pasdeb
          end do
        end do

        jfich = j0 - 1
        do initl = 1,ninitl(ifich)
          jfich = jfich + 1
          n = ne_initl(initl,ifich) 
          Ef(1:n,jfich) = Ef(1:n,jfich) + decal_initl(initl,ifich)
        end do

        Close(2)
      end do

      Energphot = .false.

      jfich = 0
      do ifich = 1,nfich
        do initl = 1,ninitl(ifich)
          jfich = jfich + 1
          if( Ef(ne_initl(initl,ifich),jfich) <= Eseuil(ifich) .or.
     &        Eseuil(ifich) <= Ef(ne_initl(initl,ifich),jfich)
     &                                          - Ef(1,jfich) ) cycle
          Energphot = .true.
          Ef(1:ne_initl(initl,ifich),jfich)
     &       = Ef(1:ne_initl(initl,ifich),jfich) - Eseuil(ifich)
        end do
      end do

      Emin = Ef(1,1)
      Emax = Ef(ne_initl(1,1),1)

      jfich = 0
      do ifich = 1,nfich
        do initl = 1,ninitl(ifich)
          jfich = jfich + 1
          Emin = min( Ef(1,jfich), Emax )
          Emax = min( Ef(ne_initl(initl,ifich),jfich), Emax )
        end do
      end do

      nfich_tot = sum( ninitl(1:nfich) )
      allocate( fichdone(nfich_tot) )

      nes = 10000

      do i = 1,2

        fichdone(:) = .false.

        jfich = 0
        boucle_0: do ifich = 1,nfich
          do initl = 1,ninitl(ifich)
            jfich = jfich + 1
            if( abs( Ef(1,jfich) - Emin ) < 1.e-10_db ) exit boucle_0 
          end do
        end do boucle_0
        ifichref = ifich
        jfichref = jfich
        initlref = initl
        fichdone(jfichref) = .true.

        je = 0
        do ie = 1,nes
          je = je + 1
          E = Ef(je,jfichref)
          if( i == 2 ) Es(ie) = E
          if( E > Emax - 1.e-10_db ) exit

          if( je == ne_initl(initlref,ifichref) ) then 
            jfich = 0
            boucle_1: do ifich = 1,nfich
              do initl = 1,ninitl(ifich)
                jfich = jfich + 1

                if( Ef(ne_initl(initl,ifich),jfich) > E+1.e-10_db ) then
                  do je = ne_initl(initl,ifich),1,-1
                    if( Ef(je,jfich) < E - 1.e-10_db ) exit
                  end do
                  ifichref = ifich 
                  jfichref = jfich
                  initlref = initl
                  fichdone(jfichref) = .true.
                  exit boucle_1
                endif  
              end do
            end do boucle_1
          endif
  
          jfich = 0
          boucle_2: do ifich = 1,nfich
            do initl = 1,ninitl(ifich)
              jfich = jfich + 1
              if( fichdone(jfich) ) cycle
              if( E > Ef(nsup(initl,ifich)+1,jfich) + 1.e-10_db ) then
                do je = nsup(initl,ifich)+1,ne_initl(initl,ifich)
                  if( Ef(je,jfich) > E + 1.e-10_db ) exit
                end do 
                je = je - 1
                ifichref = ifich 
                jfichref = jfich
                initlref = initl
                fichdone(jfichref) = .true.
                cycle
              endif
            end do
          end do boucle_2

        end do

        if( i == 1 ) then
          if( E > Emax + 1.e-10_db ) then
            nes = ie - 1
          else
            nes = ie
          endif
          allocate( Es(nes) )
        endif

      end do

      deallocate( fichdone )
      deallocate( Ef )

      if( nes == 0 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,147) 
        end do
        stop
      endif

! La preparation est terminee.

      ifich = 0
      do jfich = 1,mfich
        if( run_done(jfich) ) cycle
        ifich = ifich + 1

        if( Fermip ) then
          EFermi_orig = Efermip(ifich)
        else
          EFermi_orig = En_fermi(ifich)
        endif

        do initl = 1,ninitl(ifich)

          if( n_selec_core /= 0 ) then
            do i = 1,n_selec_core
              if( initl == num_core(i) ) exit
            end do
            if( i > n_selec_core ) cycle 
          endif

          if( decferm ) then  
            EFermi = EFermi_orig ! Fermi ne suit pas le decalage de la gamme 
          else
            EFermi = EFermi_orig + decal_initl(initl,ifich)
          endif
          nenerg = ne_initl(initl,ifich) 

          if( .not. scan_true ) then
            allocate( Adafs(nenerg,1,npldafs) )
            Adafs(:,:,:) = (0._db, 0._db)
            if( Cor_abs ) then
              allocate( mu(nenerg,1,npldafs,2) )
              mu(:,:,:,:) = (0._db, 0._db)
            endif
          endif 
          allocate( Xanes(nenerg,nxan) )
          Xanes(:,:) = 0._db
 
! -- Lecture -----------------------------------------------------------

          open(2, file = fichin(ifich), status='old', iostat=istat) 
          if( istat /= 0 ) call write_open_error(fichin(ifich),istat,1)
 
          read(2,*)
         
          n = nnombre(2,Length_line)
          if( n > 0 ) then
            read(2,*) ( fr(ipl), fi(ipl), ipl = 1,npldafs )
            read(2,*) ( fr(ipl), fi(ipl), ipl = 1,npldafs)
            dph(:) = cmplx( fr(:), fi(:),db )
            dpht(:) = dpht(:) + Pds(ifich) * dph(:)
	      if( Cor_abs ) then 
              read(2,*)  natomsym, axyz(:), angxyz(:),
     &                  ( hkl_dafs(:,ipl), ipl = 1,npldafs)
              axyz(:) = axyz(:) / bohr
            endif 
          elseif( Tenseur ) then
            f0(:) = ( 0._db, 0._db ) 
            dph(:) = ( 0._db, 0._db )
            dph(1) = ( 1._db, 0._db )
            if( Tenseur_car ) then
              dph(4) = ( 1._db, 0._db )
              dph(6) = ( 1._db, 0._db )
            endif 
            dpht(:) = dpht(:) + Pds(ifich) * dph(:)
          endif

          read(2,*)

          allocate( Energ(nenerg) )

          if( Dafs ) then

            nphi(:) = 1

            if( Cor_abs ) allocate( mua_r(npldafs,2) )
            if( Full_self_abs ) allocate( mua_i(npldafs,2) )

            do ie = nsup(initl,ifich)+1,nenerg
              if( Tenseur .and. .not. Magn ) then
                Read(2,*) Energ(ie), ((fr(ipl), ipl = 1,npldafs),
     &                                                   j = 1,initl)
                if( .not. scan_true ) Adafs(ie,1,:)
     &                                = cmplx(fr(:), 0._db,db)
              elseif( nxan > 0 ) then
                if( Full_self_abs ) then
                  Read(2,*) Energ(ie),(( Xanes(ie,ipl), ipl = 1,nxan ),  
     &               ( fr(ipl), fi(ipl), (mua_r(ipl,i), mua_i(ipl,i), 
     &               i = 1,2), ipl = 1,npldafs ), j = 1,initl) 
                elseif( self_abs ) then
                  Read(2,*) Energ(ie),(( Xanes(ie,ipl), ipl = 1,nxan ),  
     &               ( fr(ipl), fi(ipl), (mua_r(ipl,i), 
     &               i = 1,2), ipl = 1,npldafs ), j = 1,initl) 
                else
                  Read(2,*) Energ(ie), (( Xanes(ie,ipl), ipl = 1,nxan ),  
     &            ( fr(ipl), fi(ipl), ipl = 1,npldafs ),j = 1,initl)
                endif
              else
                if( Full_self_abs ) then
                  Read(2,*) Energ(ie),  
     &              (( fr(ipl), fi(ipl), ( mua_r(ipl,i), mua_i(ipl,i), 
     &              i = 1,2), ipl = 1,npldafs ), j = 1,initl)
                elseif( self_abs ) then
                  Read(2,*) Energ(ie),  
     &              (( fr(ipl), fi(ipl), ( mua_r(ipl,i), 
     &              i = 1,2 ), ipl = 1,npldafs ), j = 1,initl)
                else
                  Read(2,*) Energ(ie), ( ( fr(ipl),fi(ipl),
     &                     ipl = 1,npldafs), j = 1,initl )
                endif
              endif
              if( .not. scan_true ) then
                Adafs(ie,1,:) = cmplx( fr(:), fi(:), db)
                if( Full_self_abs ) mu(ie,1,:,:)
     &                    = cmplx( mua_r(:,:), mua_i(:,:), db)
                if( Self_abs ) mu(ie,1,:,:)
     &                    = cmplx( mua_r(:,:), 0._db, db)
              endif
            end do

            if( Cor_abs ) deallocate( mua_r)
            if( Full_self_abs ) deallocate( mua_i )

          else

            do ie = nsup(initl,ifich)+1,nenerg
              Read(2,*) Energ(ie), (( Xanes(ie,ipl), ipl = 1,nxan ),
     &                                   j = 1,initl)
            end do

          endif

          close(2)

          nphim = 1
          if( scan_true ) then
            open(2, file = fichscanin(ifich), status='old',iostat=istat) 
            if( istat /= 0 )
     &        call write_open_error(fichscanin(ifich),istat,1)
            n = nnombre(2,Length_line)
            do ipl = 1,npldafs
              read(2,*) nphi(ipl)
              nphim = max( nphim, nphi(ipl) ) 
            end do

            if( ifich == 1 .and. ( initl == 1 .or. initl ==num_core(1)))
     &        then
              allocate( angle(nphim) )
              allocate( f0scan(nphim,npldafs) )
              allocate( phdtscan(nphim,npldafs) )
            endif

            allocate( Adafs(nenerg,nphim,npldafs) )
            if( Cor_abs ) allocate( mu(nenerg,nphim,npldafs,2) )

            Adafs(:,:,:) = ( 0._db, 0._db )
            if( Cor_abs ) mu(:,:,:,:) = ( 0._db, 0._db )
            do ie = nsup(initl,ifich)+1,nenerg
              read(2,*) mots
              do ipl = 1,npldafs
                if( Dafs_bio ) then
                  read(2,*) hkl_dafs(:,ipl)
                else
                  read(2,*) mots
                endif
                do i = 1,nphi(ipl)
                  if( Self_abs ) then
                    read(2,*) angle(i),
     &                  ( a, b, a3, a4, j = 1,initl ),
     &                  ( c, d, c, d, j = initl+1,ninitl(ifich) ),
     &                        a1, b1, a2, b2 
                    mu(ie,i,ipl,1) = cmplx( a3, 0._db, db )
                    mu(ie,i,ipl,2) = cmplx( a4, 0._db, db )
                  elseif( Full_self_abs ) then
                    read(2,*) angle(i),
     &                  ( a, b, a3, b3, a4, b4, j = 1,initl ),
     &                  ( c, d, c, d, c, d, j = initl+1,ninitl(ifich) ),
     &                        a1, b1, a2, b2 
                    mu(ie,i,ipl,1) = cmplx( a3, b3, db )
                    mu(ie,i,ipl,2) = cmplx( a4, b4, db )
                  else
                    if( Dafs_bio ) then
                      read(2,*) ( a, b, j = 1,initl ),
     &                       ( c, d, j = initl+1,ninitl(ifich) ),
     &                        a1, b1, a2, b2
                    else
                      read(2,*) angle(i), ( a, b, j = 1,initl ),
     &                       ( c, d, j = initl+1,ninitl(ifich) ),
     &                        a1, b1, a2, b2
                    endif
                  endif
                  Adafs(ie,i,ipl) = cmplx( a, b, db)
                  f0scan(i,ipl) = cmplx( a1, b1, db)
                  phdtscan(i,ipl) = cmplx( a2, b2, db)
                end do
              end do
            end do
            close(2)

          endif

          Energ(:) = Energ(:) / rydb

          do ie = nsup(initl,ifich),1,-1
            Energ(ie) = Energ(ie+1) - pasdeb
          end do

          if( Energphot ) Energ(1:nenerg)
     &                    = Energ(1:nenerg) - Eseuil(ifich)
          Energ(1:nenerg) = Energ(1:nenerg) + decal_initl(initl,ifich)

          if( ifich == 1 .and. ( initl == 1 .or. initl == num_core(1) ))
     &      then
            allocate( p1f(nes) )
            allocate( p2f(nes) )
            allocate( indf(nes) )
          endif
          do ie = 1,nes
            E = Es(ie)
            E = max(E ,Energ(1) )
            E = min(E ,Energ(nenerg) )
            do i = 2,nenerg
              if( Energ(i) > E - 1.e-10_db ) exit
            end do
            p1 = ( E - Energ(i-1) ) / ( Energ(i) - Energ(i-1) )
            p2 = 1 - p1
            p1f(ie) = Pds(ifich) * p1
            p2f(ie) = Pds(ifich) * p2
            indf(ie) = i
          end do

          if( ifich == 1 .and. ( initl == 1 .or. initl == num_core(1) ))
     &      then
            allocate( Xs(nes,nxan) )
            Xs(:,:) = 0._db
            if( Dafs ) then
              allocate( As(nes,nphim,npldafs) )
              As(:,:,:) = (0._db,0._db)
              if( Cor_abs ) then
                allocate( mus(nes,nphim,npldafs,2) )
                mus(:,:,:,:) = (0._db,0._db)
              endif
            endif
          endif

          if( Conv_done .or. Skip_run(ifich) ) goto 1040

          nenerge = nenerg + njp

          if( Cor_abs ) then
            allocate( Ad(nenerg,nphim,3*npldafs) )
          else
            allocate( Ad(nenerg,nphim,npldafs) )
          endif
          allocate( Xa(nenerge,nxan) )

          if( seah .or. Arc ) then
            nelor = nenerg
            allocate( Elor(nelor) )
            allocate( betalor(nelor) )
            Elor(:) = Energ(:)
          endif

          if( Photo_emission .or. Green_int ) then 
            betalor(:) = 0._db
          elseif( seah ) then 
            call seahdench(asea,Efermi,Gamma_max,nelor,Elor,betalor)
          elseif( Arc ) then
            call gammarc(Ecent,Elarg,Gamma_max,Efermi,nelor,Elor,
     &                   betalor)
          endif

          if( .not. Gamma_hole_imp ) then
            if( ninitl(ifich) == 1 .or. ninit1 == ninitl(ifich) 
     &         .or. initl <= ninit1 ) then
              js = jseuil
            else
              js = jseuil + 1
            endif
            call tab_width(Eseuil(1),Gamma_h,js,nseuil,numat)
          elseif( ngamh == 1 .or. ninitl(ifich) == 1 ) then
            Gamma_h = Gamma_hole(1)
          elseif( ngamh == ninitl(ifich) ) then 
            Gamma_h = Gamma_hole(initl)
          elseif( initl <= ninit1 ) then
            Gamma_h = Gamma_hole(1)
          else
            Gamma_h = Gamma_hole(2)
          endif
          betalor(:) = betalor(:) + Gamma_h

          if( Check_conv ) then
            ipr1 = 3
          else
            ipr1 = 6
          endif
          ipr2 = 6
          do ipr = ipr1,ipr2,3
            if( ifich == 1 .and. 
     &          ( initl == 1 .or. initl == num_core(1) ) )  then
              if( seah ) then
                Write(ipr,'(/A)') ' Seah-Dench model'
                Write(ipr,150) Gamma_max*rydb, asea
              elseif( Arc ) then
                Write(ipr,'(/A)') ' Arctangent model'
                Write(ipr,160) Gamma_max*rydb, Ecent*rydb, Elarg*rydb
              endif
            endif
            if( nfich == 1 .and. ninitl(ifich) == 1 ) then
              Write(ipr,170) Gamma_h*rydb, Efermi*rydb 
            elseif( nfich == 1 ) then
              Write(ipr,172) Gamma_h*rydb, Efermi*rydb, initl 
            elseif( ninitl(ifich) == 1 ) then
              Write(ipr,174) Gamma_h*rydb, Efermi*rydb, ifich
            else
              Write(ipr,176) Gamma_h*rydb, Efermi*rydb, ifich, initl 
            endif 
          end do

          if( ifich == 1 .and. ( initl == 1 .or. initl == num_core(1) )
     &      .and. Gamma_max > eps10 ) then
            de_obj = ( Elor(nelor) - Elor(1) ) / 15
            E_obj = Elor(1) - 0.00001_db 
            do ie = 1,nelor
              Gamm = betalor(ie) - Gamma_h
              if( Gamm > 0._db ) exit
            end do
            n = ie - 1
            do ie = 1,nelor
              if( ie > 1 .and. ie < n ) cycle
              E = Elor(ie) - V0muf
              Gamm = betalor(ie) - Gamma_h 
              if( Gamm > 0._db ) then
                alambda = sqrt( 2 / ( sqrt( E**2 + Gamm**2 ) - E ) ) 
              else
                alambda = 0._db
              endif
              do ipr = ipr1,ipr2,3
                if( ie == 1 ) Write(ipr,178)  
                if( ipr == 6 ) then
                  if( Elor(ie) < E_obj .and. ie /= nelor ) cycle
                  E_obj = E_obj + de_obj
                endif
                if( alambda * bohr < 10000000 ) then 
                  Write(ipr,180) Elor(ie) * Rydb, betalor(ie) * Rydb, 
     &                           alambda * bohr
                else
                  Write(ipr,185) Elor(ie) * Rydb, betalor(ie) * Rydb, 
     &                           alambda * bohr
                endif
              end do
            end do
          endif

          allocate( dampl(nenerg) )
          allocate( bb(nenerge) )
          allocate( Ephoton(nenerge) )
          allocate( e1(nenerge) )
          allocate( e2(nenerge) )
          Ephoton(:) = 0._db
          allocate( lori(nenerge) )
          allocate( lorix(nenerge) )
          if( Dafs ) allocate( lorr(nenerge) )

          call cflor(bb,betalor,Efermi,Ephoton,Elor,Energ,ie1,ie2,nef,
     &                   nelor,nenerg,nenerge,Photo_emission)
          Ephoton(:) = Ephoton(:) + Eseuil(ifich)

          Extrap = .false.
          if( Dafs ) then
            if( .not. no_extrap ) then
              do ipl = 1,npldafs
                if( abs( dph(ipl) ) > 1.e-10_db ) then
                  Extrap = .true.
                  exit
                endif
              end do
            endif
          endif

          fprime_atom = ( icheck > 2 ) .and. ifich == 1 .and. 
     &          ( initl == 1 .or. initl == num_core(1) ) 
          if( Extrap .or. fprime_atom .or. Cor_abs ) then
            call extrapat(bb(nenerg),
     &          dampl,Energ,Eseuil(ifich),Extrap,fpp0,
     &          fprime_atom,numat,nenerg)
             dampl(:) = dampl(:) / ninitl(ifich)
             fpp0 = fpp0 / ninitl(ifich)
          endif

! Convolution par la lorentzienne

          if( Dafs ) Ad(1:nenerg,:,1:npldafs)
     &                             = Adafs(1:nenerg,:,1:npldafs)
          if( Cor_abs ) then
            do i = 1,2
	        do ipl = 1,npldafs
                Ad(1:nenerg,:,i*npldafs+ipl) = mu(1:nenerg,:,ipl,i)
              end do
            end do
          endif

          do ipl = 1,nxan
            Xa(1:nenerg,ipl) = Xanes(1:nenerg,ipl)
! Extrapolation
            do ie = nenerg+1,nenerge 
              Xa(ie,ipl) = Xa(nenerg,ipl)
            end do
          end do

          gamma = .false.
          do ie = nef,nenerge
            if( abs( bb(ie) ) > 1.e-10_db ) then
              gamma = .true.
              exit
            endif
          end do

          ne2 = ie2

          do ie = ie1,ie2
            if( ie == ie1 ) then
              if( Photo_emission ) then
                e1(ie) = 1.5 * Ephoton(1) - 0.5 * Ephoton(2)
              else
                e1(ie) = Efermi + Eseuil(ifich)
                if( ie == 1 ) then
                  e1m = 1.5 * Ephoton(1) - 0.5 * Ephoton(2)
                  e1(ie) = max( e1(ie), e1m )
                endif
              endif
            else
              e1(ie) = 0.5 * ( Ephoton(ie) + Ephoton(ie-1) )
            endif
            if( ie == ne2 ) then
              e2(ie) = 1.5 * Ephoton(ne2) - 0.5 * Ephoton(ne2-1)
            else 
              e2(ie) = 0.5 * ( Ephoton(ie) + Ephoton(ie+1) )
            endif
          end do

          if( Photo_emission ) then
            nen2 = ie2
          else
            nen2 = nenerg
          endif

          if( gamma .and. .not. Green_int ) then

            do ie = 1,nenerg
  
              do j = ie1,ne2

                if( .not. Gamma_var ) then 
                  bba = bb(ie)
                  if( bba >= 0._db ) then
                    bbb = max( bba, 1.E-08_db )
                  else
                    bbb = bba
                  endif
                  de2 = ( e2(j) - Ephoton(ie) ) / bbb
                  de1 = ( e1(j) - Ephoton(ie) ) / bbb
                  lorix(j) = atan( de1 ) - atan( de2 )
                endif
                if( Dafs .or. Gamma_var ) then 
                  bba = bb(j)
                  if( bba >= 0._db ) then
                    bbb = max( bba, 1.E-08_db )
                  else
                    bbb = bba
                  endif
                  de2 = ( e2(j) - Ephoton(ie) ) / bbb
                  de1 = ( e1(j) - Ephoton(ie) ) / bbb
                  lori(j) = atan( de1 ) - atan( de2 )
                  if( Gamma_var ) lorix(j) = lori(j)
                  if( Dafs .and. j <= nenerg )
     &              lorr(j) = 0.5 * log( (1 + de1**2) / (1 + de2**2) )
                endif
              end do

              do ipl = 1,nxan
                Xanes(ie,ipl) = - sum( lorix(ie1:ne2)
     &                        * Xa(ie1:ne2,ipl) ) / pi
                write(98,*) Xanes(ie,ipl)  
              end do

              if( .not. Green_plus ) lori(ie1:nen2) = - lori(ie1:nen2) 

              if( .not. Green_int ) then

                do ipl = 1,npldafs
                  do ip = 1,nphi(ipl)

                    if( Tenseur ) then
    
                      Adafs(ie,ip,ipl) = sum( 
     &                      cmplx( lorr(ie1:nen2), lori(ie1:nen2),db ) 
     &                     * real( Ad(ie1:nen2,ip,ipl),db ) )

                    else

                      Adafs(ie,ip,ipl) = sum( 
     &                        cmplx(lorr(ie1:nen2),lori(ie1:nen2),db)
     &                      * Ad(ie1:nen2,ip,ipl) )

! Ici la partie imaginaire devient le Kramers Koenig de la partie reelle
                      if( Cor_abs ) then
                        do i = 1,2
                          jpl = ipl + i*npldafs
                          mu(ie,ip,ipl,i) = sum( 
     &                       cmplx( lorr(ie1:nen2), lori(ie1:nen2), db)
     &                      * Ad(ie1:nen2,ip,jpl) ) / pi
                        end do
                      endif
  
                    endif

                  end do
                end do

              endif

            end do

          else

            if( Dafs ) then

              do ie = 1,nenerg

                do j = ie1,nen2
                  de2 = e2(j) - Ephoton(ie)
                  de1 = e1(j) - Ephoton(ie)
                  if( j == ie ) then
                    if( abs( de2 + de1 ) < 1.e-10_db ) then
                      lorr(j) = 0._db     
                    elseif( abs( de1 ) < 1.e-10_db ) then
                      lorr(j) = - log( 4._db )     
                    else
                      lorr(j) = log( - de1 / de2 )
                    endif
                  else
                    lorr(j) = log( de1 / de2 ) 
                  endif
                end do

                do ipl = 1,npldafs
                  do ip = 1,nphi(ipl)

                    if( Tenseur .and. .not. Green_int ) then
                      Adafs(ie,ip,ipl) = sum( lorr(ie1:nen2)
     &                   * real( Ad(ie1:nen2,ip,ipl),db ) ) 
                      if( ie >= ie1 ) Adafs(ie,ip,ipl) =Adafs(ie,ip,ipl)
     &                       - img * pi * real( Ad(ie,ip,ipl),db )
                    else
                      if( .not. Green_int ) then
                        Adafs(ie,ip,ipl) =
     &                       sum( lorr(ie1:nen2) * Ad(ie1:nen2,ip,ipl) ) 
                        if( ie >= ie1 ) Adafs(ie,ip,ipl) =
     &                      Adafs(ie,ip,ipl) + img * pi * Ad(ie,ip,ipl)
                      endif
! Ici pour correct d'abs, le deuxi�me xanes devient le Kramers Koenig
! du premier xanes
                      if( Cor_abs ) then
                        do i = 1,2
                          jpl = ipl + i*npldafs
                          mu(ie,ip,ipl,i) = 
     &                        sum( lorr(ie1:nen2) * Ad(ie1:nen2,ip,jpl))
                          if( ie >= ie1 ) mu(ie,ip,ipl,i)
     &                                      = mu(ie,ip,ipl,i)
     &                                      - img * pi * Ad(ie,ip,jpl)
                          mu(ie,ip,ipl,i) = mu(ie,ip,ipl,i) / pi 
                        end do
                      endif
                    endif
                  end do
                end do
              end do
              
            endif

            if( Photo_emission ) then
              do ipl = 1,nxan
                Xanes(nef+1:nenerg,ipl) = 0._db
              end do
            elseif( .not. Green_int ) then
              do ipl = 1,nxan
                Xanes(1:nef-1,ipl) = 0._db
              end do
            endif

          endif 

! On passe en convention cristallo avec f" positif
! (equivaut a Green_moins)
          if( Dafs .and. Green_plus ) then
            Adafs(:,:,:) = conjg( Adafs(:,:,:) )
            if( Cor_abs ) mu(:,:,:,:) = Conjg( mu(:,:,:,:) )
          endif

! On veut f' et f" en nombre d'electrons
          if( Tenseur ) then
            do ie = 1,nenerg            
              Eph = Energ(ie) + Eseuil(ifich)
              ct_nelec = conv_mbarn_nelec(Eph)
              adafs(ie,:,:) = ct_nelec * adafs(ie,:,:)
            end do
          endif

          if( Extrap ) then
            do ipl = 1,npldafs
              if( Signal_Sph .and. ipl > 2 ) exit
              do ip = 1,nphi(ipl)
                if( ip == 1 ) then
                  cf = dph(ipl)
                else
                  cf = phdtscan(ip,ipl)
                endif
                if( Green_int ) then
                  do ie = 1,nenerg
                    Adafs(ie,ip,ipl) = Adafs(ie,ip,ipl) + img * fpp0 *cf
                  end do
                else
                  do ie = 1,nenerg
                    Adafs(ie,ip,ipl) = Adafs(ie,ip,ipl)
     &                            + ( dampl(ie) + img * fpp0 ) * cf
                  end do
                endif
              end do 
            end do 
          endif

! On passe en convention cristallo avec f" positif
! (equivaut a Green_moins)
          if( Green_plus .and. Bormann ) then
            Adafs(:,:,:) = conjg( Adafs(:,:,:) )
          endif

          if( Cor_abs ) then
            Volume_maille = Cal_Volume_maille(axyz,angxyz)
            do ipl = 1,npldafs
              if( Full_self_abs .and.
     &              ( mod(ipl,4) == 2 .or. mod(ipl,4) == 3 ) ) cycle
              do ie = 1,nenerg
                fac = natomsym * 100 / ( Volume_maille * bohr**3
     &                         * Conv_mbarn_nelec(Ephoton(ie)) * pi )
                mu(ie,:,ipl,:) = mu(ie,:,ipl,:) + fac * dampl(ie)
              end do
            end do
! La partie reelle est l'absorption
            mu(:,:,:,:) = img * Conjg( mu(:,:,:,:) )
          endif

          if( ncal > 1 ) then
            iscr = 100 + jfich
            rewind( iscr )
            do ie = 1,nenerg
              write(iscr,*) Xanes(ie,1:nxan)
            end do
            if( Dafs ) then
              do ie = 1,nenerg
                write(iscr,*) ( ( Adafs(ie,ip,ipl), ip = 1,nphi(ipl) ),
     &                            ipl = 1,npldafs )
              end do
            endif
          endif

          deallocate( Ad )
          deallocate( bb )
          deallocate( dampl )
          deallocate( Ephoton )
          deallocate( e1 )
          deallocate( e2 )
          deallocate( lori )
          deallocate( lorix )
          if( Dafs ) deallocate( lorr )
          deallocate( Xa )
          if( seah .or. Arc ) then
            deallocate( Elor )
            deallocate( betalor )
          endif

 1040     if( Conv_done .or. Skip_run(ifich)  ) then
            iscr = 100 + jfich
            rewind( iscr )
            do ie = 1,nenerg
              read(iscr,*) Xanes(ie,1:nxan)
            end do
            if( Dafs ) then
              do ie = 1,nenerg
                read(iscr,*) ( ( Adafs(ie,ip,ipl), ip = 1,nphi(ipl) ),
     &                            ipl = 1,npldafs )
              end do
            endif
          endif

          do ie = 1,nes
            i = indf(ie)
            if( Dafs) As(ie,:,:) = As(ie,:,:) + p2f(ie) * Adafs(i-1,:,:)
     &                                        + p1f(ie) * Adafs(i,:,:)
            if( Cor_abs ) mus(ie,:,:,:) = mus(ie,:,:,:)
     &                                  + p2f(ie) * mu(i-1,:,:,:)
     &                                  + p1f(ie) * mu(i,:,:,:)
            Xs(ie,:) = Xs(ie,:) + p2f(ie) * Xanes(i-1,:)
     &                          + p1f(ie) * Xanes(i,:)
          end do

          deallocate( Adafs )
          if( Cor_abs ) deallocate( mu )
          deallocate( Xanes )
          deallocate( Energ )

        end do ! fin boucle initl

      end do ! fin boucle fichier

      deallocate( run_done )
      if( .not. ( seah .or. Arc ) ) then
        deallocate( Elor )
        deallocate( betalor )
      endif
      deallocate( Skip_run )

      if( .not. Forbidden .and. Dafs .and. .not. Tenseur ) then
        do ipl = 1,npldafs
          do ip = 1,nphi(ipl)
            if( nphi(ipl) == 1 ) then
              As(:,ip,ipl) = S0_2 * As(:,ip,ipl) + f0(ipl)
            else
              As(:,ip,ipl) = S0_2 * As(:,ip,ipl) + f0scan(ip,ipl)
            endif
          end do 
        end do 
      endif

! On calcule en fait epsilon - 1 = -4*pi*Densita_atom*r0/k2
! On est en unites atomiques r0 --> r0 / a0 = alfa^2
! k = 1/2 alfa E
! ct = -4*pi*Densite_atom * alfa^2  * 4 / ( alfa^2 * E^2 )
!    = -16*pi * Densite_atom / E^2 
      if( Bormann ) then
        do ie = 1,nes            
          Eph = Es(ie) + Esmin
          ct_epsilon = - 16 * pi * Densite_atom / Eph**2 
          As(ie,:,:) = ct_epsilon * As(ie,:,:)
        end do
      endif
      
      if( Cor_abs ) then
! mus est deja en micrometre -1
	  do ipl = 1,npldafs
          if( Full_self_abs .and. 
     &            ( mod(ipl,4) == 2 .or. mod(ipl,4) == 3 ) ) cycle
          mus(:,:,ipl,:) = mus(:,:,ipl,:) + fpp_avantseuil
	  end do
      endif

      if( abs( S0_2 - 1._db ) > eps10 ) Xs(:,:) = S0_2 * Xs(:,:)
 
! Convolution par une gaussienne
      if( abs(deltar) > 1.e-10_db .or. abs(vibration) > 1.e-10_db ) then

        allocate( Yr(nes) )
        allocate( Yi(nes) )

        do ipl = 1,nxan
          Yr(:) = Xs(:,ipl)  
          call gaussi(deltar,Efermi,Es,nes,vibration,Yr)
          Xs(:,ipl) = Yr(:)  
        end do

        do ipl = 1,npldafs
          do ip = 1,nphi(ipl)
            Yr(:) = real( As(:,ip,ipl), db )  
            call gaussi(deltar,Efermi,Es,nes,vibration,Yr)
            Yi(:) = aimag( As(:,ip,ipl) )  
            call gaussi(deltar,Efermi,Es,nes,vibration,Yi)
            As(:,ip,ipl) = cmplx( Yr(:), Yi(:),db )  

            if( Cor_abs ) then
              do i = 1,2
                Yr(:) = real( mus(:,ip,ipl,i), db )  
                call gaussi(deltar,Efermi,Es,nes,vibration,Yr)
                Yi(:) = aimag( mus(:,ip,ipl,i) )  
                call gaussi(deltar,Efermi,Es,nes,vibration,Yi)
                mus(:,ip,ipl,i) = cmplx( Yr(:), Yi(:),db )  
	        end do
            endif

          end do
        end do

        deallocate( Yr )
        deallocate( Yi )

      endif

      if( Cor_abs ) then

        allocate( Icor(nes,nphim,npldafs) )
        allocate( Icirc(nes,nphim,npldafs) )
        allocate( Icirccor(nes,nphim,npldafs) )
        allocate( Idcor(nes,nphim,npldafs) )
        allocate( Icircdcor(nes,nphim,npldafs) )

        call Corr_abs(angxyz,As,axyz,Circular,d_dead,Double_cor,
     &      Eseuil(1),fpp_avantseuil,hkl_dafs,hkl_S,Icor,Icirc,Icirccor,
     &      Icircdcor,Idcor,mus,nes,nphi,nphim,npldafs,Self_abs)

      endif

!---- Ecriture -------------------------------------------------

      if( Tenseur ) then
        n_col = 2*npldafs
      elseif( Bormann ) then
        n_col = nxan + 2*npldafs
      else
        n_col = nxan + npldafs
        if( Dafs_bio ) then
          if( icheck > 1 ) then
             n_col = n_col + 4*npldafs
             nw = 5
          else
             nw = 1
          endif
        else  
          if( fprim ) n_col = n_col + 2*npldafs
          if( Self_abs ) n_col = n_col + 3*npldafs
          if( Full_self_abs ) n_col = n_col + 5*npldafs
          if( Double_cor ) n_col = n_col + npldafs
          if( Circular ) n_col = n_col + 2*npldafs
          if( Circular .and. Double_cor ) n_col = n_col + npldafs
        endif
      endif

      allocate( nom_col(n_col) )

      call col_name(Bormann,Circular,Cor_abs,Dafs_bio,Double_cor,
     &    fichin(1),
     &    fprim,Full_self_abs,hkl_dafs,Length_line,Length_word,n_col,
     &    nom_col,npldafs,npldafs_b,nxan,Self_abs,Tenseur)

      allocate( Tens(n_col) )

      if( fprim ) then
        do ipl = 1,npldafs
          if( abs( dpht(ipl) ) > 1.e-10_db ) cycle 
          dpht(ipl) = (1._db,0._db)
        end do
      endif

      if( Energphot ) Es(:) = Es(:) + Esmin 

      do ie = 1,nes

        do ipl = 1,nxan
          Tens(ipl) = Xs(ie,ipl) 
        end do

        jpl = nxan

        if( Dafs_bio ) then

          do ipl = 1,npldafs
            Ip_sig = real(As(ie,1,ipl),db)**2 + aimag(As(ie,1,ipl))**2
            Ip_pi = real(As(ie,2,ipl),db)**2 + aimag(As(ie,2,ipl))**2
            Im_sig = real(As(ie,3,ipl),db)**2 + aimag(As(ie,3,ipl))**2
            Im_pi = real(As(ie,4,ipl),db)**2 + aimag(As(ie,4,ipl))**2
            do i = 1,nw
              jpl = jpl + 1
              select case(i)
                case(1)
                  Tens(jpl) = sqrt( Ip_sig + Ip_pi )
     &                       - sqrt( Im_sig + Im_pi )
                case(2)
                  Tens(jpl) = Ip_sig
                case(3)
                  Tens(jpl) = Ip_pi
                case(4)
                  Tens(jpl) = Im_sig
                case(5)
                  Tens(jpl) = Im_pi
              end select
            end do
          end do

        else

          do ipl = 1,npldafs
            if( .not. ( Tenseur .or. Bormann ) ) then
              jpl = jpl + 1
              Tens(jpl) = abs( As(ie,1,ipl) )**2
            endif
            if( fprim .or. Tenseur .or. Bormann ) then
              jpl = jpl + 1
              Tens(jpl) = real( As(ie,1,ipl), db )
              jpl = jpl + 1
              Tens(jpl) = aimag( As(ie,1,ipl) )
              if( Tenseur .or. Bormann ) cycle 
            endif
            if( Cor_abs ) then
              jpl = jpl + 1
              Tens(jpl) = Icor(ie,1,ipl)
            endif
            if( Double_cor ) then
              jpl = jpl + 1
              Tens(jpl) = Idcor(ie,1,ipl)
            endif
            if( Cor_abs ) then
              do i = 1,2
                jpl = jpl + 1
                Tens(jpl) = real( mus(ie,1,ipl,i), db )
                if( Full_self_abs ) then
                  jpl = jpl + 1
                  Tens(jpl) = aimag( mus(ie,1,ipl,i) )
                endif
              end do 
            endif
            if( Circular .and. mod(ipl,4) == 0 ) then
              do i = 1,4
                jpl = jpl + 1
                Tens(jpl) = Icirc(ie,1,ipl-4+i)
                jpl = jpl + 1
                Tens(jpl) = Icirccor(ie,1,ipl-4+i)
                if( Double_cor ) then
                  jpl = jpl + 1
                  Tens(jpl) = Icircdcor(ie,1,ipl-4+i)
                endif
              end do
            endif
          end do

        endif

        if( Fit_cal ) then
! Commande l'ecriture dans un fichier temporaire
          n = -iscratchconv  
        else
          n = 0
        endif
        call write_out(rdum,rdum,Densite_atom,0._db,Efermi,Es(ie),
     &                 0._db,0._db,.false.,idum,ie,Length_word,
     &                 0,n_col,jpl,0,1,convolution_out,nom_col,
     &                 1,0,0,0,0,n,cdum,
     &                 cdum,Tens,v0muf,.false.,0)

      end do

      deallocate( Tens )

      if( Scan_true .and. .not. Dafs_bio ) then
 
        Open(7, file = fichscanout)

! On enleve le '_0'
        do i = 1,n_col
          nomab = adjustl( nom_col(i) )
          if( nomab(1:1) /= 'I' ) cycle
          l = len_trim(nomab)
          nomab(l-1:l) = '  ' 
          nom_col(i) = nomab
        end do

        do ie = 1,nes
          Write(7,*)
          Write(7,200) Es(ie) * rydb 
          jpl = nxan

          do ipl = 1,npldafs

            jpl = jpl + 1
            Write(7,*) nom_col(jpl)
            do ip = 1,nphi(ipl)
              Write(7,210) angle(ip),
     &              real(As(ie,ip,ipl),db)**2 + aimag(As(ie,ip,ipl))**2
            end do

            if( fprim ) jpl = jpl + 2

            if( Cor_abs ) then
              jpl = jpl + 1
              Write(7,*) nom_col(jpl)
              do ip = 1,nphi(ipl)
                Write(7,210) angle(ip), Icor(ie,ip,ipl)
              end do
              if( Double_cor  ) then
                jpl = jpl + 1
                Write(7,*) nom_col(jpl)
                do ip = 1,nphi(ipl)
                  Write(7,210) angle(ip), Idcor(ie,ip,ipl)
                end do
	        endif
              jpl = jpl + 2  ! pour sauter les mu
              if( Full_self_abs ) jpl = jpl + 2  ! pour sauter les mu
	      endif

            if( Circular .and. mod(ipl,4) == 0 ) then
              do i = 1,4
                jpl = jpl + 1
                Write(7,*) nom_col(jpl)
                do ip = 1,nphi(ipl)
                  Write(7,210) angle(ip), Icirc(ie,ip,ipl-4+i)
                end do
                jpl = jpl + 1
                Write(7,*) nom_col(jpl)
                do ip = 1,nphi(ipl)
                  Write(7,210) angle(ip), Icirccor(ie,ip,ipl-4+i)
                end do
                if( Double_cor  ) then
                  jpl = jpl + 1
                  Write(7,*) nom_col(jpl)
                  do ip = 1,nphi(ipl)
                    Write(7,210) angle(ip), Icircdcor(ie,ip,ipl-4+i)
                  end do
	          endif
	        end do
            endif

          end do
        end do

      endif

      if( Scan_true .and. .not. Dafs_bio ) Close(7)

      if( Dafs) deallocate( As )
      deallocate( decal_initl )
      if( Cor_abs ) deallocate( Icirc, Icirccor, Icircdcor, Icor, Idcor)
      deallocate( indf )
      deallocate( En_fermi ) 
      deallocate( Eph1 ) 
      deallocate( Epsii ) 
      deallocate( Eseuil ) 
      deallocate( fichin )
      deallocate( Es ) 
      deallocate( Efermip )
      deallocate( fichscanin )
	if( Cor_abs ) deallocate( mus )
      deallocate( ne ) 
      deallocate( ninitl )
      deallocate( nom_col )
      deallocate( nsup )
      deallocate( ne_initl )
      deallocate( Pds )
      deallocate( p1f )
      deallocate( p2f )
      deallocate( Xs )
      if( Dafs ) then
        deallocate( dph )
        deallocate( dpht )
        deallocate( fr ); deallocate( fi ) 
        deallocate( f0 )
        deallocate( nphi )
        if( Cor_abs .or. Dafs_bio ) deallocate( hkl_dafs ) 
        if( scan_true ) then
          deallocate( angle )   
          deallocate( phdtscan )   
          deallocate( f0scan )   
        endif 
      endif 

      if( ncal > 1 .and. ical == ncal+1 ) then
        do jfich = 1,mfich
          iscr = 100 + jfich
          Close( iscr )
        end do
      endif

      if( Another_one ) goto 1000

      return

 9999 call write_err_form(itape1,keyword)

      return
  110 format(///' Line not understood in the indata file :'//1x,A//,
     &        ' If it is supposed to be a keyword check the spelling.',/
     &        ' Tabulations are forbidden.'/
     &        ' If it is a line containing numbers, check:'/
     &        '   - How many numbers are supposed to be there.'/
     &        '   - Numbers must be separated by spaces.'/
     &      '   - Are there extra characters (points, tabulations...'//)
  115 format(/' ---- Convolution ------',100('-'))
  120 format(///' The output file has the same name than one of the ',
     &          ' input files.',/
     &          ' This last will be overwritten !',//
     &          ' Are you sure you want to continue ? (y/n) :')
  130 format(///'     Input files are different,'/
     &          '  file 1 : nxan =',i2,', npldafs =',i2,/
     &          '  file ',i2,' : nxan =',i2,', npldafs =',i2)
  140 format(//' For the convolution, the number of energy must be',
     &         ' greater than one !'/,
     &         ' It is not the case in the file:'//A//)
  145 format(///' Error under the keyword Par_',a6,'in the indata file:'
     &       //' The wanted file is the number',i3,' !',/
     &         ' There are only',i3,' files in the job !'//)
  147 format(///' Taking into account the energy shifts, there is no',
     &       ' overlap between',/
     &   ' the energy ranges of the different absorbing atoms !',/
     &   ' No summation and convolution are possible !' //)
  150 format('    Gamma_max  =',f7.2,',  Aseah =',f7.2)
  160 format('    Gamma_max  =',f7.2,',  Ecent =',f7.2,', Elarg =',f7.2)
  170 format('    Gamma_hole =',f7.2,', Efermi =',f7.2,' eV')
  172 format('    Gamma_hole =',f7.2,', Efermi + Shift =',f7.2,
     &       ', initl =',i3)
  174 format('    Gamma_hole =',f7.2,', Efermi + Shift =',f7.2,', site',
     &       i3)
  176 format('    Gamma_hole =',f7.2,', Efermi + Shift =',f7.2,', site',
     &       i3,', initl =',i3)
  178 format(/'      E_(eV)    Width_(eV) lambda_(A)')
  180 format(3f12.3)
  185 format(2f12.3,1p,e12.4)
  200 format(f10.3,1p,240e13.5)
  210 format(f7.1,1p,3e13.5)
      end

!***********************************************************************

      subroutine col_name(Bormann,Circular,Cor_abs,Dafs_bio,Double_cor,
     &    fichin,
     &    fprim,Full_self_abs,hkl_dafs,Length_line,Length_word,n_col,
     &    nom_col,npldafs,npldafs_b,nxan,Self_abs,Tenseur)

      use declarations
      implicit none

      integer:: i, ii, ipl, j, k, l, istat, Length_line, Length_word,
     &  n, n_col, n_col_in, nnombre, npldafs, npldafs_b, nxan

      character(len=132):: fichin
      character(len=Length_word):: nomab
      character(len=length_line):: motl
      character(len=Length_word), dimension(n_col):: nom_col
      character(len=Length_word), dimension(:), allocatable:: nom_col_in

      integer, dimension(3,npldafs_b):: hkl_dafs

      logical:: Bormann, Circular, Cor_abs, Dafs_bio, Double_cor, fprim,
     &  Full_self_abs, Self_abs, Tenseur

      open(2, file = fichin, status='old', iostat=istat) 
      if( istat /= 0 ) call write_open_error(fichin,istat,1)
 
      do i = 1,5
        n = nnombre(2,Length_line)
        if( n == 0 ) exit
        read(2,*)
      end do

      read(2,'(A)') motl

      if( Dafs_bio ) then
        n_col_in = 1 + nxan
      else
        n_col_in = 1 + nxan + 2 * npldafs
        if( Self_abs ) n_col_in = n_col_in + 2*npldafs
        if( Full_self_abs ) n_col_in = n_col_in + 4*npldafs
      endif
      
      allocate( nom_col_in(n_col_in) )

      call Extract_word(Length_line,Length_word,motl,nom_col_in,
     &                  n_col_in,n_col_in)

      i = 0
      j = 1   ! Energy

! Boucle sur les xanes. On supprime le suffixe sur l'etat initial.
      do ipl = 1,nxan
        i = i + 1
        j = j + 1
        nomab = nom_col_in(j)
        l = len_trim( nomab )
        if( nomab(l-1:l-1) == '_' ) nomab(l-1:l) = '  '
        if( nomab(l-2:l-1) == '_L' ) nomab(l-2:l) = '   '
        if( nomab(l-2:l-1) == '_M' ) nomab(l-2:l) = '   '
        if( nomab(l-2:l-1) == '_N' ) nomab(l-2:l) = '   '
        if( nomab(l-2:l-1) == '_O' ) nomab(l-2:l) = '   '
        nom_col(i) = nomab
      end do

      if( Dafs_bio ) then

        if( n_col == nxan + npldafs ) then
          n = 1
        else
          n = 5
        endif

        do ipl = 1,npldafs
          do j = 1,n
            i = i + 1
            nomab = '('
            do k = 1,3
              if( j == 4 .or. j == 5 ) then
                call ad_number(-hkl_dafs(k,ipl),nomab,Length_word)
              else
                call ad_number(hkl_dafs(k,ipl),nomab,Length_word)
              endif
              l = len_trim(nomab) + 1
              if( k == 3 ) then
                nomab(l:l) = ')'
              else
                nomab(l:l) = ','
              endif
            end do
            if( j == 2 .or. j == 4 ) then
              nomab(l+1:l+1) = 's'
            elseif( j == 3 .or. j == 5 ) then
              nomab(l+1:l+1) = 'p'
            endif
            call center_word(nomab,Length_word)
            nom_col(i) = nomab
          end do          
        end do

      else

        do ipl = 1,npldafs

          if( Tenseur ) then

            i = i + 1
            j = j + 1
            nomab = nom_col_in(j)
            l = len_trim( nomab )
            nomab(l+1:l+1) = 'p'
            nom_col(i) = nomab

            i = i + 1
            nomab(l+1:l+1) = 's'
            nom_col(i) = nomab

            cycle

          elseif( Bormann ) then

            i = i + 1
            j = j + 2
            nomab = nom_col_in(j)
            nomab(1:1) = 'e'
            nom_col(i) = nomab

            i = i + 1
            l = min( len_trim( nomab ) + 1, Length_word )
	      nomab(l:l) = 'i'
            nom_col(i) = nomab

            cycle

          endif

          i = i + 1
          j = j + 2
          nomab = nom_col_in(j)
          nomab(1:1) = 'I'
          nom_col(i) = nomab

          if( fprim ) then
            i = i + 1
            nomab = nom_col_in(j)
            l = len_trim( nomab )
            do k = 1,l-1
              nomab(k:k) = nomab(k+1:k+1)
            end do
	      nomab(l:l) = 'p'
            nom_col(i) = nomab

            i = i + 1
	      nomab(l:l) = 's'
            nom_col(i) = nomab
          endif

          if( Cor_abs ) then
            i = i + 1
            nomab = nom_col_in(j)
            l = len_trim( nomab )
            do k = min(l+1,Length_word),2,-1
              nomab(k:k) = nomab(k-1:k-1)
            end do
            nomab(1:2) = 'Ic'
            nom_col(i) = nomab
          endif

          if( Double_cor ) then
            i = i + 1
            nomab = nom_col(i-1) 
            nomab(2:2) = 'd'
            nom_col(i) = nomab
          endif

          if( Cor_abs ) then
            do ii = 1,4
	        if( .not. Full_self_abs .and. ii > 2 ) exit
              i = i + 1
              j = j + 1
              nom_col(i) = nom_col_in(j)
            end do
          endif

          if( Circular .and. mod(ipl,4) == 0 ) then
            do ii = 1,4

              i = i + 1
              nomab = nom_col_in(j-4)
              nomab(1:1) = 'I'
              l = len_trim( nomab )
              do k = 2,l-1
                if( nomab(k:k+1) == 'pp' ) exit
              end do
              select case(ii)
                case(1)
                  nomab(k:k+1) = 'rs'
                case(2)
                  nomab(k:k+1) = 'rp'
                case(3)
                  nomab(k:k+1) = 'ls'
                case(4)
                  nomab(k:k+1) = 'lp'
              end select
              nom_col(i) = nomab

              i = i + 1
              do k = min(l+1,Length_word),2,-1
                nomab(k:k) = nomab(k-1:k-1)
              end do
              nomab(1:2) = 'Ic'
              nom_col(i) = nomab

              if( Double_cor ) then
                i = i + 1
                nomab(2:2) = 'd'
                nom_col(i) = nomab
              endif

            end do
          endif

        end do 

      endif

      return
      end

!***********************************************************************

! Vient en gros de M. O. Krause et J. H. Oliver,
! J. Phys. Chem. Ref. Data 8, 329 (1979)
! Seuil L2 modifie 5/11/2010 pour Fe Co Ni

      subroutine tab_width(Eseuil,Gamma_hole,jseuil,nseuil,numat)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      select case(nseuil)

        case(1)   ! Seuil K

          Es = Eseuil * Rydb

          Gamma_hole = 0.13007 + 0.00011544 * Es
     &               + 9.2262E-09 * Es**2 - 3.2344E-14 * Es**3
        case(2)

          select case(jseuil)

            case(1,0)   ! Seuil L1   (0 pour les anciennes versions)
              if( numat <= 20 ) then
                p = ( numat - 1 ) / 19.
                Gamma_hole = (1 - p ) * 0.07 + p * 4.
              elseif( numat <= 35 ) then
                p = ( numat - 20 ) / 15.
                Gamma_hole = (1 - p ) * 4. + p * 7.
              elseif( numat <= 50 ) then
                p = ( numat - 35 ) / 15.
                Gamma_hole = (1 - p ) * 7. + p * 4.
              elseif( numat <= 75 ) then
                p = ( numat - 50 ) / 25.
                Gamma_hole = (1 - p ) * 4. + p * 8.
              else
                p = ( numat - 75 ) / 20.
                Gamma_hole = (1 - p ) * 8. + p * 19.
              endif

            case(2)   ! Seuil L2
              if( numat <= 25 ) then
!                p = ( numat - 1 ) / 25.
!                Gamma_hole = (1 - p ) * 0.001 + p * 1.7
                Gamma_hole = 0.8_db           ! Modif 18 / 4 / 2011
              elseif( numat == 26 ) then      ! Modif 5 / 11 / 2010
                Gamma_hole = 1._db            ! 1.70  (valeurs precedantes)
              elseif( numat == 27 ) then
                Gamma_hole = 0.8_db           ! 1.32
              elseif( numat == 28 ) then
                Gamma_hole = 0.6_db           ! 1.34
              elseif( numat <= 31 ) then
                p = ( numat - 26 ) / 5.
                Gamma_hole = (1 - p ) * 1.7 + p * 0.8
              elseif( numat <= 60 ) then
                p = ( numat - 31 ) / 29.
                Gamma_hole = (1 - p ) * 0.8 + p * 3.5
              elseif( numat <= 80 ) then
                p = ( numat - 60 ) / 20.
                Gamma_hole = (1 - p ) * 3.5 + p * 5.
              else
                p = ( numat - 80 ) / 15.
                Gamma_hole = (1 - p ) * 5. + p * 10.
              endif

! Augmentation pour tenir compte du canal de desexcitation supplementaire
              Gamma_hole = Gamma_hole + 0.5_db

            case(3)   ! Seuil L3
              if( numat <= 25 ) then
!                p = ( numat - 1 ) / 25.
!                Gamma_hole = (1 - p ) * 0.001 + p * 1.7
                Gamma_hole = 0.8_db           ! Modif 18 / 4 / 2011
              elseif( numat == 26 ) then      ! Modif 5 / 11 / 2010
                Gamma_hole = 1._db            ! 1.70  (valeurs precedantes)
              elseif( numat == 27 ) then
                Gamma_hole = 0.8_db           ! 1.32
              elseif( numat == 28 ) then
                Gamma_hole = 0.6_db           ! 1.34
              elseif( numat <= 31 ) then
                p = ( numat - 26 ) / 5.
                Gamma_hole = (1 - p ) * 1.7 + p * 0.8
              elseif( numat <= 60 ) then
                p = ( numat - 31 ) / 29.
                Gamma_hole = (1 - p ) * 0.8 + p * 3.5
              elseif( numat <= 80 ) then
                p = ( numat - 60 ) / 20.
                Gamma_hole = (1 - p ) * 3.5 + p * 5.
              else
                p = ( numat - 80 ) / 15.
                Gamma_hole = (1 - p ) * 5. + p * 10.
              endif

          end select

! Pour les seuils M, N, O, les valeurs ne sont pas optimisees
        case default

          select case(jseuil)

            case(1,2,3)

              Gamma_hole = 0.2_db

            case(4,5,6,7)

              Gamma_hole = 0.1_db

        end select

      end select

      Gamma_hole = Gamma_hole / Rydb

      return
      end

!***********************************************************************

      subroutine seahdench(A,Efermi,Gamma_max,nelor,Elor,betalor)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db):: Elor(nelor), betalor(nelor), lambda

      do ie = 1,nelor
        E = Elor(ie) - Efermi
        Ep = E
        if( E > 0._db .and. Ep > 0._db .and.
     &                      ( A > 0._db .or. Gamma_max > 0._db ) ) then
          lambda = 0._db
          if( A > 0._db ) lambda = lambda + 1 / ( A * sqrt(Ep) )
          if( Gamma_max > 0._db ) lambda = lambda + sqrt(Ep) / Gamma_max
          betalor(ie) = sqrt(E) / lambda
        else
          betalor(ie) = 0.
        endif
      end do

      return
      end

!***********************************************************************

      subroutine gammarc(Ecent,Elarg,Gamma_max,Efermi,nelor,
     &                   Elor,betalor)
     
      use declarations
      implicit none

      integer ie, nelor

      real(kind=db):: E, Ec, Ecent, Efermi, El, Elarg, Gamma_max, p
      real(kind=db), dimension(nelor):: Elor, betalor
 
      Ec = max( Ecent, 1.E-10_db ) 
      El = max( Elarg, 1.E-10_db ) 
      p = ( pi / 3 ) * Gamma_max / El

      do ie = 1,nelor
        E = Elor(ie) - Efermi
        if ( E <= 0._db ) then
          betalor(ie) = 0._db
        else
          betalor(ie) = Gamma_max
     &                * ( 0.5 + atan( p*(E/Ec - (Ec/E)**2)) / pi)
        endif 
      end do

      return
      end

!***********************************************************************

! Convolution par une gaussienne

      subroutine gaussi(deltar,Efermi,Energ,nenerg,vibration,Y)

      use declarations    
      implicit real(kind=db) (a-h,o-z)

! nj : Points au dela de la gamme pour diminuer les effets de bord de la
! convolution  
      parameter(nj=10,n1m=1-nj)
      real(kind=db) de(n1m:nenerg+nj), Ef(n1m:nenerg+nj), Energ(nenerg), 
     &             gaus(n1m:nenerg+nj), Y(nenerg), 
     &             Xa(n1m:nenerg+nj)

      Ef(1:nenerg) = Energ(1:nenerg)
      Xa(1:nenerg) = Y(1:nenerg)

! Creation des points au dela de la gamme 
      ne1 = 1 - nj
!      ne2 = nenerg + nj
      def = Ef(2) - Ef(1)
      do ie = 0,ne1,-1 
        Ef(ie) = Ef(ie+1) - def
        Xa(ie) = Xa(1)
      end do
      ne2 = nenerg

      de(ne1) = Ef(ne1+1) - Ef(ne1) 
      do ie = ne1+1,ne2-1
        de(ie) = 0.5 * ( Ef(ie+1) - Ef(ie-1) )
      end do 
      de(ne2) = ( Ef(ne2) - Ef(ne2-1) ) 

!      fnorm = 1 / ( b * sqrt( 2 * pi ) )

      if( abs(deltar) < 1.e-10_db .and. abs(vibration) < 1.e-10_db )
     &  return

      do ie = 1,nenerg

        vib = 2 * vibration * ( Ef(ie) - Efermi + 0.5_db )
        vib = max( 0._db, vib )
        b = deltar + vib
        if( abs(b) < 1.e-10_db ) then
          Y(ie) = Xa(ie)
          cycle
        endif
             
        gaus(:) = 0._db
        Pdt = 0._db
        do je = ne1,ne2
          n = max( int( 10 * de(je) / b ), 1 )
          pas = de(je) / ( n + 1 )
          if ( je == ne1 ) then
            E = Ef(je) - 0.5 * ( Ef(je+1) -  Ef(je) )
          else
            E = Ef(je) - 0.5 * ( Ef(je) -  Ef(je-1) )
          endif 
          do i = 1,n
            E = E + pas
            if( ( E < Ef(je) .and. je /= ne1 ) .or. je == ne2 ) then
              p = ( E - Ef(je-1) ) / ( Ef(je) - Ef(je-1) )
              Yint = (1 - p ) * Xa(je-1) + p * Xa(je)
            else
              p = ( E - Ef(je) ) / ( Ef(je+1) - Ef(je) )
              Yint = (1 - p ) * Xa(je) + p * Xa(je+1)
            endif
            fac = -0.5 * ( ( E - Ef(ie) ) / b )**2
            if( fac > -600._db ) then
              efac = exp( fac )
              gaus(je) = gaus(je) + efac * Yint
              Pdt = Pdt + efac * de(je) / n
            endif
          end do 
          gaus(je) = ( gaus(je) / n ) * de(je)
        end do
        Y(ie) = sum( gaus(ne1:ne2) ) / Pdt

      end do

      return
      end

!***********************************************************************

! Calcul des coefficients de la lorentzienne

      subroutine cflor(bb,betalor,Efermi,Eph,Elor,Energ,ie1,ie2,nef,
     &                 nelor,nenerg,nenerge,Photo_emission)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      logical:: Photo_emission

      real(kind=db):: bb(nenerge), betalor(nelor), Elor(nelor),  
     &               Energ(nenerg), Eph(nenerge)

      Eph(1:nenerg) = Energ(1:nenerg)

! Creation des points au dela de la gamme 
      def = Eph(2) - Eph(1)
      def = Eph(nenerg) - Eph(nenerg-1)
      do ie = nenerg+1,nenerge 
        Eph(ie) = Eph(ie-1) + def
      end do

! Les etats en dessous de Fermi sont occupes
      do ie = 1,nenerge
        if( Eph(ie) > Efermi - 1.e-10_db ) exit
      end do
      nef = ie

      if( Photo_emission ) nef = max(1, nef-1)

      if( Photo_emission ) then
        ie1 = 1
        ie2 = nef
      else
        ie1 = nef
        ie2 = nenerge
      endif

      do ie = 1,nenerge
        if( Eph(ie) <= Elor(1) ) then
          bb(ie) = betalor(1) 
        elseif( Eph(ie) >= Elor(nelor) ) then
          bb(ie) = betalor(nelor)
        else 
          do j = 2,nelor
            if( Elor(j) >= Eph(ie) ) exit
          end do
          p = ( Eph(ie) - Elor(j-1) ) / ( Elor(j) - Elor(j-1) )
          bb(ie) = p * betalor(j) + ( 1 - p ) * betalor(j-1)
        endif
        bb(ie) = 0.5_db * bb(ie) ! bb est en fait Gamma / 2
      end do

      return
      end

!***********************************************************************

! Calcul de l'extrapolation atomique

      subroutine extrapat(bb_nenerg,dampl,
     &                    Energ,Eseuil,Extrap,fpp0,
     &                    fprime_atom,numat,nenerg)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      parameter(nim=10000)

      character(len=2):: Chemical_symbol
      complex(kind=db):: dampl(nenerg)
      logical:: Extrap, fprime_atom

      real(kind=db):: ei(nim), ei1(nim), ei2(nim), Energ(nenerg),
     &               fppn(nim)


      if( fprime_atom ) then
        ephoton = max( Eseuil - 200 / rydb, 10 / rydb )
        dde = 0.5 /rydb
        ephm = Eseuil + 5000. / rydb
        ephi = Eseuil + 500. / rydb
        write(3,110) Chemical_symbol(numat)
        write(3,'(/A)') '   Energy_eV       f_prim      f_second'
        do i = 1,nim
          if( ephoton < ephi ) then
            ephoton = ephoton + dde
          else
            ephoton = ephoton + 10 * dde
          endif
          if( ephoton > ephm ) exit
          call fprime(numat,ephoton,fpp,fp)
          write(3,'(3f13.6)') ephoton*rydb, fp, fpp
        end do
        if( .not. Extrap ) return
      endif

      ephoton = Eseuil - 2 / rydb
      call fprime(numat,ephoton,fpp0,fp)

      ephoton = Eseuil + Energ(nenerg)
      ephm = Eseuil + 10000._db / rydb
      if( nenerg > 1 ) then
        dde = Energ(nenerg) - Energ(nenerg-1)
      else
        dde = 2 / rydb
      endif

      alfa = 1.02_db
      f = 0.5_db * alfa
      do i = 1,nim
        dde = alfa * dde
        ephoton = ephoton + dde 
        if( ephoton > ephm ) exit
        call fprime(numat,ephoton,fpp,fp)
        ei(i) = ephoton
        ei2(i) = ei(i) + f * dde
        ei1(i) = ei(i) - 0.5_db * dde
        fppn(i) = - ( fpp - fpp0 ) / pi 
      end do
      ni = min(i-1,nim)

      bb2 = bb_nenerg**2

      do ie = 1,nenerg
        e1 = Energ(ie) + Eseuil
        if( abs(bb_nenerg) < 1.e-10_db ) then
          dampl(ie) = sum( fppn(1:ni)
     &    * log( ( ei2(1:ni) - e1 ) / (ei1(1:ni) - e1 ) ) )
        else
          dampl(ie) = sum( fppn(1:ni)
     &          * ( 0.5*log( ( ( ei2(1:ni) - e1 )**2 + bb2 ) 
     &                     / ( ( ei1(1:ni) - e1 )**2 + bb2 ) ) 
     &            - img * ( atan( ( ei2(1:ni) - e1 ) / bb_nenerg )
     &                    - atan( ( ei1(1:ni) - e1 ) / bb_nenerg ) ) ) )  
        endif
      end do

      return
  110 format(/1x,a2,' atomic scattering scattering amplitude')
      end

!***********************************************************************

      subroutine extract_word(length,length_word,mot,word,n_word,n_dim)

      character(len=length):: mot
      character(len=length_word):: mota
      character(len=length_word), dimension(n_dim):: word

      i_word = 0
      i = 0
      do ii = 1,length
        i = i + 1
        if( i > length ) exit
        if( mot(i:i) == ' ' ) cycle
        i_word = i_word + 1
        mota = ' '
        do j = i,length          
          if( mot(j:j) == ' ' ) exit
          mota(j-i+1:j-i+1) = mot(j:j)
        end do
        i = j
        word(i_word) = mota
        if( i_word == n_word ) exit          
      end do

      return
      end

!***********************************************************************

! Calcul de l'intensite corrige d'absorption en tenant compte de la
! birefringence

      subroutine Corr_abs(angxyz,As,axyz,Circular,d_dead,Double_cor,
     &        Eseuil,fpp_avantseuil,hkl_dafs,hkl_S,Icor,Icirc,Icirccor,
     &        Icircdcor,Idcor,mus,nes,nphi,nphim,npldafs,Self_abs)

      use declarations
      implicit none

      integer:: i, ie, inda, indp, ip, ipl, ipl_pas, j, k, nes, nphim,
     &          npldafs
      integer, dimension(3):: hkl, hkl_S
      integer, dimension(npldafs):: nphi
      integer, dimension(3,npldafs):: hkl_dafs

      logical:: Check, Circular, Double_cor, Self_abs

      complex(kind = db):: dmu_i, dmu_s, Expo, fps, fpp, fsp, fss,  
     &     mu_i, mu_pp_i, mu_pp_s, mu_ps_i, mu_ps_s, mu_s,
     &     mu_sp_i, mu_sp_s, mu_ss_i, mu_ss_s, mups_i, mups_s, musp_i,
     &     musp_s, t_i, t_s, tau_i, tau_s, Trace 

      complex(kind=db), dimension(2):: Ex_i, Ex_s
      complex(kind=db), dimension(4):: Ex
      complex(kind=db), dimension(2,2):: F, M1, M2, M3, Mat_pol
      complex(kind=db), dimension(2,2,2):: Tr_i, Tr_s
      complex(kind=db), dimension(2,2,4):: ATFT, TFT
      complex(kind=db), dimension(nes,nphim,npldafs):: As
      complex(kind=db), dimension(nes,nphim,npldafs,2):: mus

      real(kind=db):: alfa_S, Attenuation, corr_ref, Correction, csi,
     &       css, d_dead, Eseuil, eta, fpp_avantseuil, mum, mum_0, 
     &       P1, P2, P3, Resultat, theta_B

      real(kind=db), dimension(3):: angxyz, axyz
      real(kind=db), dimension(2,2):: Analyseur
      real(kind=db), dimension(nes,nphim,npldafs):: Icirc, Icirccor,
     &                                       Icircdcor, Icor, Idcor

      Check = .false.

      if( Self_abs ) then
        ipl_pas = 1
      else
        ipl_pas = 4
      endif

      do ipl = 1,npldafs,ipl_pas

        hkl(:) = hkl_dafs(:,ipl)
        call Angle_Bragg(alfa_S,angxyz,axyz,Eseuil,hkl,hkl_S,theta_B)

        csi = 1 / sin( theta_B + alfa_S )
        css = 1 / sin( theta_B - alfa_S )
        mum_0 = fpp_avantseuil * ( csi + css )
        corr_ref = mum_0 / exp( - mum_0 * d_dead )

        do ip = 1,nphi(ipl)

          do ie = 1,nes

            if( Self_abs .or. Double_cor ) mum = 
     &                            real( mus(ie,ip,ipl,1), db )*csi
     &                          + real( mus(ie,ip,ipl,2), db )*css
            if( Double_cor ) Correction = mum_0 / mum

            if( self_abs ) then
              Attenuation = exp( - mum * d_dead ) * corr_ref
              Icor(ie,ip,ipl) = ( Attenuation / mum ) 
     &                        * abs( As(ie,ip,ipl) )**2 
              if( Double_cor ) Idcor(ie,ip,ipl) = Icor(ie,ip,ipl)
     &                                          * Correction 
              cycle
            endif

            mu_ss_i = mus(ie,ip,ipl,1)
            mu_ss_s = mus(ie,ip,ipl,2)
            mu_sp_i = mus(ie,ip,ipl+1,1)
            mu_sp_s = mus(ie,ip,ipl+1,2)
            mu_ps_i = mus(ie,ip,ipl+2,1)
            mu_ps_s = mus(ie,ip,ipl+2,2)
            mu_pp_i = mus(ie,ip,ipl+3,1)
            mu_pp_s = mus(ie,ip,ipl+3,2)

            fss = As(ie,ip,ipl)
            fsp = As(ie,ip,ipl+1)
            fps = As(ie,ip,ipl+2)
            fpp = As(ie,ip,ipl+3)

            mu_i = ( mu_ss_i + mu_pp_i ) * csi / 4
            mu_s = ( mu_ss_s + mu_pp_s ) * css / 4

            mum = 0.5_db*real( mu_ss_i + mu_pp_i, db )*csi
     &          + 0.5_db*real( mu_ss_s + mu_pp_s, db )*css

            Attenuation = exp( - mum * d_dead ) * corr_ref

            tau_i = sqrt( ( mu_pp_i - mu_ss_i )**2 + 4*mu_sp_i*mu_ps_i )
            tau_s = sqrt( ( mu_pp_s - mu_ss_s )**2 + 4*mu_sp_s*mu_ps_s )

            t_i = csi * tau_i / 4
            t_s = css * tau_s / 4

            if( abs(tau_i) > eps10 ) then
              mups_i = 2*mu_ps_i / tau_i 
              musp_i = 2*mu_sp_i / tau_i 
              dmu_i = ( mu_pp_i - mu_ss_i ) / tau_i
            else
              mups_i = (0._db, 0._db) 
              musp_i = (0._db, 0._db) 
              dmu_i = (0._db, 0._db)
            endif 

            if( abs(tau_s) > eps10 ) then
              mups_s = 2*mu_ps_s / tau_s
              musp_s = 2*mu_sp_s / tau_s
              dmu_s = ( mu_pp_s - mu_ss_s ) / tau_s
            else
              mups_s = (0._db, 0._db)
              musp_s = (0._db, 0._db)
              dmu_s = (0._db, 0._db)
            endif 

! Matrices transmittance

            Tr_i(1,1,1) = 1 + dmu_i; Tr_i(1,2,1) = - musp_i 
            Tr_i(2,1,1) = - mups_i;  Tr_i(2,2,1) = 1 - dmu_i 

            Tr_i(1,1,2) = 1 - dmu_i; Tr_i(1,2,2) = musp_i 
            Tr_i(2,1,2) = mups_i;    Tr_i(2,2,2) = 1 + dmu_i 

            Tr_i(:,:,:) = Tr_i(:,:,:) / 2 

            Tr_s(1,1,1) = 1 + dmu_s; Tr_s(1,2,1) = - musp_s 
            Tr_s(2,1,1) = - mups_s;  Tr_s(2,2,1) = 1 - dmu_s 

            Tr_s(1,1,2) = 1 - dmu_s; Tr_s(1,2,2) = musp_s 
            Tr_s(2,1,2) = mups_s;    Tr_s(2,2,2) = 1 + dmu_s 

            Tr_s(:,:,:) = Tr_s(:,:,:) / 2 

! Exposant des exponentielles

            Ex_i(1) = - mu_i + t_i;  Ex_i(2) = - mu_i - t_i 

            Ex_s(1) = - mu_s + t_s;  Ex_s(2) = - mu_s - t_s 

! Amplitude de diffusion

            F(1,1) = fss; F(1,2) = fps
            F(2,1) = fsp; F(2,2) = fpp

! Produit
            k = 0
            do i = 1,2
              do j = 1,2
                k = k + 1
                M1(:,:) = Tr_s(:,:,i)
                M2(:,:) = Tr_i(:,:,j)
                M3 = Matmul( M1, Matmul( F, M2 ) )
                TFT(:,:,k) = M3
                Ex(k) = Ex_s(i) + Ex_i(j) 
              end do
            end do

            if( check ) then
              Write(6,'(A)') ' tau_i tau_s'
	        Write(6,'(1p,8e13.5)') tau_i, tau_s, t_i, t_s 
              Write(6,'(A)') ' Tr_i'
              do i = 1,2
	          Write(6,'(1p,8e13.5)') Tr_i(i,:,1), Tr_i(i,:,2) 
              end do
              Write(6,'(A)') ' Tr_s'
              do i = 1,2
	          Write(6,'(1p,8e13.5)') Tr_s(i,:,1), Tr_s(i,:,2) 
              end do
              Write(6,'(A)') ' Ex_i'
	          Write(6,'(1p,8e13.5)') Ex_i(:) 
              Write(6,'(A)') ' Ex_s'
	          Write(6,'(1p,8e13.5)') Ex_s(:) 
              Write(6,'(A)') ' F'
              do i = 1,2
	          Write(6,'(1p,8e13.5)') F(i,:)
              end do
              do k = 1,4
	          Write(6,'(a9,i2)') ' TFT, k =', k
	          Write(6,'(1p,8e13.5)') TFT(1,:,k), Ex(k)
	          Write(6,'(1p,8e13.5)') TFT(2,:,k)
              end do
            endif

            do inda = 1,2   ! Boucle sur l'analyseur

              if( inda == 1 ) then
                eta = 0._db   ! Polarisation analyseur sigma
              else
                eta = pi / 2  ! Polarisation analyseur pi
              endif
              Analyseur(1,1) = Cos( eta ); Analyseur(1,2) = - Sin( eta ) 
              Analyseur(2,1) = 0._db;      Analyseur(2,2) = 0._db 

              if( Check ) then
	          Write(6,'(/A)') ' Analyseur'
	          Write(6,'(2f8.5)') Analyseur(1,:)
	          Write(6,'(2f8.5)') Analyseur(2,:)
              endif

              do i = 1,4
                M1(:,:) = TFT(:,:,i)
                M2 = Matmul( Analyseur, M1 )
                ATFT(:,:,i) = M2
              end do

              do indp = 1,4   ! Boucle sur la polarisation

                if( .not. circular .and. indp > 2 ) exit

                select case(indp)
                  case(1)  ! Polarisation sigma
                    P1 = 1._db; P2 = 0._db; P3 = 0._db
                  case(2)  ! Polarisation pi
                    P1 = -1._db; P2 = 0._db; P3 = 0._db
                  case(3)  ! Polarisation circulaire droite
                    P1 = 0._db; P2 = 0._db; P3 = 1._db
                  case(4)  ! Polarisation circulaire gauche
                    P1 = 0._db; P2 = 0._db; P3 = -1._db
                end select

                Mat_pol(1,1) = 1._db + P1;  Mat_pol(1,2) = P2 - img*P3
                Mat_pol(2,1) = P2 + img*P3; Mat_pol(2,2) = 1._db - P1

                Mat_pol = Mat_pol / 2

                Resultat = 0._db
                do i = 1,4
                  do j = 1,4
                    M1(:,:) = ATFT(:,:,i)
                    M2(:,:) = ATFT(:,:,j)
                    M2 = Transpose( Conjg( M2 ) )
                    M3 = Matmul( M1, Matmul( Mat_pol, M2 ) )
                    Trace = M3(1,1) + M3(2,2)
                    Expo = Ex(i) + Conjg( Ex(j) )
! Signe moins car integrale de 0 a l'infini = - Primitive(0)
                    Resultat = Resultat - Real ( Trace / Expo, db )   
                  end do
                end do
                Resultat = Resultat * Attenuation

                if( check ) then
	            Write(6,'(/A)') ' Polar'
	            Write(6,'(4f8.3)') Mat_pol(1,:)
	            Write(6,'(4f8.3)') Mat_pol(2,:)
                  Write(6,'(/A)') ' Attenuation, Resultat'
	            Write(6,'(1p,5e13.5)') Attenuation, Resultat
                endif
                       
                if( inda == 1 ) then

                  select case(indp)
                    case(1)  ! sigma - sigma
                      Icor(ie,ip,ipl) = Resultat
                      if( Double_cor ) Idcor(ie,ip,ipl) = Resultat
     &                     * Corr_ref
     &                     / Real( mu_ss_i*csi + mu_ss_s*css, db) 
                    case(2)  ! pi - sigma
                      Icor(ie,ip,ipl+2) = Resultat
                      if( Double_cor ) Idcor(ie,ip,ipl+2) = Resultat
     &                     * Corr_ref
     &                     / Real( mu_pp_i*csi + mu_ss_s*css, db) 
                    case(3)  ! droite - sigma
                      Icirc(ie,ip,ipl) = aimag( Conjg(fps)*fss )
     &                  + Real( fps*Conjg(fps) + fss*Conjg(fss), db )/2
                      Icirccor(ie,ip,ipl) = Resultat
                      if( Double_cor ) Icircdcor(ie,ip,ipl) = Resultat
     &                    * Corr_ref
     &                    / Real( 0.5_db*( mu_ss_i + mu_pp_i )*csi
     &                                               + mu_ss_s*css, db) 
                    case(4)  ! gauche - sigma
                      Icirc(ie,ip,ipl+2) = - aimag( Conjg(fps)*fss )
     &                  + Real( fps*Conjg(fps) + fss*Conjg(fss), db )/2
                      Icirccor(ie,ip,ipl+2) = Resultat
                      if( Double_cor ) Icircdcor(ie,ip,ipl+2) = Resultat
     &                    * Corr_ref
     &                    / Real( 0.5_db*( mu_ss_i + mu_pp_i )*csi
     &                                               + mu_ss_s*css, db) 
                  end select

                else

                  select case(indp)
                    case(1)  ! sigma - pi
                      Icor(ie,ip,ipl+1) = Resultat
                      if( Double_cor ) Idcor(ie,ip,ipl+1) = Resultat
     &                    * Corr_ref
     &                    / Real( mu_ss_i*csi + mu_pp_s*css, db) 
                    case(2)  ! pi - pi
                      Icor(ie,ip,ipl+3) = Resultat
                      if( Double_cor ) Idcor(ie,ip,ipl+3) = Resultat
     &                    * Corr_ref
     &                    / Real( mu_pp_i*csi + mu_pp_s*css, db) 
                    case(3)  ! droite - pi
                      Icirc(ie,ip,ipl+1) = aimag( Conjg(fpp)*fsp )
     &                  + Real( fpp*Conjg(fpp) + fsp*Conjg(fsp), db )/2
                      Icirccor(ie,ip,ipl+1) = Resultat
                      if( Double_cor ) Icircdcor(ie,ip,ipl+1) = Resultat
     &                    * Corr_ref
     &                    / Real( 0.5_db*( mu_ss_i + mu_pp_i )*csi
     &                                               + mu_pp_s*css, db) 
                    case(4)  ! gauche - pi
                      Icirc(ie,ip,ipl+3) = - aimag( Conjg(fpp)*fsp )
     &                  + Real( fpp*Conjg(fpp) + fsp*Conjg(fsp), db )/2
                      Icirccor(ie,ip,ipl+3) = Resultat
                      if( Double_cor ) Icircdcor(ie,ip,ipl+3) = Resultat
     &                    * Corr_ref
     &                    / Real( 0.5_db*( mu_ss_i + mu_pp_i )*csi
     &                                               + mu_pp_s*css, db) 
                  end select
                endif

              end do ! Boucle Polarisation
            end do ! Boucle analyseur

          end do
        end do
      end do

      return
      end

!***********************************************************************

      subroutine Angle_Bragg(alfa_S,angxyz,axyz,Eseuil,hkl,hkl_S,
     &                       theta_B)

      use declarations
      implicit none

      integer:: i, j
      integer, dimension(3):: hkl, hkl_S

      real(kind=db):: alfa_S, detmat, dhkl, Eseuil, fac, konde, qkn,
     &                qkn_S, theta_B, vol
      real(kind=db), dimension(3):: angxyz, axyz, cosdir, hklred, 
     &                              qk, qk_S, vx, vy, vz, wx, wy, wz
      real(kind=db), dimension(0:3):: det 
      real(kind=db), dimension(3,3):: Cubmat, mat

      character(len=5):: Struct

      call cal_cubmat(angxyz,Cubmat,struct)

      konde = 0.5_db * alfa_sf * Eseuil

      cosdir(:) = cos( angxyz(:) * pi / 180 )

! Base du reseau directe, exprimee dans une base orthonormee
      vx(:) = Cubmat(:,1) * axyz(:)
      vy(:) = Cubmat(:,2) * axyz(:)
      vz(:) = Cubmat(:,3) * axyz(:)

! wx, wy, wz : base du reseau reciproque exprimee dans la meme base
      call prodvec(wx,vy,vz)

      vol = sum( wx(:) * vx(:) )
      wx(:) = wx(:) / vol
      call prodvec(wy,vz,vx)
      wy(:) = wy(:) / vol
      call prodvec(wz,vx,vy)
      wz(:) = wz(:) / vol

      qk_S(:) = hkl_S(1) * wx(:) + hkl_S(2) * wy(:) + hkl_S(3) * wz(:)
      qkn_S = sqrt( sum( qk_S(:)**2 ) )

      if( abs(qkn_S) > eps10 ) qk_S(:) = qk_S(:) / qkn

      if( hkl(1) == 0 .and. hkl(2) == 0 .and. hkl(3) == 0 ) then

        dhkl = 0._db
        theta_B = 0._db

      else
 
        hklred(:) = hkl(:) / axyz(:)
        do i = 0,3
          do j = 1,3
            mat(j,j) = 1._db
          end do
          mat(1,2) = cosdir(3)
          mat(1,3) = cosdir(2)
          mat(2,3) = cosdir(1)
          mat(2,1) = mat(1,2)
          mat(3,1) = mat(1,3)
          mat(3,2) = mat(2,3)
          if( i > 0 ) mat(:,i) = hklred(:)
          det(i) = detmat(mat)
        end do
! Distance interplan
        dhkl = sqrt( det(0) / sum( hklred(1:3) * det(1:3) ) )
        fac = abs( pi / ( konde * dhkl ) )

        if( fac > 1._db ) then
          call write_error
          write(6,120) hkl(:)
          stop
        endif

        theta_B = asin( pi / ( konde * dhkl ) )

      endif

! Calcul angle avec plan de surface

      if( qkn_S > eps10 ) then

        qk(:) = hkl(1) * wx(:) + hkl(2) * wy(:) + hkl(3) * wz(:)
        qkn = sqrt( sum( qk(:)**2 ) )
        if( abs(qkn) > eps10 ) then
          qk(:) = qk(:) / qkn
          alfa_S = acos( sum( qk(:)*qk_S(:) ) ) 
        else
          alfa_S = 0._db
        endif
      else
        alfa_S = 0._db
      endif

      return
  120 format(//' The reflection number',i3,' : (h,k,l) = (',3i3,')',
     &           ' does not exist at this energy !'//)
      end

