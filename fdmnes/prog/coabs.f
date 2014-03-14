! FDMNES subroutine

! Calculate the absorption cross sections and the RXS amplitudes

      subroutine write_coabs(Allsite,angxyz,axyz,Base_spin,
     &          Core_resolved,Dafs,Dafs_bio,Densite_atom,E_cut,
     &          E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Energ,Energphot,
     &          Extract,Epsii,Eseuil,Final_tddft,fpp_avantseuil,
     &          Full_self_abs,Green_int,Green_plus,hkl_dafs,
     &          iabsorig,icheck,ie,ie_computer,
     &          Int_tens,isigpi,isymeq,jseuil,length_word,ltypcal,M1M1,
     &          Moyenne,mpinodee,n_multi_run,n_oo,natomsym,nbseuil,
     &          ncolm,ncolr,ncolt,nenerg,ninit1,ninitlr,nomabs,           
     &          nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs,
     &          nphim,nplr,
     &          nplrm,nseuil,nspin,numat_abs,nxanout,pdp,phdafs,
     &          phdf0t,phdt,pol,poldafse,poldafss,Rot_int,sec_atom,
     &          secdd_a,secdd_m_a,secdq_a,secdq_m_a,secdo_a,secdo_m_a,
     &          secmd_a,secmd_m_a,secmm_a,secmm_m_a,secoo_a,secoo_m_a,
     &          secqq_a,secqq_m_a,self_abs,Spinorbite_p,Taux_eq,V0muf,
     &          vecdafse,vecdafss,vec,Volume_maille,Xan_atom)

      use declarations
      implicit none

      integer:: he, hs, i, ia, iabsorig, ib, ic1, ic2, icheck, icn1,
     &  icn2, id, ie,
     &  ie_computer, ig, ii, ind_mu, initlr, ip, ipl, ipldafs, iseuil,
     &  isym, ixandafs, j, j1, je, jhe, jhs, jpl, js, jseuil, ke, ks, l,
     &  length_word, ll, long, mpinodee, n_dim, n_tens, n_tens_dd, 
     &  n_tens_qq, n_tens_t, n_tens_max, n_multi_run, n_oo, n1, n2, na,
     &  natomsym, nb, nbseuil, nc, nccm, ncolm, n_tens_dq, ncolr,  
     &  ncolt, nenerg, ninit1, ninitlr, nl, np, npldafs, 
     &  nphim, nplt, nplr, nplrm, nseuil, nspin, numat_abs, nxanout, nw 

      parameter( n_tens_dd=9, n_tens_dq=15, n_tens_qq=25,
     &           n_tens_t = n_tens_dd + n_tens_dq + n_tens_qq,
     &           n_tens_max = 8 + 2 * n_tens_t + 2 * n_tens_dq ) 

      character(len=length_word):: nomab
      character(len=length_word), dimension(ncolm):: nomabs
      character(len=length_word), dimension(ncolm*ninitlr):: title
      character(len=13), dimension(nplrm):: ltypcal
      character(len=132) nomfich, nomfich_cal_convt, nomfich_s,
     &                   nomfichdafst, nomficht
      character(len=2310) mot

      complex(kind=db):: cf, ph, ph_m, sec
      complex(kind=db), dimension(3):: plae, plas, uae, uas
      complex(kind=db), dimension(8*ninitlr):: compnum
      complex(kind=db), dimension(3,3):: secdd, secmd, secmm
      complex(kind=db), dimension(3,3,3):: secdq
      complex(kind=db), dimension(3,3,3,3):: secdo, secqq
      complex(kind=db), dimension(3,3,3,3,3,3):: Mat6
      complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia,
     &         secddia_m, secmdia, secmdia_m, secmmia, secmmia_m
      complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia,
     &                                                    secdqia_m
      complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secdoia,
     &                                   secdoia_m, secqqia, secqqia_m
      complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:natomsym)::
     &                                              secooia, secooia_m  
      complex(kind=db), dimension(3,3,ninitlr,0:mpinodee-1):: secdd_a,  
     &                secdd_m_a, secmd_a, secmd_m_a, secmm_a, secmm_m_a
      complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodee-1):: secdq_a,
     &                                                       secdq_m_a
      complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodee-1)::
     &                        secdo_a, secdo_m_a, secqq_a, secqq_m_a
      complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodee-1)::
     &                                            secoo_a, secoo_m_a
      complex(kind=db), dimension(npldafs):: phdtem, phdf0t1, phdt1
      complex(kind=db), dimension(3,nplrm):: pol
      complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt 
      complex(kind=db), dimension(natomsym,npldafs):: phdafs 
      complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss 
     
      complex(kind=db), dimension(:,:,:,:), allocatable :: ampldafs,
     &  ampldafsdd, ampldafsdo, ampldafsdq, ampldafsmd,
     &  ampldafsmm, ampldafsoo, ampldafsqq
      complex(kind=db), dimension(:,:,:,:,:), allocatable :: mu,
     &  mudd, mudo, mudq, mumd, mumm, muoo, muqq

      integer, dimension(0):: idum
      integer, dimension(natomsym):: isymeq
      integer, dimension(npldafs):: nphi_dafs
      integer, dimension(3,npldafs):: hkl_dafs
      integer, dimension(npldafs,2):: isigpi

      logical Allsite, Base_spin, Cartesian_tensor, Cor_abs, 
     &  Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, E_vec, Dafs, 
     &  Dafs_bio, Final_tddft, Energphot, Extract, Full_self_abs,
     &  Green_int, Green_int_mag, Green_plus, idafs, M1M1, magn_sens,     
     &  Moyenne, mu_cal, Self_abs,   
     &  Spherical_tensor, Spinorbite_p, Tens_comp, xan_atom
 
      real(kind=db):: ang, c_micro, c_milli, cst, conv_mbarn_nelec,
     &  dang, Densite_atom, E_cut, 
     &  eph2, Ephoton, Ephseuil, fpp_avantseuil, V0muf, Volume_maille
 
      real(kind=db), dimension(0):: rdum
      real(kind=db), dimension(ninitlr) :: ct_nelec, Epsii
      real(kind=db), dimension(nbseuil):: Eseuil
      real(kind=db), dimension(3):: angxyz, axyz, voae, voas
      real(kind=db), dimension(3,3):: matopsym, Rot_int
      real(kind=db), dimension(nenerg) :: Energ
      real(kind=db), dimension(ninitlr) :: sec_atom
      real(kind=db), dimension(3,nplrm) :: vec
      real(kind=db), dimension(nplrm,2) :: pdp
      real(kind=db), dimension(ncolm*ninitlr) :: tens
      real(kind=db), dimension(natomsym) :: Taux_eq
      real(kind=db), dimension(n_tens_max,ninitlr,0:natomsym):: Int_tens
      real(kind=db), dimension(ncolr,ninitlr,0:natomsym):: secabs, 
     &      secabsdd, secabsdq, secabsdo, secabsmd, secabsmm, secabsoo,
     &      secabsqq
      real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss 
     
      common/cartesian/ cartesian_tensor 
      common/spheric/ spherical_tensor 

      if( icheck > 0 ) write(3,110)

      if( E1E1 ) secabsdd(:,:,:) = ( 0._db, 0._db )
      if( E2E2 ) secabsqq(:,:,:) = ( 0._db, 0._db ) 
      if( E1E3 ) secabsdo(:,:,:) = ( 0._db, 0._db )
      if( E3E3 ) secabsoo(:,:,:) = ( 0._db, 0._db )
      if( M1M1 ) secabsmm(:,:,:) = ( 0._db, 0._db )
      if( E1M1 ) secabsmd(:,:,:) = ( 0._db, 0._db )
      if( E1E2 ) secabsdq(:,:,:) = ( 0._db, 0._db )

      if( ( jseuil > 1 .and. nspin == 2 ) .or. Spinorbite_p .or. M1M1
     &   .or. E1M1 ) then
        magn_sens = .true.
      else
        magn_sens = .false.
      endif

      if( Green_int .and. Magn_sens ) then
        Green_int_mag = .true.
      else
        Green_int_mag = .false.
      endif

      Cor_abs = Full_self_abs .or. Self_abs

      if( Allsite ) then
        na = natomsym
        nb = natomsym
      else
        na = 0
        nb = 1
      endif

      if( dafs ) then

        allocate( ampldafs(npldafs,nphim,ninitlr,0:natomsym) )
        if( E1E1 )
     &    allocate( ampldafsdd(npldafs,nphim,ninitlr,0:natomsym) )
        if( E1E2 )
     &    allocate( ampldafsdq(npldafs,nphim,ninitlr,0:natomsym) )
        if( E2E2 )
     &    allocate( ampldafsqq(npldafs,nphim,ninitlr,0:natomsym) )
        if( E1E3 )
     &    allocate( ampldafsdo(npldafs,nphim,ninitlr,0:natomsym) )
        if( E3E3 )
     &    allocate( ampldafsoo(npldafs,nphim,ninitlr,0:natomsym) )
        if( M1M1 )
     &    allocate( ampldafsmm(npldafs,nphim,ninitlr,0:natomsym) )
        if( E1M1 )
     &    allocate( ampldafsmd(npldafs,nphim,ninitlr,0:natomsym) )
        if( Cor_abs ) then
          allocate( mu(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( E1E1 )
     &      allocate( mudd(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( E1E2 )
     &      allocate( mudq(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( E2E2 )
     &      allocate( muqq(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( E1E3 )
     &      allocate( mudo(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( E3E3 )
     &      allocate( muoo(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( M1M1 )
     &      allocate( mumm(npldafs,nphim,2,ninitlr,0:natomsym) )
          if( E1M1 )
     &      allocate( mumd(npldafs,nphim,2,ninitlr,0:natomsym) )
        endif

      endif        

! Correction du terme magnetique en cas de green_moins.
      if( .not. Green_plus .and. .not. Green_int ) then
        if( E1E1 ) secdd_a(:,:,:,ie_computer)
     &                      = conjg( secdd_a(:,:,:,ie_computer) )
! Comme dans convolution, on prend le complexe conjugue, le cas 
! Green_plus est a prendre avant le img facteur du vecteur d'onde dans
! l'operateur quadrupolaire.
        if( E1E2 ) secdq_a(:,:,:,:,ie_computer)
     &                      = conjg( secdq_a(:,:,:,:,ie_computer) ) 
        if( E2E2 ) secqq_a(:,:,:,:,:,ie_computer)
     &                      = conjg( secqq_a(:,:,:,:,:,ie_computer) ) 
        if( E1E3 ) secdo_a(:,:,:,:,:,ie_computer)
     &                      = conjg( secdo_a(:,:,:,:,:,ie_computer) ) 
        if( E3E3 ) secoo_a(:,:,:,:,:,ie_computer)
     &                      = conjg( secoo_a(:,:,:,:,:,ie_computer) ) 
        if( E1M1 ) secmd_a(:,:,:,ie_computer)
     &                      = conjg( secmd_a(:,:,:,ie_computer) )
        if( M1M1 ) secmm_a(:,:,:,ie_computer)
     &                      = conjg( secmm_a(:,:,:,ie_computer) )
      endif

      do initlr = 1,ninitlr       ! ----------> Boucle sur les seuils ou les etats initiaux

        if( Core_resolved .and. .not. Final_tddft ) then
          if( initlr <= ninit1 ) then  
            iseuil = 1
          else
            iseuil = min(2, nbseuil)
          endif
        elseif( Final_tddft ) then
          iseuil = min(2, nbseuil) 
        else
          iseuil = initlr
        endif

        Ephseuil = Energ(ie)
        Ephoton = Ephseuil + Eseuil(iseuil)
! Pour les seuils de tres basse Energie:
        Ephoton = max(0.001_db/Rydb, Ephoton)
        if( Energphot ) Ephseuil = Ephoton

        ct_nelec(initlr) = conv_mbarn_nelec(Ephoton) 
        eph2 = 0.5 * Ephoton**2
! Pour avoir les tenseurs et sections efficace en Megabarn
        cst = eph2 / ct_nelec(initlr)

! Les tenseurs sont convertis en megabarn
        if( .not. Extract ) then
          if( xan_atom ) sec_atom(initlr) = sec_atom(initlr) * cst 
          if( E1E1 ) secdd_a(:,:,initlr,ie_computer)
     &                   = secdd_a(:,:,initlr,ie_computer) * cst
          if( E1E2 ) secdq_a(:,:,:,initlr,ie_computer)
     &                   = secdq_a(:,:,:,initlr,ie_computer) * cst 
          if( E2E2 ) secqq_a(:,:,:,:,initlr,ie_computer)
     &                   = secqq_a(:,:,:,:,initlr,ie_computer) * cst 
          if( E1E3 ) secdo_a(:,:,:,:,initlr,ie_computer)
     &                   = secdo_a(:,:,:,:,initlr,ie_computer) * cst 
          if( E3E3 ) secoo_a(:,:,:,:,initlr,ie_computer)
     &                   = secoo_a(:,:,:,:,initlr,ie_computer) * cst 
          if( E1M1 ) secmd_a(:,:,initlr,ie_computer)
     &            = - secmd_a(:,:,initlr,ie_computer) * cst
          if( M1M1 ) secmm_a(:,:,initlr,ie_computer)
     &            = secmm_a(:,:,initlr,ie_computer) * cst
          if( Green_int_mag ) then
            if( E1E1 ) secdd_m_a(:,:,initlr,ie_computer)
     &                   = secdd_m_a(:,:,initlr,ie_computer) * cst
            if( E1E2 ) secdq_m_a(:,:,:,initlr,ie_computer)
     &                   = secdq_m_a(:,:,:,initlr,ie_computer) * cst 
            if( E2E2 ) secqq_m_a(:,:,:,:,initlr,ie_computer)
     &                   = secqq_m_a(:,:,:,:,initlr,ie_computer) * cst 
            if( E1E3 ) secdo_m_a(:,:,:,:,initlr,ie_computer)
     &                   = secdo_m_a(:,:,:,:,initlr,ie_computer) * cst 
            if( E3E3 ) secoo_m_a(:,:,:,:,initlr,ie_computer)
     &                   = secoo_m_a(:,:,:,:,initlr,ie_computer) * cst 
            if( E1M1 ) secmd_m_a(:,:,initlr,ie_computer)
     &           = - secmd_m_a(:,:,initlr,ie_computer) * cst
            if( M1M1 ) secmm_m_a(:,:,initlr,ie_computer)
     &          = secmm_m_a(:,:,initlr,ie_computer) * cst
          endif
        endif

        do ia = 1,natomsym

          isym = abs( isymeq(ia) )
          call opsym(isym,matopsym)
          if( base_spin ) then
            matopsym = matmul( matopsym, rot_int )
            matopsym = matmul( transpose(rot_int), matopsym )
          endif

          if( E1E1 ) then
            secdd(:,:) = secdd_a(:,:,initlr,ie_computer)
            if( isym /= 1 ) call rot_tensor_2( secdd, matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                            secdd(:,:) = conjg( secdd(:,:) )
            secddia(:,:,initlr,ia) = secdd(:,:)
            if( Green_int_mag ) then
              secdd(:,:) = secdd_m_a(:,:,initlr,ie_computer)
              if( isym /= 1 ) call rot_tensor_2( secdd, matopsym )
              if( isymeq(ia) < 0 ) secdd(:,:) = - secdd(:,:)
              secddia_m(:,:,initlr,ia) = secdd(:,:)
            endif
          endif

          if( E1E2 ) then
            secdq(:,:,:) = secdq_a(:,:,:,initlr,ie_computer)
            if( isym /= 1 ) call rot_tensor_3( secdq, matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                            secdq(:,:,:) = conjg( secdq(:,:,:) )
            secdqia(:,:,:,initlr,ia) = secdq(:,:,:)
            if( Green_int_mag ) then
              secdq(:,:,:) = secdq_m_a(:,:,:,initlr,ie_computer)
              if( isym /= 1 ) call rot_tensor_3( secdq, matopsym )
              if( isymeq(ia) < 0 ) secdq(:,:,:) = - secdq(:,:,:)
              secdqia_m(:,:,:,initlr,ia) = secdq(:,:,:)
            endif
          endif

          if( E2E2 ) then
            secqq(:,:,:,:) = secqq_a(:,:,:,:,initlr,ie_computer)
            if( isym /= 1 ) call rot_tensor_4( secqq, matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                          secqq(:,:,:,:) = conjg( secqq(:,:,:,:) )
            secqqia(:,:,:,:,initlr,ia) = secqq(:,:,:,:)
            if( Green_int_mag ) then
              secqq(:,:,:,:) = secqq_m_a(:,:,:,:,initlr,ie_computer)
              if( isym /= 1 ) call rot_tensor_4( secqq, matopsym )
              if( isymeq(ia) < 0 ) secqq(:,:,:,:) = - secqq(:,:,:,:)
              secqqia_m(:,:,:,:,initlr,ia) = secqq(:,:,:,:)
            endif
         endif

          if( E1E3 ) then
            secdo(:,:,:,:) = secdo_a(:,:,:,:,initlr,ie_computer)
            if( isym /= 1 ) call rot_tensor_4( secdo, matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                          secdo(:,:,:,:) = conjg( secdo(:,:,:,:) )
            secdoia(:,:,:,:,initlr,ia) = secdo(:,:,:,:)
            if( Green_int_mag ) then
              secdo(:,:,:,:) = secdo_m_a(:,:,:,:,initlr,ie_computer)
              if( isym /= 1 ) call rot_tensor_4( secdo, matopsym )
              if( isymeq(ia) < 0 ) secdo(:,:,:,:) = - secdo(:,:,:,:)
              secdoia_m(:,:,:,:,initlr,ia) = secdo(:,:,:,:)
            endif
          endif

          if( E3E3 ) then
            jhe = 0
            do je = 1,3
              do he = 1,3
                jhe = jhe + 1
                jhs = 0
                do js = 1,3
                  do hs = 1,3
                    jhs = jhs + 1
                    Mat6(:,je,he,:,js,hs)
     &                = secoo_a(:,jhe,:,jhs,initlr,ie_computer) 
                  end do
                end do
              end do
            end do
            if( isym /= 1 ) call rot_tensor_6( Mat6, Matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                                        Mat6 = conjg( Mat6 )
            jhe = 0
            do je = 1,3
              do he = 1,3
                jhe = jhe + 1
                jhs = 0
                do js = 1,3
                  do hs = 1,3
                    jhs = jhs + 1
                    secooia(:,jhe,:,jhs,initlr,ia)
     &                                          = Mat6(:,je,he,:,js,hs) 
                  end do
                end do
              end do
            end do
            if( Green_int_mag ) then
              jhe = 0
              do je = 1,3
                do he = 1,3
                  jhe = jhe + 1
                  jhs = 0
                  do js = 1,3
                    do hs = 1,3
                      jhs = jhs + 1
                      Mat6(:,je,he,:,js,hs)
     &                  = secoo_m_a(:,jhe,:,jhs,initlr,ie_computer) 
                    end do
                  end do
                end do
              end do
              if( isym /= 1 ) call rot_tensor_6( Mat6, Matopsym )
              if( isymeq(ia) < 0 .and. .not. Green_int )
     &                                          Mat6 = conjg( Mat6 )
              jhe = 0
              do je = 1,3
                do he = 1,3
                  jhe = jhe + 1
                  jhs = 0
                  do js = 1,3
                    do hs = 1,3
                      jhs = jhs + 1
                      secooia_m(:,jhe,:,jhs,initlr,ia)
     &                                          = Mat6(:,je,he,:,js,hs) 
                    end do
                  end do
                end do
              end do
            endif
          endif

          if( E1M1 ) then
            secmd(:,:) = secmd_a(:,:,initlr,ie_computer)
            if( isym /= 1 ) call rot_tensor_2( secmd, matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                               secmd(:,:) = conjg( secmd(:,:) )
            secmdia(:,:,initlr,ia) = secmd(:,:)
            if( Green_int_mag ) then
              secmd(:,:) = secmd_m_a(:,:,initlr,ie_computer)
              if( isym /= 1 ) call rot_tensor_2( secmd, matopsym )
              if( isymeq(ia) < 0 ) secmd(:,:) = - secmd(:,:)
              secmdia_m(:,:,initlr,ia) = secmd(:,:)
            endif
          endif

          if( M1M1 ) then
            secmm(:,:) = secmm_a(:,:,initlr,ie_computer)
            if( isym /= 1 ) call rot_tensor_2( secmm, matopsym )
            if( isymeq(ia) < 0 .and. .not. Green_int )
     &                                secmm(:,:) = conjg( secmm(:,:) )
            secmmia(:,:,initlr,ia) = secmm(:,:)
            if( Green_int_mag ) then
              secmm(:,:) = secmm_m_a(:,:,initlr,ie_computer)
              if( isym /= 1 ) call rot_tensor_2( secmm, matopsym )
              if( isymeq(ia) < 0 ) secmm(:,:) = - secmm(:,:)
              secmmia_m(:,:,initlr,ia) = secmm(:,:)
            endif
          endif

        end do

      end do

      E_vec = E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. E1M1 .or. M1M1

      if( dafs ) then
        phdt1(:) = phdt(:,1)
        phdf0t1(:) = phdf0t(:,1)
      endif

      jpl = 0

      if( Cor_abs ) then
        nplt = nplr + 3*npldafs
      else
        nplt = nplr + npldafs
      endif

      do ixandafs = 1,2

        do ipl = 1,nplt

          mu_cal = .false.
          idafs = .false.

          if( ipl > nplr ) then
            if( Cor_abs ) then
              ipldafs = ( ipl - nplr + 2 ) / 3
              ind_mu = mod(ipl - nplr + 2, 3) 
              if( ind_mu == 0 ) then
                idafs = .true.
              else
                mu_cal = .true.
              endif
            else
              ipldafs = ipl - nplr
              idafs = .true.
            endif
          else
            if( ixandafs == 2 ) cycle
            idafs = .false.
            jpl = jpl + 1
            ipldafs = 0
          endif

          if( idafs .and. ixandafs == 1 ) cycle
          if( .not. idafs .and. ixandafs == 2 ) cycle

          tens_comp = magn_sens .or. Green_int .or. idafs
          
          if( idafs .and. ipldafs > 1 ) then
            if( ( hkl_dafs(1,ipldafs) == hkl_dafs(1,ipldafs-1) ) 
     &        .and. ( hkl_dafs(2,ipldafs) == hkl_dafs(2,ipldafs-1) )
     &        .and. ( hkl_dafs(3,ipldafs) == hkl_dafs(3,ipldafs-1) ) )
     &        goto 1010
          endif
          if( ipl > 1 .and. .not. idafs ) goto 1010
          
          secddia(:,:,:,0) = (0._db,0._db)
          secdqia(:,:,:,:,0) = (0._db,0._db)
          secdqia_m(:,:,:,:,0) = (0._db,0._db)
          secqqia(:,:,:,:,:,0) = (0._db,0._db)
          secdoia(:,:,:,:,:,0) = (0._db,0._db)
          secooia(:,:,:,:,:,0) = (0._db,0._db)
          secmdia(:,:,:,0) = (0._db,0._db)
          secmdia_m(:,:,:,0) = (0._db,0._db)
          secmmia(:,:,:,0) = (0._db,0._db)
          if( Green_int_mag ) then
            secddia_m(:,:,:,0) = (0._db,0._db)
            secqqia_m(:,:,:,:,:,0) = (0._db,0._db)
            secdoia_m(:,:,:,:,:,0) = (0._db,0._db)
            secooia_m(:,:,:,:,:,0) = (0._db,0._db)
            secmmia_m(:,:,:,0) = (0._db,0._db)
          endif

          do ia = 1,natomsym

            if( idafs ) then
              if( Green_plus ) then
! Le exp(iQr) est converti. On recupere le complexe conjugue dans convolution.
                ph = conjg( phdafs(ia,ipldafs) )
              else
                ph = phdafs(ia,ipldafs)
              endif
            else
              ph = (1._db, 0._db) * Taux_eq(ia)
            endif
            ph_m = img * ph
                    
            if( E1E1 ) secddia(:,:,:,0) = secddia(:,:,:,0) 
     &                                   + ph * secddia(:,:,:,ia)

            if( E1E2 ) then
              if( Green_int ) then
                secdqia(:,:,:,:,0) = secdqia(:,:,:,:,0) 
     &                        + ph * secdqia(:,:,:,:,ia)
              else
                secdqia(:,:,:,:,0) = secdqia(:,:,:,:,0) 
     &                        + ph * real( secdqia(:,:,:,:,ia), db)
                if( magn_sens ) secdqia_m(:,:,:,:,0)
     &                   = secdqia_m(:,:,:,:,0)
     &                        + ph_m * aimag( secdqia(:,:,:,:,ia) )
              endif
            endif

            if( E2E2 ) secqqia(:,:,:,:,:,0) = secqqia(:,:,:,:,:,0)
     &                                      + ph * secqqia(:,:,:,:,:,ia)

            if( E1E3 ) secdoia(:,:,:,:,:,0) = secdoia(:,:,:,:,:,0)
     &                                      + ph * secdoia(:,:,:,:,:,ia)

            if( E3E3 ) secooia(:,:,:,:,:,0) = secooia(:,:,:,:,:,0)
     &                                      + ph * secooia(:,:,:,:,:,ia)

            if( E1M1 ) then
              if( Green_int ) then
                secmdia(:,:,:,0) = secmdia(:,:,:,0)
     &                           + ph * secmdia(:,:,:,ia)
              else
                secmdia(:,:,:,0) = secmdia(:,:,:,0)
     &                     + ph * real( secmdia(:,:,:,ia), db )

                if( magn_sens ) secmdia_m(:,:,:,0)
     &                     = secmdia_m(:,:,:,0)
     &                     + ph_m * aimag( secmdia(:,:,:,ia) )
              endif
            endif

            if( M1M1 ) secmmia(:,:,:,0) = secmmia(:,:,:,0)
     &                                  + ph * secmmia(:,:,:,ia)

            if( Green_int_mag ) then
              if( E1E1 ) secddia_m(:,:,:,0) = secddia_m(:,:,:,0) 
     &                                      + ph * secddia_m(:,:,:,ia)
              if( E1E2 ) secdqia_m(:,:,:,:,0) = secdqia_m(:,:,:,:,0)
     &                                      + ph * secdqia_m(:,:,:,:,ia)
              if( E2E2 ) secqqia_m(:,:,:,:,:,0) = secqqia_m(:,:,:,:,:,0)
     &                                    + ph * secqqia_m(:,:,:,:,:,ia)

              if( E1E3 ) secdoia_m(:,:,:,:,:,0) = secdoia_m(:,:,:,:,:,0)
     &                                    + ph * secdoia_m(:,:,:,:,:,ia)
              if( E3E3 ) secooia_m(:,:,:,:,:,0) = secooia_m(:,:,:,:,:,0)
     &                                    + ph * secooia_m(:,:,:,:,:,ia)
              if( E1M1 ) secmdia_m(:,:,:,0) = secmdia_m(:,:,:,0)
     &                                    + ph * secmdia_m(:,:,:,ia)
              if( M1M1 ) secmmia_m(:,:,:,0) = secmmia_m(:,:,:,0)
     &                                    + ph * secmmia_m(:,:,:,ia)
            endif

          end do

 1010     continue

          if( cartesian_tensor .and. ( ipl == 1 .or. idafs ) ) then
            do ib = 0,nb
              if( natomsym == 1 .and. ib > 0 ) cycle
              if( ib == 0 ) then
                ia = 1
              elseif( ib == 1 ) then
                ia = 0
              else
                ia = ib
              endif
              if( ia /= 0 .and. ipl > 1 ) cycle
              call write_cartesian_tensor(Densite_atom,E_cut,E1E2,E2E2,
     &             Ephseuil,
     &             Epsii,Eseuil(nbseuil),ia,ie,ipldafs,jseuil,
     &             length_word,M1M1,magn_sens,natomsym,ninit1,ninitlr,
     &             nomfich_s,nseuil,numat_abs,secddia,secdqia,
     &             secdqia_m,secqqia,secmdia,tens_comp,V0muf,
     &             Core_resolved)
            end do
          endif
             
          if( spherical_tensor .and. ( ipl == 1 .or. idafs ) ) then
            do ib = 0,nb
              if( natomsym == 1 .and. ib > 0 ) cycle
              if( ib == 0 ) then
                ia = 1
              elseif( ib == 1 ) then
                ia = 0
              else
                ia = ib
              endif
              if( ia /= 0 .and. ipl > 1 ) cycle
              if( idafs ) then 
                plae(:) = poldafse(:,ipldafs,1)
                plas(:) = poldafss(:,ipldafs,1)
                voae(:) = vecdafse(:,ipldafs,1)
                voas(:) = vecdafss(:,ipldafs,1)
              endif
              call spherical_tensor_cal(ct_nelec,Densite_atom,E_cut,
     &          E1E1,E1E2,E2E2,
     &          Energ,Ephseuil,Epsii,Eseuil(nbseuil),ia,icheck,ie,
     &          Int_tens,
     &          ipl,ipldafs,jseuil,Length_word,magn_sens,moyenne,
     &          natomsym,ncolm,nenerg,ninitlr,nomfich_s,npldafs,nplr,
     &          nplrm,nplt,nseuil,numat_abs,pdp,phdf0t1,phdt1,plae,pol,
     &          plas,secddia,secdqia,secdqia_m,secqqia,v0muf,voae,
     &          vec,voas)
            end do
          endif

          if( ( icheck > 0 .and. ipl == 1 ) .or.
     &        ( icheck > 1 .and. idafs ) ) then

            do initlr = 1,ninitlr

              if( Core_resolved .and. .not. Final_tddft ) then
                if( initlr <= ninit1 ) then  
                  iseuil = 1
                else
                  iseuil = min(2, nbseuil)
                endif
              elseif( Final_tddft ) then
                iseuil = min(2, nbseuil) 
              else
                iseuil = initlr
              endif

              if( Final_tddft ) then
                if( nbseuil == 2 ) then
                  write(3,120) achar(nseuil+74)//achar(jseuil+iseuil+46)
     &                                         //achar(jseuil+iseuil+47)
                else
                  write(3,120) achar(nseuil+74)//achar(jseuil+iseuil+46)
                endif
              elseif( nseuil == 0 ) then  ! optic
                write(3,120) 'Opt'
              else
                if( Core_resolved ) then
                  write(3,130) initlr, ninitlr
                elseif( ninitlr > 1 ) then
                  write(3,135) initlr, ninitlr
                endif
              endif
            
              do ib = 0,nb
                if( natomsym == 1 .and. ib > 0 ) cycle
                if( ib == 0 ) then
                  ia = 1
                elseif( ib == 1 ) then
                  ia = 0
                else
                  ia = ib
                endif
                if( ia /= 0 .and. ipl > 1 ) cycle
                if( ipl > 1 .and. .not. idafs ) cycle
                if( idafs .and. icheck < 2 ) cycle
                if( ipl > 1 ) write(3,140) ipldafs
                if( E1E1 ) then
                  if( ia == 1 ) then
                    if( Green_int_mag ) then
                      write(3,141)
                    elseif( Green_int ) then
                      write(3,142)
                    else
                      write(3,143)
                    endif
                  elseif( ia == 0 ) then
                    if( Green_int_mag ) then
                      write(3,144)
                    elseif( Green_int ) then
                      write(3,145)
                    else
                      write(3,146)
                    endif
                  else
                    if( Green_int_mag ) then
                      write(3,147) ia
                    elseif( Green_int ) then
                      write(3,148) ia
                    else
                      write(3,149) ia
                    endif
                  endif
                  do ke = 1,3
                    if( Green_int_mag ) then
                      write(3,150) secddia(ke,:,initlr,ia),
     &                               secddia_m(ke,:,initlr,ia)
                    elseif( tens_comp ) then
                      write(3,150) secddia(ke,:,initlr,ia)
                    else
                      write(3,150) real( secddia(ke,:,initlr,ia) )
                    endif
                  end do
                endif
                if( E1E2 ) then
                  do ke = 1,3
                    if( ia == 1 ) then
                      if( Green_int_mag ) then
                        write(3,160) ke
                      elseif( Green_int ) then
                        write(3,161) ke
                      else
                        write(3,162) ke
                      endif
                    elseif( ia == 0 ) then
                      if( Green_int_mag ) then
                        write(3,163) ke
                      elseif( Green_int ) then
                        write(3,164) ke
                      elseif( magn_sens ) then
                        write(3,165) ke
                      else
                        write(3,166) ke
                      endif
                    else
                      if( Green_int_mag ) then
                        write(3,167) ia, ke
                      elseif( Green_int ) then
                        write(3,168) ia, ke
                      else
                        write(3,169) ia, ke
                      endif
                    endif
                    do ks = 1,3
                      if( ( magn_sens .and. ia == 0 )
     &                               .or. Green_int_mag ) then
                        write(3,150) secdqia(ke,ks,:,initlr,ia),
     &                               secdqia_m(ke,ks,:,initlr,ia)
                      elseif( tens_comp ) then
                        write(3,150) secdqia(ke,ks,:,initlr,ia)
                      else
                        write(3,150) real(secdqia(ke,ks,:,initlr,ia),db)
                      endif
                    end do
                  end do
                endif
                if( E2E2 ) then
                  do js = 1,3
                    do ks = 1,3
                      if( ia == 1 ) then
                        if( Green_int_mag ) then
                          write(3,171) ks, js
                        elseif( Green_int ) then
                          write(3,172) ks, js
                        else
                          write(3,173) ks, js
                        endif
                      elseif( ia == 0 ) then
                        if( Green_int_mag ) then
                          write(3,174) ks, js
                        elseif( Green_int ) then
                          write(3,175) ks, js
                        else
                          write(3,176) ks, js
                        endif
                      else
                        if( Green_int_mag ) then
                          write(3,177) ia, ks, js
                        elseif( Green_int ) then
                          write(3,178) ia, ks, js
                        else
                          write(3,179) ia, ks, js
                        endif
                      endif
                      do ke = 1,3
                        if( Green_int_mag ) then
                          write(3,150) secqqia(ke,1:3,ks,js,initlr,ia),
     &                                 secqqia_m(ke,1:3,ks,js,initlr,ia)
                        elseif( tens_comp ) then
                          write(3,150) secqqia(ke,1:3,ks,js,initlr,ia)
                        else
                          write(3,150)
     &                          real( secqqia(ke,1:3,ks,js,initlr,ia) )
                        endif
                      end do
                    end do
                  end do
                endif
                if( E1E3 ) then
                  do ke = 1,3
                    do ks = 1,3
                      if( ia == 1 ) then
                        if( Green_int_mag ) then
                          write(3,181) ke, ks
                        elseif( Green_int ) then
                          write(3,182) ke, ks
                        else
                          write(3,183) ke, ks
                        endif
                      elseif( ia == 0 ) then
                        if( Green_int_mag ) then
                          write(3,184) ke, ks
                        elseif( Green_int ) then
                          write(3,185) ke, ks
                        else
                          write(3,186) ke, ks
                        endif
                      else
                        if( Green_int_mag ) then
                          write(3,187) ia, ks, js
                        elseif( Green_int ) then
                          write(3,188) ia, ke, ks
                        else
                          write(3,189) ia, ke, ks
                        endif
                      endif
                      do j1 = 1,3
                        if( Green_int_mag ) then
                          write(3,150) secdoia(ke,ks,j1,:,initlr,ia),
     &                                 secdoia_m(ke,ks,j1,:,initlr,ia)
                        elseif( tens_comp ) then
                          write(3,150) secdoia(ke,ks,j1,:,initlr,ia)
                        else
                          write(3,150)
     &                          real( secdoia(ke,ks,j1,:,initlr,ia) )
                        endif
                      end do
                    end do
                  end do
                endif
                if( E1M1 ) then
                  if( ia == 1 ) then
                    if( Green_int_mag ) then
                      write(3,190) 
                    elseif( Green_int ) then
                      write(3,191)
                    else
                      write(3,192)
                    endif
                  elseif( ia == 0 ) then
                    if( Green_int_mag ) then
                      write(3,193) 
                    elseif( Green_int ) then
                      write(3,194)
                    elseif( magn_sens ) then
                      write(3,195) 
                    else
                      write(3,196)
                    endif
                  else
                    if( Green_int_mag ) then
                      write(3,197) ia
                    elseif( Green_int ) then
                      write(3,198) ia
                    else
                      write(3,199) ia
                    endif
                  endif
                  do ke = 1,3
                    if( ( magn_sens .and. ia == 0 )
     &                               .or. Green_int_mag ) then
                      write(3,150) secmdia(ke,:,initlr,ia),
     &                               secmdia_m(ke,:,initlr,ia)
                    elseif( tens_comp ) then
                      write(3,150) secmdia(ke,:,initlr,ia)
                    else
                      write(3,150) real( secmdia(ke,:,initlr,ia) )
                    endif
                  end do
                endif
                if( M1M1 ) then
                  if( ia == 1 ) then
                    if( Green_int_mag ) then
                      write(3,201)
                    elseif( Green_int ) then
                      write(3,202)
                    else
                      write(3,203)
                    endif
                  elseif( ia == 0 ) then
                    if( Green_int_mag ) then
                      write(3,204)
                    elseif( Green_int ) then
                      write(3,205)
                    else
                      write(3,206)
                    endif
                  else
                    if( Green_int_mag ) then
                      write(3,207) ia
                    elseif( Green_int ) then
                      write(3,208) ia
                    else
                      write(3,209) ia
                    endif
                  endif
                  do ke = 1,3
                    if( Green_int_mag ) then
                      write(3,150) secmmia(ke,:,initlr,ia),
     &                               secmmia_m(ke,:,initlr,ia)
                    elseif( tens_comp ) then
                      write(3,150) secmmia(ke,:,initlr,ia)
                    else
                      write(3,150) real( secmmia(ke,:,initlr,ia) )
                    endif
                  end do
                endif

                if( E3E3 ) then
                  do hs = 1,3
                    do js = 1,3
                      jhs = 3 * ( js - 1 ) + hs 
                      if( ia == 1 ) then
                        if( Green_int_mag ) then
                          write(3,211) (ks, js, hs, ks = 1,3)
                        elseif( Green_int ) then
                          write(3,212) (ks, js, hs, ks = 1,3)
                        else
                          write(3,213) (ks, js, hs, ks = 1,3)
                        endif
                      elseif( ia == 0 ) then
                        if( Green_int_mag ) then
                          write(3,214) (ks, js, hs, ks = 1,3)
                        elseif( Green_int ) then
                          write(3,215) (ks, js, hs, ks = 1,3)
                        else
                          write(3,216) (ks, js, hs, ks = 1,3)
                        endif
                      else
                        if( Green_int_mag ) then
                          write(3,217) ia, (ks, js, hs, ks = 1,3)
                        elseif( Green_int ) then
                          write(3,218) ia, (ks, js, hs, ks = 1,3)
                        else
                          write(3,219) ia, (ks, js, hs, ks = 1,3)
                        endif
                      endif
                      do he = 1,3
                        do je = 1,3
                          jhe = 3 * ( je - 1 ) + he
                          if( Green_int_mag ) then
                            write(3,150)
     &                        ( secooia(1:3,jhe,ks,jhs,initlr,ia),
     &                          secooia_m(1:3,jhe,ks,jhs,initlr,ia),
     &                                                        ks = 1,3)
                          elseif( tens_comp ) then
                            write(3,150)
     &                        ( secooia(1:3,jhe,ks,jhs,initlr,ia),
     &                                                        ks = 1,3)
                          else
                            write(3,150)
     &                      ( real( secooia(1:3,jhe,ks,jhs,initlr,ia) ), 
     &                                                        ks = 1,3)
                          endif
                        end do
                      end do
                    end do
                  end do
                endif

              end do
            end do

          endif

          if( idafs .or. mu_cal ) then
            np = nphi_dafs(ipldafs)
          else
            np = 1
          endif

          do ip = 1,np
        
            if( .not. ( idafs .or. mu_cal ) ) then

              plae(:) = pol(:,ipl)
              plas(:) = pol(:,ipl)
              if( E_vec ) voae(:) = vec(:,ipl)
              if( E_vec ) voas(:) = vec(:,ipl)

            elseif( idafs ) then

              plae(:) = poldafse(:,ipldafs,ip)
              plas(:) = poldafss(:,ipldafs,ip)
              if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
              if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)

            else ! calcul mu

              if( Full_self_abs ) then

                select case(mod(ipldafs,4))
                  case(1,0)
                    if( ind_mu == 1 ) then   ! entrant
                      plae(:) = poldafse(:,ipldafs,ip) 
                      plas(:) = poldafse(:,ipldafs,ip)
                      if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)  
                      if( E_vec ) voas(:) = vecdafse(:,ipldafs,ip)  
                    else
                      plae(:) = poldafss(:,ipldafs,ip) 
                      plas(:) = poldafss(:,ipldafs,ip) 
                      if( E_vec ) voae(:) = vecdafss(:,ipldafs,ip)  
                      if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)
                    endif  
                  case(2)
                    if( ind_mu == 1 ) then   ! entrant
                      plae(:) = poldafse(:,ipldafs,ip) 
                      plas(:) = poldafse(:,ipldafs+1,ip)
                      if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)  
                      if( E_vec ) voas(:) = vecdafse(:,ipldafs+1,ip)  
                    else
                      plae(:) = poldafss(:,ipldafs+1,ip) 
                      plas(:) = poldafss(:,ipldafs+2,ip) 
                      if( E_vec ) voae(:) = vecdafss(:,ipldafs+1,ip)  
                      if( E_vec ) voas(:) = vecdafss(:,ipldafs+2,ip)
                    endif  
                  case(3)
                    if( ind_mu == 1 ) then   ! entrant
                      plae(:) = poldafse(:,ipldafs,ip) 
                      plas(:) = poldafse(:,ipldafs-1,ip)
                      if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)  
                      if( E_vec ) voas(:) = vecdafse(:,ipldafs-1,ip)  
                    else
                      plae(:) = poldafss(:,ipldafs+1,ip) 
                      plas(:) = poldafss(:,ipldafs,ip) 
                      if( E_vec ) voae(:) = vecdafss(:,ipldafs+1,ip)  
                      if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)
                    endif  
                end select

              else  ! Self_abs

                if( ind_mu == 1 ) then
                  plae(:) = poldafse(:,ipldafs,ip)
                  plas(:) = poldafse(:,ipldafs,ip)
                  if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
                  if( E_vec ) voas(:) = vecdafse(:,ipldafs,ip)
                else
                  plae(:) = poldafss(:,ipldafs,ip)
                  plas(:) = poldafss(:,ipldafs,ip)
                  if( E_vec ) voae(:) = vecdafss(:,ipldafs,ip)
                  if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)
                endif

              endif

            endif

            if( Green_plus .and. ( idafs .or. mu_cal ) ) then
! Dans convolution on reprend le complexe conjugue de l'ensemble
! polarisation x diffusion
              plae(:) = conjg( plae(:) )
              plas(:) = conjg( plas(:) )
            endif

            if( M1M1 .or. E1M1 ) then
              uae(1) = voae(2) * plae(3) - voae(3) * plae(2) 
              uae(2) = voae(3) * plae(1) - voae(1) * plae(3) 
              uae(3) = voae(1) * plae(2) - voae(2) * plae(1) 
              uas(1) = voas(2) * plas(3) - voas(3) * plas(2) 
              uas(2) = voas(3) * plas(1) - voas(1) * plas(3) 
              uas(3) = voas(1) * plas(2) - voas(2) * plas(1) 
            endif

            do ia = 0,na
              do initlr = 1,ninitlr

                if( E1E1 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    sec = sec + plae(ke) 
     &                  * sum( conjg(plas(:)) * secddia(:,ke,initlr,ia))
                  end do
                  if( Green_int_mag ) then
                    do ke = 1,3
                      sec = sec + plae(ke) * sum( conjg(plas(:)) 
     &                          * secddia_m(:,ke,initlr,ia) )
                    end do
                  endif 
! Il manque un facteur pi qui a deja ete pris en compte dans le calcul
! du tenseur dans la routine Tens_ab.
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsdd(ipldafs,ip,initlr,ia)=sec
                  elseif( mu_cal ) then
                    mudd(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsdd(jpl,initlr,ia) = real( sec,db )
                  endif
                endif

                if( E1E2 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    do ks = 1,3
                      if( ia == 0 ) then
                        sec = sec + conjg( plas(ks) ) * plae(ke)
     &                      * sum( voae(:) * secdqia(ks,ke,:,initlr,ia)
     &                         - voas(:) * secdqia(ke,ks,:,initlr,ia) )
                        if( magn_sens )     
     &                    sec = sec + conjg( plas(ks) ) * plae(ke)
     &                     * sum( voae(:)*secdqia_m(ks,ke,:,initlr,ia)
     &                        + voas(:)*secdqia_m(ke,ks,:,initlr,ia) )
                      else
                        sec = sec + conjg( plas(ks) ) * plae(ke)
     &                      * sum( voae(:) * secdqia(ks,ke,:,initlr,ia)
     &                      - voas(:)*conjg(secdqia(ke,ks,:,initlr,ia)))
                      endif
                    end do
                  end do

                  if( Green_plus .and. idafs ) then
                    sec = - img * sec  
                  else
                    sec = img * sec      ! C'est ici qu'on recupere le img
                  endif
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsdq(ipldafs,ip,initlr,ia) = sec
                  elseif( mu_cal ) then
                    mudq(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsdq(jpl,initlr,ia) = real( sec,db )
                  endif
                endif

                if( E2E2 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    do je = 1,3
                      do ks = 1,3
                        sec = sec
     &                      + conjg( plas(ks) ) * plae(ke) * voae(je)
     &                      * sum(voas(:)*secqqia(ks,:,ke,je,initlr,ia))
                        if( .not. Green_int_mag ) cycle
                        sec = sec
     &                      + conjg( plas(ks) ) * plae(ke) * voae(je)
     &                    * sum(voas(:)*secqqia_m(ks,:,ke,je,initlr,ia))
                      end do
                    end do
                  end do
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsqq(ipldafs,ip,initlr,ia) = sec
                  elseif( mu_cal ) then
                    muqq(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsqq(jpl,initlr,ia) = real( sec,db )
                  endif
                endif

                if( E1E3 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    do ks = 1,3
                      do j1 = 1,3
                        sec = sec + conjg( plas(ks) ) * plae(ke)
     &                      * ( voae(j1) * sum( voae(:)
     &                      * secdoia(ks,ke,j1,:,initlr,ia) )
     &                      + voas(j1) * sum( voas(:)
     &                      * conjg( secdoia(ke,ks,j1,:,initlr,ia) ) ) )
                        if( .not. Green_int_mag ) cycle
                        sec = sec + conjg( plas(ks) ) * plae(ke)
     &                      * ( voae(j1) * sum( voae(:)
     &                      * secdoia_m(ks,ke,j1,:,initlr,ia) )
     &                      + voas(j1) * sum( voas(:)
     &                      * conjg( secdoia_m(ke,ks,j1,:,initlr,ia) )))
                      end do
                    end do
                  end do
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsdo(ipldafs,ip,initlr,ia)=sec
                  elseif( mu_cal ) then
                    mudo(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsdo(jpl,initlr,ia) = real( sec, db )
                  endif
                endif

                if( E3E3 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    do je = 1,3
                      do he = 1,3
                        jhe = 3 * ( je - 1 ) + he
                        do ks = 1,3
                          do js = 1,3
                            do hs = 1,3
                              jhs = 3 * ( js - 1 ) + hs
                              sec = sec
     &                          +conjg( plas(ks) ) * voas(js) * voas(hs)
     &                          * plae(ke) * voae(je) * voae(he)
     &                          * secooia(ks,jhs,ke,jhe,initlr,ia)
                              if( .not. Green_int_mag ) cycle
                              sec = sec
     &                          +conjg( plas(ks) ) * voas(js) * voas(hs)
     &                          * plae(ke) * voae(je) * voae(he)
     &                          * secooia_m(ks,jhs,ke,jhe,initlr,ia)
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsoo(ipldafs,ip,initlr,ia) = sec
                  elseif( mu_cal ) then
                    muoo(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsoo(jpl,initlr,ia) = real( sec,db )
                  endif
                endif

                if( M1M1 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    sec = sec + uae(ke) 
     &                  * sum( conjg(uas(:)) * secmmia(:,ke,initlr,ia) )
                    if( .not. Green_int_mag ) cycle
                    sec = sec + uae(ke) 
     &                  * sum( conjg(uas(:))*secmmia_m(:,ke,initlr,ia) )
                  end do
! Il manque un facteur pi qui a deja ete pris en compte dans le calcul
! du tenseur dans la routine Tens_ab.
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsmm(ipldafs,ip,initlr,ia)=sec
                  elseif( mu_cal ) then
                    mumm(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsmm(jpl,initlr,ia) = real( sec, db )
                  endif
                endif

                if( E1M1 ) then
                  sec = (0._db,0._db)
                  do ke = 1,3
                    if( ia == 0 ) then
                      sec = sec + conjg( uas(ke) ) 
     &                    * sum( plae(:) * secmdia(:,ke,initlr,ia) )
     &                    + uae(ke) 
     &                    * sum( conjg(plas(:))*secmdia(:,ke,initlr,ia))
                      if( magn_sens )     
     &                  sec = sec + conjg( uas(ke) ) 
     &                    * sum( plae(:) * secmdia_m(:,ke,initlr,ia) )
     &                    - uae(ke) * sum( conjg( plas(:) ) 
     &                    * secmdia_m(:,ke,initlr,ia) )
                    else
                      sec = sec + conjg( uas(ke) ) 
     &                    * sum( plae(:) * secmdia(:,ke,initlr,ia) )
     &                    + uae(ke) 
     &                    * conjg( sum(plas(:)*secmdia(:,ke,initlr,ia)))
                      if( .not. Green_int_mag ) cycle
                      sec = sec + conjg( uas(ke) ) 
     &                    * sum( plae(:) * secmdia_m(:,ke,initlr,ia) )
     &                    + uae(ke) 
     &                  * conjg( sum(plas(:)*secmdia_m(:,ke,initlr,ia)))
                    endif
                  end do
! Il manque un facteur pi qui a deja ete pris en compte dans le calcul
! du tenseur dans la routine Tens_ab.
                  if( idafs ) then
                    if( Green_int ) sec = pi * img * conjg( sec )
                    ampldafsmd(ipldafs,ip,initlr,ia)=sec
                  elseif( mu_cal ) then
                    mumd(ipldafs,ip,ind_mu,initlr,ia) = sec
                  else
                    secabsmd(jpl,initlr,ia) = real( sec, db )
                  endif
                endif

              end do   ! Fin de la boucle sur initlr
            end do   ! Fin de la boucle sur les atomes

            if( ipl > nplr ) cycle
            if( ltypcal(ipl) == 'xanes circ d' ) jpl = jpl + 1

          end do

        end do   ! Fin de la boucle sur les polarisation
      end do   

      if( Moyenne ) then
        if( Xan_atom ) then
          i = ncolr - 1
        else
          i = ncolr
        endif
        ipl = 0
        do j = 1,ncolr
          if( ipl > 1 ) then
            if( ltypcal(ipl) == 'xanes circ d' ) cycle
          endif
          ipl = ipl + 1
          if( E1E1 ) secabsdd(i,:,:) = secabsdd(i,:,:)
     &                               + pdp(ipl,1) * secabsdd(j,:,:)
          if( E2E2 ) secabsqq(i,:,:) = secabsqq(i,:,:) 
     &                               + pdp(ipl,2) * secabsqq(j,:,:)
          if( E1E3 ) secabsdo(i,:,:) = secabsdo(i,:,:)
     &                               + pdp(ipl,1) * secabsdo(j,:,:)
          if( E3E3 ) secabsoo(i,:,:) = secabsoo(i,:,:)
     &                               + pdp(ipl,1) * secabsoo(j,:,:)
          if( M1M1 ) secabsmm(i,:,:) = secabsmm(i,:,:)
     &                               + pdp(ipl,1) * secabsmm(j,:,:)
          if( E1M1 ) secabsmd(i,:,:) = secabsmd(i,:,:)
     &                               + pdp(ipl,1) * secabsmd(j,:,:)
          if( ipl == nplr ) exit
        end do
        if( E1E2 ) secabsdq(i,:,:) = (0._db,0._db)
      endif

      jpl = 0
      do ipl = 1,nplr
        jpl = jpl + 1
        if( ltypcal(ipl) /= 'xanes circ d' ) cycle
        jpl = jpl + 1
        ig = jpl - 1
        id = jpl - 2
        if( E1E1 ) secabsdd(jpl,:,:) = secabsdd(id,:,:)
     &                               - secabsdd(ig,:,:)
        if( E1E2 ) secabsdq(jpl,:,:) = secabsdq(id,:,:)
     &                               - secabsdq(ig,:,:)
        if( E2E2 ) secabsqq(jpl,:,:) = secabsqq(id,:,:)
     &                               - secabsqq(ig,:,:)
        if( E1E3 ) secabsdo(jpl,:,:) = secabsdo(id,:,:)
     &                               - secabsdo(ig,:,:)
        if( E3E3 ) secabsoo(jpl,:,:) = secabsoo(id,:,:)
     &                               - secabsoo(ig,:,:)
        if( M1M1 ) secabsmm(jpl,:,:) = secabsmm(id,:,:)
     &                               - secabsmm(ig,:,:)
        if( E1M1 ) secabsmd(jpl,:,:) = secabsmd(id,:,:)
     &                               - secabsmd(ig,:,:)
      end do

      if( xan_atom ) then
        secabsdd(ncolr,:,0) = sec_atom(:) * natomsym
        do ia = 1,na
          secabsdd(ncolr,:,ia) = sec_atom(:)
        end do
      endif

      secabs(:,:,:) = 0._db
      if( E1E1 ) secabs(:,:,:) = secabs(:,:,:) + secabsdd(:,:,:)
      if( E1E2 ) secabs(:,:,:) = secabs(:,:,:) + secabsdq(:,:,:)
      if( E2E2 ) secabs(:,:,:) = secabs(:,:,:) + secabsqq(:,:,:)
      if( E1E3 ) secabs(:,:,:) = secabs(:,:,:) + secabsdo(:,:,:)
      if( E3E3 ) secabs(:,:,:) = secabs(:,:,:) + secabsoo(:,:,:)
      if( M1M1 ) secabs(:,:,:) = secabs(:,:,:) + secabsmm(:,:,:)
      if( E1M1 ) secabs(:,:,:) = secabs(:,:,:) + secabsmd(:,:,:)

      if( Dafs ) then
        ampldafs(:,:,:,:) = (0._db,0._db)
        if( E1E1 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsdd(:,:,:,:)
        if( E1E2 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsdq(:,:,:,:)
        if( E2E2 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsqq(:,:,:,:)
        if( E1E3 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsdo(:,:,:,:)
        if( E3E3 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsoo(:,:,:,:)
        if( M1M1 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsmm(:,:,:,:)
        if( E1M1 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:)
     &                               + ampldafsmd(:,:,:,:)
      endif

      if( Cor_abs ) then
        mu(:,:,:,:,:) = (0._db,0._db)
        if( E1E1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mudd(:,:,:,:,:)
        if( E1E2 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mudq(:,:,:,:,:)
        if( E2E2 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + muqq(:,:,:,:,:)
        if( E1E3 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mudo(:,:,:,:,:)
        if( E3E3 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + muoo(:,:,:,:,:)
        if( M1M1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mumm(:,:,:,:,:)
        if( E1M1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mumd(:,:,:,:,:)
      endif

! Conversion en nombre d'electrons
      if( dafs ) then
        do initlr = 1,ninitlr
          ampldafs(:,:,initlr,:) = ct_nelec(initlr)
     &                           * ampldafs(:,:,initlr,:)
          if( E1E1 ) ampldafsdd(:,:,initlr,:) = ct_nelec(initlr)
     &                                        * ampldafsdd(:,:,initlr,:)
          if( E1E2 ) ampldafsdq(:,:,initlr,:) = ct_nelec(initlr) 
     &                                        * ampldafsdq(:,:,initlr,:)
          if( E2E2 ) ampldafsqq(:,:,initlr,:) = ct_nelec(initlr) 
     &                                        * ampldafsqq(:,:,initlr,:)
          if( M1M1 ) ampldafsmm(:,:,initlr,:) = ct_nelec(initlr) 
     &                                        * ampldafsmm(:,:,initlr,:)
          if( E1M1 ) ampldafsmd(:,:,initlr,:) = ct_nelec(initlr) 
     &                                        * ampldafsmd(:,:,initlr,:)
        end do
      endif

! Conversion en micrometres
      if( Cor_abs ) then
        c_micro = 100 / ( Volume_maille * bohr**3 )
        mu(:,:,:,:,:) = c_micro * mu(:,:,:,:,:)
        if( E1E1 ) mudd(:,:,:,:,:) = c_micro * mudd(:,:,:,:,:)
        if( E1E2 ) mudq(:,:,:,:,:) = c_micro * mudq(:,:,:,:,:)
        if( E2E2 ) muqq(:,:,:,:,:) = c_micro * muqq(:,:,:,:,:)
        if( E1E3 ) mudo(:,:,:,:,:) = c_micro * mudo(:,:,:,:,:)
        if( E3E3 ) muoo(:,:,:,:,:) = c_micro * muoo(:,:,:,:,:)
        if( M1M1 ) mumm(:,:,:,:,:) = c_micro * mumm(:,:,:,:,:)
        if( E1M1 ) mumd(:,:,:,:,:) = c_micro * mumd(:,:,:,:,:)
      endif
      if( nseuil == 0 ) then ! cas de l'optique: en millimetre^-1
        c_milli = 100000 / ( Volume_maille * bohr**3 )
        secabs(:,:,:) = c_milli * secabs(:,:,:)
        if( E1E1 ) secabsdd(:,:,:) = c_milli * secabsdd(:,:,:)
        if( E1E2 ) secabsdq(:,:,:) = c_milli * secabsdq(:,:,:)
        if( E2E2 ) secabsqq(:,:,:) = c_milli * secabsqq(:,:,:)
        if( E1E3 ) secabsdo(:,:,:) = c_milli * secabsdo(:,:,:)
        if( E3E3 ) secabsoo(:,:,:) = c_milli * secabsoo(:,:,:)
        if( M1M1 ) secabsmm(:,:,:) = c_milli * secabsmm(:,:,:)
        if( E1M1 ) secabsmd(:,:,:) = c_milli * secabsmd(:,:,:)
      endif

      if( icheck > 0 ) then
        do ia = 0,na
          if( ia == 0 ) write(3,283) ct_nelec(:) * pi
          if( ia == 0 ) then
            write(3,285)
          else
            write(3,290) ia
          endif
          do ipl = 1,ncolt
            nomab = nomabs(ipl)
            call center_word(nomab,Length_word)
            nomabs(ipl) = nomab
          end do
          nccm = 36
          nl = 1 + ( ncolr - 1 ) / nccm

          do initlr = 1,ninitlr

            if( ninitlr > 1 ) write(3,295) initlr

            do i = 1,nl
              ic1 = 1 + ( i - 1 ) * nccm
              ic2 = min( i * nccm, ncolr )
              write(3,300) nomabs(ic1:ic2)
              write(3,310) Ephseuil*rydb, secabs(ic1:ic2,initlr,ia)
              if( E1E1 .and. ( E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. 
     &        E3E3 .or. E1M1) ) write(3,320) secabsdd(ic1:ic2,initlr,ia)
              if( E1E2 .and. ( E1E1 .or. E2E2 .or. E1E3 .or. M1M1 .or.  
     &        E3E3 .or. E1M1) ) write(3,330) secabsdq(ic1:ic2,initlr,ia)
              if( E2E2 .and. ( E1E1 .or. E2E2 .or. E1E3 .or. M1M1 .or. 
     &        E3E3 .or. E1M1) ) write(3,340) secabsqq(ic1:ic2,initlr,ia)
              if( E1E3 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. M1M1 .or. 
     &        E3E3 .or. E1M1) ) write(3,350) secabsdo(ic1:ic2,initlr,ia)
              if( E3E3 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. M1M1 .or. 
     &        E1E3 .or. E1M1) ) write(3,351) secabsoo(ic1:ic2,initlr,ia)
              if( M1M1 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. E1E3 .or. 
     &        E3E3 .or. E1M1) ) write(3,352) secabsmm(ic1:ic2,initlr,ia)
              if( E1M1 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. E1E3 .or. 
     &        E3E3 .or. M1M1) ) write(3,354) secabsmd(ic1:ic2,initlr,ia)
            end do
            if( dafs ) then
              if( self_abs ) then
                nc = 4
              elseif( Full_self_abs ) then
                nc = 6
              else
                nc = 2
              endif
              nl = 1 + ( nc * npldafs - 1 ) / nccm
              do i = 1,nl
                icn1 = 1 + ( i - 1 ) * nccm
                icn2 = min( i * nccm, nc * npldafs )
                ic1 = 1 + ( i - 1 ) * (nccm/nc)
                ic2 = min( i * (nccm/nc), npldafs )
                write(3,360) nomabs(ncolr+icn1:ncolr+icn2)
                if( self_abs ) then
                  write(3,370) (ampldafs(j,1,initlr,ia),
     &                 Real(mu(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,320) (ampldafsdd(j,1,initlr,ia),
     &                 Real(mudd(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,330) (ampldafsdq(j,1,initlr,ia),
     &                 Real(mudq(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,340) (ampldafsqq(j,1,initlr,ia),
     &                 Real(muqq(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,350) (ampldafsdo(j,1,initlr,ia),
     &                 Real(mudo(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,350) (ampldafsoo(j,1,initlr,ia),
     &                 Real(mudo(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. E1E3 .or. E1M1 ) )
     &              write(3,352) (ampldafsmm(j,1,initlr,ia),
     &                 Real(mumm(j,1,:,initlr,ia)), j = ic1,ic2)
                  if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. E1E3 .or. M1M1 ) )
     &              write(3,354) (ampldafsmd(j,1,initlr,ia),
     &                 Real(mumd(j,1,:,initlr,ia)), j = ic1,ic2)
                elseif(Full_self_abs ) then
                  write(3,370) (ampldafs(j,1,initlr,ia),
     &                 mu(j,1,:,initlr,ia), j = ic1,ic2)
                  if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,320) (ampldafsdd(j,1,initlr,ia),
     &                 mudd(j,1,:,initlr,ia), j = ic1,ic2)
                  if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,330) (ampldafsdq(j,1,initlr,ia),
     &                 mudq(j,1,:,initlr,ia), j = ic1,ic2)
                  if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,340) (ampldafsqq(j,1,initlr,ia),
     &                 muqq(j,1,:,initlr,ia), j = ic1,ic2)
                  if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,350) (ampldafsdo(j,1,initlr,ia),
     &                 mudo(j,1,:,initlr,ia), j = ic1,ic2)
                  if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,350) (ampldafsoo(j,1,initlr,ia),
     &                 mudo(j,1,:,initlr,ia), j = ic1,ic2)
                  if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. E1E3 .or. E1M1 ) )
     &              write(3,352) (ampldafsmm(j,1,initlr,ia),
     &                 mumm(j,1,:,initlr,ia), j = ic1,ic2)
                  if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. E1E3 .or. M1M1 ) )
     &              write(3,354) (ampldafsmd(j,1,initlr,ia),
     &                 mumd(j,1,:,initlr,ia), j = ic1,ic2)
                else
                  write(3,370) ampldafs(ic1:ic2,1,initlr,ia)
                  if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,320) ampldafsdd(ic1:ic2,1,initlr,ia)
                  if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,330) ampldafsdq(ic1:ic2,1,initlr,ia)
                  if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,340) ampldafsqq(ic1:ic2,1,initlr,ia)
                  if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,350) ampldafsdo(ic1:ic2,1,initlr,ia)
                  if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3
     &                               .or. M1M1 .or. E1M1 ) )
     &              write(3,350) ampldafsoo(ic1:ic2,1,initlr,ia)
                  if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. E1E3 .or. E1M1 ) )
     &              write(3,352) ampldafsmm(ic1:ic2,1,initlr,ia)
                  if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3
     &                               .or. E1E3 .or. M1M1 ) )
     &              write(3,354) ampldafsmd(ic1:ic2,1,initlr,ia)
                endif
              end do
            endif
          end do
        end do
      endif

      if( ie == 1 ) write(6,392) nenerg
      n_dim = ncolm * ninitlr

      do ia = 0,na

        nomficht = nomfich
        nomfichdafst = nomfich

        if( ia > 0 ) then
          long = len_trim(nomficht)
          nomficht(long+1:long+5) = '_atom'
          call ad_number(ia,nomficht,132)
          nomfichdafst(long+1:long+5) = '_atom'
          call ad_number(ia,nomfichdafst,132)
        endif
        long = len_trim(nomficht)

        if( Final_tddft .and. .not. Extract ) then
          nomficht(long+1:long+6) = '_tddft'
          nomfichdafst(long+1:long+11) = '_tddft_scan'
        else
          nomfichdafst(long+1:long+5) = '_scan'
        end if

        if( n_multi_run > 1 ) then
          l = len_trim(nomficht)
          nomficht(l+1:l+1) = '_'
          call ad_number(iabsorig,nomficht,132)
          l = len_trim(nomfichdafst)
          nomfichdafst(l+1:l+1) = '_'
          call ad_number(iabsorig,nomfichdafst,132)
        endif
        l = len_trim(nomficht)
        nomficht(l+1:l+4) = '.txt'
        l = len_trim(nomfichdafst)
        nomfichdafst(l+1:l+4) = '.txt'

        if( ie == 1 ) nomfich_cal_convt = nomficht

        n_tens = 0

        do initlr = 1,ninitlr

          n1 = ( initlr - 1 ) * ( ncolt-nxanout+1 ) + 1
          n2 = initlr * ( ncolt-nxanout+1 )
          title(n1:n2) = nomabs(nxanout:ncolt)

          if( ninitlr > 1 ) then
            do i = n1, n2
              nomab = title(i)
              ll = len_trim( nomab )
              if( ll > length_word - 3 ) cycle 
              if( nomab(ll:ll) /= '>' .and. n2-n1+1 /= ncolr ) cycle
              ll = ll + 1
              nomab(ll:ll) = '_'
              if( .not. Core_resolved ) then
                ll = ll + 1
                nomab(ll:ll) = achar(nseuil+74)
                ll = ll + 1
                nomab(ll:ll) = achar(jseuil+initlr+47)
              else
                call ad_number(initlr,nomab,length_word)
              endif
              call center_word(nomab,Length_word)
              title(i) = nomab 
            end do
          endif

          ipl = n_tens + ncolr - nxanout + 1
          Tens(n_tens+1:ipl) = secabs(nxanout:ncolr,initlr,ia) 
          do ipldafs = 1,npldafs
            if( ia == 0 ) then
              cf = ampldafs(ipldafs,1,initlr,ia) 
            else
              if( Green_plus ) then
! Le exp(iQr) est converti. On recupere le complexe conjugue dans convolution.
                cf = conjg( phdafs(ia,ipldafs) )
     &                                  * ampldafs(ipldafs,1,initlr,ia)
              else
                cf = phdafs(ia,ipldafs) * ampldafs(ipldafs,1,initlr,ia)
              endif
            endif
            ipl = ipl + 1
            Tens(ipl) = real( cf,db )
            ipl = ipl + 1
            Tens(ipl) = aimag( cf )
            if( self_abs ) then
              do i = 1,2
                ipl = ipl + 1
                Tens(ipl) = real( mu(ipldafs,1,i,initlr,ia), db )
              end do
            elseif( Full_self_abs ) then
              do i = 1,2
                ipl = ipl + 1
                Tens(ipl) = real( mu(ipldafs,1,i,initlr,ia), db )
                ipl = ipl + 1
                Tens(ipl) = aimag( mu(ipldafs,1,i,initlr,ia) )
              end do
            endif
          end do

          n_tens = ipl

        end do

        if( ia == 0 ) then
          phdtem(:) = phdt(:,1)
        else
          phdtem(:) = phdafs(ia,:)
        endif

! Ecriture dans le fichier
        if( Full_self_abs .or. Self_abs ) then
          call write_out(angxyz,axyz,Densite_atom,fpp_avantseuil,E_cut,
     &          Ephseuil,
     &          Epsii,Eseuil(nbseuil),Green_int,hkl_dafs,ie,Length_word,
     &          jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,
     &          npldafs,npldafs,3,npldafs,nseuil,numat_abs,phdtem,
     &          phdf0t1,tens,v0muf,Core_resolved,natomsym)
        else
          call write_out(rdum,rdum,Densite_atom,fpp_avantseuil,E_cut,
     &          Ephseuil,
     &          Epsii,Eseuil(nbseuil),Green_int,idum,ie,Length_word,
     &          jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,
     &          npldafs,npldafs,0,0,nseuil,numat_abs,phdtem,
     &          phdf0t1,tens,v0muf,Core_resolved,natomsym)
        endif

! Ecriture a l'ecran
        if( ia == 0 )
     &    call write_out(rdum,rdum,Densite_atom,fpp_avantseuil,E_cut,
     &          Ephseuil,
     &          Epsii,Eseuil(nbseuil),Green_int,idum,ie,Length_word,
     &          jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,
     &          npldafs,npldafs,0,0,nseuil,-1,phdtem,
     &          phdf0t1,tens,v0muf,Core_resolved,natomsym)

        if( Dafs .and. nphim > 1 ) then
          if( ie == 1 ) then
            open(7, file = nomfichdafst)
            do ipl = 1,npldafs
              write(7,400) nphi_dafs(ipl)
            end do
          else
            open(7, file = nomfichdafst, position='append')
          endif
          mot = ' '
          i = 6 
          do initlr = 1,ninitlr
            do ii = 1,3
              select case(ii)
                case(1)
                  nomab = 'Amplitude_'
                case(2) 
                  if( .not. ( Full_self_abs .or. Self_abs ) ) cycle
                  if( Full_self_abs ) then
                    nomab = 'mu_in_'
                  else
                    nomab = 'mu_io_'
                  endif 
                case(3) 
                  if( .not. Full_self_abs ) cycle
                  nomab = 'mu_ou_' 
              end select
              if( .not. Core_resolved ) then
                nomab(11:11) = achar(nseuil+74)
                nomab(12:12) = achar(jseuil+initlr+47)
              else
                call ad_number(initlr,nomab,length_word)
              endif
              l = len_trim(nomab)
              mot(i:i+l-1) = nomab(1:l)
              i = i+26
            end do
          end do
          i = i - 9
          mot(i+1:i+26) = '   Non-resonant Amplitude '
          i = i + 26
          mot(i+1:i+26) = ' e_s.e_i * Somme_exp(iQR) '
          i = i + 26

          l = len_trim(mot)
          write(7,405) Ephseuil*rydb, mot(1:l)

          do ipl = 1,npldafs

            if( Dafs_bio ) then
              write(7,407) hkl_dafs(:,ipl)
            else
              write(7,410) hkl_dafs(:,ipl), isigpi(ipl,:)
              dang = 360._db / nphi_dafs(ipl)
            endif

            do ip = 1,nphi_dafs(ipl)
              ang = ( ip - 1 ) * dang
              if( ia == 0 ) then
                cf = (1._db,0._db)
              else
                cf = phdafs(ia,ipl)
                if( Green_plus ) cf = conjg( cf )
              endif
              nw = 0
              do initlr = 1,ninitlr
                nw = nw + 1
                compnum(nw) = cf * ampldafs(ipl,ip,initlr,ia)
                if( Full_self_abs ) then
                  do ind_mu = 1,2
                    nw = nw + 1
                    compnum(nw) = mu(ipl,ip,ind_mu,initlr,ia)
                  end do
                elseif( Self_abs ) then
                  nw = nw + 1
                  compnum(nw) = Cmplx( Real(mu(ipl,ip,1,initlr,ia), db),
     &                               Real(mu(ipl,ip,2,initlr,ia),db),db)
                endif
              end do
              nw = nw + 1
              compnum(nw) = phdf0t(ipl,ip)
              nw = nw + 1
              compnum(nw) = phdt(ipl,ip)
              if( Dafs_bio ) then
                write(7,420) compnum(1:nw)
              else
                write(7,430) ang, compnum(1:nw)
              endif
            end do
          end do
          close(7)
        endif

      end do  ! fin boucle sur atomes

      if( dafs ) then
        deallocate( ampldafs )
        if( E1E1 ) deallocate( ampldafsdd )
        if( E1E2 ) deallocate( ampldafsdq )
        if( E2E2 ) deallocate( ampldafsqq )
        if( E1E3 ) deallocate( ampldafsdo )
        if( E3E3 ) deallocate( ampldafsoo )
        if( M1M1 ) deallocate( ampldafsmm )
        if( E1M1 ) deallocate( ampldafsmd )
        if( Cor_abs ) then
          deallocate( mu )
          if( E1E1 ) deallocate( mudd )
          if( E1E2 ) deallocate( mudq )
          if( E2E2 ) deallocate( muqq )
          if( E1E3 ) deallocate( mudo )
          if( E3E3 ) deallocate( muoo )
          if( M1M1 ) deallocate( mumm )
          if( E1M1 ) deallocate( mumd )
        endif
      endif

      return
  110 format(/' ---- Coabs --------',100('-'))
  120 format(/' Threshold ',a3)
  130 format(/'   Initial state =',i3,' on',i3)
  135 format(/'            Edge =',i3,' on',i3)
  140 format(/8x,' RXS polarization number',i3)
  141 format(/' Tensor_dd(ke,ks), Prototypical atom,',4x,
     &        ' Green integral, Non magnetic part',82x,'Magnetic part')
  142 format(/' Tensor_dd(ke,ks), Prototypical atom,',4x,
     &        ' Green integral')
  143 format(/' Tensor_dd(ke,ks), Prototypical atom')
  144 format(/' Crystal Tensor_dd(ke,ks),',
     &        ' Green integral, Non magnetic part',95x,'Magnetic part')
  145 format(/' Crystal Tensor_dd(ke,ks), Green integral')
  146 format(/' Crystal Tensor_dd(ke,ks)')
  147 format(/' Atom ',i3,' Tensor_dd(ke,ks),',
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  148 format(/' Atom ',i3,' Tensor_dd(ke,ks),',
     &        ' Green integral')
  149 format(/' Atom ',i3,' Tensor_dd(ke,ks)')
  150 format(1p,12(6e18.10,2x))
  160 format(/' Tensor_dq(',i1,',ks,j2), Prototypical atom,',2x,
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  161 format(/' Tensor_dq(',i1,',ks,j2), Prototypical atom,',2x,
     &        ' Green integral')
  162 format(/' Tensor_dq(',i1,',ks,j2), Prototypical atom')
  163 format(/' Crystal Tensor_dq(',i1,',ks,j2),',
     &        ' Green integral, Non magnetic part',89x,'Magnetic part')
  164 format(/' Crystal Tensor_dq(',i1,',ks,j2),',
     &        ' Green integral')
  165 format(/' Crystal Tensor_dq(',i1,',ks,j2),',
     &        '                 Non magnetic part',89x,'Magnetic part')
  166 format(/' Crystal Tensor_dq(',i1,',ks,j2)')
  167 format(/' Atom ',i3,' Tensor_dq(ke,ks),',
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  168 format(/' Atom ',i3,' Tensor_dq(',i1,',ks,j2),  Green integral')
  169 format(/' Atom ',i3,' Tensor_dq(',i1,',ks,j2)')
  171 format(/' Tensor_qq(ke,je,',i1,',',i1,'), Prototypical atom,',
     &        ' Green integral, Non magnetic part',78x,'Magnetic part')
  172 format(/' Tensor_qq(ke,je,',i1,',',i1,'), Prototypical atom,',
     &        ' Green integral')
  173 format(/' Tensor_qq(ke,je,',i1,',',i1,'), Prototypical atom')
  174 format(/' Crystal Tensor_qq(ke,je,',i1,',',i1,'),',
     &        ' Green integral, Non magnetic part',87x,'Magnetic part')
  175 format(/' Crystal Tensor_qq(ke,je,',i1,',',i1,'),',
     &        ' Green integral')
  176 format(/' Crystal Tensor_qq(ke,je,',i1,',',i1,')')
  177 format(/' Atom ',i3,' Tensor_qq(ke,je,',i1,',',i1,'),',
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  178 format(/' Atom ',i3,' Tensor_qq(ke,je,',i1,',',i1,'),',
     &        ' Green integral')
  179 format(/' Atom ',i3,' Tensor_qq(ke,je,',i1,',',i1,')')

  181 format(/' Tensor_do(',i1,',',i1,',j1,j2), Prototypical atom,',
     &        ' Green integral, Non magnetic part',78x,'Magnetic part')
  182 format(/' Tensor_do(',i1,',',i1,',j1,j2), Prototypical atom,',
     &        ' Green integral')
  183 format(/' Tensor_do(',i1,',',i1,',j1,j2), Prototypical atom')
  184 format(/' Crystal Tensor_do(',i1,',',i1,',j1,j2),',
     &        ' Green integral, Non magnetic part',87x,'Magnetic part')
  185 format(/' Crystal Tensor_do(',i1,',',i1,',j1,j2),',
     &        ' Green integral')
  186 format(/' Crystal Tensor_do(',i1,',',i1,',j1,j2)')
  187 format(/' Atom ',i3,' Tensor_do(',i1,',',i1,',j1,j2),',
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  188 format(/' Atom ',i3,' Tensor_do(',i1,',',i1,',j1,j2),',
     &        ' Green integral')
  189 format(/' Atom ',i3,' Tensor_do(',i1,',',i1,',j1,j2)')

  190 format(/' Tensor_md(ke,ks), Prototypical atom,',4x,
     &        ' Green integral, Non magnetic part',82x,'Magnetic part')
  191 format(/' Tensor_md(ke,ks), Prototypical atom,',4x,
     &        ' Green integral')
  192 format(/' Tensor_md(ke,ks), Prototypical atom')
  193 format(/' Crystal Tensor_md(ke,ks),',
     &        ' Green integral, Non magnetic part',89x,'Magnetic part')
  194 format(/' Crystal Tensor_md(ke,ks),',
     &        ' Green integral')
  195 format(/' Crystal Tensor_md(ke,ks),',
     &        '                 Non magnetic part',89x,'Magnetic part')
  196 format(/' Crystal Tensor_md(ke,ks)')
  197 format(/' Atom ',i3,' Tensor_md(ke,ks),',
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  198 format(/' Atom ',i3,' Tensor_md(ke,ks),  Green integral')
  199 format(/' Atom ',i3,' Tensor_md(ke,ks)')
  201 format(/' Tensor_mm(ke,ks), Prototypical atom,',4x,
     &        ' Green integral, Non magnetic part',82x,'Magnetic part')
  202 format(/' Tensor_mm(ke,ks), Prototypical atom,',4x,
     &        ' Green integral')
  203 format(/' Tensor_mm(ke,ks), Prototypical atom')
  204 format(/' Crystal Tensor_mm(ke,ks),',
     &        ' Green integral, Non magnetic part',91x,'Magnetic part')
  205 format(/' Crystal Tensor_mm(ke,ks), Green integral')
  206 format(/' Crystal Tensor_mm(ke,ks)')
  207 format(/' Atom ',i3,' Tensor_mm(ke,ks),',
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  208 format(/' Atom ',i3,' Tensor_mm(ke,ks),',
     &        ' Green integral')
  209 format(/' Atom ',i3,' Tensor_mm(ke,ks)')
  211 format(/3(' Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Prototypical atom,',
     &        ' Green integral, Non magnetic part',78x,'Magnetic part')
  212 format(/3(' Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Prototypical atom,',
     &        ' Green integral')
  213 format(/3(' Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Prototypical atom')
  214 format(/3(' Crystal Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Green integral, Non magnetic part',87x,'Magnetic part')
  215 format(/3(' Crystal Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Green integral')
  216 format(/3(' Crystal Tensor_oo(ke,je,he',3(',',i1),')',36x))
  217 format(/' Atom ',i3,3(' Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Green integral, Non magnetic part',80x,'Magnetic part')
  218 format(/' Atom ',i3,3(' Tensor_oo(ke,je,he',3(',',i1),')',36x),
     &        ' Green integral')
  219 format(/' Atom ',i3,3(' Tensor_oo(ke,je,he',3(',',i1),')',36x))

  283 format(/' Conversion factor (numb. of electron/Mbarn) =',10f10.5)
  285 format(/'   Total signal')
  290 format(/'   Signal atom',i3)
  295 format(/'   Core state or edge',i3)
  300 format(/4x,'Energy',320a13)
  310 format(f10.3,1p,320e13.5)
  320 format('     E1-E1',1p,320e13.5)
  330 format('     E1-E2',1p,320e13.5)
  340 format('     E2-E2',1p,320e13.5)
  350 format('     E1-E3',1p,320e13.5)
  351 format('     E3-E3',1p,320e13.5)
  352 format('     M1-M1',1p,320e13.5)
  354 format('     M1-E1',1p,320e13.5)
  360 format(/10x,1p,320a13)
  370 format('  Ampldafs',1p,320e13.5)
  392 format(' Number of Energies =',i5)
  400 format(i5,4x,' = Number of angles')
  405 format(f10.3,A)
  407 format(3i5,' = (h,k,l)')
  410 format(' (h,k,l) = ',3i3,', sigpi =',2i3)
  420 format(7x,1p,320e13.5)
  430 format(f7.1,1p,320e13.5)
      end

!***********************************************************************

! Calcul du facteur de conversion Mbarn --> nbr. d'elec (divise par pi)

      function conv_mbarn_nelec(Ephoton)

      use declarations
      implicit real(kind=db) (a-h,o-z)

! Calcul de la constante multiplicative.
  ! ptrans = S02 fixe a 1.
      ptrans = 1
  ! alfa_sf = e*e/(2*epsilon0*h*c) est la constante de structure fine.
      cst = quatre_pi * pi * ptrans * alfa_sf * Ephoton
  ! pour avoir le resultat en megabarn (10E-18 cm2)
      cst = 100 * bohr**2 * cst

      eph2 = 0.5 * Ephoton**2
! Constante multiplication pour avoir le resultat en nombre d'electron
      conv_mbarn_nelec = eph2 / cst

      return
      end

!***********************************************************************

      subroutine write_cartesian_tensor(Densite_atom,E_cut,E1E2,E2E2,
     &               Ephseuil,
     &               Epsii,Eseuil,ia,ie,ipldafs,jseuil,Length_word,
     &               M1M1,magn_sens,natomsym,ninit1,ninitlr,nomfich_s,
     &               nseuil,numat_abs,secddia,secdqia,secdqia_m,
     &               secqqia,secmdia,tens_comp,v0muf,Core_resolved)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      parameter( n_dim=10*(168+12) )

      character(len=132) nomficht, nomfich_s
      character(len=Length_word) mot
      character(len=Length_word), dimension(n_dim):: nomtens

      complex(kind=db), dimension(1):: cdum
      complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia,
     &                                                      secmdia
      complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia,
     &                                             secdqia_m
      complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia

      integer, dimension(0):: idum

      logical Core_resolved, E1E2, E2E2, M1M1, magn_sens, tens_comp

      real(kind=db), dimension(0):: rdum
      real(kind=db), dimension(ninitlr):: Epsii
      real(kind=db), dimension(n_dim):: Tens

      nomficht = nomfich_s
      long = len_trim(nomficht)
      if( ia == 0 ) then
        nomficht(long+1:long+9) = '_car_xtal'
        if( ipldafs > 0 ) then
          nomficht(long+10:long+13) = '_rxs'
          call ad_number(ipldafs,nomficht,132)
        endif
      else
        nomficht(long+1:long+9) = '_car_atom'
        call ad_number(ia,nomficht,132)
      endif
      long = len_trim(nomficht)
      nomficht(long+1:long+4) = '.txt'

      index = 0

      do initlr = 1,ninitlr

        mot = ' '
        mot(1:2) = 'D_'
        do i = 1,3
          mot(3:3) = achar(i+119) 
          do j = i,3
            index = index + 1
            mot(4:4) = achar(j+119) 
            if( tens_comp ) mot(5:6) = '_r'
            Tens(index) = real( secddia(i,j,initlr,ia),db ) 
            nomtens(index) = mot
            if( tens_comp ) then
              index = index + 1 
              mot(6:6) = 'i'
              nomtens(index) = mot
              Tens(index) = aimag( secddia(i,j,initlr,ia) ) 
            endif 
          end do
        end do
        if( E1E2 ) then
          mot(1:2) = 'I_'
          do i = 1,3
            mot(3:3) = achar(i+119) 
            do j = 1,3
              mot(4:4) = achar(j+119) 
              do k = j,3
                mot(5:5) = achar(k+119)
                index = index + 1 
                if( tens_comp ) mot(6:7) = '_r'
                nomtens(index) = mot
                Tens(index) = real( secdqia(i,j,k,initlr,ia),db ) 
                if( tens_comp ) then
                  index = index + 1 
                  mot(7:7) = 'i'
                  nomtens(index) = mot
                  Tens(index) = aimag( secdqia(i,j,k,initlr,ia) ) 
                endif
                if( ia == 0 .and. magn_sens ) then
                  index = index+1
                  mot(7:11) = 'r_mag'
                  nomtens(index) = mot
                  Tens(index) = real( secdqia_m(i,j,k,initlr,ia),db ) 
                  index = index+1
                  mot(7:7) = 'i'
                  nomtens(index) = mot
                  mot(7:11) = '    '
                  Tens(index) = aimag( secdqia_m(i,j,k,initlr,ia) ) 
                endif 
              end do
            end do
          end do
        endif
        if( E2E2 ) then
          mot = ' '
          mot(1:2) = 'Q_'
          do i = 1,3
            mot(3:3) = achar(i+119)
            do j = i,3
              mot(4:4) = achar(j+119)
              do k = 1,3
                mot(5:5) = achar(k+119)
                do l = k,3
                  mot(6:6) = achar(l+119)
                  index = index + 1
                  if( tens_comp ) mot(7:8) = '_r'
                  nomtens(index) = mot
                  Tens(index) = real( secqqia(i,j,k,l,initlr,ia),db ) 
                  if( tens_comp ) then
                    index = index + 1 
                    mot(8:8) = 'i'
                    nomtens(index) = mot
                    Tens(index) = aimag( secqqia(i,j,k,l,initlr,ia) ) 
                  endif 
                end do
              end do
            end do
          end do
          do ijk = 1,6
            index = index + 1
            select case(ijk)
              case(1)
                i = 1; j = 2; k = 1; l = 2 
              case(2)
                i = 1; j = 2; k = 1; l = 3 
              case(3)
                i = 1; j = 3; k = 1; l = 3 
              case(4)
                i = 1; j = 3; k = 2; l = 2 
              case(5)
                i = 1; j = 3; k = 2; l = 3 
              case(6)
                i = 2; j = 3; k = 2; l = 3 
            end select
            mot(3:3) = achar(i+119)
            mot(4:4) = achar(j+119)
            mot(5:5) = achar(k+119)
            mot(6:6) = achar(l+119)
            if( tens_comp ) mot(7:8) = '_r'
            nomtens(index) = mot
            Tens(index) = real( secqqia(i,j,k,l,initlr,ia),db ) 
            if( tens_comp ) then
              index = index + 1 
              mot(8:8) = 'i'
              nomtens(index) = mot
              Tens(index) = aimag( secqqia(i,j,k,l,initlr,ia) ) 
            endif 
          end do
        endif
        if( M1M1 ) then
          mot = ' '
          mot(1:2) = 'M_'
          do i = 1,3
            mot(3:3) = achar(i+119) 
            do j = i,3
              index = index + 1
              mot(4:4) = achar(j+119) 
              if( tens_comp ) mot(5:6) = '_r'
              Tens(index) = real( secmdia(i,j,initlr,ia),db ) 
              nomtens(index) = mot
              if( tens_comp ) then
                index = index + 1 
                mot(6:6) = 'i'
                nomtens(index) = mot
                Tens(index) = aimag( secmdia(i,j,initlr,ia) ) 
              endif 
            end do
          end do
        endif
        n_tens = index

      end do

      call write_out(rdum,rdum,Densite_atom,0._db,E_cut,Ephseuil,
     &         Epsii,Eseuil,.false.,idum,ie,Length_word,
     &         jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,nomtens,
     &         1,0,0,0,nseuil,numat_abs,cdum,
     &         cdum,tens,v0muf,Core_resolved,0)

      return
      end

!***********************************************************************

      subroutine write_out(angxyz,axyz,Densite_atom,fpp_avantseuil,
     &          E_cut,Ephseuil,
     &          Epsii,Eseuil,Green_int,hkl_dafs,ie,Length_word,
     &          jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,
     &          np,npp,nppa,npps,nseuil,numat,ph1,
     &          ph2,tens,v0muf,Core_resolved,natomsym)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      parameter(n_tens_max = 10000)

      integer, dimension(3,npps):: hkl_dafs
      character(len=132) nomficht
      character(len=Length_word):: mot
      character(len=Length_word), dimension(n_dim):: title
      character(len=10+(n_tens-2*npp)*Length_word):: dummy

      complex(kind=db), dimension(np):: ph1, ph2

      logical Core_resolved, Green_int

      real(kind=db), dimension(nppa):: angxyz, axyz
      real(kind=db), dimension(ninitlr):: Epsii
      real(kind=db), dimension(n_dim):: Tens

      dummy = ' '

      if( numat == -1 ) then
        ipr = 6
      elseif( numat < -1 ) then
        ipr = - numat 
      else
        ipr = 4
      endif

      if( ipr == 6 ) then
        n = min(n_tens,4)
      else
        n = n_tens
      endif

      if( n_tens > n_tens_max ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,105) n_tens
        end do
        stop
      endif

      if( ie == 1 ) then
        if( numat >= 0 ) then
          open(ipr, file = nomficht)
        elseif( numat < -1 ) then
          open(ipr, status='scratch')
        endif
        if( numat > 0 ) then
          if( Core_resolved ) then
            icor = ninit1
          else
            icor = 1
          endif
          if( Green_int ) icor = - icor

          if( ninitlr == 1 ) then
            write(ipr,110) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          elseif( ninitlr == 2 ) then
            write(ipr,111) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          elseif( ninitlr == 4 ) then
            write(ipr,112) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          elseif( ninitlr == 6 ) then
            write(ipr,113) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          elseif( ninitlr == 8 ) then
            write(ipr,114) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          elseif( ninitlr == 10 ) then
            write(ipr,115) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          else     ! 14
            write(ipr,116) Eseuil*rydb, numat, nseuil, 
     &        jseuil, fpp_avantseuil, v0muf*rydb,
     &        E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom
          endif
        endif
         
        if( npp > 0 .and. numat >= 0 ) then
          select case(Length_word)
            case(11)
              write(ipr,121) dummy, ph2(1:npp)
              write(ipr,121) dummy, ph1(1:npp)
            case(12)
              write(ipr,122) dummy, ph2(1:npp)
              write(ipr,122) dummy, ph1(1:npp)
            case(13)
              write(ipr,123) dummy, ph2(1:npp)
              write(ipr,123) dummy, ph1(1:npp)
            case(14)
              write(ipr,124) dummy, ph2(1:npp)
              write(ipr,124) dummy, ph1(1:npp)
            case(15)
              write(ipr,125) dummy, ph2(1:npp)
              write(ipr,125) dummy, ph1(1:npp)
            case(16)
              write(ipr,126) dummy, ph2(1:npp)
              write(ipr,126) dummy, ph1(1:npp)
            case(17)
              write(ipr,127) dummy, ph2(1:npp)
              write(ipr,127) dummy, ph1(1:npp)
            case default
              call write_error
              do ipr = 6,9,3
                write(ipr,130) Length_word
              end do
              stop
          end select
          if( nppa > 0 ) write(ipr,135) natomsym, axyz(:)*bohr,
     &                         angxyz(:), ( hkl_dafs(:,i), i = 1,npp )
        endif
        do i = 1,n
          mot = title(i) 
          call center_word( mot, Length_word )
          title(i) = mot 
        end do
        write(ipr,140) title(1:n)
      elseif( numat >= 0 ) then
        open(ipr, file = nomficht, position='append')
      endif

      select case(Length_word)
        case(11)
          write(ipr,151) Ephseuil*rydb, Tens(1:n)
        case(12)
          write(ipr,152) Ephseuil*rydb, Tens(1:n)
        case(13)
          write(ipr,153) Ephseuil*rydb, Tens(1:n)
        case(14)
          write(ipr,154) Ephseuil*rydb, Tens(1:n)
        case(15)
          write(ipr,155) Ephseuil*rydb, Tens(1:n)
        case(16)
          write(ipr,156) Ephseuil*rydb, Tens(1:n)
        case(17)
          write(ipr,157) Ephseuil*rydb, Tens(1:n)
        case default
          call write_error
          do ipr = 6,9,3
            write(ipr,130) Length_word
          end do
          stop
      end select

      if( numat >= 0 ) close(ipr)

      return
  105 format(//' The number of column to be written is ',i5,//
     &  ' This is greater than the maximum possible value',
     &  ' given in the routine write_out in the file coabs.f !',/
     &  ' To change that you must modify the formats 121 up to',
     &  ' 157 in this routine,'/,
     &  ' and increase to the same value the parameter n_tens_max.',/
     &  '  Then you compile again.'//)
  110 format(f10.3,i5,2i3,3f12.5,2i3,f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii, Atom_density')
  111 format(f10.3,i5,2i3,3f12.5,2i3,2f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii(1), Epsii(2), Atom_density')
  112 format(f10.3,i5,2i3,3f12.5,2i3,4f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii(1...4), Atom_density')
  113 format(f10.3,i5,2i3,3f12.5,2i3,6f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii(1...6), Atom_density')
  114 format(f10.3,i5,2i3,3f12.5,2i3,8f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii(1...8), Atom_density')
  115 format(f10.3,i5,2i3,3f12.5,2i3,10f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii(1...10), Atom_density')
  116 format(f10.3,i5,2i3,3f12.5,2i3,14f13.3,f13.9,
     &' = E_edge, Z, n_edge, j_edge,',
     &' Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1,',
     &' Epsii(1...14), Atom_density')
  121 format(A,1p,10000e11.3)
  122 format(A,1p,10000e12.4)
  123 format(A,1p,10000e13.5)
  124 format(A,1p,10000e14.6)
  125 format(A,1p,10000e15.7)
  126 format(A,1p,10000e16.8)
  127 format(A,1p,10000e17.9)
  130 format(//' Length_word =',i3,
     &         ' This parameter must be set between 11 and 17 !'//)
  135 format(i4,3f10.5,3x,3f10.5,10000(14x,3i4))
  140 format('    Energy',10000A)
  151 format(f10.3,1p,10000e11.3)
  152 format(f10.3,1p,10000e12.4)
  153 format(f10.3,1p,10000e13.5)
  154 format(f10.3,1p,10000e14.6)
  155 format(f10.3,1p,10000e15.7)
  156 format(f10.3,1p,10000e16.8)
  157 format(f10.3,1p,10000e17.9)
      end

!***********************************************************************

      subroutine center_word( mot, Length_word )

      character(len=*) mot
      character(len=Length_word) mot2

      mot2 = ' '
      mot = adjustl( mot )
      l = len_trim( mot )
      lshift = ( Length_word - l + 1 ) / 2
      lshift = max( 0, lshift )
      lm = min( l, Length_word ) 
      mot2(1+lshift:lm+lshift) = mot(1:lm)
      mot = mot2

      return
      end

!***********************************************************************

      subroutine ad_number(ib,nomfich,Length)

      character(len=Length) nomfich

      l = len_trim(nomfich)

      if( ib < 0 ) then
        l = l + 1
        if( l <= Length ) nomfich(l:l) = '-' 
      endif

      i = abs(ib)

      do iu = 2,10
        if( i / 10**(iu-1) < 1 ) exit
      end do
      iumax = iu - 1

      do iu = iumax,1,-1

        ipuis = 10**(iu-1)

        in = i / ipuis

! S'il n'y a pas la place on ecrit que les derniers chiffres
        if( l + iu <= Length ) then
          l = l + 1
          nomfich(l:l) = achar(in+48)
        elseif( l + iu == Length-1 .and. nomfich(l:l) == '_' ) then
          nomfich(l:l) = achar(in+48)
        endif

        i = i - ipuis * in

      end do

      return
      end

!***********************************************************************

      subroutine spherical_tensor_cal(ct_nelec,Densite_atom,E_cut,
     &          E1E1,E1E2,E2E2,
     &          Energ,Ephseuil,Epsii,Eseuil,ia,icheck,ie,Int_tens,
     &          kpl,ipldafs,jseuil,Length_word,magn_sens,moyenne,
     &          natomsym,ncolm,nenerg,ninitlr,nomfich_s,npldafs,nplr,
     &          nplrm,nplt,nseuil,numat_abs,pdp,phdf0t,phdt,plae,pol,
     &          plas,secddia,secdqia,secdqia_m,secqqia,v0muf,voae,
     &          vec,voas)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      parameter( n_tens_dd=9, n_tens_dq=15, n_tens_qq=25,
     &           n_tens_t = n_tens_dd + n_tens_dq + n_tens_qq,
     &           n_tens_max = 8 + 2 * n_tens_t + 2 * n_tens_dq ) 

      character(len=132) nomfich_s, nomficht

      complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
      complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq,
     &                                        Sph_tensor_dq_m
      complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq
      complex(kind=db), dimension(3,3):: secdd
      complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia
      complex(kind=db), dimension(3,3,3):: secdq
      complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia,
     &                                                    secdqia_m
      complex(kind=db), dimension(3):: plae, plas
      complex(kind=db), dimension(3,3,3,3):: secqq
      complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia
      complex(kind=db), dimension(3,nplrm):: pol
      complex(kind=db), dimension(0:nplt,n_tens_dd):: Tensor_pol_dd
      complex(kind=db), dimension(0:nplt,2*n_tens_dq):: Tensor_pol_dq
      complex(kind=db), dimension(0:nplt,n_tens_qq):: Tensor_pol_qq
      complex(kind=db), dimension(npldafs):: phdf0t, phdt 

      logical E1E1, E1E2, E2E2,
     &        magn_sens, moyenne, writbav, writout 

      real(kind=db), dimension(3):: voae, voas
      real(kind=db), dimension(ninitlr) :: ct_nelec, Epsii
      real(kind=db), dimension(nenerg) :: Energ
      real(kind=db), dimension(3,nplrm):: vec
      real(kind=db), dimension(ncolm,2) :: pdp
      real(kind=db), dimension(n_tens_max,ninitlr,0:natomsym):: Int_tens

      if( kpl == 1 ) then
        ipl1 = 1
        ipl2 = nplr
      else
        ipl1 = ipldafs
        ipl2 = ipldafs
      endif

      do ipl = ipl1,ipl2

        if( kpl == 1 ) then
          plae(:) = pol(:,ipl)
          plas(:) = pol(:,ipl)
          voae(:) = vec(:,ipl)
          voas(:) = vec(:,ipl)
        endif

        if( E1E1 ) call Tensor_pol_dd_cal(ipl,n_tens_dd,nplt,
     &                             plae,plas,Tensor_pol_dd)

        if( E1E2 ) call Tensor_pol_dq_cal(ipl,n_tens_dq,nplt,
     &                             plae,plas,voae,voas,Tensor_pol_dq)

        if( E2E2 ) call Tensor_pol_qq_cal(ipl,n_tens_qq,nplt,
     &                             plae,plas,voae,voas,Tensor_pol_qq)
      end do

      if( kpl == 1 .and. moyenne ) then
        ipl0 = 0
        if( E1E1 ) then
          do i = 1,n_tens_dd 
            Tensor_pol_dd(0,i) = sum( pdp(ipl1:ipl2,1)
     &                              * Tensor_pol_dd(ipl1:ipl2,i) )
          end do
        endif
        if( E1E2 ) then
          do i = 1,2*n_tens_dq 
            Tensor_pol_dq(0,i) = (0._db,0._db)
          end do
        endif
        if( E2E2 ) then 
          do i = 1,n_tens_qq 
            Tensor_pol_qq(0,i) = sum( pdp(ipl1:ipl2,2)
     &                              * Tensor_pol_qq(ipl1:ipl2,i) )
          end do
        endif
      else
        ipl0 = ipl1
      endif

      if( icheck > 1 .and. ie == 1 ) then 
        if( E1E1 ) then 
          if( ipldafs > 0 ) then
            write(3,110) 'Dipole-dipole', ipldafs
          else
            if( ipl0 == 0 ) then
              write(3,120) 'Dipole-dipole'
            else
              write(3,125) 'Dipole-dipole'
            endif
          endif 
          do i = 1,n_tens_dd 
            write(3,130) i, Tensor_pol_dd(ipl0:ipl2,i)
          end do
        endif

        if( E1E2 ) then 
          if( ipldafs > 0 ) then
            write(3,110) 'Dipole-quadrupole', ipldafs
          else
            if( ipl0 == 0 ) then
              write(3,120) 'Dipole-quadrupole'
            else
              write(3,125) 'Dipole-quadrupole'
            endif
          endif 
          do i = 1,2*n_tens_dq 
            write(3,130) i, Tensor_pol_dq(ipl0:ipl2,i)
          end do
        endif

        if( E2E2 ) then 
          if( ipldafs > 0 ) then
            write(3,110) 'Quadrupole-quadrupole', ipldafs
          else
            if( ipl0 == 0 ) then
              write(3,120) 'Quadrupole-quadrupole'
            else
              write(3,125) 'Quadrupole-quadrupole'
            endif
          endif 
          do i = 1,n_tens_qq 
            write(3,130) i, Tensor_pol_qq(ipl0:ipl2,i)
          end do
        endif
      endif

      writout = ia == 0 .or. kpl == 1
      writbav = icheck > 1 .and. writout

      if( writbav ) then
        if( ia == 0 ) then
          if( kpl == 1 ) then
            write(3,142)
          else
            write(3,143) ipldafs
          endif
        else 
          write(3,144) ia
        endif
      endif

      do initlr = 1,ninitlr

        if( writbav ) write(3,145) initlr 

        if( E1E1 ) then
          secdd(:,:) = secddia(:,:,initlr,ia)
          call Sph_tensor_dd_cal(n_tens_dd,secdd,Sph_tensor_dd)
          Sph_tensor_dd(:) = ct_nelec(initlr) * Sph_Tensor_dd(:)
          if( writbav ) then
            write(3,148) 
            if( ipl0 <= 1 ) then
              write(3,150) Real( Sph_tensor_dd(:) )
            else
              write(3,155) Sph_tensor_dd(:)
            endif 
          endif
        endif

        if( E1E2 ) then
          secdq(:,:,:) = secdqia(:,:,:,initlr,ia)
          call Sph_tensor_dq_cal(n_tens_dq,secdq,Sph_tensor_dq)
          Sph_tensor_dq(:) = ct_nelec(initlr) * Sph_Tensor_dq(:)
          if( ia == 0 .and. magn_sens ) then
            secdq(:,:,:) = secdqia_m(:,:,:,initlr,ia)
            call Sph_tensor_dq_cal(n_tens_dq,secdq,Sph_tensor_dq_m)
            Sph_tensor_dq_m(:) = ct_nelec(initlr)
     &                           * Sph_Tensor_dq_m(:)
          endif
          if( writbav ) then
            write(3,158) 
            if( ipl0 <= 1 .and. ia == 0 .and. magn_sens ) then
              write(3,165) ( Real( Sph_tensor_dq(i) ), 
     &                       Real( Sph_tensor_dq_m(i) ), i = 1,3 )
              write(3,175) ( Real( Sph_tensor_dq(i) ),
     &                       Real( Sph_tensor_dq_m(i) ), i = 4,8 ) 
              write(3,185) ( Real( Sph_tensor_dq(i) ),
     &                       Real( Sph_tensor_dq_m(i) ), i = 9,15 )
            elseif( ia == 0 .and. magn_sens ) then
              write(3,166) ( Sph_tensor_dq(i), Sph_tensor_dq_m(i),
     &                       i = 1,3 )
              write(3,176) ( Sph_tensor_dq(i), Sph_tensor_dq_m(i),
     &                       i = 4,8 ) 
              write(3,186) ( Sph_tensor_dq(i), Sph_tensor_dq_m(i),
     &                       i = 9,15 )
            elseif( ipl0 <= 1 ) then
              write(3,160) Real( Sph_tensor_dq(1:3) )
              write(3,170) Real( Sph_tensor_dq(4:8) ) 
              write(3,180) Real( Sph_tensor_dq(9:15) )
            else
              write(3,162) Sph_tensor_dq(1:3)
              write(3,172) Sph_tensor_dq(4:8) 
              write(3,182) Sph_tensor_dq(9:15)
            endif
          endif
        endif

        if( E2E2 ) then
          secqq(:,:,:,:) = secqqia(:,:,:,:,initlr,ia)
          call Sph_tensor_qq_cal(n_tens_qq,secqq,Sph_tensor_qq)
          Sph_tensor_qq(:) = ct_nelec(initlr) * Sph_Tensor_qq(:)
          if( writbav ) then
            write(3,188) 
            if( ipl0 <= 1 ) then
              write(3,190) Real( Sph_tensor_qq(1:9) ) 
              write(3,200) Real( Sph_tensor_qq(10:n_tens_qq) )
            else
              write(3,210) Sph_tensor_qq(1:9) 
              write(3,220) Sph_tensor_qq(10:n_tens_qq)
            endif
          endif
        endif

        nomficht = nomfich_s
        if( ninitlr > 1 ) then
          long = len_trim(nomficht)
          nomficht(long+1:long+4) = '_g'
          long = long + 2
          call ad_number(initlr,nomficht,132)
        endif

        call write_phys(ct_nelec(initlr),Densite_atom,E_cut,E1E1,E1E2,
     &    E2E2,Energ,Ephseuil,Epsii,
     &    Eseuil,ia,ie,initlr,Int_tens,ipl0,ipl2,ipldafs,Length_word,
     &    jseuil,magn_sens,n_tens_dd,n_tens_dq,n_tens_max,
     &    n_tens_qq,n_tens_t,natomsym,nenerg,ninitlr,
     &    nomficht,npldafs,nplt,nseuil,numat_abs,phdf0t,phdt,
     &    Sph_tensor_dd,Sph_tensor_dq,Sph_tensor_dq_m,
     &    Sph_tensor_qq,Tensor_pol_dd,Tensor_pol_dq,Tensor_pol_qq,
     &    v0muf,writout)

      end do

      return
  110 format(/1x,A,' polarisation tensor for the RXS:',
     &       /'  i ',5x,' Reflection number',i3)
  120 format(/1x,A,' polarisation tensor for the xanes:',
     &     /'  i        <xanes>           ipl = 0, 1, 2...')
  125 format(/1x,A,' polarisation tensor for the xanes:',
     &     /'  i        ipl = 0, 1, 2...')
  130 format(i3,25(1x,2f9.5))
  142 format(/' Spherical tensors for the unit cell')
  143 format(/' Spherical tensors for the unit cell for the RXS',
     &  ' (num. of electron), reflection number',i3)
  144 format(/' Spherical tensors (numb. of electron) for the atom',
     &        ' number :',i3)
  145 format(/'    initlr =',i3)
  148 format(/' Dipole-dipole spherical tensor (numb. of',
     &        ' electron):')
  150 format(/1p,
     & ' rank 0, non-magnetic scalar :',/
     & '    D(00)                      =',e13.5,//
     & ' rank 1, magnetic dipole :',/
     & '    D(10)                 = lz =',e13.5,/
     & '   (D(11)-D(1-1))/sqrt(2) =-lx =',e13.5,/
     & ' -i(D(11)+D(1-1))/sqrt(2) = ly =',e13.5,//
     & ' rank 2, non-magnetic quadrupole :',/
     & '    D(20)                      =',e13.5,/
     & '   (D(21)-D(2-1))/sqrt(2)      =',e13.5,/
     & ' -i(D(21)+D(2-1))/sqrt(2)      =',e13.5,/
     & ' -i(D(22)-D(2-2))/sqrt(2)      =',e13.5,/
     & '   (D(22)+D(2-2))/sqrt(2)      =',e13.5)
  155 format(1p,
     & ' rank 0, non-magnetic scalar :',/
     & '    D(00)                      =',2e13.5,//
     & ' rank 1, magnetic dipole :',/
     & '    D(10)                 = lz =',2e13.5,/
     & '   (D(11)-D(1-1))/sqrt(2) =-lx =',2e13.5,/
     & ' -i(D(11)+D(1-1))/sqrt(2) = ly =',2e13.5,//
     & ' rank 2, non-magnetic quadrupole :',/
     & '    D(20)                      =',2e13.5,/
     & '   (D(21)-D(2-1))/sqrt(2)      =',2e13.5,/
     & ' -i(D(21)+D(2-1))/sqrt(2)      =',2e13.5,/
     & ' -i(D(22)-D(2-2))/sqrt(2)      =',2e13.5,/
     & '   (D(22)+D(2-2))/sqrt(2)      =',2e13.5)
  158 format(/' Dipole-quadrupole spherical tensor (numb. of',
     &        ' electron):')
  160 format(/42x,'non-magnetic ',/ ' rank 1 :',1p,/
     & '    I(10)                          = nz =',e13.5,/
     & '   (I(11)-I(1-1))/sqrt(2)          = nx =',e13.5,/
     & ' -i(I(11)+I(1-1))/sqrt(2)          = ny =',e13.5)
  162 format(/42x,'non-magnetic ',/ ' rank 1 :',1p,/
     & '    I(10)                          = nz =',2e13.5,/
     & '   (I(11)-I(1-1))/sqrt(2)          = nx =',2e13.5,/
     & ' -i(I(11)+I(1-1))/sqrt(2)          = ny =',2e13.5)
  165 format(/37x,'non-magnetic ',25x,'magnetic',/' rank 1 :',1p,/
     & '    I(10)                          = nz =',e13.5,
     &                              '              Toroiz =',e13.5,/
     & '   (I(11)-I(1-1))/sqrt(2)          = nx =',e13.5,
     &                              '              Toroix =',e13.5,/
     & ' -i(I(11)+I(1-1))/sqrt(2)          = ny =',e13.5,
     &                              '              Toroiy =',e13.5)
  166 format(/49x,'non-magnetic ',37x,'magnetic',/' rank 1 :',1p,/
     & '    I(10)                          = nz =',2e13.5,
     &                              '              Toroiz =',2e13.5,/
     & '   (I(11)-I(1-1))/sqrt(2)          = nx =',2e13.5,
     &                              '              Toroix =',2e13.5,/
     & ' -i(I(11)+I(1-1))/sqrt(2)          = ny =',2e13.5,
     &                              '              Toroiy =',2e13.5)
  170 format(/' rank 2 :',1p,/
     & '  -iI(20)                  =  lz*Toroiz =',e13.5,/
     & ' -i(I(21)-I(2-1))/sqrt(2)  = (l,Toroi)2 =',e13.5,/
     & '   (I(21)+I(2-1))/sqrt(2)  = (l,Toroi)2 =',e13.5,/
     & '   (I(22)-I(2-2))/sqrt(2)  = (l,Toroi)2 =',e13.5,/
     & ' -i(I(22)+I(2-2))/sqrt(2)  = (l,Toroi)2 =',e13.5)
  172 format(/' rank 2 :',1p,/
     & '  -iI(20)                  =  lz*Toroiz =',2e13.5,/
     & ' -i(I(21)-I(2-1))/sqrt(2)  = (l,Toroi)2 =',2e13.5,/
     & '   (I(21)+I(2-1))/sqrt(2)  = (l,Toroi)2 =',2e13.5,/
     & '   (I(22)-I(2-2))/sqrt(2)  = (l,Toroi)2 =',2e13.5,/
     & ' -i(I(22)+I(2-2))/sqrt(2)  = (l,Toroi)2 =',2e13.5)
  175 format(/' rank 2 :',1p,/
     & '  -iI(20)                  =  lz*Toroiz =',e13.5,
     &                              '               nz*lz =',e13.5,/
     & ' -i(I(21)-I(2-1))/sqrt(2)  = (l,Toroi)2 =',e13.5,
     &                              '              (n,l)2 =',e13.5,/
     & '   (I(21)+I(2-1))/sqrt(2)  = (l,Toroi)2 =',e13.5,
     &                              '              (n,l)2 =',e13.5,/
     & '   (I(22)-I(2-2))/sqrt(2)  = (l,Toroi)2 =',e13.5,
     &                              '              (n,l)2 =',e13.5,/
     & ' -i(I(22)+I(2-2))/sqrt(2)  = (l,Toroi)2 =',e13.5,
     &                              '              (n,l)2 =',e13.5)
  176 format(/' rank 2 :',1p,/
     & '  -iI(20)                  =  lz*Toroiz =',2e13.5,
     &                              '               nz*lz =',2e13.5,/
     & ' -i(I(21)-I(2-1))/sqrt(2)  = (l,Toroi)2 =',2e13.5,
     &                              '              (n,l)2 =',2e13.5,/
     & '   (I(21)+I(2-1))/sqrt(2)  = (l,Toroi)2 =',2e13.5,
     &                              '              (n,l)2 =',2e13.5,/
     & '   (I(22)-I(2-2))/sqrt(2)  = (l,Toroi)2 =',2e13.5,
     &                              '              (n,l)2 =',2e13.5,/
     & ' -i(I(22)+I(2-2))/sqrt(2)  = (l,Toroi)2 =',2e13.5,
     &                              '              (n,l)2 =',2e13.5)
  180 format(/' rank 3 :',1p,/
     & '    I(30)              = nz*(3lz2 - l2) =',e13.5,/
     & '   (I(31)-I(3-1))/sqrt(2) = (n,(l,l)2)3 =',e13.5,/
     & ' -i(I(31)+I(3-1))/sqrt(2) = (n,(l,l)2)3 =',e13.5,/
     & ' -i(I(32)-I(3-2))/sqrt(2) = (n,(l,l)2)3 =',e13.5,/
     & '   (I(32)+I(3-2))/sqrt(2) = (n,(l,l)2)3 =',e13.5,/
     & '   (I(33)-I(3-3))/sqrt(2) = (n,(l,l)2)3 =',e13.5,/
     & ' -i(I(33)+I(3-3))/sqrt(2) = (n,(l,l)2)3 =',e13.5)
  182 format(/' rank 3 :',1p,/
     & '    I(30)              = nz*(3lz2 - l2) =',2e13.5,/
     & '   (I(31)-I(3-1))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,/
     & ' -i(I(31)+I(3-1))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,/
     & ' -i(I(32)-I(3-2))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,/
     & '   (I(32)+I(3-2))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,/
     & '   (I(33)-I(3-3))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,/
     & ' -i(I(33)+I(3-3))/sqrt(2) = (n,(l,l)2)3 =',2e13.5)
  185 format(/' rank 3 :',1p,/
     & '    I(30)              = nz*(3lz2 - l2) =',e13.5,
     &                              '  Toroiz*(3lz2 - l2) =',e13.5,/
     & '   (I(31)-I(3-1))/sqrt(2) = (n,(l,l)2)3 =',e13.5,
     &                              '     (Toroi,(l,l)2)3 =',e13.5,/
     & ' -i(I(31)+I(3-1))/sqrt(2) = (n,(l,l)2)3 =',e13.5,
     &                              '     (Toroi,(l,l)2)3 =',e13.5,/
     & ' -i(I(32)-I(3-2))/sqrt(2) = (n,(l,l)2)3 =',e13.5,
     &                              '     (Toroi,(l,l)2)3 =',e13.5,/
     & '   (I(32)+I(3-2))/sqrt(2) = (n,(l,l)2)3 =',e13.5,
     &                              '     (Toroi,(l,l)2)3 =',e13.5,/
     & '   (I(33)-I(3-3))/sqrt(2) = (n,(l,l)2)3 =',e13.5,
     &                              '     (Toroi,(l,l)2)3 =',e13.5,/
     & ' -i(I(33)+I(3-3))/sqrt(2) = (n,(l,l)2)3 =',e13.5,
     &                              '     (Toroi,(l,l)2)3 =',e13.5)
  186 format(/' rank 3 :',1p,/
     & '    I(30)              = nz*(3lz2 - l2) =',2e13.5,
     &                              '  Toroiz*(3lz2 - l2) =',2e13.5,/
     & '   (I(31)-I(3-1))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,
     &                              '     (Toroi,(l,l)2)3 =',2e13.5,/
     & ' -i(I(31)+I(3-1))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,
     &                              '     (Toroi,(l,l)2)3 =',2e13.5,/
     & ' -i(I(32)-I(3-2))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,
     &                              '     (Toroi,(l,l)2)3 =',2e13.5,/
     & '   (I(32)+I(3-2))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,
     &                              '     (Toroi,(l,l)2)3 =',2e13.5,/
     & '   (I(33)-I(3-3))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,
     &                              '     (Toroi,(l,l)2)3 =',2e13.5,/
     & ' -i(I(33)+I(3-3))/sqrt(2) = (n,(l,l)2)3 =',2e13.5,
     &                              '     (Toroi,(l,l)2)3 =',2e13.5)
  188 format(/' Quadrupole-quadrupole spherical tensor (numb. of',
     &        ' electron) :')
  190 format(/' rank 0, scalar :',1p,/
     & '    Q(00)                      =',e13.5,//
     & ' rank 1, magnetic dipole :',/
     & '    Q(10)                 = lz =',e13.5,/
     & '   (Q(11)-Q(1-1))/sqrt(2) =-lx =',e13.5,/
     & ' -i(Q(11)+Q(1-1))/sqrt(2) = ly =',e13.5,//
     & ' Rank 2, non-magnetic quadrupole :',/
     & '    Q(20)                      =',e13.5,/
     & '   (Q(21)-Q(2-1))/sqrt(2)      =',e13.5,/
     & ' -i(Q(21)+Q(2-1))/sqrt(2)      =',e13.5,/
     & ' -i(Q(22)-Q(2-2))/sqrt(2)      =',e13.5,/
     & '   (Q(22)+Q(2-2))/sqrt(2)      =',e13.5,/)
  200 format(' Rank 3, magnetic octupole :',/
     & '    Q(30)                      =',e13.5,/
     & '   (Q(31)-Q(3-1))/sqrt(2)      =',e13.5,/
     & ' -i(Q(31)+Q(3-1))/sqrt(2)      =',e13.5,/
     & ' -i(Q(32)-Q(3-2))/sqrt(2)      =',e13.5,/
     & '   (Q(32)+Q(3-2))/sqrt(2)      =',e13.5,/
     & ' -i(Q(33)-Q(3-3))/sqrt(2)      =',e13.5,/
     & ' -i(Q(33)+Q(3-3))/sqrt(2)      =',e13.5,//
     & ' Rank 4, non-magnetic hexadecapole :',/
     & '    Q(40)                      =',e13.5,/
     & '   (Q(41)-Q(4-1))/sqrt(2)      =',e13.5,/
     & ' -i(Q(41)+Q(4-1))/sqrt(2)      =',e13.5,/
     & ' -i(Q(42)-Q(4-2))/sqrt(2)      =',e13.5,/
     & '   (Q(42)+Q(4-2))/sqrt(2)      =',e13.5,/
     & '   (Q(43)-Q(4-3))/sqrt(2)      =',e13.5,/
     & ' -i(Q(43)+Q(4-3))/sqrt(2)      =',e13.5,/
     & ' -i(Q(44)-Q(4-4))/sqrt(2)      =',e13.5,/
     & '   (Q(44)+Q(4-4))/sqrt(2)      =',e13.5)
  210 format(/1p,
     & ' rank 0, scalar :',/
     & '    Q(00)                      =',2e13.5,//
     & ' rank 1, magnetic dipole :',/
     & '    Q(10)                 = lz =',2e13.5,/
     & '   (Q(11)-Q(1-1))/sqrt(2) =-lx =',2e13.5,/
     & ' -i(Q(11)+Q(1-1))/sqrt(2) = ly =',2e13.5,//
     & ' Rank 2, non-magnetic quadrupole :',/
     & '    Q(20)                      =',2e13.5,/
     & '   (Q(21)-Q(2-1))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(21)+Q(2-1))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(22)-Q(2-2))/sqrt(2)      =',2e13.5,/
     & '   (Q(22)+Q(2-2))/sqrt(2)      =',2e13.5,/)
  220 format(
     & ' Rank 3, magnetic octupole :',/
     & '    Q(30)                      =',2e13.5,/
     & '   (Q(31)-Q(3-1))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(31)+Q(3-1))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(32)-Q(3-2))/sqrt(2)      =',2e13.5,/
     & '   (Q(32)+Q(3-2))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(33)-Q(3-3))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(33)+Q(3-3))/sqrt(2)      =',2e13.5,//
     & ' Rank 4, non-magnetic hexadecapole :',/
     & '    Q(40)                      =',2e13.5,/
     & '   (Q(41)-Q(4-1))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(41)+Q(4-1))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(42)-Q(4-2))/sqrt(2)      =',2e13.5,/
     & '   (Q(42)+Q(4-2))/sqrt(2)      =',2e13.5,/
     & '   (Q(43)-Q(4-3))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(43)+Q(4-3))/sqrt(2)      =',2e13.5,/
     & ' -i(Q(44)-Q(4-4))/sqrt(2)      =',2e13.5,/
     & '   (Q(44)+Q(4-4))/sqrt(2)      =',2e13.5)
      end

!***********************************************************************

      subroutine Sph_tensor_dd_cal(n_tens_dd,secdd,Sph_tensor_dd)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
      complex(kind=db), dimension(3,3):: secdd

! Tenseur 0

      Sph_tensor_dd(1) = ( 1 / sqrt(3._db) )
     &                 * ( secdd(1,1) + secdd(2,2) + secdd(3,3) )

! Tenseur 1
! Les composantes de ce tenseur sont en cas de seuil K : -lx, ly et lz.
      fac = 1 / sqrt( 2._db )

! Multiplie par - img, equivalent a prendre la partie imaginaire quand
! il n'y a pas de multiplication par le terme de Bragg.
! lz = D(01)
      Sph_tensor_dd(2) = - img * fac * ( secdd(1,2) - secdd(2,1) )

! -lx = (1/sqrt(2))*(D(11)-D(-11))
      Sph_tensor_dd(3) = img * fac * ( secdd(2,3) - secdd(3,2) )

! ly = (-i/sqrt(2))*(D(-11)+D(11))
      Sph_tensor_dd(4) = - img * fac * ( secdd(1,3) - secdd(3,1) )

! Tenseur 2

! D02 = (1/sqrt(6))*(2*Dzz-Dxx-Dyy)
      fac = 1 / sqrt( 6._db )
      Sph_tensor_dd(5) = fac
     &                 * ( 2*secdd(3,3) - secdd(1,1)  - secdd(2,2) )

      fac = 1 / sqrt( 2._db ) 

! (1/sqrt(2))*(D(12) - D(-12))
      Sph_tensor_dd(6) = - fac * ( secdd(1,3) + secdd(3,1) )

! (-i/sqrt(2))*(D(12) + D(-12))
      Sph_tensor_dd(7) = - fac * ( secdd(2,3) + secdd(3,2) )

! (-i/sqrt(2))*(D(22) - D(-22))
      Sph_tensor_dd(8) = fac * ( secdd(1,2) + secdd(2,1) )

! (1/sqrt(2))*(D(22) + D(-22))
      Sph_tensor_dd(9) = fac * ( secdd(1,1) - secdd(2,2) )

      return
      end

!***********************************************************************

      subroutine Sph_tensor_dq_cal(n_tens_dq,secdq,Sph_tensor_dq)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db), dimension(3,3,3):: secdq
      complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq

! Tenseur 1

  ! I(10)
      Sph_tensor_dq(1) = - ( 1 / sqrt(15._db) )
     &   * ( 3 * secdq(1,1,3) + 3 * secdq(2,2,3)
     &     + 2 * secdq(3,3,3) - secdq(3,1,1) - secdq(3,2,2) )

      fac = 2 / sqrt(60._db)

! (1/sqrt(2))*(I(11) - I(1-1))
      Sph_tensor_dq(2) = fac
     &   * ( 2 * secdq(1,1,1) - secdq(1,2,2) - secdq(1,3,3)
     &     + 3 * ( secdq(2,1,2) + secdq(3,1,3) ) )

! (-i/sqrt(2))*(I(11) + I(-1-1))
      Sph_tensor_dq(3) = fac
     &   * ( 2 * secdq(2,2,2) - secdq(2,1,1) - secdq(2,3,3)
     &     + 3 * ( secdq(1,1,2) + secdq(3,2,3) ) )

! Tenseur 2

  ! -i*I(20)
      Sph_tensor_dq(4) = secdq(1,2,3) - secdq(2,1,3)

      fac = 1 / sqrt(3._db)

  ! (-i/sqrt(2))*(I(21) - I(2-1))
      Sph_tensor_dq(5) = fac
     &   * ( secdq(2,1,1) - secdq(2,3,3) - secdq(1,1,2) + secdq(3,2,3) )

  ! (1/sqrt(2))*(I(21) + I(2-1))
      Sph_tensor_dq(6) = fac
     &   * ( secdq(1,2,2) - secdq(1,3,3) - secdq(2,1,2) + secdq(3,1,3) )

      Sph_tensor_dq(7) = fac
     &   * ( secdq(1,1,3) - secdq(2,2,3) - secdq(3,1,1) + secdq(3,2,2) )

  ! (-i/sqrt(2))*(I(22) + I(2-2))
      Sph_tensor_dq(8) = fac
     &   * ( secdq(1,2,3) + secdq(2,1,3) - secdq(3,1,2) - secdq(3,2,1) )

! Tenseur 3

      Sph_tensor_dq(9) = 1 / sqrt(10._db)
     &           * ( 2 * secdq(3,3,3) - 2 * secdq(1,1,3) - secdq(3,1,1) 
     &             - 2 * secdq(2,2,3) - secdq(3,2,2) )

      fac = 1 / sqrt(60._db)

      Sph_tensor_dq(10) = fac
     &           * ( 3 * secdq(1,1,1) + secdq(1,2,2) + 2 * secdq(2,1,2) 
     &             - 4 * secdq(1,3,3) - 8 * secdq(3,1,3) )

  ! (-i/sqrt(2))*(I(31) + I(3-1))
      Sph_tensor_dq(11) = fac
     &           * ( 3 * secdq(2,2,2) + secdq(2,1,1) + 2 * secdq(1,2,1) 
     &             - 4 * secdq(2,3,3) - 8 * secdq(3,2,3) )

      fac = 1 / sqrt(6._db)

  ! (-i/sqrt(2))*(I(32) - I(3-2))
      Sph_tensor_dq(12) = 2 * fac
     &                  * ( secdq(1,2,3) + secdq(2,1,3) + secdq(3,1,2) ) 

      Sph_tensor_dq(13) = fac
     &                  * ( 2 * secdq(1,1,3) - 2 * secdq(2,2,3)
     &                        + secdq(3,1,1) - secdq(3,2,2) ) 

      fac = 0.5_db

      Sph_tensor_dq(14) = fac
     &             * ( secdq(1,2,2) + 2 * secdq(2,1,2) - secdq(1,1,1) ) 

  ! (-i/sqrt(2))*(I(33) + I(3-3))
      Sph_tensor_dq(15) = - fac
     &             * ( secdq(2,1,1) + 2 * secdq(1,2,1) - secdq(2,2,2) ) 

      return
      end

!***********************************************************************

      subroutine Sph_tensor_qq_cal(n_tens_qq,secqq,Sph_tensor_qq)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db), dimension(3,3,3,3):: secqq
      complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq

! La multiplication par - img correspond a la partie imaginaire quand le
! tenseur n'est pas multiplie par le terme de Bragg.

! Tenseur 0, scalaire, signal isotropique quadrupolaire

      fac = 2 / sqrt(45._db) 
  ! Q(00)
      Sph_tensor_qq(1) = fac
     &      * ( 3 * ( secqq(1,3,1,3) + secqq(2,3,2,3) + secqq(1,2,1,2) )
     &              + secqq(1,1,1,1) + secqq(2,2,2,2) + secqq(3,3,3,3)
     &      - 0.5 * ( secqq(1,1,2,2) + secqq(1,1,3,3) + secqq(3,3,2,2)  
     &            + secqq(2,2,1,1) + secqq(3,3,1,1) + secqq(2,2,3,3) ) )  

! Tenseur 1, vecteur, lz, lx, ly, magnetique

      fac = 2 / sqrt(10._db) 

  ! Q(10)
      Sph_tensor_qq(2) = img * fac
     &           * ( secqq(2,1,1,1) - secqq(1,2,2,2) + secqq(2,3,1,3)
     &             - secqq(1,1,2,1) + secqq(2,2,1,2) - secqq(1,3,2,3) )

  ! (1/sqrt(2)) * ( Q(11) - Q(1-1) )
      Sph_tensor_qq(3) = - img * fac
     &           * ( secqq(3,2,2,2) - secqq(2,3,3,3) + secqq(3,1,2,1)
     &             - secqq(2,2,3,2) + secqq(3,3,2,3) - secqq(2,1,3,1) )

  ! -i * (1/sqrt(2)) * ( Q(11) + Q(1-1) )
      Sph_tensor_qq(4) = - img * fac
     &           * ( secqq(1,3,3,3) - secqq(3,1,1,1) + secqq(1,2,3,2)
     &             - secqq(3,3,1,3) + secqq(1,1,3,1) - secqq(3,2,1,2) )

! Tenseur 2 : quadrupole non magnetique

      fac = 2 / ( 3 * sqrt(14._db) ) 
  ! Q(20)
      Sph_tensor_qq(5) = fac
     &       * ( 6 * secqq(1,2,1,2) - 3 * secqq(1,3,1,3)
     &         - 3 * secqq(2,3,2,3) + secqq(1,1,1,1) + secqq(2,2,2,2) 
     &         - 2 * secqq(1,1,2,2) - 2 * secqq(2,2,1,1)
     &         - 2 * secqq(3,3,3,3) 
     &         + secqq(1,1,3,3) + secqq(2,2,3,3) 
     &         + secqq(3,3,1,1) + secqq(3,3,2,2) ) 

      fac = 2 / sqrt(42._db) 

  ! (1/sqrt(2)) * ( Q(21) - Q(2-1) )
      Sph_tensor_qq(6) = fac
     &                 * ( 3 * secqq(2,3,1,2) + secqq(3,3,1,3) 
     &                   + secqq(1,1,1,3) - 2 * secqq(2,2,1,3)
     &                   + 3 * secqq(1,2,2,3) + secqq(1,3,3,3) 
     &                   + secqq(1,3,1,1) - 2 * secqq(1,3,2,2) )

  ! (-i/sqrt(2)) * ( Q(21) + Q(2-1) )
      Sph_tensor_qq(7) = fac
     &                 * ( 3 * secqq(1,3,1,2) + secqq(3,3,2,3) 
     &                       + secqq(2,2,2,3) - 2 * secqq(1,1,2,3)
     &                   + 3 * secqq(1,2,1,3) + secqq(2,3,3,3) 
     &                       + secqq(2,3,2,2) - 2 * secqq(2,3,1,1) )

      fac = 2 / sqrt( 42._db ) 

  ! (-i/sqrt(2)) * ( Q(22) - Q(2-2) )
      Sph_tensor_qq(8) = fac
     &         * ( 2 * secqq(3,3,1,2) - secqq(1,1,1,2) - secqq(2,2,2,1) 
     &           - 3 * secqq(1,3,2,3)
     &           + 2 * secqq(1,2,3,3) - secqq(1,2,1,1) - secqq(2,1,2,2) 
     &           - 3 * secqq(2,3,1,3) )

  ! (1/sqrt(2)) * ( Q(22) + Q(2-2) )
      Sph_tensor_qq(9) = fac
     &         * ( secqq(3,3,1,1) - secqq(3,3,2,2)
     &           + secqq(1,1,3,3) - secqq(2,2,3,3) 
     &           + secqq(2,2,2,2) - secqq(1,1,1,1)
     &           + 3 * secqq(2,3,2,3) - 3 * secqq(1,3,1,3) )

! Tenseur 3 : octupole magnetique

      fac = 1 / sqrt(10._db) 

  ! Q(30)
      Sph_tensor_qq(10) = img * fac
     &       * ( secqq(1,2,1,1) - secqq(2,1,2,2) + 4 * secqq(1,3,2,3)
     &         - secqq(1,1,1,2) + secqq(2,2,2,1) - 4 * secqq(2,3,1,3) )

      fac = 0.5_db / sqrt(15._db)

  ! (1/sqrt(2)) * ( Q(31) - Q(3-1) )
      Sph_tensor_qq(11) = img *  fac
     &                  * ( 6 * secqq(1,2,1,3) + 4 * secqq(3,3,2,3)
     &                    - 5 * secqq(1,1,2,3) + secqq(2,2,2,3)
     &                    - 6 * secqq(1,3,1,2) - 4 * secqq(2,3,3,3)
     &                    + 5 * secqq(2,3,1,1) - secqq(2,3,2,2) )

  ! (-i/sqrt(2)) * ( Q(31) + Q(3-1) )
      Sph_tensor_qq(12) = - img * fac
     &                  * ( 6 * secqq(1,2,2,3) + 4 * secqq(3,3,1,3)
     &                    - 5 * secqq(2,2,1,3) + secqq(1,1,1,3)
     &                    - 6 * secqq(2,3,1,2) - 4 * secqq(1,3,3,3)
     &                    + 5 * secqq(1,3,2,2) - secqq(1,3,1,1) )

      fac = 1 / sqrt(6._db)

  ! (-i/sqrt(2)) * ( Q(32) - Q(3-2) )
      Sph_tensor_qq(13) = - img * fac
     &           * ( secqq(1,1,3,3) + secqq(2,2,1,1) + secqq(3,3,2,2)
     &             - secqq(3,3,1,1) - secqq(1,1,2,2) - secqq(2,2,3,3) )

  ! (1/sqrt(2)) * ( Q(32) + Q(3-2) )
      Sph_tensor_qq(14) = img * fac
     &       * ( 2 * secqq(1,2,3,3) - secqq(1,2,1,1) - secqq(1,2,2,2)
     &         - 2 * secqq(3,3,1,2) + secqq(1,1,1,2) + secqq(2,2,1,2) )

      fac = 0.5_db

  ! (1/sqrt(2)) * ( Q(33) - Q(3-3) )
      Sph_tensor_qq(15) = img * fac
     &       * ( secqq(2,3,1,1) - secqq(2,3,2,2) + 2 * secqq(1,3,1,2)
     &         - secqq(1,1,2,3) + secqq(2,2,2,3) - 2 * secqq(1,2,1,3) )

  ! (-i/sqrt(2)) * ( Q(33) + Q(3-3) )
      Sph_tensor_qq(16) = - img * fac
     &       * ( secqq(1,3,1,1) - secqq(1,3,2,2) + 2 * secqq(2,1,3,2)
     &         - secqq(1,1,1,3) + secqq(2,2,1,3) - 2 * secqq(3,2,2,1) )

! Tenseur 4 : hexadecapole non magnetique

      fac = 1 / ( 2. * sqrt(70._db) ) 

  ! Q(40)
      Sph_tensor_qq(17) = fac
     &   * ( 3 * secqq(1,1,1,1) + 3 * secqq(2,2,2,2)
     &         + 8 * secqq(3,3,3,3) + secqq(1,1,2,2) + secqq(2,2,1,1)
     &         - 4 * secqq(3,3,1,1) - 4 * secqq(3,3,2,2)
     &         - 4 * secqq(1,1,3,3) - 4 * secqq(2,2,3,3)
     &         + 4 * secqq(1,2,1,2)
     &        - 16 * secqq(1,3,1,3) - 16 * secqq(2,3,2,3) )

      fac = 0.5_db / sqrt(7._db) 

  ! (1/sqrt(2)) * ( Q(41) - Q(4-1) )
      Sph_tensor_qq(18) = fac
     &      * ( 3 * secqq(1,3,1,1) + secqq(1,3,2,2) - 4 * secqq(3,3,1,3) 
     &        + 2 * secqq(2,3,1,2)
     &        + 3 * secqq(1,1,1,3) + secqq(2,2,1,3) - 4 * secqq(1,3,3,3) 
     &        + 2 * secqq(1,2,2,3) )

  ! (-i/sqrt(2)) * ( Q(41) + Q(4-1) )
      Sph_tensor_qq(19) = fac
     &      * ( 3 * secqq(2,3,2,2) + secqq(2,3,1,1) - 4 * secqq(2,3,3,3) 
     &        + 2 * secqq(1,3,1,2)
     &        + 3 * secqq(2,2,2,3) + secqq(1,1,2,3) - 4 * secqq(3,3,2,3) 
     &        + 2 * secqq(1,2,1,3) )

      fac = 0.5_db / sqrt(14._db) 

  ! (-i/sqrt(2)) * ( Q(42) - Q(4-2) )
      Sph_tensor_qq(20) = 2 * fac
     &      * ( 2 * secqq(3,3,1,2) - secqq(1,1,1,2) - secqq(2,2,2,1) 
     &        + 4 * secqq(1,3,2,3)
     &        + 2 * secqq(1,2,3,3) - secqq(1,2,1,1) - secqq(2,1,2,2) 
     &        + 4 * secqq(2,3,1,3) )

  ! (1/sqrt(2)) * ( Q(42) + Q(4-2) )
      Sph_tensor_qq(21) = fac
     &    * ( 2 * secqq(3,3,1,1) - 2 * secqq(3,3,2,2)  
     &      + 2 * secqq(1,1,3,3) - 2 * secqq(2,2,3,3) 
     &      - 2 * secqq(1,1,1,1) + 2 * secqq(2,2,2,2) 
     &      + 8 * secqq(1,3,1,3) - 8 * secqq(2,3,2,3) )

      fac = 0.5_db 

  ! (1/sqrt(2)) * ( Q(43) - Q(4-3) )
      Sph_tensor_qq(22) = fac
     &       * ( secqq(1,3,2,2) - secqq(1,3,1,1) + 2 * secqq(2,3,1,2) 
     &         + secqq(2,2,1,3) - secqq(1,1,1,3) + 2 * secqq(1,2,2,3) ) 

  ! (-i/sqrt(2)) * ( Q(43) + Q(4-3) )
      Sph_tensor_qq(23) = fac
     &       * ( secqq(2,3,2,2) - secqq(2,3,1,1) - 2 * secqq(1,3,1,2)
     &         + secqq(2,2,2,3) - secqq(1,1,2,3) - 2 * secqq(1,2,1,3) ) 

      fac = 1 / sqrt(2._db) 

  ! (-i/sqrt(2)) * ( Q(44) - Q(4-4) )
      Sph_tensor_qq(24) = fac * ( secqq(1,2,1,1) - secqq(1,2,2,2)
     &                          + secqq(1,1,1,2) - secqq(2,2,1,2) ) 

  ! (1/sqrt(2)) * ( Q(44) + Q(4-4) )
      Sph_tensor_qq(25) = 0.5_db * fac
     &         * ( secqq(1,1,1,1) + secqq(2,2,2,2) - secqq(1,1,2,2)
     &           - secqq(2,2,1,1) - 4 * secqq(1,2,1,2) ) 

      return
      end

!***********************************************************************

      subroutine Tensor_pol_dd_cal(ipl,n_tens_dd,nplt,pe,ps,
     &                             Tensor_pol_dd)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
      complex(kind=db), dimension(3):: pe, ps
      complex(kind=db), dimension(n_tens_dd):: Tens
      complex(kind=db), dimension(0:nplt,n_tens_dd):: Tensor_pol_dd
 
      Px = Pe(1); Qx = conjg( Ps(1) )
      Py = Pe(2); Qy = conjg( Ps(2) )
      Pz = Pe(3); Qz = conjg( Ps(3) )

      fac = 1 / sqrt(3._db)

      Tens(1) = fac * ( Qx*Px + Qy*Py + Qz*Pz ) 

      fac = 1 / sqrt(2._db)

      Tens(2) = - img * fac * ( Qx*Py - Qy*Px ) 

      Tens(3) = img * fac * ( Qy*Pz - Qz*Py ) 

      Tens(4) = fac * ( Qx*Pz - Qz*Px ) 

      fac = 1 / sqrt(6._db)

      Tens(5) = fac * ( 2*Qz*Pz - Qx*Px - Qy*Py ) 

      fac = 1._db / sqrt( 2._db )

      Tens(6) = - fac * ( Qx*Pz + Qz*Px ) 

      Tens(7) = - fac * img * ( Qy*Pz + Qz*Py ) 

      Tens(8) = fac * img * ( Qx*Py + Qy*Px ) 

      Tens(9) = fac * ( Qx*Px - Qy*Py )  

      Tensor_pol_dd(ipl,:) = Tens(:)

      return
      end

!***********************************************************************

      subroutine Tensor_pol_dq_cal(ipl,n_tens_dq,nplt,pe,ps,ve,vs,
     &                             Tensor_pol_dq)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
      complex(kind=db), dimension(3):: pe, ps
      complex(kind=db), dimension(2*n_tens_dq):: Tens
      complex(kind=db), dimension(0:nplt,2*n_tens_dq):: Tensor_pol_dq
 
      real(kind=db), dimension(3):: ve, vs

      Px = Pe(1); Qx = conjg( Ps(1) )
      Py = Pe(2); Qy = conjg( Ps(2) )
      Pz = Pe(3); Qz = conjg( Ps(3) )
      
      Vx = Ve(1); Wx = Vs(1)
      Vy = Ve(2); Wy = Vs(2)
      Vz = Ve(3); Wz = Vs(3)

      Tens(:) = 0._db

      j = 0

      do i = 1,2

        if( i == 2 ) then
          Wx = - Wx; Wy = - Wy; Wz = - Wz
          j = n_tens_dq
        endif

  ! T(10)
        Tens(j+1) = - ( 1 / sqrt(15._db) ) 
     &     * ( ( 1.5*Qx*Px + 1.5*Qy*Py + 2*Qz*Pz ) * ( Vz - Wz )
     &       + ( 1.5*Qx*Pz - Qz*Px ) * Vx - ( 1.5*Qz*Px - Qx*Pz ) * Wx 
     &       + ( 1.5*Qy*Pz - Qz*Py ) * Vy - ( 1.5*Qz*Py - Qy*Pz ) * Wy ) 

      fac = 1 / sqrt( 60._db )

  ! (T(11)-T(1-1))/sqrt(2)
        Tens(j+2) = fac 
     &     * ( ( 4*Qx*Px + 3*Qy*Py + 3*Qz*Pz ) * ( Vx - Wx )
     &       + ( 3*Qy*Px - 2*Qx*Py ) * Vy - ( 3*Qx*Py - 2*Qy*Px ) * Wy
     &       + ( 3*Qz*Px - 2*Qx*Pz ) * Vz - ( 3*Qx*Pz - 2*Qz*Px ) * Wz )

  ! (T(11)+T(1-1))/sqrt(2)
        Tens(j+3) = fac * img 
     &     * ( ( 4*Qy*Py + 3*Qx*Px + 3*Qz*Pz ) * ( Vy - Wy )
     &       + ( 3*Qx*Py - 2*Qy*Px ) * Vx - ( 3*Qy*Px - 2*Qx*Py ) * Wx
     &       + ( 3*Qz*Py - 2*Qy*Pz ) * Vz - ( 3*Qy*Pz - 2*Qz*Py ) * Wz )

  ! T(20)
       Tens(j+4) = 0.5_db * img 
     &        * ( ( Qx*Py - Qy*Px ) * ( Vz + Wz )
     &            - Qy*Pz*Vx + Qz*Py*Wx + Qx*Pz*Vy - Qz*Px*Wy ) 

        fac = 1 / sqrt( 12._db )

  ! (T(21)-T(2-1))/sqrt(2)
        Tens(j+5) = fac * img 
     &        * ( ( Qz*Pz - Qx*Px ) * ( Vy - Wy )
     &          + ( Qz*Py - 2*Qy*Pz ) * Vz - ( Qy*Pz - 2*Qz*Py ) * Wz
     &          + ( 2*Qy*Px - Qx*Py ) * Vx - ( 2*Qx*Py - Qy*Px ) * Wx )

  ! (T(21)+T(2-1))/sqrt(2)
        Tens(j+6) = fac 
     &        * ( ( Qz*Pz - Qy*Py ) * ( Vx - Wx )
     &          + ( Qz*Px - 2*Qx*Pz ) * Vz - ( Qx*Pz - 2*Qz*Px ) * Wz
     &          + ( 2*Qx*Py - Qy*Px ) * Vy - ( 2*Qy*Px - Qx*Py ) * Wy )

  ! (T(22)-T(2-2))/sqrt(2)
        Tens(j+7) = fac 
     &        * ( ( Qx*Px - Qy*Py ) * ( Vz - Wz )
     &          + ( Qx*Pz - 2*Qz*Px ) * Vx - ( Qz*Px - 2*Qx*Pz ) * Wx
     &          + ( 2*Qz*Py - Qy*Pz ) * Vy - ( 2*Qy*Pz - Qz*Py ) * Wy )

  ! (T(22)+T(2-2))/sqrt(2)
        Tens(j+8) = fac * img 
     &        * ( ( Qx*Py + Qy*Px ) * ( Vz - Wz )
     &          + ( Qx*Pz - 2*Qz*Px ) * Vy - ( Qz*Px - 2*Qx*Pz ) * Wy
     &          - ( 2*Qz*Py - Qy*Pz ) * Vx + ( 2*Qy*Pz - Qz*Py ) * Wx )

  ! T(30)
        Tens(j+9) = ( 1 / sqrt(10._db) )  
     &        * ( ( 2*Qz*Pz - Qx*Px - Qy*Py ) * ( Vz - Wz )
     &                    - ( Qx*Pz + Qz*Px ) * ( Vx - Wx ) 
     &                    - ( Qy*Pz + Qz*Py ) * ( Vy - Wy ) ) 

        fac = 1 / sqrt( 60._db )

  ! (T(31)-T(3-1))/sqrt(2)
        Tens(j+10) = fac 
     &        * ( ( 3*Qx*Px + Qy*Py - 4*Qz*Pz ) * ( Vx - Wx )
     &                      + ( Qx*Py + Qy*Px ) * ( Vy - Wy ) 
     &                  - 4 * ( Qx*Pz + Qz*Px ) * ( Vz - Wz ) ) 

  ! (T(31)+T(3-1))/sqrt(2)
        Tens(j+11) = fac * img 
     &        * ( ( 3*Qy*Py + Qx*Px - 4*Qz*Pz ) * ( Vy - Wy )
     &                      + ( Qx*Py + Qy*Px ) * ( Vx - Wx ) 
     &                  - 4 * ( Qy*Pz + Qz*Py ) * ( Vz - Wz ) ) 

        fac = 1 / sqrt( 6._db )

  ! (T(32)-T(3-2))/sqrt(2)
        Tens(j+12) = fac * img * ( ( Qx*Py + Qy*Px ) * ( Vz - Wz )
     &                           + ( Qx*Pz + Qz*Px ) * ( Vy - Wy ) 
     &                           + ( Qy*Pz + Qz*Py ) * ( Vx - Wx ) ) 

  ! (T(32)+T(3-2))/sqrt(2)
        Tens(j+13) = fac  * ( ( Qx*Px - Qy*Py ) * ( Vz - Wz )
     &                      + ( Qx*Pz + Qz*Px ) * ( Vx - Wx ) 
     &                      - ( Qy*Pz + Qz*Py ) * ( Vy - Wy ) )

  ! (T(33)-T(3-3))/sqrt(2)
        Tens(j+14) = 0.5_db * ( ( Qy*Py - Qx*Px ) * ( Vx - Wx )
     &                       + ( Qy*Px + Qx*Py ) * ( Vy - Wy ) ) 

  ! (T(33)+T(3-3))/sqrt(2)
        Tens(j+15) = 0.5_db * img * ( ( Qy*Py - Qx*Px ) * ( Vy - Wy )
     &                             - ( Qy*Px + Qx*Py ) * ( Vx - Wx ) ) 

      end do

      Tensor_pol_dq(ipl,:) = Tens(:) 

      return
      end

!***********************************************************************

      subroutine Tensor_pol_qq_cal(ipl,n_tens_qq,nplt,pe,ps,ve,vs,
     &                             Tensor_pol_qq)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
      complex(kind=db), dimension(3):: pe, ps
      complex(kind=db), dimension(n_tens_qq):: Tens
      complex(kind=db), dimension(0:nplt,n_tens_qq):: Tensor_pol_qq
 
      real(kind=db), dimension(3):: ve, vs

      Px = Pe(1); Qx = conjg( Ps(1) ) 
      Py = Pe(2); Qy = conjg( Ps(2) )
      Pz = Pe(3); Qz = conjg( Ps(3) )
      
      Vx = Ve(1); Wx = Vs(1)
      Vy = Ve(2); Wy = Vs(2)
      Vz = Ve(3); Wz = Vs(3)

      fac = 1 / ( 3 * sqrt( 5._db ) )

! Scalaire

      Tens(1) = 1.5 * fac 
     &   * ( Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx 
     &     + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy 
     &     + Qx*Wy*Px*Vy + Qx*Wy*Py*Vx + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx )
      Tens(1) = Tens(1) + 2 * fac 
     &   * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy + Qz*Wz*Pz*Vz )
     &        - fac 
     &   * ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qx*Wx*Pz*Vz + Qz*Wz*Px*Vx 
     &     + Qy*Wy*Pz*Vz + Qz*Wz*Py*Vy ) 

! Dipole magnetique

      fac = 1 / sqrt(10._db)

      Tens(2) = fac * img 
     &   * ( Qx*Wy*Px*Vx + Qy*Wx*Px*Vx - Qx*Wy*Py*Vy - Qy*Wx*Py*Vy 
     &     - Qx*Wx*Px*Vy - Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx ) 
     &        + 0.5 * fac * img 
     &   * ( Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx 
     &     - Qx*Wz*Py*Vz - Qx*Wz*Pz*Vy - Qz*Wx*Py*Vz - Qz*Wx*Pz*Vy ) 

      Tens(3) = 0.5 * fac * img 
     &   * ( Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx 
     &     - Qx*Wz*Px*Vy - Qx*Wz*Py*Vx - Qz*Wx*Px*Vy - Qz*Wx*Py*Vx ) 
     &        - fac * img 
     &   * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy - Qy*Wy*Py*Vz - Qy*Wy*Pz*Vy 
     &     - Qy*Wz*Pz*Vz - Qz*Wy*Pz*Vz + Qy*Wz*Py*Vy + Qz*Wy*Py*Vy ) 

      Tens(4) = 0.5 * fac 
     &   * ( Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy 
     &     - Qy*Wz*Px*Vy - Qy*Wz*Py*Vx - Qz*Wy*Px*Vy - Qz*Wy*Py*Vx ) 
     &        - fac 
     &   * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx - Qx*Wx*Px*Vz - Qx*Wx*Pz*Vx 
     &     - Qx*Wz*Pz*Vz - Qz*Wx*Pz*Vz + Qx*Wz*Px*Vx + Qz*Wx*Px*Vx ) 

! Quadrupole non-magnetique

      fac = 1 / ( 3 * sqrt(14._db) )

      Tens(5) = 3 * fac 
     &   * ( Qx*Wy*Px*Vy + Qx*Wy*Py*Vx + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx ) 
     &        - 1.5 * fac 
     &   * ( Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx 
     &     + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy ) 
     &        + 2 * fac 
     &   * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy + Qx*Wx*Pz*Vz + Qz*Wz*Px*Vx 
     &     + Qy*Wy*Pz*Vz + Qz*Wz*Py*Vy) 
     &        - 4 * fac 
     &   * ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qz*Wz*Pz*Vz ) 

      fac = 1 / sqrt( 42._db )

      Tens(6) = 1.5 * fac 
     &   * ( Qy*Wz*Px*Vy + Qz*Wy*Px*Vy + Qy*Wz*Py*Vx + Qz*Wy*Py*Vx 
     &     + Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy ) 
     &        + fac 
     &   * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx + Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx 
     &     + Qx*Wz*Pz*Vz + Qz*Wx*Pz*Vz + Qx*Wz*Px*Vx + Qz*Wx*Px*Vx ) 
     &        - 2 * fac 
     &   * ( Qy*Wy*Px*Vz + Qy*Wy*Pz*Vx + Qx*Wz*Py*Vy + Qz*Wx*Py*Vy )

      Tens(7) = 1.5 * fac * img 
     &   * ( Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx 
     &     + Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx ) 
     &        + fac * img 
     &   * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy + Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy 
     &     + Qy*Wz*Pz*Vz + Qz*Wy*Pz*Vz + Qy*Wz*Py*Vy + Qz*Wy*Py*Vy ) 
     &        - 2 * fac * img 
     &   * ( Qx*Wx*Py*Vz + Qx*Wx*Pz*Vy + Qy*Wz*Px*Vx + Qz*Wy*Px*Vx )

      fac = 1 / sqrt( 42._db )

      Tens(8) = - fac * img 
     &   * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx 
     &     + Qx*Wy*Px*Vx + Qy*Wx*Px*Vx + Qx*Wy*Py*Vy + Qy*Wx*Py*Vy ) 
     &        + 2 * fac * img
     &   * ( Qz*Wz*Px*Vy + Qz*Wz*Py*Vx + Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz )   
     &        - 1. 5 * fac * img 
     &   * ( Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy 
     &     + Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx ) 

      Tens(9) = 2 * fac 
     &   * ( Qz*Wz*Px*Vx - Qz*Wz*Py*Vy + Qx*Wx*Pz*Vz - Qy*Wy*Pz*Vz 
     &     + Qy*Wy*Py*Vy - Qx*Wx*Px*Vx ) 
     &        + 1. 5 * fac 
     &   * ( Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy 
     &     - Qx*Wz*Px*Vz - Qx*Wz*Pz*Vx - Qz*Wx*Px*Vz - Qz*Wx*Pz*Vx ) 

! Octupole magnetique

      fac = 1 / sqrt( 10._db )

      Tens(10) = 0.5 * fac * img 
     &   * ( Qx*Wy*Px*Vx - Qx*Wy*Py*Vy + Qy*Wx*Px*Vx - Qy*Wx*Py*Vy 
     &     - Qx*Wx*Px*Vy + Qy*Wy*Px*Vy - Qx*Wx*Py*Vx + Qy*Wy*Py*Vx )
     &         + fac * img 
     &   * ( Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy 
     &     - Qy*Wz*Px*Vz - Qz*Wy*Px*Vz - Qy*Wz*Pz*Vx - Qz*Wy*Pz*Vx ) 

      fac = 1 / ( 4 * sqrt( 15._db ) )

      Tens(11) = 3 * fac * img 
     &   * ( Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx 
     &     - Qx*Wz*Px*Vy - Qx*Wz*Py*Vx - Qz*Wx*Px*Vy - Qz*Wx*Py*Vx ) 
     &         + 4 * fac * img 
     &   * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy - Qy*Wz*Pz*Vz - Qz*Wy*Pz*Vz ) 
     &         + 5 * fac * img 
     &   * ( Qy*Wz*Px*Vx + Qz*Wy*Px*Vx - Qx*Wx*Py*Vz - Qx*Wx*Pz*Vy ) 
     &         + fac * img 
     &   * ( Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy - Qy*Wz*Py*Vy - Qz*Wy*Py*Vy ) 

      Tens(12) = 3 * fac 
     &   * ( Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy 
     &     - Qy*Wz*Px*Vy - Qy*Wz*Py*Vx - Qz*Wy*Px*Vy - Qz*Wy*Py*Vx ) 
     &         + 4 * fac 
     &   * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx - Qx*Wz*Pz*Vz - Qz*Wx*Pz*Vz ) 
     &         + 5 * fac 
     &   * ( Qx*Wz*Py*Vy + Qz*Wx*Py*Vy - Qy*Wy*Px*Vz - Qy*Wy*Pz*Vx ) 
     &         + fac 
     &   * ( Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx - Qx*Wz*Px*Vx - Qz*Wx*Px*Vx ) 

      fac = 1 / sqrt( 6._db )

      Tens(13) = fac 
     &   * ( Qx*Wx*Pz*Vz - Qz*Wz*Px*Vx + Qy*Wy*Px*Vx - Qx*Wx*Py*Vy 
     &     + Qz*Wz*Py*Vy - Qy*Wy*Pz*Vz ) 

      Tens(14) = fac * img 
     &   * ( Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz - Qz*Wz*Px*Vy - Qz*Wz*Py*Vx ) 
     &         + 0.5 * fac * img 
     &   * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx 
     &     - Qx*Wy*Px*Vx - Qy*Wx*Px*Vx - Qx*Wy*Py*Vy - Qy*Wx*Py*Vy ) 

      fac = 0.25_db

      Tens(15) = fac * img 
     &   * ( Qy*Wz*Px*Vx - Qy*Wz*Py*Vy + Qz*Wy*Px*Vx - Qz*Wy*Py*Vy 
     &     + Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx 
     &     - Qx*Wx*Py*Vz + Qy*Wy*Py*Vz - Qx*Wx*Pz*Vy + Qy*Wy*Pz*Vy 
     &     - Qx*Wy*Px*Vz - Qy*Wx*Px*Vz - Qx*Wy*Pz*Vx - Qy*Wx*Pz*Vx ) 

      Tens(16) = fac 
     &   * ( Qx*Wz*Px*Vx - Qx*Wz*Py*Vy + Qz*Wx*Px*Vx - Qz*Wx*Py*Vy 
     &     - Qy*Wz*Px*Vy - Qy*Wz*Py*Vx - Qz*Wy*Px*Vy - Qz*Wy*Py*Vx 
     &     - Qx*Wx*Px*Vz + Qy*Wy*Px*Vz - Qx*Wx*Pz*Vx + Qy*Wy*Pz*Vx 
     &     + Qx*Wy*Py*Vz + Qy*Wx*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Pz*Vy ) 

! Hexadecapole non-magnetique

      fac = 1 / ( 2 * sqrt( 70._db ) ) 

      Tens(17) = 3 * fac 
     &   * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy ) 
     &         + fac
     &   * ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qx*Wy*Px*Vy + Qx*Wy*Py*Vx 
     &     + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx ) 
     &         + 8 * fac
     &   * Qz*Wz*Pz*Vz 
     &         - 4 * fac
     &   * ( Qz*Wz*Px*Vx + Qz*Wz*Py*Vy + Qx*Wx*Pz*Vz + Qy*Wy*Pz*Vz 
     &     + Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx 
     &     + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy ) 

      fac = 1 / ( 4 * sqrt( 7._db ) ) 

      Tens(18) = 3 * fac 
     &   * ( Qx*Wz*Px*Vx + Qz*Wx*Px*Vx + Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx ) 
     &         + fac
     &   * ( Qx*Wz*Py*Vy + Qz*Wx*Py*Vy + Qy*Wy*Px*Vz + Qy*Wy*Pz*Vx  
     &     + Qy*Wz*Px*Vy + Qy*Wz*Py*Vx + Qz*Wy*Px*Vy + Qz*Wy*Py*Vx 
     &     + Qy*Wx*Pz*Vy + Qx*Wy*Pz*Vy + Qx*Wy*Py*Vz + Qy*Wx*Py*Vz ) 
     &         - 4 * fac
     &   * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx + Qx*Wz*Pz*Vz + Qz*Wx*Pz*Vz ) 

      Tens(19) = 3 * fac * img 
     &   * ( Qy*Wz*Py*Vy + Qz*Wy*Py*Vy + Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy ) 
     &         + fac * img
     &   * ( Qy*Wz*Px*Vx + Qz*Wy*Px*Vx + Qx*Wx*Py*Vz + Qx*Wx*Pz*Vy  
     &     + Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx 
     &     + Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx ) 
     &         - 4 * fac * img
     &   * ( Qy*Wz*Pz*Vz + Qz*Wy*Pz*Vz + Qz*Wz*Pz*Vy + Qz*Wz*Py*Vz ) 

      fac = 1 / sqrt( 14._db ) 

      Tens(20) = fac * img 
     &   * ( Qz*Wz*Px*Vy + Qz*Wz*Py*Vx + Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz  
     &     + Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy 
     &     + Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx ) 
     &         - 0.5 * fac * img
     &   * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx  
     &     + Qx*Wy*Px*Vx + Qy*Wx*Px*Vx + Qx*Wy*Py*Vy + Qy*Wx*Py*Vy ) 

      Tens(21) = fac 
     &   * ( Qz*Wz*Px*Vx - Qz*Wz*Py*Vy - Qx*Wx*Px*Vx + Qy*Wy*Py*Vy  
     &     + Qx*Wx*Pz*Vz - Qy*Wy*Pz*Vz + Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx 
     &     + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx - Qy*Wz*Py*Vz - Qy*Wz*Pz*Vy
     &     - Qz*Wy*Py*Vz - Qz*Wy*Pz*Vy ) 

      fac = 0.25_db 

      Tens(22) = - fac 
     &   * ( Qx*Wz*Px*Vx - Qx*Wz*Py*Vy + Qz*Wx*Px*Vx - Qz*Wx*Py*Vy  
     &     - Qy*Wz*Py*Vx - Qy*Wz*Px*Vy - Qz*Wy*Py*Vx - Qz*Wy*Px*Vy 
     &     - Qx*Wy*Pz*Vy + Qx*Wx*Px*Vz - Qy*Wy*Px*Vz + Qx*Wx*Pz*Vx
     &     - Qy*Wy*Pz*Vx - Qy*Wx*Py*Vz - Qx*Wy*Py*Vz - Qy*Wx*Pz*Vy ) 

      Tens(23) = - fac * img 
     &   * ( Qy*Wz*Px*Vx - Qy*Wz*Py*Vy + Qz*Wy*Px*Vx - Qz*Wy*Py*Vy  
     &     + Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx 
     &     + Qy*Wx*Pz*Vx + Qx*Wx*Py*Vz - Qy*Wy*Py*Vz + Qx*Wx*Pz*Vy
     &     - Qy*Wy*Pz*Vy + Qx*Wy*Px*Vz + Qy*Wx*Px*Vz + Qx*Wy*Pz*Vx ) 

      fac = 1 / ( 2 * sqrt( 2._db ) ) 

      Tens(24) = fac * img 
     &   * ( Qx*Wy*Px*Vx - Qx*Wy*Py*Vy + Qy*Wx*Px*Vx - Qy*Wx*Py*Vy  
     &     + Qx*Wx*Px*Vy + Qx*Wx*Py*Vx - Qy*Wy*Px*Vy - Qy*Wy*Py*Vx ) 

      Tens(25) = fac 
     &   * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy - Qx*Wx*Py*Vy - Qy*Wy*Px*Vx  
     &     - Qx*Wy*Px*Vy - Qx*Wy*Py*Vx - Qy*Wx*Px*Vy - Qy*Wx*Py*Vx ) 

      Tensor_pol_qq(ipl,:) = Tens(:) 

      return
      end

!***********************************************************************

! Ecriture des fonctions physiques

      subroutine write_phys(ct_nelec,Densite_atom,E_cut,E1E1,E1E2,E2E2,
     &      Energ,Ephseuil,Epsii,
     &      Eseuil,ia,ie,initlr,Int_tens,ipl0,ipl2,ipldafs,Length_word,
     &      jseuil,magn_sens,n_tens_dd,n_tens_dq,n_tens_max,
     &      n_tens_qq,n_tens_t,natomsym,nenerg,ninitlr,
     &      nomfich1,npldafs,nplt,nseuil,numat_abs,phdf0t,phdt,
     &      Sph_tensor_dd,Sph_tensor_dq,Sph_tensor_dq_m,
     &      Sph_tensor_qq,Tensor_pol_dd,Tensor_pol_dq,Tensor_pol_qq,
     &      v0muf,writout)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      character(len=Length_word) mot
      character(len=132) nomfich1, nomficht
      character(len=8), dimension(49):: nomtens
      character(len=Length_word), dimension(n_tens_max):: nomten

      complex(kind=db):: cf, cg, Ten, Ten_m
      complex(kind=db), dimension(1):: cdum
      complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
      complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq,
     &                                        Sph_tensor_dq_m
      complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq
      complex(kind=db), dimension(0:nplt,n_tens_dd):: Tensor_pol_dd
      complex(kind=db), dimension(0:nplt,2*n_tens_dq):: Tensor_pol_dq
      complex(kind=db), dimension(0:nplt,n_tens_qq):: Tensor_pol_qq
      complex(kind=db), dimension(npldafs):: phdf0t, phdt
      complex(kind=db), dimension(n_tens_max):: ph0, phtem, Resul

      integer, dimension(0):: idum

      logical E1E1, E1E2, E2E2,
     &        magn_sens, polarise, spherical_signal, writout 

      real(kind=db), dimension(0):: rdum
      real(kind=db), dimension(n_tens_max,ninitlr,0:natomsym):: Int_tens
      real(kind=db), dimension(n_tens_max):: Int_tenst, Tens
      real(kind=db), dimension(nenerg):: Energ
      real(kind=db), dimension(ninitlr):: Epsii

      common/polarise/ polarise
      common/spherical_signal/ spherical_signal

      data nomtens/ 
     &  '  D(00) ','  lz_dd ',' -lx_dd ','  ly_dd ','  D(20) ',
     &  '  D(21)d','-iD(21)s','-iD(22)d','  D(22)s',
     &  '  I(10) ','  I(11)d','-iI(11)s','-iI(20) ','-iI(21)d',
     &  '  I(21)s','  I(22)d','-iI(22)s','  I(30) ','  I(31)d',
     &  '  I(31)s','-iI(32)d',' I(32)s ','  I(33)d','-iI(33)s',
     &  '  Q(00) ','  lz_qq ',' -lx_qq ','  ly_qq ','  Q(20) ',
     &  '  Q(21)d','-iQ(21)s','-iQ(22)d','  Q(22)s','  Q(30) ',
     &  '  Q(31)d','-iQ(31)s','-iQ(32)d','  Q(32)s','-iQ(33)d',
     &  '-iQ(33)s','  Q(40) ','  Q(41)d','  Q(41)s','-iQ(42)d',
     &  '  Q(42)s','  Q(43)d','-iQ(43)s','-iQ(44)d','  Q(44)s'/

      if( writout ) then

        j = 0
        do itens = 1,n_tens_t

          if( itens <= n_tens_dd ) then

            if( .not. E1E1 ) cycle
            i = itens
            if( i == 1 ) then
! On divise la premiere composante du tenseur spherique par rac(3)
! (premiere composante du tenseur de polarisation) pour obtenir le
! terme de diffusion isotrope
              Ten = Sph_tensor_dd(i) / sqrt(3._db)
            elseif( i >= 2 .and. i <= 4 ) then
! On divise les composantes 2, 3 et 4 du tenseur spherique par rac(2)
! (composantes du tenseur de polarisation) pour obtenir le
! moment magnetique
              Ten = Sph_tensor_dd(i) / sqrt(2._db)
            else
              Ten = Sph_tensor_dd(i)
            endif

          elseif( itens > n_tens_dd 
     &                    .and. itens <= n_tens_dd + n_tens_dq ) then

            if( .not. E1E2 ) cycle
            i = itens - n_tens_dd
            Ten = Sph_tensor_dq(i)
            if( ia == 0 .and. magn_sens ) Ten_m = Sph_tensor_dq_m(i)

          elseif(  itens > n_tens_dd + n_tens_dq) then

            if( .not. E2E2 ) cycle
            i = itens - n_tens_dd - n_tens_dq
            if( i == 1 ) then
! On divise la premiere composante du tenseur spherique par rac(3)
! (premiere composante du tenseur de polarisation) pour obtenir le
! terme de diffusion isotrope
              Ten = Sph_tensor_qq(i) / sqrt(3._db)
            elseif( i >= 2 .and. i <= 4 ) then
! On divise les composantes 2, 3 et 4 du tenseur spherique par rac(2)
! (composantes du tenseur de polarisation) pour obtenir le
! moment magnetique
              Ten = Sph_tensor_qq(i) / sqrt(2._db)
            else
              Ten = Sph_tensor_qq(i)
            endif

          endif

          j = j + 1
          Tens(j) = real( Ten,db ) 
          if( ipldafs > 0 .or. ( magn_sens .and. itens > n_tens_dd 
     &                    .and. itens <= n_tens_dd + n_tens_dq ) ) then
            mot = nomtens(itens) 
            mot = adjustl( mot )
            l = len_trim( mot )
            mot(l+1:l+2) = '_r'
            nomten(j) = mot
            j = j + 1
            Tens(j) = aimag( Ten )
            mot(l+2:l+2) = 'i'
            nomten(j) = mot
          else 
            nomten(j) = nomtens(itens)
          endif
          if( ia == 0 .and. magn_sens .and. itens > n_tens_dd 
     &                    .and. itens <= n_tens_dd + n_tens_dq ) then
            j = j + 1
            Tens(j) = real( Ten_m,db ) 
            mot(l+2:l+3) = 'rm' 
            nomten(j) = mot
            j = j + 1
            Tens(j) = aimag( Ten_m ) 
            mot(l+2:l+3) = 'im' 
            nomten(j) = mot
          endif

        end do

        n_tens = j

        nomficht = nomfich1
        long = len_trim(nomficht)
        nomficht(long+1:long+4) = '_sph'
        if( ia > 0 ) then
          nomficht(long+5:long+9) = '_atom'
          call ad_number(ia,nomficht,132)
        else
          nomficht(long+5:long+9) = '_xtal'
        endif
        if( ipldafs > 0 ) then
          long = len_trim(nomficht)
          nomficht(long+1:long+4) = '_rxs'
          call ad_number(ipldafs,nomficht,132)
        endif
        long = len_trim(nomficht)
        nomficht(long+1:long+4) = '.txt'

        call write_out(rdum,rdum,Densite_atom,0._db,E_cut,Ephseuil,
     &          Epsii,Eseuil,.false.,idum,ie,Length_word,
     &          jseuil,n_tens_max,n_tens,0,ninitlr,nomficht,nomten,
     &          1,0,0,0,nseuil,numat_abs,cdum,
     &          cdum,Tens,v0muf,.false.,0)

      endif

! Integrale
      if( nenerg > 1 .and. ipldafs == 0 .and. writout ) then
        nomficht = nomfich1
        long = len_trim(nomficht)
        nomficht(long+1:long+4) = '_sph'
        if( ia > 0 ) then
          nomficht(long+5:long+9) = '_atom'
          call ad_number(ia,nomficht,132)
        else
          nomficht(long+5:long+9) = '_xtal'
        endif
        long = len_trim(nomficht)
        nomficht(long+1:long+8) = '_int.txt'

        if( ie == 1 ) then
          do i = 1,n_tens
            mot = nomten(i)
            mot = adjustl( mot )
            long = len_trim( mot )
            if( long < Length_word - 2 ) then
              mot(1:Length_word) = 'I_' // mot(1:Length_word-2)
            elseif( long == Length_word - 2 ) then
              mot(1:Length_word-1) = 'I' // mot(1:Length_word-1)
            endif 
            nomten(i) = mot
          end do
          de = Energ(2) - Energ(1) 
          Int_tens(1:n_tens,initlr,ia) = de * Tens(1:n_tens)
        else
          if( ie == nenerg ) then
            de = Energ(ie) - Energ(ie-1) 
          else 
            de = 0.5 * ( Energ(ie+1) -  Energ(ie-1) ) 
          endif
          Int_tens(1:n_tens,initlr,ia) = Int_tens(1:n_tens,initlr,ia)
     &                          + de * Tens(1:n_tens)
        endif

        Int_tenst(1:n_tens) = Int_tens(1:n_tens,initlr,ia) 
        call write_out(rdum,rdum,Densite_atom,0._db,E_cut,Ephseuil,
     &          Epsii,Eseuil,.false.,idum,ie,Length_word,
     &          jseuil,n_tens_max,n_tens,0,ninitlr,nomficht,nomten,
     &          1,0,0,0,nseuil,numat_abs,cdum,
     &          cdum,Int_tenst,v0muf,.false.,0)

      endif

      if( .not. spherical_signal ) return

! Calcul des tenseurs appliques aux reflexions RXS et au xanes

      if( E1E1 .and. E1E2 .and. E2E2 ) then
        j0 = 4
      elseif( ( E1E1 .and. E1E2 ) .or. ( E1E1 .and. E2E2 )
     &   .or. ( E1E2 .and. E2E2 ) ) then
        j0 = 3
      else
        j0 = 1
      endif
      jdd = 0; jdq = 0; jqq = 0
      if( E1E1 ) jdd = min(j0,2)
      if( E1E2 ) then
        if( E2E2 ) then      
          jdq = j0 - 1
        else
          jdq = j0
        endif
      endif
      if( E2E2 ) jqq = j0

      if( polarise ) then
        iplf = ipl2
      else
        iplf = ipl0
      endif

      do ipl = ipl0,iplf

        Resul(:) = (0._db,0._db)

        j = j0
        if( ipldafs > 0 ) then
          jj = 2 * j
        else
          jj = j
        endif  
        do itens = 1,n_tens_t

          if( itens <= n_tens_dd ) then
            if( .not. E1E1 ) cycle
            i = itens
            j = j + 1
            Resul(j) = Tensor_pol_dd(ipl,i) * Sph_tensor_dd(i)

! Devant le produit, il faut mettre un signe -1 devant les tenseurs
! impairs
            if( i >= 2 .and. i <= 4 ) Resul(j) = - Resul(j)   
! On recupere le img omis dans les tenseurs imaginaires
            if( i == 4 .or. i == 7 .or. i == 8 ) Resul(j) = img*Resul(j)
! On recupere le -1 devant les tenseurs differences
            if( i == 3 .or. i == 6 .or. i == 8 ) Resul(j) = - Resul(j)
! On recupere le (-1)**m
            if( i == 3 .or. i == 4 .or. i == 6 .or. i == 7)
     &                                         Resul(j) = - Resul(j)

            Resul(jdd) = Resul(jdd) + Resul(j)

          elseif( itens > n_tens_dd 
     &                   .and. itens <= n_tens_dd + n_tens_dq ) then

            if( .not. E1E2 ) cycle
            i = itens - n_tens_dd
            do is = 1,2
              if( is == 2 .and. .not. ( ia == 0 .and. magn_sens ) ) exit  
              j = j + 1
! Les parties reelles et imaginaires sont a considerer avant la
! multiplication par le img eventuel.
              if( is == 1 ) then
                Resul(j) = Sph_tensor_dq(i) * Tensor_pol_dq(ipl,i) 
              else
                Resul(j) = Sph_tensor_dq_m(i)
     &                   * Tensor_pol_dq(ipl,i+n_tens_dq)
              endif

! Devant le produit, il faut mettre un signe -1 devant les tenseurs
! impairs
              if( i >= 4 .and. i <= 8 ) Resul(j) = - Resul(j)
! On recupere le img omis dans les tenseurs
              if( i == 3 .or. i == 4 .or. i == 5 .or. i == 8 
     &         .or. i == 11 .or. i == 12 .or. i == 15 ) 
     &                                  Resul(j) = img * Resul(j)

! On multiplie par -1 devant les tenseurs differences
              if( i == 2 .or. i == 5 .or. i == 7 .or. i == 10
     &          .or. i == 12 .or. i == 14 ) Resul(j) = - Resul(j)
! On multiplie par (-1)**m
              if( i == 2 .or. i == 3 .or. i == 5 .or. i == 6
     &         .or. i == 10 .or. i == 11 .or. i == 14 .or. i == 15 )
     &                                      Resul(j) = - Resul(j)
! On recupere le img exterieur au tenseur propre au dipole-quadrupole
              Resul(j) = img * Resul(j) 

              Resul(jdq) = Resul(jdq) + Resul(j)

            end do

          elseif( itens > n_tens_dd + n_tens_dq ) then

            if( .not. E2E2 ) cycle
            i = itens - n_tens_dd - n_tens_dq
            j = j + 1

            Resul(j) = Tensor_pol_qq(ipl,i) * Sph_tensor_qq(i)
! Devant le produit, il faut mettre un signe -1 devant les tenseurs
! impairs
            if( ( i >= 2 .and. i <= 4 ) .or. ( i >= 10 .and. i <= 16 ) )
     &          Resul(j) = - Resul(j)
! On recupere le img omis dans les tenseurs imaginaires
            if( i == 4 .or. i == 7 .or. i == 8 .or. i == 12 
     &         .or. i == 13 .or. i == 16 .or. i == 19 .or. i == 20  
     &          .or. i == 23 .or. i == 24 ) Resul(j) = img * Resul(j)
! On multiplie par -1 devant les tenseurs differences
            if( i == 3 .or. i == 6 .or. i == 8 .or. i == 11 .or. i == 13
     &         .or. i == 15 .or. i == 18 .or. i == 20 .or. i == 22
     &         .or. i == 24 ) Resul(j) = - Resul(j)
! On multiplie par (-1)**m
            if( i == 3 .or. i == 4 .or. i == 6 .or. i == 7
     &         .or. i == 11 .or. i == 12 .or. i == 15  .or. i == 16
     &         .or. i == 18 .or. i == 19 .or. i == 22  .or. i == 23 )
     &                       Resul(j) = - Resul(j)

            Resul(jqq) = Resul(jqq) + Resul(j)

          endif

          mot = nomtens(itens)
          l = len_trim( mot )
          jj = jj + 1
          if( itens > n_tens_dd 
     &                   .and. itens <= n_tens_dd + n_tens_dq
     &          .and. ia == 0 .and. magn_sens ) then
            if( ipldafs > 0 ) then
              mot(l+1:l+2) = '_r'
              nomten(jj) = mot
              tens(jj) = real( Resul(j-1),db )
              jj = jj + 1
              mot(l+1:l+2) = '_i'
              nomten(jj) = mot
              tens(jj) = aimag( Resul(j-1) )
              jj = jj + 1
              mot(l+1:l+3) = '_mr'
              nomten(jj) = mot
              tens(jj) = real( Resul(j),db )
              jj = jj + 1
              mot(l+1:l+3) = '_mi'
              nomten(jj) = mot
              tens(jj) = aimag( Resul(j) )
            else
              nomten(jj) = mot
              tens(jj) = real( Resul(j-1),db )
              jj = jj + 1
              mot(l+1:l+2) = '_m'
              nomten(jj) = mot
              tens(jj) = real( Resul(j),db )
            endif
          else
            if( ipldafs > 0 ) then
              mot(l+1:l+2) = '_r'
              nomten(jj) = mot
              tens(jj) = real( Resul(j),db )
              jj = jj + 1
              mot(l+1:l+2) = '_i'
              nomten(jj) = mot
              tens(jj) = aimag( Resul(j) )
            else
              nomten(jj) = mot
              tens(jj) = real( Resul(j),db )
            endif
          endif

        end do

        n_tens = jj

        if( j0 > 1 ) Resul(1) = sum( Resul(2:j0) )

        i = 0
        j = 0
        do it = 1,j0
          i = i + 1
          if( it == jdd ) then
            nomten(i) = 'Sum_dd'
          elseif( it == jdq ) then
            nomten(i) = 'Sum_dq'
          elseif( it == jqq ) then
            nomten(i) = 'Sum_qq'
          else
            nomten(i) = 'Sum_tot'
          endif
          j = j + 1
          Tens(i) = real( Resul(j),db )
          if( ipldafs > 0 ) then
            mot = nomten(i)
            l = len_trim( mot )
            mot(l+1:l+2) = '_r'
            nomten(i) = mot
            i = i + 1
            mot(l+1:l+2) = '_i'
            nomten(i) = mot
            Tens(i) = aimag( Resul(j) )
          endif
        end do

! On donne le Xanes en Megabarns
        if( ipldafs == 0 ) then
          do i = 1,n_tens
             if( abs(Tens(i)) > 1e-20_db ) Tens(i) = Tens(i) / ct_nelec
           end do
        endif

! Ecriture des tenseurs appliques aux reflexions RXS et au xanes

        nomficht = nomfich1

        long = len_trim(nomficht)
        nomficht(long+1:long+11) = '_sph_signal'
        if( ia > 0 ) then
          nomficht(long+12:long+16) = '_atom'
          call ad_number(ia,nomficht,132)
        endif
        long = len_trim(nomficht)
        if( ipldafs > 0 ) then
          nomficht(long+1:long+4) = '_rxs'
          call ad_number(ipldafs,nomficht,132)
        elseif( ipl == 0 ) then
          nomficht(long+1:long+4) = '_xan'
        else
          nomficht(long+1:long+4) = '_pol'
          call ad_number(ipl,nomficht,132)
        endif
        long = len_trim(nomficht)
        nomficht(long+1:long+4) = '.txt'

        if( ipldafs > 0 ) then
          if( ia == 0 ) then
            cf = phdt(ipldafs)
            cg = phdf0t(ipldafs) 
          else
            cf = (1._db,0._db)
            cg = (0._db,0._db)
          endif 
          phtem(:) = cf
          ph0(:) = cg 
          n_tens2 = n_tens / 2
          call write_out(rdum,rdum,Densite_atom,0._db,E_cut,Ephseuil,
     &          Epsii,Eseuil,.false.,idum,ie,Length_word,
     &          jseuil,n_tens_max,n_tens,0,ninitlr,nomficht,nomten,
     &          n_tens_max,n_tens2,0,0,nseuil,numat_abs,phtem,
     &          ph0,Tens,v0muf,.false.,0)
        else
          call write_out(rdum,rdum,Densite_atom,0._db,E_cut,Ephseuil,
     &           Epsii,Eseuil,.false.,idum,ie,Length_word,
     &           jseuil,n_tens_max,n_tens,0,ninitlr,nomficht,nomten,
     &           1,0,0,0,nseuil,numat_abs,cdum,
     &           cdum,Tens,v0muf,.false.,0)
        endif

      end do

      return
      end

