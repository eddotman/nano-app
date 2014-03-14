! FDMNES subroutine

! Sous-ensemble de routines qui servent a la procedure tddft. 
! Pour chaque energie, on calcule les valeurs du noyau K.
! A la fin de la boucle sur les energies, on calcule chi_0 et chi.

!***********************************************************************

! nbseuil : nr de seuil
! ninitl  : nbr d'etat initiaux
! ninitlu = ninitl  en Core_resolved
!         = nbseuil sinon
! ninitlr = ninitlu en DFT
!         = 1       en TDDFT
! ninitlt = 1       en DFT
!         = ninitlu en TDDFT
! ninitlv = nbseuil en DFT
!         = ninitlu en TDDFT
  
      subroutine main_tddft(alfpot,angxyz,Allsite,
     &      Atom_nonsph,Atomic_scr,axyz,Base_spin,BSE,coef_g,
     &      Core_resolved,Dafs,Dafs_bio,Delta_edge,Delta_Eseuil,
     &      Densite_atom,Dipmag,
     &      dv_ex_nex,dv0bdcF,Dyn_eg,Dyn_g,E_cut,E_cut_imp,E_Fermi_man,
     &      E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Ecent,Ef,Efato,Elarg,
     &      Energ_s,Energphot,Epsii_a,Extract,Epsii_moy,Eseuil,Estart,
     &      fpp_avantseuil,Full_atom,Full_potential,Full_self_abs,
     &      Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_int,
     &      Green_plus,hkl_dafs,Hubb_a,Hubb_d,icheck,iaabsi,
     &      iabsorig,iapot,iaprotoi,imag_taull,imoy,
     &      imoy_out,imparite,iprabs,iprabs_nonexc,is_g,isigpi,isymeq,
     &      itabs,itypepr,jseuil,Kern_fac,ldip,length_word,
     &      lmaxabs_t,lmaxat0,lmaxfree,lmoins1,loct,lplus1,lqua,
     &      lseuil,ltypcal,m_g,m_hubb,M1M1,Magnetic,
     &      mix_repr,Moyenne,mpinodes,mpirank,msymdd,msymddi,msymdq, 
     &      msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,
     &      multi_run,n_atom_0,
     &      n_atom_ind,n_atom_proto,n_multi_run,n_oo,n_tens_max,natome,
     &      natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ngamh,ngrph,
     &      ninit1,ninitl,ninitlu,nlmamax,nlmamax_u,nlmsa0,
     &      nlmsam,nomabs,nomfich,nomfich_cal_tddft_conv,
     &      nomfich_s,nomfich_tddft_data,Nonexc_g,
     &      nphi_dafs,nphim,npldafs,nplr,nplrm,npoint,npsom,nptmoy,
     &      nptmoy_out,nr,nrato,nrm,nseuil,nspin,nspino,ntype,
     &      numat,nxanout,Octupole,pdp,phdafs,phdf0t,
     &      phdt,pol,poidsov,poidsov_out,poldafse,poldafss,psii,
     &      Quadrupole,r,rato,Recup_tddft_data,Relativiste,
     &      Repres_comp,rho,rhoato_abs,rhons,Rmtg,Rmtsd,
     &      rof0,rot_atom_abs,Rot_int,RPALF,rs,rsato,rsato_abs,
     &      Self_abs,Solsing_only,Spinorbite,Taux_eq,Tddft_mix,
     &      Tddft_so,tpt,V_intmax,V_hubb,V0muf,Vcato,vec,
     &      vecdafse,vecdafss,Vh,Vhns,Volume_maille,Vxc,Vxcato,Workf,
     &      xyz,Ylm_comp)

      use declarations
      implicit none
      include 'mpif.h'

      integer:: i, iaabsi, iabsorig, icheck_s, ie, ie_computer,
     &  imparite, index_e, initl, ip_max, ip0, iprabs_nonexc, iprabs,
     &  isp, itabs, je, jseuil, l,
     &  length_word, lmax, lmax_probe, lmaxabs_t, lmaxat0,
     &  lseuil, m_hubb, mpinodes, mpirank, multi_run, n_atom_0,
     &  n_atom_ind, n_atom_proto, n_multi_run, n_oo, n_tens_max, n_vr_0,
     &  n_vr_ind, natome, natomsym, nbseuil, ncolm,  
     &  ncolr, ncolt, ndim1, ndim2, nenerg, nenerg_s, ngamh, nge,
     &  ngrph, ninit1, ninitl, ninitlr,
     &  ninitlt, ninitlu, ninitlv, nlm_probe, nlma, nlma2,
     &  nlmam, nlmam2, nlmam_u, nlmamax, nlmamax_u, 
     &  nlmsam, nphim, npldafs, nplr, nplrm, npoint,
     &  npsom, nptmoy, nptmoy_out, nr, nrm, ns_dipmag, nseuil, nspin,
     &  nspin_t, nspino, nspino_t, ntype, numat, nxanout
      integer, dimension(2):: mix_repr     
      integer, dimension(30):: icheck
      integer, dimension(npoint):: imoy, imoy_out
      integer, dimension(natome):: iaprotoi
      integer, dimension(natomsym):: isymeq
      integer, dimension(0:ntype):: nrato
      integer, dimension(0:n_atom_proto):: iapot, itypepr
      integer, dimension(ninitl):: is_g
      integer, dimension(npldafs):: nphi_dafs
      integer, dimension(ninitl,2):: m_g
      integer, dimension(npldafs,2):: isigpi
      integer, dimension(3,npldafs):: hkl_dafs
      integer, dimension(natome,ngrph), intent(in):: nlmsa0
      integer, dimension(nlmsam,natome,ngrph):: lato, mato
      integer, dimension(3):: ldip
      integer, dimension(3,3):: lqua, msymdd, msymddi
      integer, dimension(3,3,3):: loct, msymdq, msymdqi
      integer, dimension(3,3,3,3)::  msymdo, msymdoi, msymqq, msymqqi 
      integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

      parameter(ninitlr=1)

      character(len=132):: nomfich, nomfich_s, nomfich_tddft_data,
     &                     nomfich_cal_convt, nomfich_cal_tddft_conv
      character(len=13), dimension(nplrm):: ltypcal
      character(len=length_word), dimension(ncolm):: nomabs

      complex(kind=db), dimension(3,nplrm):: pol
      complex(kind=db), dimension(natomsym,npldafs):: phdafs 
      complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt 
      complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss 
      complex(kind=db), dimension(nenerg_s,nlmamax,nspin,
     &                                         nspino,nbseuil):: rof0
      complex(kind=db), dimension(3,3,ninitlr,0:mpinodes-1):: secdd,
     &                        secdd_m, secmd, secmd_m, secmm, secmm_m
      complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodes-1):: secdq,
     &                                                        secdq_m
      complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodes-1):: secdo,
     &                                          secdo_m, secqq, secqq_m
      complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1):: 
     &                                              secoo, secoo_m
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                          V_hubb
      complex(kind=db), dimension(:,:), allocatable:: Trans, V
      complex(kind=db), dimension(:,:,:), allocatable:: V_hubb_t
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: rof0_e
      complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: Chi_0
      complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Chi

      logical:: Allsite, Atom_nonsph, Atomic_scr, Base_spin, 
     &  BSE, Core_resolved, Dafs, Dafs_bio, Dipmag, Dyn_eg, Dyn_g,
     &  E_Fermi_man, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Energphot,
     &  Extract, Final_optic, Final_tddft, Full_atom, Full_potential,
     &  Full_self_abs,
     &  Gamma_hole_imp, Gamma_tddft, Green_i, Green_int, Green_plus,
     &  Hubb_a, Hubb_d, lmaxfree, lmoins1, lplus1, M1M1, Magnetic, 
     &  Moyenne, Nonexc_g, Octupole, Quadrupole, Radial_comp, 
     &  Recup_tddft_data, Relativiste, RPALF, Self_abs, Solsing,
     &  Solsing_only, Spinorbite, Spinorbite_t, Tddft_mix, Tddft_so, 
     &  V0muf, Xan_atom, Ylm_comp, Ylm_compa 
      logical, dimension(ngrph):: Repres_comp

      real(kind=sg) time
      
      real(kind=db):: alfpot, Delta_edge, Delta_Eseuil, Densite_atom,
     &   E_cut,
     &   E_cut_imp, E_cut_tddft, Ec_min, Ecent, Ecmax, EFermi_min,
     &   Elarg, Energ_t, Enervide, Epsii_moy, Estart,
     &   fpp_avantseuil, Gamma_max, Kern_fac, Rmtsd, V_intmax,
     &   Volume_maille, Workf
      real(kind=db), dimension(3):: angxyz, axyz
      real(kind=db), dimension(6):: tp
      real(kind=db), dimension(11):: Gamma_hole 
      real(kind=db), dimension(12):: tpt 
      real(kind=db), dimension(nspin):: dv0bdcF, V0bd, V0bd_out
      real(kind=db), dimension(npoint):: Ef, poidsov, poidsov_out, rs,
     &                                   Vh, Vhns, rhons
      real(kind=db), dimension(nenerg_s):: Energ_s
      real(kind=db), dimension(nbseuil):: Eseuil
      real(kind=db), dimension(ninitlu):: Epsii_a
      real(kind=db), dimension(ninitlr):: Epsii, sec_atom
      real(kind=db), dimension(nr):: r, rsato_abs
      real(kind=db), dimension(natomsym) :: Taux_eq
      real(kind=db), dimension(nrm):: dv_ex_nex
      real(kind=db), dimension(0:n_atom_proto):: Rmtg 
      real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int
      real(kind=db), dimension(nplrm,2) :: pdp
      real(kind=db), dimension(3,nplrm) :: vec
      real(kind=db), dimension(ninitl,2):: coef_g
      real(kind=db), dimension(nrm,nbseuil):: psii
      real(kind=db), dimension(nr,nspin):: rhoato_abs
      real(kind=db), dimension(4,npsom):: xyz
      real(kind=db), dimension(n_tens_max,0:natomsym):: Int_tens
      real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: Efato,
     &                                                rsato, Vcato
      real(kind=db), dimension(npoint,nspin):: rho, Vxc
      real(kind=db), dimension(0:nrm,0:ntype):: rato
      real(kind=db), dimension(nr,2,2):: fxc
      real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss 
      real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: Vxcato

      real(kind=db), dimension(nenerg_s,nlmamax_u,nspin,nlmamax_u,nspin)
     &                                                    :: imag_taull
      real(kind=db), dimension(:), allocatable:: Decal_initl,  
     &                               Ecinetic_e, Eimag, Energ, V0bd_t      
      real(kind=db), dimension(:,:), allocatable:: Ecinetic, V0bd_e,
     &                                             Vrato_t, Vr
      real(kind=db), dimension(:,:,:), allocatable:: Vrato, Vrato_e 
      real(kind=db), dimension(:,:,:,:,:,:), allocatable:: zet
      real(kind=db), dimension(:,:,:,:,:,:,:), allocatable::  Kern

      if( icheck(1) > 0 ) write(3,100)
 
      Final_optic = .false.
      Final_tddft = .true.
      sec_atom(:) = 0._db

      if( Tddft_so ) then
        nspino_t = 2
        nspin_t = 2
        Spinorbite_t = .true.
      else
        nspin_t = nspin
        nspino_t = nspino
        Spinorbite_t = Spinorbite
      endif

      n_vr_0 = 1
      n_vr_ind = ninitlu
      ninitlt = ninitlu
      ninitlv = ninitlu

      allocate( Ecinetic(nspin_t,ninitlt) )
      allocate( Ecinetic_e(nspin_t) )
      allocate( V0bd_e(nspin_t,ninitlt) )    

      allocate( Decal_initl(ninitlt) )    

! Epsii_moy est la moyenne pour le 2eme seuil (si 2 seuils)
! par exemple L3 )
      Decal_initl(:) = Epsii_a(:) - Epsii_moy

      Epsii(1) = Epsii_moy

      if( Gamma_tddft ) then
        Green_i = .true.
      else
        Green_i = Green_int
      endif

      Ylm_compa = .true.  ! pour routine tenseur

      if( Hubb_a .and. .not. Ylm_comp .and. Ylm_compa ) then
        l = m_hubb
        allocate( V(-l:l,-l:l) )
        allocate( Trans(-l:l,-l:l) )
        Call cal_trans_lh(l,Trans)
        do isp = 1,nspin
          V(-l:l,-l:l) = V_hubb(-l:l,-l:l,isp) 
!          call trans_rc(l,V)
          V = Matmul( Transpose(Trans), Matmul( V, Conjg(Trans) )) 
          V_hubb(-l:l,-l:l,isp) = V(-l:l,-l:l)  
        end do
        deallocate( Trans, V )
      endif
      allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin_t) )
      do isp = 1,nspin_t
        V_hubb_t(:,:,isp) = V_hubb(:,:,min(isp,nspin)) 
      end do

      if( Recup_tddft_data ) call Read_tddft_data(E_cut,Energ_s,
     &                 imag_taull,mpirank,multi_run,
     &                 n_multi_run,nbseuil,nenerg_s,nlmamax,
     &                 nlmamax_u,nomfich_tddft_data,nspin,nspino,rof0)

! Elaboration de la grille etendue pour la Tddft
      call dim_grille_tddft(Energ_s,Delta_Eseuil,Estart,nbseuil,
     &                            nenerg,nenerg_s)

      allocate( Energ(nenerg) )
      allocate( Eimag(nenerg) )

      call grille_tddft(Energ,Energ_s,Delta_Eseuil,Estart,
     &                        icheck(22),nbseuil,nenerg,nenerg_s)

! On garde l'energie imaginaire nulle
      Eimag(:) = 0._db             
      Solsing = .false.

! Calcul de Chi_0
      allocate( Chi_0(nenerg,nlmamax_u*nspino,nlmamax_u*nspino,
     &                                     nspin,nspin,ninitlt) ) 

      if( E_Fermi_man ) then
        E_cut_tddft = E_cut_imp
      else
        E_cut_tddft = E_cut
      endif

      call Chi_0_int(Chi_0,Core_resolved,Decal_initl,Delta_edge,
     &      E_cut_tddft,Ecent,EFermi_min,Elarg,Energ,Energ_s,
     &      Eseuil(nbseuil),Gamma_hole,Gamma_hole_imp,Gamma_max,
     &      Gamma_tddft,icheck(23),imag_taull,imparite,jseuil,lmaxabs_t,
     &      lseuil,mpinodes,mpirank,nbseuil,nenerg,nenerg_s,
     &      ngamh,ninit1,ninitlt,nlmamax,nlmamax_u,
     &      nspin,nseuil,nspino,numat,rof0,Spinorbite)

! Valeur eventuellement decalee vers le bas pour ne rien couper a la
! convolution.
      E_cut = EFermi_min

      Xan_atom = .false.

! Calcul du noyau xc
      if( .not. RPALF ) call fxcorr(alfpot,fxc,icheck(24),Magnetic,nr,
     &                              nspin,r,rhoato_abs,rsato_abs)

      if( Octupole ) then
        lmax_probe = lseuil + 3
      elseif( Quadrupole ) then
        lmax_probe = lseuil + 2
      else
        lmax_probe = lseuil + 1
      endif
      nlm_probe = (lmax_probe + 1)**2
      if( Dipmag ) then
        ns_dipmag = 4
      else
        ns_dipmag = 1
      endif  
      if( Dipmag ) then
        ip0 = 0
      else
        ip0 = 1
      endif
      if( Octupole ) then
        ip_max = 3
      else
        ip_max = 2
      endif

      allocate( Chi(nlm_probe*nspino_t,nlm_probe*nspino_t,2,2,
     &                  ns_dipmag,ninitl,ninitl) )
      allocate( Vrato(0:nrm,nspin,ninitlt) ) 
      allocate( Vrato_e(nr,nspin_t,ninitlt) ) 

! Boucle sur l'energie
      nge = ( nenerg - 1 ) / mpinodes + 1

      index_e = 0

      boucle_energ: do je = 1,nge

        if( mpirank == 0 ) then
          call CPU_TIME(time)
          tp(1) = real(time,db)
        endif

        ie = ( je - 1 ) * mpinodes + mpirank + 1

        if( ie > nenerg ) goto 1010
 
        index_e = index_e + 1

        icheck_s = max( icheck(16), icheck(22) )
        if( icheck_s == 1 ) icheck_s = 0

        do initl = 1,ninitlt

! Calcul du potentiel dans l'etat excite
          allocate( Vr(npoint,nspin) )
          Enervide = Energ(ie) - Workf  - Decal_initl(initl)
                  
          call potex(Atom_nonsph,axyz,alfpot,.true.,
     &        dv0bdcF,Ef,Efato,Energ(ie),Enervide,Final_tddft,Full_atom,
     &        iaabsi,iapot,iaprotoi,icheck_s,imoy,imoy_out,initl,
     &        iprabs_nonexc,iprabs,itabs,
     &        itypepr,Magnetic,n_atom_0,n_atom_ind,
     &        n_atom_proto,n_vr_0,n_vr_ind,natome,
     &        Nonexc_g,npoint,npsom,nptmoy,
     &        nptmoy_out,nrato,nrm,nspin,ntype,poidsov,
     &        poidsov_out,rato,rho,rhons,Rmtg,rs,rsato,V_intmax,Vcato,
     &        Vh,Vhns,Vr,Vxc,Vrato,Vxcato,V0bd,V0bd_out,xyz)

!Vr ne sert pas en tddft
          deallocate( Vr )

          if( nspin == nspin_t ) then
            V0bd_e(:,initl) = V0bd(:) 
            Vrato_e(1:nr,:,initl) = Vrato(1:nr,:,initl)
          else
            do isp = 1,nspin_t
              V0bd_e(isp,initl) = V0bd(1) 
              Vrato_e(1:nr,isp,initl) = Vrato(1:nr,1,initl)
            end do
          endif
          Ecinetic(:,initl) = Enervide - V0bd_e(:,initl)

        end do

        Ecmax = 0._db
        do initl = 1,ninitlt
          do isp = 1,nspin_t
! L'energie cinetique minimale est prise egale a celle du niveau de
! Fermi. 
            Ec_min = E_cut_tddft - Workf - V0bd_e(isp,initl)
            Ec_min = max( Ec_min, 0.01_db / Rydb ) 
            Ecinetic(isp,initl) = max( Ecinetic(isp,initl), Ec_min )
          end do 
          do isp = 1,nspin_t 
            Ecmax = max( Ecmax, Ecinetic(isp,initl) )
          end do    
        end do    

! En Tddft sur 2 seuils, le premier seuil a plus haute energie peut
! avoir un lmax plus grand (et inutile).
        call clmax(Ecmax,Rmtg(iprabs),lmaxat0,lmax,numat,lmaxfree)
        lmax = min(lmax,lmaxabs_t)
        nlmam = ( lmax + 1 )**2
        nlmam_u = 0
        do l = 0,lmax
          if( (imparite /= mod(l,2)) .and. (imparite /= 2) ) cycle
          nlmam_u = nlmam_u + 2 * l + 1
        end do
        if( Full_potential .or. Hubb_a .and. .not. Hubb_d ) then
          nlmam2 = nlmam
        else
          nlmam2 = 1
        endif

        Radial_comp = Eimag(ie) > eps10 .or. ( Hubb_a .and. Ylm_comp )

        allocate( zet(nr,nlmam,nlmam2,nspino_t,nspin_t,ninitlt) )  
        zet(:,:,:,:,:,:) = 0._db
            
        if( mpirank == 0 ) then
          call CPU_TIME(time)
          tp(2) = real(time,db)
        endif

        icheck_s = max( icheck(18), icheck(22) )

        do initl = 1,ninitlt

          allocate( Vrato_t(nr,nspin_t) )
          allocate( V0bd_t(nspin_t) )

          Ecinetic_e(:) = Ecinetic(:,initl)

          Energ_t = Energ(ie) - Decal_initl(initl)
          Enervide = Energ(ie) - Workf  - Decal_initl(initl)
          Enervide = max( Enervide, E_cut_tddft - Workf )
          Vrato_t(1:nr,:) = Vrato_e(1:nr,:,initl)
          V0bd_t(:) = V0bd_e(:,initl)

          call Radial_wave(Ecinetic_e,Eimag(ie),Energ_t,Enervide,
     &          Full_potential,Hubb_a,Hubb_d,icheck_s,initl,
     &          lmax,m_hubb,ninitlt,nlmam,nlmam2,nr,nspin_t,nspino_t,
     &          numat,r,Radial_comp,Rmtg,Rmtsd,Spinorbite_t,V_hubb_t,
     &          V_intmax,V0bd_t,Vrato_t,zet)

          deallocate( Vrato_t )
          deallocate( V0bd_t )

        end do

        if( mpirank == 0 ) then
          call CPU_TIME(time)
          tp(3) = real(time,db)
        endif

        allocate( Kern(nlmam_u*nspino_t,nlmam_u*nspino_t,2,2,
     &                       ns_dipmag,ninitl,ninitl) )

! Calcul du noyau
        call kernel(Atomic_scr,BSE,coef_g,Core_resolved,Dipmag,
     &              dv_ex_nex,Dyn_eg,Dyn_g,fxc,Kern,Kern_fac,
     &              icheck(24),ie,
     &              imparite,lmax,lseuil,m_g,nbseuil,ninit1,
     &              ninitl,ninitlt,nlmam,nlmam2,nlmam_u,nr,nrm,
     &              ns_dipmag,nspin_t,nspino_t,
     &              psii,r,Rmtsd,RPALF,Spinorbite_t,zet)

        deallocate( zet )

        call Cal_Chi(Chi,Chi_0,coef_g,Core_resolved,Dipmag,
     &            Energ(ie),iaabsi,icheck(25),ie,imparite,Kern,
     &            lato,lmax_probe,
     &            lmax,mato,mix_repr,natome,nenerg,ninit1,ninitl,
     &            ninitlt,ngrph,nlm_probe,nlmam_u,nlmamax_u,nlmsa0,
     &            nlmsam,nomfich,ns_dipmag,nspin,nspino,nspino_t,
     &            Repres_comp,Spinorbite,Tddft_mix,Tddft_so,Ylm_comp)

        deallocate( Kern )

        if( mpirank == 0 ) then
          call CPU_TIME(time)
          tp(4) = real(time,db)
        endif

! Calcul des tenseurs cartesiens

        nlma = nlm_probe
        if( Hubb_a .or. Full_potential ) then
          nlma2 = nlma
        else
          nlma2 = 1
        endif

        ndim1 = ns_dipmag
        ndim2 = ninitl

        allocate( rof0_e(nenerg,nlmamax,nspin_t,nspino_t,nbseuil) )

        icheck_s = max( icheck(22), icheck(20) )

        call tenseur_car(Base_spin,coef_g,
     &            Core_resolved,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Ecinetic,
     &            Eimag(ie),Energ(ie),Enervide,Eseuil,Final_optic,
     &            Final_tddft,
     &            Full_potential,Green_i,Green_plus,Hubb_a,Hubb_d,
     &            icheck_s,ie,ip_max,ip0,is_g,lmax_probe,ldip,
     &            lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb,M1M1,
     &            mpinodes,mpirank,msymdd,msymddi,msymdq, 
     &            msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,
     &            n_oo,nbseuil,ndim1,ndim2,nenerg,ninit1,
     &            ninitl,ninitlr,ninitlt,ninitlv,nlma,nlma2,nlmamax,
     &            nr,nrm,nspin_t,nspino_t,numat,psii,r,
     &            Relativiste,Rmtg,Rmtsd,
     &            rof0_e,rot_atom_abs,Rot_int,secdd,secdd_m,secdo,
     &            secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,
     &            secoo,secoo_m,secqq,secqq_m,Solsing,Solsing_only,
     &            Spinorbite_t,Chi,.true.,V_hubb_t,V_intmax,V0bd_e,
     &            Vrato_e,Ylm_compa)

       deallocate( rof0_e )

 1010   continue

        if( mpirank == 0 ) then        
          call CPU_TIME(time)
          tp(5) = real(time,db)
        endif

        if( mpinodes > 1 ) then

          call MPI_RECV_all(E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,M1M1,
     &                 mpinodes,mpirank,n_oo,ninitlr,secdd,secdo,   
     &                 secdq,secmd,secmm,secoo,secqq)      
          if( Green_i ) call MPI_RECV_all(E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,
     &               M1M1,mpinodes,mpirank,n_oo,ninitlr,secdd_m,secdo_m,   
     &               secdq_m,secmd_m,secmm_m,secoo_m,secqq_m)      

        endif

        if( mpirank /= 0 ) cycle

        do ie_computer = 0,mpinodes-1

          ie = ( je - 1 ) * mpinodes + ie_computer + 1

          if( ie > nenerg ) exit      

          call write_coabs(Allsite,angxyz,axyz,Base_spin,
     &          Core_resolved,Dafs,Dafs_bio,Densite_atom,E_cut,
     &          E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Energ,Energphot,
     &          Extract,Epsii,Eseuil,Final_tddft,fpp_avantseuil,
     &          Full_self_abs,Green_i,Green_plus,hkl_dafs,
     &          iabsorig,icheck(21),ie,ie_computer,
     &          Int_tens,isigpi,isymeq,jseuil,length_word,ltypcal,M1M1,
     &          Moyenne,mpinodes,n_multi_run,n_oo,natomsym,nbseuil,
     &          ncolm,ncolr,ncolt,nenerg,ninit1,ninitlr,nomabs,
     &          nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs,
     &          nphim,nplr,
     &          nplrm,nseuil,nspin_t,numat,nxanout,pdp,phdafs,
     &          phdf0t,phdt,pol,poldafse,poldafss,Rot_int,sec_atom,
     &          secdd,secdd_m,secdq,secdq_m,secdo,secdo_m,
     &          secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,
     &          secqq,secqq_m,Self_abs,Spinorbite_t,Taux_eq,V0muf,
     &          vecdafse,vecdafss,vec,Volume_maille,Xan_atom)

          if( ie == 1 ) nomfich_cal_tddft_conv = nomfich_cal_convt
              
        end do
                
        if( mpirank == 0 ) then        
          call CPU_TIME(time)
          tp(6) = real(time,db)
          do i = 5,9
            tpt(i) = tpt(i) + tp(i-3) - tp(i-4)
          end do
        endif

      end do boucle_energ   ! Fin de la boucle sur l'energie.

      deallocate( Chi_0 )
      deallocate( Chi )    

      deallocate( Decal_initl )
      deallocate( Ecinetic )
      deallocate( Ecinetic_e )
      deallocate( Energ )
      deallocate( Eimag )
      deallocate( V0bd_e )    
      deallocate( V_hubb_t )
      deallocate( Vrato, Vrato_e )    

      return
  100 format(/1x,120('-')//,' Cycle Tddft')
      end

!***********************************************************************

! Sous-programme qui definit les dimensions de la nouvelle gamme 
! d'energie pour le calcul tddft. Sert lorsqu'on a deux seuils qui se
! recouvrent

      subroutine dim_grille_tddft(Energ_s,Delta_Eseuil,Estart,nbseuil,
     &                            nenerg,nenerg_s)

      use declarations
      implicit none
      
      integer,intent(in):: nbseuil, nenerg_s
      integer,intent(out):: nenerg

      real(kind=db),intent(in):: Delta_Eseuil, Estart
      real(kind=db),dimension(nenerg_s),intent(in):: Energ_s
      
      integer ie, je0

      real(kind=db):: Pasdeb
 
      if( nbseuil == 1 ) then

        nenerg = nenerg_s 

      else

        do ie = 1,nenerg_s
! Fonctionne meme si on est en energie de photon
          if( Energ_s(ie) - Energ_s(1) > Delta_Eseuil - eps10 ) exit 
        end do
        je0 = ie - 1
        do ie = 1,nenerg_s
          if(  Delta_Eseuil + Energ_s(ie) > Energ_s(nenerg_s) ) exit
        end do
        nenerg = je0 + ie

      end if

      if( Energ_s(1) > Estart - 1.e-10_db ) then
        pasdeb = 0.5_db / rydb
        nenerg = nenerg + nint( ( Energ_s(1) - Estart ) / pasdeb )
      endif

      return
      end

!***********************************************************************

! Sous-programme qui definit une nouvelle gamme d'energie pour le calcul 
! tddft; elle sert lorsque deux seuils se recouvrent

      subroutine grille_tddft(Energ,Energ_s,Delta_Eseuil,Estart,icheck,
     &                        nbseuil,nenerg,nenerg_s)

      use declarations
      implicit none
      
      integer,intent(in):: icheck, nbseuil, nenerg, nenerg_s

      real(kind=db),intent(in):: Delta_Eseuil, Estart
      real(kind=db),dimension(nenerg_s),intent(in):: Energ_s

      real(kind=db),dimension(nenerg),intent(out):: Energ
    
      integer ie, je0, n

      real(kind=db):: Pasdeb

      if( Energ_s(1) > Estart - 1.e-10_db ) then
        pasdeb = 0.5_db / rydb
        n = nint( ( Energ_s(1) - Estart ) / pasdeb )
        do ie = 1,n
          Energ(ie) = Energ_s(1) + ( ie - n - 1 ) * Pasdeb
        end do
      else
        n = 0 
      endif

      if( nbseuil == 1 ) then

        Energ(n+1:nenerg_s) = Energ_s(1:nenerg_s) 

      else

        do ie = 1,nenerg_s
          Energ(n+ie) = Energ_s(ie)
! Fonctionne meme si on est en energie de photon
          if( Energ_s(ie) - Energ_s(1) > Delta_Eseuil - eps10 ) exit 
        end do
        je0 = n + ie - 1
        do ie = 1,nenerg_s
          Energ( je0 + ie ) = Delta_Eseuil + Energ_s(ie)
          if(  Delta_Eseuil + Energ_s(ie) > Energ_s(nenerg_s) ) exit
        end do         
      end if

      if( icheck > 2 ) then
        write(3,100)
        write(3,110)
        write(3,120) Energ(:)*rydb
      end if

      return
  100 format(/' ---- Grille_tddft -------',100('-'))
  110 format(/'The energy grid for the tddft calculation:')
  120 format(5f13.7)
      end

!***********************************************************************

! Calcul de Chi_0
   
      subroutine Chi_0_int(Chi_0,Core_resolved,Decal_initl,Delta_edge,
     &     EFermi,Ecent,EFermi_min,Elarg,Energ,Energ_s,
     &     Eseuil,Gamma_hole,Gamma_hole_imp,Gamma_max,
     &     Gamma_tddft,icheck,imag_taull,imparite,jseuil,lmaxabs_t,
     &     lseuil,mpinodes,mpirank,nbseuil,nenerg,nenerg_s,
     &     ngamh,ninit1,ninitlt,nlmamax,nlmamax_u,
     &     nspin,nseuil,nspino,numat,rof0,Spinorbite)

      use declarations
      implicit none
      include 'mpif.h'

! Declarations des donnees d'entree:
      integer, intent(in):: icheck, imparite, jseuil, lmaxabs_t, lseuil,  
     &  mpinodes, mpirank, nbseuil, nenerg, nenerg_s, ngamh, ninit1, 
     &  ninitlt, nlmamax, nlmamax_u, nseuil, nspin, nspino, numat

      logical, intent(in):: Core_resolved, Gamma_hole_imp, Gamma_tddft

      real(kind=db), intent(in):: Delta_edge, Ecent, EFermi, Elarg,
     &                            Eseuil, Gamma_max  
      real(kind=db), dimension(ninitlt), intent(in):: Decal_initl
      real(kind=db), dimension(nenerg), intent(in):: Energ
      real(kind=db), dimension(nenerg_s), intent(in):: Energ_s

      character(len=2) ch2

      integer ie, ie_saut, ief, iem, initl, ipr, iseuil, iso, iso1,  
     &  iso2, isp, isp1, isp2, iss1, iss2, je, jnitl, js, l, l1, l2, lm,    
     &  lm1, lm1g, lm2, lm2g, lmv1, lmv2, m, m1, m2, mv1,
     &  mv2, nenerge, nie, njp     

      complex(kind=db),dimension(nenerg_s,nlmamax,nspin,nspino,
     &                                                   nbseuil):: rof0
      complex(kind=db),dimension(nenerg,nlmamax_u*nspino,
     &                     nlmamax_u*nspino,nspin,nspin,ninitlt):: Chi_0 

      logical:: rofp_use, saut, Spinorbite
      logical, dimension(nlmamax,nspin,nspino,nbseuil):: lms_exist

      real(kind=db) alfa, Chi_0_r, Chi_0_i, dch, dde, def,
     &        delta, EFermi_i, EFermi_min, Ephm, Ephoton, fp, fpp0,
     &        f_interp1, f_interp2, f_interp3, num, param, pasmin, t1,
     &        t2, x, x1, x2, x3, x4, y, y1, y2, y3, y4 
      real(kind=db), dimension(10):: Gamma_hole  
      real(kind=db),dimension(nenerg_s,nlmamax_u,nspin,nlmamax_u,nspin)
     &                                                    :: imag_taull

      real(kind=db),dimension(:),allocatable:: e1, e2, gamma        
      real(kind=db),dimension(:),allocatable:: Energe
      real(kind=db),dimension(:),allocatable:: fpp
      real(kind=db),dimension(:,:,:,:,:,:,:),allocatable:: fppn

      Chi_0(:,:,:,:,:,:) = (0._db,0._db)

      if( .not. Gamma_tddft )  then
        pasmin = 1.0_db / rydb
! le pas minimum
        do ie = 1,nenerg_s-1
          pasmin = min( pasmin, Energ_s(ie+1)-Energ_s(ie) )
        end do
        param = 0.5_db * pasmin
      end if

      if( icheck > 0 ) write(3,100)

      Ephm = Eseuil + 10000._db / rydb
      if( nenerg_s > 1 ) then
        def = Energ_s(nenerg_s) - Energ_s(nenerg_s-1)
      else
        def = 2 / rydb
      endif

      alfa = 1.02_db
      dde = def

      Ephoton = Eseuil + Energ_s(nenerg_s)
      do ie = 1,10000
        dde = alfa * dde
        Ephoton = Ephoton + dde 
        if( Ephoton > Ephm ) exit
      end do
      njp = ie - 1

! On definit les gammes etendues en energie 
      nenerge = nenerg_s + njp
      allocate( e1(nenerge) )
      allocate( e2(nenerge) )
      allocate( Energe(nenerge) )
      allocate( gamma(nenerge) )
      allocate( fppn(nenerg_s-1:nenerge,nlmamax_u,nspino,nlmamax_u,
     &                nspino,nspin,nspin) )
      allocate( fpp(nenerg_s-1:nenerge) )

      Ephoton = Eseuil - 2 / rydb
      call fprime(numat,Ephoton,fpp0,fp)  ! fpp0 = avantseuil de f"

      Energe(1:nenerg_s) = Energ_s(1:nenerg_s) 
      dde = def
      do ie = nenerg_s-1,nenerge 
        if( ie > nenerg_s ) then
          dde = alfa * dde
          Energe(ie) = Energe(ie-1) + dde
        endif

        Ephoton = Energe(ie) + Eseuil
! fpp(ie) = f'' atomique pour une energie du spectre etendu
        call fprime(numat,Ephoton,fpp(ie),fp)
        fpp(ie) = fpp(ie) - fpp0
      end do                
      
      EFermi_min = EFermi
      do initl = 1,ninitlt
        delta = decal_initl(initl) 
        if( abs(Delta_edge) > eps10 .and. 
     &      ( ( Core_resolved .and. initl <= ninit1 ) .or.
     &        ( .not. Core_resolved .and. initl < ninitlt ) ) )
     &                         delta = delta + Delta_edge
        EFermi_min = min( EFermi_min, EFermi + delta )
      end do 

! Orbitales existantes (cas du spin-orbite)
      lms_exist(:,:,:,:) = abs( rof0(nenerg_s,:,:,:,:) ) > eps10

      nie = ( ninitlt - 1 ) / mpinodes + 1

      boucle_seuil: do jnitl = 1,nie

        initl = ( jnitl - 1 ) * mpinodes + mpirank + 1

        if( initl > ninitlt ) goto 1010

        if( mpinodes == 1 ) then
          write(6,102) initl, ninitlt
        else
          write(6,103) initl, ninitlt, mpirank
        endif

        if( Core_resolved ) then
          if( initl > ninit1 ) then
            iseuil = 2
          else
            iseuil = 1
          endif
        else
          iseuil = initl
        endif

! Energe est la gamme interne shiftee
        delta = decal_initl(initl)
        if( abs(Delta_edge) > eps10 .and. 
     &      ( ( Core_resolved .and. initl <= ninit1 ) .or.
     &        ( .not. Core_resolved .and. initl < ninitlt ) ) )
     &                         delta = delta + Delta_edge

        EFermi_i = EFermi + delta  

        Energe(:) = Energe(:) + delta

! Elargissement
        Gamma(:) = 0._db
        if( Gamma_tddft )  then
        
          if( Gamma_max > eps10 ) call gammarc(Ecent,Elarg,Gamma_max, 
     &                                    EFermi_i,nenerge,Energe,Gamma)

          if( Gamma_hole_imp ) then
            if( ngamh == 1 ) then
              Gamma(:) = Gamma(:) + Gamma_hole(1)
            elseif( ngamh == ninitlt ) then
              Gamma(:) = Gamma(:) + Gamma_hole(initl)
            elseif( initl <= ninit1 ) then
              Gamma(:) = Gamma(:) + Gamma_hole(1)
            else
              Gamma(:) = Gamma(:) + Gamma_hole(2)
            endif
          else
            js = jseuil + iseuil - 1
            call tab_width(Eseuil,Gamma_hole(1),js,nseuil,numat)
            Gamma(:) = Gamma(:) + Gamma_hole(1)
          end if
          if( icheck > 0 ) write(3,105) initl, iseuil,
     &                                   Gamma_hole(1) * rydb

 ! On prend la mi-largeur
          Gamma(:) = Gamma(:) / 2

        elseif( initl == 1 ) then

          if( icheck > 0 ) write(3,106 )

        end if

! Extrapolation couvrant la gamme etendue d'energie

        fppn(:,:,:,:,:,:,:) = 0._db

        lm1 = 0
        do l1 = 0,lmaxabs_t
          if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
          do m1 = -l1,l1
            lm1 = lm1 + 1
                
            lm2 = 0
            do l2 = 0,lmaxabs_t
              if( mod(l2,2) /= imparite .and. imparite /= 2 ) cycle
              do m2 = -l2,l2
                lm2 = lm2 + 1

                if( l1 == lseuil+1 .and. lm1 == lm2 )  then

                  saut = .false.
                  ie_saut = -100
                  fppn(nenerg_s-1,lm1,:,lm2,:,:,:) = fpp(nenerg_s-1)                
                  fppn(nenerg_s,lm1,:,lm2,:,:,:) = fpp(nenerg_s)                

                  big_loop: do ie = nenerg_s+1,nenerge

                    if( .not. saut ) then
                      fppn(ie,lm1,:,lm2,:,:,:) = fpp(ie)                
                      if( fppn(ie,lm1,1,lm2,1,1,1)
     &                       < fppn(ie-1,lm1,1,lm2,1,1,1)) cycle
! Si f'' a un saut (<=> on tombe sur un seuil voisin) on le detecte et
! on l'efface par continuite
                      saut = .true.                          
                      ie_saut = ie
                    endif
                    if( saut) then
                      x  = Energe(ie)
                      x1 = Energe(ie_saut-2)
                      y1 = fppn(ie_saut-2,lm1,1,lm2,1,1,1)
                      x2 = Energe(ie_saut-1)
                      y2 = fppn(ie_saut-1,lm1,1,lm2,1,1,1)
                      y = f_interp1(x,x1,x2,y1,y2)
! On extrapole seulement si f''reste positif. 
                      if( y < 0._db ) exit big_loop
                      fppn(ie,lm1,:,lm2,:,:,:) = y
                    end if

                  end do big_loop

                endif

! Normalisation de la partie extrapole 
                do iso1 = 1,nspino
                  do isp1 = 1,nspin
                    do iso2 = 1,nspino
                      do isp2 = 1,nspin
                        if( Spinorbite) then
                          iss1 = iso1; iss2 = iso2
                          mv1 = m1 + isp1 - iso1 
                          mv2 = m2 + isp2 - iso2
                          if( abs(mv1) > l1 .or. abs(mv2) > l2 ) cycle 
                        else
                          if( isp1 /= isp2 ) cycle
                          iss1 = isp1; iss2 = isp2
                          mv1 = m1
                          mv2 = m2
                        endif
                        lmv1 = l1**2 + l1 + 1 + mv1
                        lmv2 = l2**2 + l2 + 1 + mv2
                        if( lms_exist(lmv1,isp1,iso1,iseuil) .and.
     &                      lms_exist(lmv2,isp2,iso2,iseuil) ) then

                          if( l1 == lseuil+1 .and. lm1 == lm2 )  then
                            num = imag_taull(nenerg_s,lm1,iss1,lm2,iss2)
     &                       * real(rof0(nenerg_s,lmv1,isp1,iso1,iseuil)
     &                       * rof0(nenerg_s,lmv2,isp2,iso2,iseuil) ) 
     &                       /fppn(nenerg_s,lm1,iso1,lm2,iso2,isp1,isp2) 

                            do ie = nenerg_s+1,nenerge 
                              fppn(ie,lm1,iso1,lm2,iso2,isp1,isp2)
     &                      = num * fppn(ie,lm1,iso1,lm2,iso2,isp1,isp2)
                            end do
                          else
                            y1 =imag_taull(nenerg_s-1,lm1,iss1,lm2,iss2)
     &                      *real(rof0(nenerg_s-1,lmv1,isp1,iso1,iseuil)
     &                      * rof0(nenerg_s-1,lmv2,isp2,iso2,iseuil) ) 
                            y2 = imag_taull(nenerg_s,lm1,iss1,lm2,iss2)
     &                       * real(rof0(nenerg_s,lmv1,isp1,iso1,iseuil)
     &                       * rof0(nenerg_s,lmv2,isp2,iso2,iseuil) ) 
                            if( y1 > y2 ) then    
! Si fppn decroit on extrapole lineairement
                              x1 = Energe(nenerg_s-1)
                              x2 = Energe(nenerg_s)  
                              do ie = nenerg_s,nenerge
                                x = Energe(ie)
                                y = f_interp1(x,x1,x2,y1,y2)
                                fppn(ie,lm1,iso1,lm2,iso2,isp1,isp2)
     &                            = max( y, 0._db )
                              end do
                            else
! Si fppn croit on prolonge par une constante
                              fppn(nenerg_s:nenerge,lm1,iso1,lm2,iso2,
     &                              isp1,isp2) = y2
                            end if
                          endif
                        else
                          fppn(:,lm1,iso1,lm2,iso2,isp1,isp2) = 0._db
                        endif

                      end do
                    end do
                  end do
                end do

              end do
            end do
          end do
        end do

! Indice du niveau de Fermi
        if( EFermi_i < Energe(1) .or. EFermi_i > Energe(nenerg_s) ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,500) EFermi*rydb, Energ_s(1)*rydb, 
     &                     Energ_s(nenerg_s)*rydb
          end do
          stop
        end if
        ief = 0
        do ie = 1,nenerg_s-1
          if( Energe(ie) > EFermi_i  ) exit
        end do
        ief = ie

! Creation des tableaux bornes des intervales:  
        do ie = 2,nenerge-1
          e1(ie) =  0.5_db * ( Energe(ie-1) + Energe(ie) )
          e2(ie) =  0.5_db * ( Energe(ie) + Energe(ie+1) )
        end do

! Bornes de l'intervale:
        e1(1) = Energe(1) - 0.5_db * ( Energe(2) - Energe(1) )
        e2(1) = 0.5_db * ( Energe(1) + Energe(2) )
        e1(nenerge) =  0.5_db * ( Energe(nenerge-1) + 
     &                                      Energe(nenerge) )
        e2(nenerge) = Energe(nenerge) + 0.5_db * ( Energe(nenerge-1) -
     &                                      Energe(nenerge) )

        e1(ief) = Efermi_i

! Boucle sur les energies des photons
        boucle_energ: do je = 1,nenerg  

          lm1 = 0
          do l1 = 0,lmaxabs_t
            if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
            do m1 = -l1,l1
              lm1 = lm1 + 1
              do iso1 = 1,nspino
                lm1g = nlmamax_u*(iso1 - 1 ) + lm1
                do isp1 = 1,nspin
                  if( Spinorbite ) then
                    mv1 = m1 + isp1 - iso1
                    if( abs(mv1) > l1 ) cycle 
                    iss1 = iso1
                  else
                    mv1 = m1
                    iss1 = isp1
                  endif
                  lmv1 = l1**2 + l1 + 1 + mv1
                  if(.not. lms_exist(lmv1,isp1,iso1,iseuil)) cycle
                  lm2 = 0
                  do l2 = 0,lmaxabs_t
                    if(mod(l2,2) /= imparite .and. imparite /= 2 ) cycle
                    do m2 = -l2,l2
                      lm2 = lm2 + 1
                      do iso2 = 1,nspino
                        lm2g = nlmamax_u*(iso2 - 1 ) + lm2
                        do isp2 = 1,nspin
                          if( Spinorbite ) then
                            mv2 = m2 + isp2 - iso2
                            if( abs(mv2) > l2 ) cycle 
                            iss2 = iso2
                          else
                            mv2 = m2
                            iss2 = isp2
                          endif
                          lmv2 = l2**2 + l2 + 1 + mv2
                          if(.not.lms_exist(lmv2,isp2,iso2,iseuil))cycle

                          if( .not. Spinorbite .and. isp1 /= isp2) cycle

                          Chi_0_r = 0._db ; Chi_0_i = 0._db   
                
                          do ie = ief,nenerge                                      

                            t1 = Energ(je) - e1(ie)  
                            t2 = Energ(je) - e2(ie) 
 
                            if( ie > nenerg_s ) then
                              dch = fppn(ie,lm1,iso1,lm2,iso2,isp1,isp2)
                            else
                              dch = imag_taull(ie,lm1,iss1,lm2,iss2)
     &                      * real( rof0(ie,lmv1,isp1,iso1,iseuil)
     &                      * rof0(ie,lmv2,isp2,iso2,iseuil) )
                            endif

                            if( Gamma_tddft ) then
                              t1 = t1 / Gamma(ie)
                              t2 = t2 / Gamma(ie)       
                              Chi_0_r = Chi_0_r - dch * 0.5_db
     &                                * log( (1+t2**2)/(1+t1**2) ) / pi
                              Chi_0_i = Chi_0_i
     &                              - dch * ( atan(t2) - atan(t1) ) / pi

                            else    ! gamma = 0

! Eviter les divergences qd les gammes se chevauchent
                              if( abs(t1) < param ) t1 =  param
                              if( abs(t2) < param ) t2 = -param
! On elimine la divergence en 1/e du f' (t1 -> -t1;
! si jamais |t1| = |t2| les contributions des deux cotes du point
! s'annulent.
! abs change le signe quand t1 et t2 sont de signes opposes.
                              if( abs(t1+t2) > eps10 ) Chi_0_r
     &                          = Chi_0_r - dch * log( abs(t2/t1) ) / pi

                              if( Energe(ie) < Energ(je) - eps10 
     &                              .or. Energ(je) < Efermi_i ) cycle
! Partie imaginaire, on interpole quand on depasse le point courant
                              if( abs( Energe(ie) - Energ(je) ) < eps10
     &                             .or. ie == 1 ) then
                                Chi_0_i = dch
                              elseif( Energe(ie-1) < Energ(je) ) then
                                if( ie-1 > nenerg_s ) then
                                  y1 = 
     &                            fppn(ie-1,lm1,iso1,lm2,iso2,isp1,isp2)
                                else
                                  y1 =imag_taull(ie-1,lm1,iss1,lm2,iss2)
     &                      *real(rof0(ie-1,lmv1,isp1,iso1,iseuil)
     &                      * rof0(ie-1,lmv2,isp2,iso2,iseuil) )
                                endif
                                x  = Energ(je)
                                x1 = Energe(ie-1)
                                x2 = Energe(ie)
                                Chi_0_i = f_interp1(x,x1,x2,y1,dch)
                              endif

                            end if

                          end do

                          Chi_0(je,lm1g,lm2g,isp1,isp2,initl) 
     &                                     = cmplx(Chi_0_r,Chi_0_i,db)

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

          do ie = 3,nenerg_s-1
            if( Energe(ie) > Energ(je) ) exit
          end do
          if( ie < ief ) then
            rofp_use = .true.
          else
            rofp_use = .false.
          endif
          lm1 = 0
          do l1 = 0,lmaxabs_t
            if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
            do m1 = -l1,l1
              lm1 = lm1 + 1
              do iso1 = 1,nspino
                lm1g = nlmamax_u * (iso1 - 1 ) + lm1
                do isp1 = 1,nspin
                  if( Spinorbite ) then
                    mv1 = m1 + isp1 - iso1
                    if( abs(mv1) > l1 ) cycle 
                    iss1 = iso1
                  else
                    mv1 = m1
                    iss1 = isp1
                  endif
                  lmv1 = l1**2 + l1 + 1 + mv1
                  if(.not. lms_exist(lmv1,isp1,iso1,iseuil)) cycle

                  lm2 = 0
                  do l2 = 0,lmaxabs_t
                    if( mod(l2,2) /= imparite .and. imparite /= 2) cycle
                    do m2 = -l2,l2
                      lm2 = lm2 + 1
                      do iso2 = 1,nspino
                        lm2g = nlmamax_u * (iso2 - 1 ) + lm2
                        do isp2 = 1,nspin
                          if( Spinorbite ) then
                            mv2 = m2 + isp2 - iso2
                            if( abs(mv2) > l2 ) cycle 
                            iss2 = iso2
                          else
                            mv2 = m2
                            iss2 = isp2
                          endif
                          lmv2 = l2**2 + l2 + 1 + mv2
                          if( .not. lms_exist(lmv2,isp2,iso2,iseuil) )
     &                                                            cycle

                          if(.not. Spinorbite .and. isp1 /= isp2) cycle

                          if( rofp_use ) then 
                            y = real( rof0(ief,lmv1,isp1,iso1,iseuil)
     &                        * rof0(ief,lmv2,isp2,iso2,iseuil) )
                          else
                            iem = max(ief,ie-2)
                            y1 = real( rof0(iem,lmv1,isp1,iso1,iseuil)
     &                         * rof0(iem,lmv2,isp2,iso2,iseuil) )
                            iem = max(ief,ie-1)
                            y2 = real( rof0(iem,lmv1,isp1,iso1,iseuil)
     &                         * rof0(iem,lmv2,isp2,iso2,iseuil) )
                            y3 = real( rof0(ie,lmv1,isp1,iso1,iseuil)
     &                          * rof0(ie,lmv2,isp2,iso2,iseuil) )
                            x1 = Energe(ie-2)
                            x2 = Energe(ie-1)
                            x3 = Energe(ie)
                            x  = Energ(je)
                            if( ie == nenerg_s ) then
                              y = f_interp2(x,x1,x2,x3,y1,y2,y3)
                            else
                              y4 = real(rof0(ie+1,lmv1,isp1,iso1,iseuil)
     &                           * rof0(ie+1,lmv2,isp2,iso2,iseuil) )
                              x4 = Energe(ie+1)
                              y = f_interp3(x,x1,x2,x3,x4,y1,y2,y3,y4)
                            endif
                          endif
                          if( abs(y) > 1.e-20_db )
     &                      Chi_0(je,lm1g,lm2g,isp1,isp2,initl)
     &                        = Chi_0(je,lm1g,lm2g,isp1,isp2,initl) / y 
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

        end do boucle_energ 

        if( icheck > 2 ) then
          if( imparite == 2 ) then
            write(3,110) initl, (((( l,ch2(m), iso, isp, isp = 1,nspin), 
     &                    iso = 1,nspino), m = -l,l), l = 0,lmaxabs_t )
          elseif( imparite == 1 ) then
            write(3,110) initl, (((( l,ch2(m), iso, isp, isp = 1,nspin), 
     &                    iso = 1,nspino), m = -l,l), l = 1,lmaxabs_t,2)
          else
            write(3,110) initl, (((( l,ch2(m), iso, isp, isp = 1,nspin), 
     &                    iso = 1,nspino), m = -l,l), l = 0,lmaxabs_t,2)
          endif 
          do ie = nenerg_s+1, nenerge
            write(3,120) Energe(ie)*rydb, 
     &       ((( fppn(ie,lm,iso,lm,iso,isp,isp), isp = 1,nspin),  
     &                   iso = 1,nspino), lm = 1,nlmamax_u)
          end do
        end if

        Energe(:) = Energe(:) - delta

 1010   continue

        if( mpinodes > 1 ) call MPI_Bcast_Chi_0(Chi_0,jnitl,mpinodes,
     &                    mpirank,nenerg,ninitlt,nlmamax_u,nspin,nspino)

      end do boucle_seuil

! Ecriture des termes diagonaux dans le fichier bavard:
      if( icheck > 1 ) then
        do initl = 1,ninitlt
          if( imparite == 2 ) then
            write(3,130) initl, (((( l,ch2(m), iso, isp, isp = 1,nspin), 
     &                    iso = 1,nspino), m = -l,l), l = 0,lmaxabs_t )
          elseif( imparite == 1 ) then
            write(3,130) initl, (((( l,ch2(m), iso, isp, isp = 1,nspin), 
     &                    iso = 1,nspino), m = -l,l), l = 1,lmaxabs_t,2)
          else
            write(3,130) initl, (((( l,ch2(m), iso, isp, isp = 1,nspin), 
     &                    iso = 1,nspino), m = -l,l), l = 0,lmaxabs_t,2)
          endif 
          do ie = 1,nenerg
            write(3,120) Energ(ie)*rydb,
     &       ((( Chi_0(ie,nlmamax_u*(iso-1)+lm,nlmamax_u*(iso-1)+lm,isp,
     &               isp,initl), isp = 1,nspin), 
     &          iso = 1,nspino), lm = 1,nlmamax_u ) 
          end do
        end do
      end if

      if( icheck > 2 ) then
        if( imparite == 2 ) then
          write(3,140) ((( l, ch2(m), isp, isp = 1,nspin),  m = -l,l), 
     &                                               l = 0,lmaxabs_t )
        elseif( imparite == 1 ) then
          write(3,140) ((( l, ch2(m), isp, isp = 1,nspin),  m = -l,l), 
     &                                               l = 1,lmaxabs_t,2)
        else
          write(3,140) ((( l, ch2(m), isp, isp = 1,nspin),  m = -l,l),
     &                                               l = 0,lmaxabs_t,2)
        endif 
        do ie = 1,nenerg_s
          write(3,120) Energ_s(ie)*rydb, 
     &                (( imag_taull(ie,lm,isp,lm,isp),
     &                   isp = 1,nspin), lm = 1,nlmamax_u ) 
        end do

        write(3,150) (((( l, ch2(m), isp, iseuil, iseuil = 1,nbseuil),  
     &                   isp = 1,nspin), m = -l,l), l = 0,lmaxabs_t )   
        do ie = 1,nenerg_s
          write(3,120) Energ_s(ie)*rydb,
     &             (( real( rof0(ie,lm,isp,min(isp,nspino),1:nbseuil)),
     &               isp = 1,nspin), lm = 1,nlmamax ) 
        end do

      end if

      deallocate( e1 )
      deallocate( e2 )
      deallocate( Energe )   
      deallocate( fpp )
      deallocate( fppn )
      deallocate( Gamma )   

      return
  100 format(/' ---- Chi_0_int -------',100('-'))
  102 format('   Chi_0 calculation, state =',i2,' on',i2)
  103 format('   Chi_0 calculation, state =',i2,' on',i2,', computer,',
     &       i3)
  105 format(/' initl =',i2,', iseuil =',i2,', Gamma_hole =',f7.3,' eV')
  106 format(/' Gamma = 0')
  110 format(/' fpp(l,m,iso,isp) for diagonal terms, initl',i2,/
     & '   Energy ',250(1x,'(',i1,',',a2,2(',',i1),')') )
  120 format(f9.3,1p,500e11.3)
  130 format(/' Chi_0(l,m,iso,isp,l,m,iso,isp) for diagonal terms,',
     & ' initl ',i2,/
     & '   Energy ',250('(',i1,',',a2,',',i1,',',i1,')',4x,'Im',6x) )
  140 format(/' imag_taull(l,m,isp,l,m,isp) for diagonal terms',/
     & '   Energy ',250(2x,'(',i1,',',a2,',',i1,') ') )
  150 format(/' Monopole rof0(l,m,isp,iso,iseuil) for iso = isp',/
     & '   Energy ',250('(',i1,',',a2,',',i1,',',i1,') ') )
  500 format(//' Error: the Fermi energy is out of the energy ',
     &         'calculation range !',/
     &          5x,'EFermi =',f10.3,' eV',/
     &          5x,'E_min  =',f10.3,' eV,  E_max  =',f10.3,' eV',//
     &         ' It is not possible in the TDDFT part !',/)
      end

!***********************************************************************

      character(len=2) function ch2(m)

      integer m

      if( m < 0 ) then
        ch2(1:1) = '-'
      elseif( m == 0 ) then
        ch2(1:1) = '.'
      else
        ch2(1:1) = '+'
      endif
      ch2(2:2) = achar( abs(m) + 48 )
 
      return
      end

!***********************************************************************

      subroutine MPI_Bcast_Chi_0(Chi_0,jnitl,mpinodes,mpirank,
     &                    nenerg,ninitlt,nlmamax_u,nspin,nspino)

      use declarations
      implicit real(kind=db) (a-h,o-z)
      include 'mpif.h'

      integer:: ie_computer, initl, jnitl, mpinodes, mpirank, ndim, 
     &          nenerg, ninitlt, nlmamax_u, nspin, nspino

      complex(kind=db),dimension(nenerg,nlmamax_u*nspino,
     &                     nlmamax_u*nspino,nspin,nspin,ninitlt):: Chi_0 

      real(kind=db),dimension(nlmamax_u*nspino,
     &                nlmamax_u*nspino,nspin,nspin):: Chi_0_i, Chi_0_r

      ndim = ( nlmamax_u * nspino * nspin )**2 

      do ie_computer = 0,mpinodes-1

        initl = ( jnitl - 1 ) * mpinodes + ie_computer + 1

        if( initl > ninitlt ) exit      

        do ie = 1,nenerg

          if( ie_computer == mpirank ) then
            Chi_0_i(:,:,:,:) = aimag( Chi_0(ie,:,:,:,:,initl) ) 
            Chi_0_r(:,:,:,:) = real( Chi_0(ie,:,:,:,:,initl), db)
          endif 

          call MPI_Bcast(Chi_0_i,ndim,MPI_REAL8,ie_computer,
     &                                           MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(Chi_0_r,ndim,MPI_REAL8,ie_computer,
     &                                           MPI_COMM_WORLD,mpierr)

          call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

          if( ie_computer /= mpirank )
     &    Chi_0(ie,:,:,:,:,initl)
     &             = cmplx( Chi_0_r(:,:,:,:), Chi_0_i(:,:,:,:),db ) 

        end do
      end do

      return
      end

!***********************************************************************

! Sous-programme qui calcule le noyau pour le calcul TDDFT

! Le developement en harmoniques spheriques du potentiel coulombien est
! pris d'apres J.Phys.A:Math. Gen. 39(2006) 8613-8630

! On le calcule dans le cas harmoniques complexe.
! En cas de calcul avec harmoniques
! reelles, Chi_0 est transforme dans la routine Cal_chi.

      subroutine kernel(Atomic_scr,BSE,coef_g,Core_resolved,Dipmag,
     &              dv_ex_nex,
     &              Dyn_eg,Dyn_g,fxc,Kern,Kern_fac,icheck,ie,
     &              imparite,lmax_t,lseuil,m_g,nbseuil,ninit1,
     &              ninitl,ninitlt,nlmam,nlmam2,nlmam_u,nr,nrm,
     &              ns_dipmag,nspin,nspino,
     &              psii,r,Rmtsd,RPALF,Spinorbite,zet)

      use declarations
      implicit none

      integer,intent(in):: icheck, ie, imparite, lmax_t, lseuil,     
     &   nbseuil, ninit1, ninitl, ninitlt, nlmam, nlmam2, nlmam_u, 
     &   nr, nrm, ns_dipmag, nspin, nspino
      integer,dimension(ninitl,2),intent(in):: m_g

      logical,intent(in):: atomic_scr, BSE, Dipmag, Dyn_eg, Dyn_g, 
     &                                                         RPALF

      real(kind=db),intent(in):: Kern_fac, Rmtsd
      real(kind=db),dimension(nrm),intent(in):: dv_ex_nex
      real(kind=db),dimension(nr),intent(in):: r
      real(kind=db),dimension(nrm,nbseuil),intent(in):: psii
      real(kind=db),dimension(nr,2,2),intent(in):: fxc
      real(kind=db),dimension(ninitl,2),intent(in):: coef_g
      real(kind=db), dimension(nr,nlmam,nlmam2,nspino,nspin,
     &                                ninitlt), intent(in):: zet

      complex(kind=db):: gaunttd
      
      integer is, is1, is2, isf1, isf2, isg1, isg12, isg2, iso1, iso2,  
     &  isp1, isp2, ist1, ist2, iz1, iz2, l0, l1, l2, lcut, lg,   
     &  lm1, lm2, lm1g, lm2g, lmp1, lmp2, lmv1, lmv2, lp1, lp2,     
     &  m0, m1, m2, mg1, mg2, mp1, mp2, mv1, mv2, pwr1, pwr2

      logical Core_resolved, Spinorbite, Ylm_comp, Zero_term

      real(kind=db):: angl, angl1, angl2, f_integr3, fac, Gaunt4Y,
     &                intrad_r, K_BSE 
      real(kind=db), dimension(nr):: f, t1, t2
      real(kind=db),dimension(nlmam_u*nspino,nlmam_u*nspino,2,2,
     &                               ns_dipmag,ninitl,ninitl):: Kern

      if( ( icheck > 1 .and. ie == 1 ) .or. icheck > 2 ) then
        write(3,100)
        write(3,110) lmax_t, nlmam_u
      endif

      Kern(:,:,:,:,:,:,:) = 0._db
      Ylm_comp = .true.

! l'etat initial est calcule en corestate:
      lg = lseuil

! L'indice 0 porte sur le developpement du potentiel coulombien, et
! l'indice g sur l'etat de coeur
! psii est reel

! le cut-off: voir les regles de selection pour les coef de gaunt

      lcut = lg + lmax_t

      do l0 = 0,lcut
        do m0 = -l0,l0
! Boucle sur les etats initiaux

          do ist1 = 1,ninitl       

            if( ist1 > ninit1 ) then
              is1 = 2
            else
              is1 = 1
            end if
            if( Core_resolved ) then
              iz1 = ist1
            else
              iz1 = is1
            endif

            do ist2 = 1,ninitl       

              if( ist2 > ninit1 ) then
                is2 = 2
              else
                is2 = 1
              end if
              if( Core_resolved ) then
                iz2 = ist2
              else
                iz2 = is2
              endif

              do isf1 = 1,2         ! Spin du premier etat final

! Pour traiter le cas des calculs non magnetiques
                isp1 = min(isf1,nspin)

                do isf2 = 1,2       ! Spin du deuxieme etat final
! noyau diagonal en spin
!                  if( isf2 /= isf1 ) cycle
                  isp2 = min(isf2,nspin)

                  lm1 = 0
                  do l1 = 0,lmax_t 
                    if( mod(l1,2) /= imparite .and. imparite /= 2) cycle
                    do m1 = -l1,l1
                      lm1 = lm1 + 1

! Si Gaunt vaut zero ce n'est pas la peine de calculer la suite
                      Zero_term = .true.
                      boucle_isg1: do isg1 = 1,2
                        if( .not. Dipmag .and. isg1 /= isf1 ) cycle
                        if( abs( coef_g(ist1,isg1) ) < eps10) cycle 
                        mg1 = m_g(ist1,isg1)
                        do iso1 = 1,nspino 
                          if( Spinorbite .and. nlmam2 == 1 ) then
                            mv1 = m1 + isf1 - iso1
                            if( abs(mv1) > l1 ) cycle
                          else
                            mv1 = m1
                          endif
                          angl1 = gaunttd(l1,mv1,l0,m0,lg,mg1,Ylm_comp)
                          if( abs( angl1 ) < eps10 ) cycle
                          Zero_term = .false.
                          exit boucle_isg1
                        end do
                      end do  boucle_isg1
                      if( Zero_term ) cycle

! boucle cas non spherique
                      lmp1 = 0
                      do lp1 = 0,lmax_t        

                        if( mod(lp1,2) /= imparite 
     &                                .and. imparite /= 2) cycle
                        do mp1 = -lp1,lp1
                          if( nlmam2 == 1 .and.
     &                      ( lp1 /= l1 .or. mp1 /= m1 ) ) cycle 
                            lmp1 = lmp1 + 1

                      do iso1 = 1,nspino 

                        if( nlmam2 == 1 ) then             
                          lm1g = nlmam_u * ( iso1 - 1) + lm1
                        else
                          lm1g = nlmam_u * ( iso1 - 1) + lmp1
                        endif 
                        if( Spinorbite .and. nlmam2 == 1 ) then
                          mv1 = m1 + isp1 - iso1
                          if( abs(mv1) > l1 ) cycle
                        else
                          mv1 = m1
                        endif
                        lmv1 = l1**2 + l1 + 1 + mv1

                        lm2 = 0
                        do l2 = 0,lmax_t        

                          if( mod(l2,2) /= imparite .and. imparite /= 2)
     &                                                             cycle
                          do m2 = -l2,l2
                            lm2 = lm2 + 1

                            Zero_term = .true.
                            boucle_isg2: do isg2 = 1,2
                              if( .not. Dipmag .and. isg2 /= isf2 )cycle
                              if( abs( coef_g(ist2,isg2) ) < eps10)cycle 
                              mg2 = m_g(ist2,isg2)
                              do iso2 = 1,nspino
                                if( Spinorbite .and. nlmam2 == 1 ) then
                                  mv2 = m2 + isf2 - iso2
                                  if( abs(mv2) > l2 ) cycle
                                else
                                  mv2 = m2
                                endif
                                angl2 = conjg( gaunttd(l2,mv2,l0,m0,
     &                                               lg,mg2,Ylm_comp) )
                                if( abs( angl2 ) < eps10 ) cycle
                                Zero_term = .false.
                                exit boucle_isg2
                              end do
                            end do boucle_isg2
                            if( Zero_term ) cycle

                            lmp2 = 0
                            do lp2 = 0,lmax_t        

                              if( mod(lp2,2) /= imparite 
     &                                      .and. imparite /= 2) cycle
                              do mp2 = -lp2,lp2
                                if( nlmam2 == 1 .and.
     &                            ( lp2 /= l2 .or. mp2 /= m2 ) ) cycle 
                                  lmp2 = lmp2 + 1

                            do iso2 = 1,nspino  

                              if( nlmam2 == 1 ) then
                                lm2g = nlmam_u * ( iso2 - 1) + lm2
                              else
                                lm2g = nlmam_u * ( iso2 - 1) + lmp2
                              endif 
                                
                              if( Spinorbite .and. nlmam2 == 1 ) then
                                mv2 = m2 + isp2 - iso2
                                if( abs(mv2) > l2 ) cycle
                              else
                                mv2 = m2
                              endif
                              lmv2 = l2**2 + l2 + 1 + mv2

                              pwr1 = l0 + 1     ! psii = r*fct
                              f(1:nr) = r(1:nr)**pwr1
     &                             * psii(1:nr,is2) 
     &                             * zet(1:nr,lmv2,lmp2,iso2,isp2,iz2)

                              call ffintegr2_r(t1,f,r,nr,1,Rmtsd)

                              pwr2 = -l0        ! psii = r*fct       
                              f(1:nr) = r(1:nr)**pwr2
     &                           * psii(1:nr,is2) 
     &                           * zet(1:nr,lmv2,lmp2,iso2,isp2,iz2)

                              call ffintegr2_r(t2,f,r,nr,-1,Rmtsd)

                              f(1:nr) = psii(1:nr,is1)
     &                           * zet(1:nr,lmv1,lmp1,iso1,isp1,iz1)
     &                                  * ( r(1:nr)**pwr2 * t1(1:nr)
     &                                    + r(1:nr)**pwr1 * t2(1:nr) )

                              intrad_r = f_integr3(r,f,1,nr,Rmtsd) 
                      ! le 2 de 2 / (r-r')
                              fac = 8 * pi * intrad_r / ( 2 * l0 + 1 )

                              do isg1 = 1,2
                                if( .not. Dipmag .and. isg1 /= isf1 )
     &                                                        cycle
                                if( abs( coef_g(ist1,isg1) ) < eps10)
     &                                                        cycle 
                                mg1 = m_g(ist1,isg1)
                                angl1 = real( gaunttd(l1,mv1,l0,m0,lg,
     &                                                mg1,Ylm_comp), db)
                                if( abs( angl1 ) < eps10 ) cycle

                                do isg2 = 1,2
                                  if( .not. Dipmag .and. isg2 /= isf2 )
     &                                                        cycle
                                  if( abs( coef_g(ist2,isg2) ) < eps10 )
     &                                                        cycle 
                                  mg2 = m_g(ist2,isg2)
                                  angl2 = real( conjg( gaunttd(l2,mv2,
     &                                    l0,m0,lg,mg2,Ylm_comp) ), db )
                                  if( abs( angl2 ) < eps10 ) cycle

                                  if( Dipmag ) then
                                    isg12 = 2 * ( isg1 - 1 ) + isg2
                                  else
                                    isg12 = 1
                                  endif

! Atomic screening
                                  if( .not. atomic_scr .or. l2 == lg + 1 
     &                                              .or. l2 == lg + 1 ) 
     &                       Kern(lm1g,lm2g,isf1,isf2,isg12,ist1,ist2)        
     &                       = Kern(lm1g,lm2g,isf1,isf2,isg12,ist1,ist2)  
     &                            + coef_g(ist1,isg1)*coef_g(ist2,isg2)       
     &                            * angl1 * angl2 * fac 

                                end do ! fin boucle isg2
                              end do ! fin boucle isg1

                            end do ! fin boucle iso2
                            end do ! fin boucle lmp2
                            end do ! fin boucle lmp2
                          end do   ! fin boucle m2
                        end do     ! fin boucle l2

                      end do ! fin boucle iso1
                      end do ! fin boucle lmp1
                      end do ! fin boucle lmp1
                    end do   ! fin boucle m1
                  end do     ! fin boucle l1

                end do   ! fin boucle isf2
              end do     ! fin boucle isf1

            end do   ! fin boucle ist2
          end do     ! fin boucle ist1

        end do  ! fin boucle m0
      end do    ! fin boucle l0

! On rajoute la partie xc en LDA:
      if( .not. RPALF ) then

        do ist1 = 1,ninitl       

          if( ist1 > ninit1 ) then
            is1 = 2
          else
            is1 = 1
          end if

          do ist2 = 1,ninitl      
 
            if( Dyn_g .and. ist1 /= ist2 ) cycle  ! Noyau dynamique, fxc
                                            ! diagonal en etats initiaux
            if( ist2 > ninit1 ) then
              is2 = 2
            else
              is2 = 1
            end if

            if( Dyn_eg .and. is1/=is2 ) cycle   ! Noyau dynamique, fxc 
                                                ! diagonal en seuils
            do isf1 = 1,2         ! Spin du premier etat final

              isp1 = min(isf1,nspin)

              do isf2 = 1,2     ! Spin du deuxieme etat final

                isp2 = min( isf2, nspin )

                lm1 = 0
                do l1 = 0,lmax_t        
                  if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
                  do m1 = -l1,l1
                    lm1 = lm1 + 1

                    lmp1 = 0
                    do lp1 = 0,lmax_t        

                      if( mod(lp1,2) /= imparite 
     &                              .and. imparite /= 2) cycle
                      do mp1 = -lp1,lp1
                        if( nlmam2 == 1 .and.
     &                    ( lp1 /= l1 .or. mp1 /= m1 ) ) cycle 
                          lmp1 = lmp1 + 1

                    do iso1 = 1,nspino 
           
                      if( nlmam2 == 1 ) then
                        lm1g = nlmam_u * ( iso1 - 1) + lm1
                      else
                        lm1g = nlmam_u * ( iso1 - 1) + lmp1
                      endif 
                      if( Spinorbite .and. nlmam2 == 1 ) then
                        mv1 = m1 + isp1 - iso1
                        if( abs(mv1) > l1 ) cycle
                      else
                        mv1 = m1
                      endif
                      lmv1 = l1**2 + l1 + 1 + mv1

                      lm2 = 0
                      do l2 = 0,lmax_t        
                        if( mod(l2,2) /= imparite .and. imparite /= 2)
     &                                                             cycle
                        do m2 = -l2,l2
                          lm2 = lm2 + 1

                          lmp2 = 0
                          do lp2 = 0,lmax_t        

                            if( mod(lp2,2) /= imparite 
     &                                    .and. imparite /= 2) cycle
                            do mp2 = -lp2,lp2
                              if( nlmam2 == 1 .and.
     &                          ( lp2 /= l2 .or. mp2 /= m2 ) ) cycle 
                                lmp2 = lmp2 + 1

                          do iso2 = 1,nspino  

                            if( nlmam2 == 1 ) then
                              lm2g = nlmam_u * ( iso2 - 1) + lm2
                            else
                              lm2g = nlmam_u * ( iso2 - 1) + lmp2
                            endif 
                              
                            if( Spinorbite .and. nlmam2 == 1 ) then
                              mv2 = m2 + isp2 - iso2
                              if( abs(mv2) > l2 ) cycle
                            else
                              mv2 = m2
                            endif
                            lmv2 = l2**2 + l2 + 1 + mv2

! psii = r*fct
                            f(1:nr) = psii(1:nr,is1) * psii(1:nr,is2) 
     &                              * fxc(1:nr,isf1,isf2)
     &                              * zet(1:nr,lmv1,lmp1,iso1,isp1,iz1)
     &                              * zet(1:nr,lmv2,lmp2,iso2,isp2,iz2)

                            fac = f_integr3(r,f,1,nr,Rmtsd)

                            do isg1 = 1,2
                              if( .not. Dipmag .and. isg1 /= isf1 )
     &                                                      cycle
                              if( abs( coef_g(ist1,isg1) ) < eps10)
     &                                                      cycle 
                              mg1 = m_g(ist1,isg1)

                              do isg2 = 1,2
                                if( .not. Dipmag .and. isg2 /= isf2 )
     &                                                      cycle
                                if( abs( coef_g(ist2,isg2) ) < eps10 )
     &                                                      cycle 
                                mg2 = m_g(ist2,isg2)

                                if( Dipmag ) then
                                  isg12 = 2 * ( isg1 - 1 ) + isg2
                                else
                                  isg12 = 1
                                endif
! Les harmoniques sont supposees complexes
! Dans Gaunt4Y, ce sont les 2 premieres harmoniques qui sont 
! complexe-conjuguees.
                                angl = Gaunt4Y(lg,mg2,l1,mv1,l2,mv2,lg,
     &                                         mg1)
! Atomic screening
                                if( .not. atomic_scr .or. l2 == lg + 1 
     &                                               .or. l2 == lg + 1 ) 
     &                       Kern(lm1g,lm2g,isf1,isf2,isg12,ist1,ist2)        
     &                       = Kern(lm1g,lm2g,isf1,isf2,isg12,ist1,ist2)  
     &                          + coef_g(ist1,isg1) * coef_g(ist2,isg2)       
     &                          * angl * fac 

                              end do ! fin boucle isg2

                            end do ! fin boucle isg1

                          end do ! fin boucle iso2
                          end do ! fin boucle lmp2
                          end do ! fin boucle lmp2
                        end do   ! fin boucle m2
                      end do     ! fin boucle l2

                    end do ! fin boucle iso1
                    end do ! fin boucle lmp1
                    end do ! fin boucle lmp1
                  end do   ! fin boucle m1
                end do     ! fin boucle l1

              end do   ! fin boucle isf2
            end do     ! fin boucle isf1

          end do   ! fin boucle ist2
        end do     ! fin boucle ist1
      end if

      if( BSE ) then

        if( ( icheck > 1 .and. ie == 1 ) .or. icheck > 2 ) write(3,120)
        do is = 1,nbseuil
          do isf1 = 1,2
            isp1 = min(isf1,nspin)
            lm1 = 0
            do l1 = 0,lmax_t        
              if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
              do m1 = -l1,l1
                lm1 = lm1 + 1

                lmp1 = 0
                do lp1 = 0,lmax_t        

                  if( mod(lp1,2) /= imparite 
     &                          .and. imparite /= 2) cycle
                  do mp1 = -lp1,lp1
                    if( nlmam2 == 1 .and.
     &                ( lp1 /= l1 .or. mp1 /= m1 ) ) cycle 
                      lmp1 = lmp1 + 1

                do iso1 = 1,nspino

                  if( nlmam2 == 1 ) then
                    lm1g = nlmam_u * ( iso1 - 1) + lm1
                  else
                    lm1g = nlmam_u * ( iso1 - 1) + lmp1
                  endif 
                              
                  if( Spinorbite .and. nlmam2 == 1 ) then
                    mv1 = m1 + isp1 - iso1
                    if( abs(mv1) > l1 ) cycle
                  else
                    mv1 = m1
                  endif
                  lmv1 = l1**2 + l1 + 1 + mv1

                  if( .not. Core_resolved ) then
                    f(1:nr) = dv_ex_nex(1:nr) 
     &              * ( zet(1:nr,lmv1,lmp1,iso1,isp1,is) * r(1:nr) )**2
                    K_BSE = f_integr3(r,f,1,nr,Rmtsd)
                  endif

                  do ist1 = 1,ninitl
                    if( ( is == 1 .and. ist1 > ninit1 ) .or.
     &                  ( is == 2 .and. ist1 <= ninit1 ) ) cycle
                    if( Core_resolved ) then
                      f(1:nr) = dv_ex_nex(1:nr) 
     &                         * ( zet(1:nr,lmv1,lmp1,iso1,isp1,ist1)
     &                           * r(1:nr) )**2
                      K_BSE = f_integr3(r,f,1,nr,Rmtsd)
                    endif
                    do isg1 = 1,2
                      if( .not. Dipmag .and. isg1 /= isf1 ) cycle
                      if( Dipmag ) then
                        isg12 = 2 * ( isg1 - 1 ) + isg1
                      else
                        isg12 = 1
                      endif

                      if( ( icheck > 1 .and. ie == 1 ) .or. icheck > 2 )
     &                  Write(3,130) l1, m1, iso1, isf1, ist1,
     &                    coef_g(ist1,isg1)**2 * K_BSE
! Atomic screening
                        if( .not. atomic_scr .or. l2 == lg + 1 
     &                                               .or. l2 == lg + 1 ) 
     &                Kern(lm1g,lm1g,isf1,isf1,isg12,ist1,ist1)        
     &                  = Kern(lm1g,lm1g,isf1,isf1,isg12,ist1,ist1)  
     &                  + coef_g(ist1,isg1)**2 * K_BSE
                    end do 
                  end do 

                end do
                end do
                end do
              end do
            end do
          end do
        end do

      endif

      if( abs( Kern_fac - 1._db ) > eps10 )
     &       Kern(:,:,:,:,:,:,:) = Kern_fac * Kern(:,:,:,:,:,:,:)

      if( ( icheck > 1 .and. ie == 1 ) .or. icheck > 2 ) then
        write(3,140) ((ist1,ist2, ist2 = 1,ninitl), ist1 = 1,ninitl)
        lm1 = 0
        do l1 = 0,lmax_t        
          if( mod(l1,2) /= imparite .and. imparite /= 2) cycle
          do m1 = -l1,l1
            lm1 = lm1 + 1
            do iso1 = 1,nspino
              lm1g = nlmam_u * ( iso1 - 1 ) + lm1
              lm2 = 0
              do l2 = 0,lmax_t        
                if( mod(l2,2) /= imparite .and. imparite /= 2) cycle
                do m2 = -l2,l2
                  lm2 = lm2 + 1
                  do iso2 = 1,nspino
                    lm2g = nlmam_u * ( iso2 - 1 ) + lm2
                    do isf1 = 1,2
                      do isf2 = 1,2
                        do isg12 = 1,ns_dipmag
                          if( sum( abs( Kern(lm1g,lm2g,isf1,isf2,isg12, 
     &                                        :,:))) < eps10 ) cycle 
                          write(3,150) l1, m1, iso1, l2, m2, iso2, isf1,
     &                      isf2, isg12,  
     &                      ((Kern(lm1g,lm2g,isf1,isf2,isg12,ist1,ist2),
     &                         ist2 = 1,ninitl), ist1 = 1,ninitl)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      return
  100 format(/' ---- Kernel -------',100('-'))
  110 format(/' lmax_t =',i2,', nlmam_u =',i3)
  120 format(/'   l   m iso isp initl     K_BSE')
  130 format(5i4,1p,e16.5)
  140 format(/' Kern(g1,g2)',/'  l1  m1 so1  l2  m2 so2 is1 is2 g12  ',
     &         100(3x,'(',i1,',',i1,')',3x))
  150 format(9i4,1x,1p,100e11.3)
      end
 
!***********************************************************************

! On calcule Chi suivant la formule : Chi = Chi0 + Chi0 * K * Chi
!                                     Chi = (1 - Chi0 * K )^(-1) * Chi0

! Chi_0 en entree est en convention cristallo
! Chi en sortie est en convention physique (compl. conj. de cristallo)

      subroutine Cal_Chi(Chi,Chi_0,coef_g,Core_resolved,Dipmag,
     &            Energ,iaabsi,icheck,ie,imparite,Kern,lato,
     &            lmax_probe,lmax_t,
     &            mato,mix_repr,natome,nenerg,ninit1,ninitl,
     &            ninitlt,ngrph,nlm_probe,nlmam_u,nlmamax_u,nlmsa0,
     &            nlmsam,nomfich,ns_dipmag,nspin,nspino,nspino_t,
     &            Repres_comp,Spinorbite,tddft_mix,Tddft_so,Ylm_comp)

      use declarations
      implicit none

      integer,intent(in):: iaabsi, ie, imparite, lmax_probe, lmax_t, 
     &  nenerg, ninit1, ninitl, ninitlt, nlmam_u, nlmamax_u,
     &  nlm_probe, ns_dipmag, nspin, nspino, nspino_t, nlmsam, 
     &  ngrph, natome

      integer, dimension(2), intent(in):: mix_repr
      integer, dimension(natome,ngrph), intent(in):: nlmsa0
      character(len=132),intent(in):: nomfich

      complex(kind=db), dimension(nlm_probe*nspino_t,nlm_probe*nspino_t,
     &                   2,2,ns_dipmag,ninitl,ninitl), intent(out):: Chi
      complex(kind=db), dimension(nenerg,nlmamax_u*nspino,
     &                     nlmamax_u*nspino,nspin,nspin,ninitlt):: Chi_0
      integer, dimension(nlmsam,natome,ngrph):: lato, mato
      complex(kind=db), dimension(:), allocatable:: V
      complex(kind=db), dimension(:,:), allocatable:: A, B, Trans

      integer i, ifac, i1, i2, icheck, idim, is, isf1, isf2, isg1,
     &        isg12, isg2, iso, iso1, iso2, isp1, isp2, ist,
     &        ist1, ist2, j, l1, l2, lm, lm1, lm1c, lm1k, lm2,
     &        lm2c, lm2k, lmv1, lmv2, lm1d, lm2d, ls, m1, m2, m1d, m2d,  
     &        ndim, ndimg, nlmam

      logical Core_resolved, Dipmag, Spinorbite, Tddft_so, Ylm_comp, 
     &         tddft_mix, Stop_job
      logical, dimension(ngrph):: Repres_comp 

      real(kind=db),intent(in):: Energ
      real(kind=db),dimension(ninitl,2),intent(in):: coef_g
      real(kind=db), dimension(nlmam_u*nspino_t,nlmam_u*nspino_t,
     &              2,2,ns_dipmag,ninitl,ninitl), intent(in):: Kern
      real(kind=db), dimension(:,:), allocatable:: K

      if( icheck > 1 ) write(3,100)
      if( icheck > 1 ) write(3,110) Energ * rydb

      Chi(:,:,:,:,:,:,:) = ( 0._db, 0._db )
      Stop_job = .false.

 1000 continue

      ndim = 0
      ndimg = 0
      do ist1 = 1,ninitl
        do isg1 = 1,2    ! Spin etat initial
          if( abs( coef_g(ist1,isg1) ) < eps10 ) cycle
          ndimg = ndimg + 1
          do isf1 = 1,2    ! Spin etat final
            if( .not. Dipmag .and. isf1 /= isg1 ) cycle
            ndim = ndim + 1
          end do
        end do
      end do

! Pour l'instant, on ne considere pas la transition dipolaire magnetique
      ndim = ndim * nlmam_u * nspino_t

! Transformation, si harmoniques reelles

      if( .not. Ylm_comp ) then

        nlmam = ( lmax_t + 1 )**2
        allocate( Trans(nlmam,nlmam) )
        allocate( A(nlmam,nlmam) )
        A(:,:) = (0._db, 0._db)
        Trans(:,:) = (0._db, 0._db)
        Call Cal_Trans(nlmam,Trans)

        do is = 1,ninitlt
          do isp1 = 1,nspin
            do isp2 = 1,nspin

              lm1 = 0
              do l1 = 0,lmax_t        
                if( mod(l1,2) /= imparite .and. imparite /=2 ) cycle
                do m1 = -l1,l1
                  lm1 = lm1 + 1
                  lmv1 = l1**2 + l1 + 1 + m1
                  lm2 = 0
                  do l2 = 0,lmax_t        
                    if( mod(l2,2) /= imparite .and. imparite /=2 ) cycle
                    do m2 = -l2,l2
                      lm2 = lm2 + 1
                      lmv2 = l2**2 + l2 + 1 + m2
! On a forcement nspino = 1
                      A(lmv1,lmv2) = Chi_0(ie,lm1,lm2,isp1,isp2,is)
                    end do
                  end do
                end do
              end do

              if( icheck > 2 ) then
                write(3,'(/A)' ) ' Chi_0'
                do i = 1,9
                  write(3,120) A(i,:)
                end do
              endif

              A = matmul( A, Trans )
              A = matmul( conjg( transpose( Trans ) ), A )

              if( icheck > 2 ) then
                write(3,'(/A)' ) ' Chi_r'
                do i = 1,nlmam
                  write(3,120) A(i,:)
                end do
              endif

              lm1 = 0
              do l1 = 0,lmax_t        
                if( mod(l1,2) /= imparite .and. imparite /=2 ) cycle
                do m1 = -l1,l1
                  lm1 = lm1 + 1
                  lmv1 = l1**2 + l1 + 1 + m1
                  lm2 = 0
                  do l2 = 0,lmax_t        
                    if( mod(l2,2) /= imparite .and. imparite /=2 ) cycle
                    do m2 = -l2,l2
                      lm2 = lm2 + 1
                      lmv2 = l2**2 + l2 + 1 + m2
                      Chi_0(ie,lm1,lm2,isp1,isp2,is) = A(lmv1,lmv2)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do

        deallocate( Trans )
        deallocate( A )

      endif

      allocate( A(ndim,ndim), B(ndim,ndim), K(ndim,ndim) )

      A(:,:) = (0._db, 0._db)
      B(:,:) = (0._db, 0._db)
      K(:,:) = 0._db

      i1 = 0
      do ist1 = 1,ninitl ! Etats initiaux en entree

        if( Core_resolved ) then
          is = ist1    
        elseif( ist1 <= ninit1 ) then
          is = 1       ! indice du seuil qui correspond pour Chi_0
        else
          is = 2
        endif

        do isg1 = 1,2    ! Spin etat initial entree

          if( abs( coef_g(ist1,isg1) ) < eps10 ) cycle

          do isf1 = 1,2    ! Spin etat final entree

            if( .not. Dipmag .and. isf1 /= isg1 ) cycle
            isp1 = min( isf1, nspin)

            do lm1 = 1,nlmam_u
              do iso1 = 1,nspino_t
                if( Tddft_so ) then
                  lm1c = lm1
                else
                  lm1c = nlmamax_u * (iso1 - 1) + lm1
                endif
                lm1k = nlmam_u * (iso1 - 1) + lm1
                i1 = i1 + 1
                i2 = 0
                do ist2 = 1,ninitl

                  do isg2 = 1,2    ! Spin etat initial sortie
                    if( abs( coef_g(ist2,isg2) ) < eps10 ) cycle

                    do isf2 = 1,2    ! Spin etat final sortie

                      if( .not. Dipmag .and. isf2 /= isg2 ) cycle
                      isp2 = min( isf2, nspin)

                      if( Dipmag ) then
                        isg12 = 2 * (isg1 - 1 ) + isg2
                      else
                        isg12 = 1
                      endif
                      do lm2 = 1,nlmam_u
                        do iso2 = 1,nspino_t
                          if( Tddft_so ) then
                            lm2c = lm2
                          else
                            lm2c = nlmamax_u * (iso2 - 1) + lm2
                          endif
                          lm2k = nlmam_u * (iso2 - 1) + lm2
                          i2 = i2 + 1
! Conjugue car convention physique
                          if( ist1 == ist2 ) then
                            if(    ( Tddft_so .and. isg1 == isg2 .and. 
     &                            isf1 == iso1 .and. isf2 == iso2 )
     &                        .or. ( .not. Tddft_so .and.
     &                          ( isg1 == isg2 .or. Spinorbite ) ) )
     &                        A(i1,i2) = Conjg(
     &                             Chi_0(ie,lm1c,lm2c,isp1,isp2,is) ) 
                          endif 

                          K(i1,i2) = Kern(lm1k,lm2k,isf1,isf2,
     &                                      isg12,ist1,ist2)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      if( icheck > 1 ) then
        write(3,'(/A)') '  Chi_0(lm, iso, isg_ist)'
        write(3,130) (((lm, iso, ist, iso = 1,nspino_t),
     &           lm = 1,nlmam_u), ist = 1,ndimg)
        do lm = 1,ndim
          write(3,140) A(lm,:)
        end do
        write(3,'(/A)') '  Kern(lm, iso, isg_ist) '
        write(3,135) (((lm, iso, ist, iso = 1,nspino_t),
     &           lm = 1,nlmam_u), ist = 1,ndimg)
        do lm = 1,ndim
          write(3,140) K(lm,:)
        end do
      endif

      do i = 1,ndim
        do j = 1,ndim
          B(i,j) = - sum( A(i,:) * K(:,j) )
        end do
      end do

      deallocate( K )

      do i = 1,ndim
        B(i,i) = 1 + B(i,i) ! B = (1 - chi0*K)
      end do

      if( icheck > 2 ) then
        write(3,'(/A)') '  (1 - Chi_0.K)(lm,iso,isg_ist)'
        write(3,130) (((lm, iso, ist, iso = 1,nspino_t),
     &           lm = 1,nlmam_u), ist = 1,ndimg)
        do lm = 1,ndim
          write(3,140) B(lm,:)
        end do
      endif

      if( Stop_job ) stop

! B = (1 - chi0*K)**(-1)
      call invcomp(ndim,B,ndim,ndim,0,Stop_job)

      if( Stop_job ) then
        deallocate( A, B )
        icheck = 3
        goto 1000
      endif

      if( icheck > 2 ) then
        write(3,'(/A)') '  (1 - Chi_0.K)**(-1)(lm,iso,isg_ist)'
        write(3,130) (((lm, iso, ist, iso = 1,nspino_t),
     &           lm = 1,nlmam_u), ist = 1,ndimg)
        do lm = 1,ndim
          write(3,140) B(lm,:)
        end do
      endif

! Chi = (1 - chi0*K)**(-1) * Chi0
! La recopie sur V fait gagner de l'espace memoire.
      allocate( V(ndim) )

      do i = 1,ndim
        V(:) = B(i,:)
        do j = 1,ndim
          B(i,j) = sum( V(:) * A(:,j) )
        end do
      end do

      deallocate( V )

      if( icheck > 1 ) then
        write(3,'(/A)') '  Chi(lm, iso, isg_ist)'
        write(3,130) (((lm, iso, ist, iso = 1,nspino_t),
     &           lm = 1,nlmam_u), ist = 1,ndimg)
        do lm = 1,ndim
          write(3,140) B(lm,:)
        end do
      endif

      if( icheck > 1 ) write(3,150)
! Remplissage de Chi pour les lm utiles (jusqu'a nlm_probe).
      i1 = 0
      do ist1 = 1,ninitl
        do isg1 = 1,2    ! Spin etat initial entree
          if( abs( coef_g(ist1,isg1) ) < eps10 ) cycle
          do isf1 = 1,2    ! Spin etat final entree
            if( .not. Dipmag .and. isf1 /= isg1 ) cycle

            do l1 = 0,lmax_t        
              if( mod(l1,2) /= imparite .and. imparite /=2 ) cycle
              do m1 = -l1,l1
                lmv1 = l1**2 + l1 + 1 + m1
                do iso1 = 1,nspino_t
                  lm1c = nlm_probe * (iso1 - 1) + lmv1
                  i1 = i1 + 1

                  i2 = 0
                  do ist2 = 1,ninitl
                    do isg2 = 1,2    ! Spin etat initial entree
                      if( abs( coef_g(ist2,isg2) ) < eps10 ) cycle
                      if( Dipmag ) then
                        isg12 = 2 * (isg1 - 1 ) + isg2
                      else
                        isg12 = 1
                      endif
                      do isf2 = 1,2    ! Spin etat final entree
                        if( .not. Dipmag .and. isf2 /= isg2 ) cycle

                        do l2 = 0,lmax_t        
                          if(mod(l2,2) /= imparite .and. imparite /=2)
     &                                                           cycle
                          do m2 = -l2,l2
                            lmv2 = l2**2 + l2 + 1 + m2
                            do iso2 = 1,nspino_t
                              lm2c = nlm_probe * (iso2 - 1) + lmv2
                              i2 = i2 + 1
                              if( l1 > lmax_probe
     &                            .or. l2 > lmax_probe ) cycle
! Oana: analyse du melange des representations
                              if( tddft_mix) then
                                do lm1 = 1, nlmsa0(iaabsi,mix_repr(1))
                                  do lm2 = 1, nlmsa0(iaabsi,mix_repr(2))
                                    if( tddft_so ) then
                                      if( l1 == lato(lm1,
     &                                        iaabsi,mix_repr(1)).and. 
     &                                    l2 == lato(lm2,
     &                                        iaabsi,mix_repr(2)).and.
     &                                    m1 - isf1 + iso1  == mato(lm1,
     &                                          iaabsi,mix_repr(1)).and.
     &                                    m2 - isf2 + iso2  == mato(lm2,
     &                                        iaabsi,mix_repr(2))) then
  ! all representations are calculated 
!!! bug car igrph pas defini
                                        if( .not. Repres_comp(1) .or.         
     &                                                  spinorbite) then 
                                          Chi(lm1c,lm2c,isf1,isf2,isg12,
     &                                             ist1,ist2) = B(i1,i2)
                                        else                                
                                          Chi(lm1c,lm2c,isf1,isf2,isg12,
     &                                             ist1,ist2) = B(i1,i2)
! the m < 0 are not contained in mato;
!  tau(l1 m1;l2 m2) = tau(l1 -m1;l2 -m2)
! for instance, starting from (2 2) I produce (-2 -2); (-2 2) and (2 -2)
                                          m1d = - m1 + isf1 - iso1 ! -m1
                                          m2d = - m2 + isf2 - iso2 ! -m2
                                          lm1d =  l1**2 + l1 + 1 + m1d                                   
                                          lm2d =  l2**2 + l2 + 1 + m2d 
                                          ifac = (-1)**(m1d+m2d)
                                          Chi(lm1d,lm2d,isf1,isf2,isg12,
     &                                      ist1,ist2) = B(i1,i2) * ifac
                                          ifac = (-1)**(m1d+m2)
                                          Chi(lm1d,lm2c,isf1,isf2,isg12,
     &                                      ist1,ist2) = B(i1,i2) * ifac
                                          ifac = (-1)**(m1+m2d)
                                          Chi(lm1c,lm2d,isf1,isf2,isg12,
     &                                      ist1,ist2) = B(i1,i2) * ifac
                                        end if
                                      end if
                                    else   ! scalar 
                                      if( l1 == lato(lm1,
     &                                          iaabsi,mix_repr(1)).and. 
     &                                    l2 == lato(lm2,
     &                                          iaabsi,mix_repr(2)).and.
     &                                    m1 == mato(lm1,
     &                                          iaabsi,mix_repr(1)).and.
     &                                    m2 == mato(lm2,
     &                                          iaabsi,mix_repr(2)))then
                                        if( .not. Repres_comp(1) .or.          
     &                                                  spinorbite) then
                                          Chi(lm1c,lm2c,isf1,isf2,isg12,
     &                                             ist1,ist2) = B(i1,i2)
                                        else
                                          lm1d =  l1**2 + l1 + 1 - m1                                   
                                          lm2d =  l2**2 + l2 + 1 - m2
                                          ifac = (-1)**(-m1-m2)
                                          Chi(lm1d,lm2d,isf1,isf2,isg12,
     &                                      ist1,ist2) = B(i1,i2) * ifac
                                          ifac = (-1)**(-m1+m2)
                                          Chi(lm1d,lm2c,isf1,isf2,isg12,
     &                                      ist1,ist2) = B(i1,i2) * ifac
                                          ifac = (-1)**(m1-m2)
                                          Chi(lm1c,lm2d,isf1,isf2,isg12,
     &                                      ist1,ist2) = B(i1,i2) * ifac
                                        end if
                                      end if
                                    end if
                                  end do
                                end do
                              else
                                Chi(lm1c,lm2c,isf1,isf2,isg12,ist1,ist2)
     &                                                      = B(i1,i2)
                              end if
                              if( icheck > 1
     &                           .and. abs( B(i1,i2) ) > eps10 )
     &                          write(3,160) l1, m1, ist1, isg1, l2,   
     &                           m2, ist2, isg2,lm1c, lm2c, B(i1,i2)
                            end do
                          end do
                        end do

                      end do
                    end do
                  end do

                end do
              end do
            end do

          end do
        end do
      end do
      
      deallocate( A, B )

      if( icheck > 1 ) call write_Chi(Chi,Energ,ninitl,icheck,ie,
     &                      nlm_probe,nomfich,ns_dipmag,nspino_t,'Chi')

      return
  100 format(/' ---- Cal_chi -------',100('-'))
  110 format(/' Energ =',f10.3,' eV')
  120 format(1p,9(1x,2e10.2))
  130 format(200(7x,'(',2(i2,','),i2,')',5x))
  135 format(200(1x,'(',2(i2,','),i2,')'))
  140 format(1p,400e11.3)
  150 format(' l1 m1 s1 g1 l2 m2 s2 g2  lm1c lm2c             chi')
  160 format(8i3,2i5,1p,2e13.5)

      end

!**********************************************************************

      subroutine write_Chi(Chi,Energ,ninitl,icheck,ie,nlm_probe,
     &                     nomfich,ns_dipmag,nspino,mot3)

      use declarations
      implicit none

      character(len=3),intent(in)::  mot3
      character(len=132),intent(in)::  nomfich

      integer,intent(in):: ninitl, icheck, ie, nlm_probe, ns_dipmag,
     &                     nspino

      complex(kind=db),dimension(nlm_probe*nspino,nlm_probe*nspino,2,2,
     &                        ns_dipmag,ninitl,ninitl),intent(in):: Chi
 
      logical diag

      real(kind=db),intent(in):: Energ

      character(len=132)::  chi_conv

      integer l, isp1, isp2, ist1, ist2, lm1, lm2

      if( icheck > 2 ) then
        diag = .false.
      else
        diag = .true.
      endif

! Nom du fichier de sortie
      chi_conv = nomfich
      l = len_trim( chi_conv )
      chi_conv(l+1:l+1) = '_'
      chi_conv(l+2:l+4) = mot3(1:3)
      chi_conv(l+5:l+8) = '.txt'

      if( ie == 1 ) then 
        open(31, file = chi_conv)
        if( diag ) then
          write(31,110) ((( mod(lm1,10), isp1, mod(ist1,10) ,
     &                    lm1 = 1,nspino*nlm_probe),             
     &                    isp1 = 1,2), ist1 = 1,ninitl)
        else
          write(31,120) (((((( mod(lm1,10), mod(lm2,10), isp1, isp2,
     &      ist1, ist2, 
     &      lm2 = 1,nspino*nlm_probe), lm1 = 1,nspino*nlm_probe),             
     &      isp2 = 1,2), isp1 = 1,2), ist2 = 1,ninitl),ist1 = 1,ninitl)
        endif
      else
        open(31, file = chi_conv, position='append')
      endif

      if( diag ) then
        write(31,130) Energ*rydb,
     &              ((( chi(lm1,lm1,isp1,isp1,1,ist1,ist1),
     &       lm1 = 1,nspino*nlm_probe),             
     &       isp1 = 1,2), ist1 = 1,ninitl)
      else
        write(31,130) Energ*rydb,
     &              (((((( chi(lm1,lm2,isp1,isp2,1,ist1,ist2),
     &      lm2 = 1,nspino*nlm_probe),lm1 = 1,nspino*nlm_probe),             
     &      isp2 = 1,2), isp1 = 1,2), ist2 = 1,ninitl), ist1 = 1,ninitl)
      endif

      Close(31)

      return
  110 format('  Energy  ',
     &        5000('  (lm=',i1,',isp=',i1,',i=',i1,')  Im'))
  120 format('  Energy  ',5000('  (',i1,',',i1,',',i1,',',i1,',',i1,
     &         ',',i1,')    Im '))
  130 format(f10.3,1p,11664e11.3)
      end

!***********************************************************************

! Fonction qui calcule le coefficient de Gaunt: 
! Int( Y(l1,m1)*Y(l2,m2)Y(l3,m3) dOmega  )
! ou Y(l2,m2) et Y(l3,m3) sont complexes et  Y(l1,m1) est soit reelle
! soit complexe

      function gaunttd(l1,m1,l2,m2,l3,m3,Ylm_comp)

      use declarations
      implicit none

      complex(kind=db):: gaunttd 

      integer, intent(in):: l1, m1, l2, m2, l3, m3
      logical, intent(in):: Ylm_comp

      real(kind=db) gr, gi, gauntcp

      if( Ylm_comp .or. m1 == 0 ) then
        gaunttd = cmplx( gauntcp(l1,m1,l2,m2,l3,m3), 0._db, db )
      else if( m1 > 0 ) then
! Gauntcp calculant Gaunt pour le complexe conjugue, on appele avec
! (-1)**m Y(l1,-m1) = Y(l1,m1)*
        gr = ( (-1)**m1 * gauntcp(l1,-m1,l2,m2,l3,m3) 
     &                  + gauntcp(l1,m1,l2,m2,l3,m3) ) / sqrt(2._db)
        gi = 0._db
        gaunttd = cmplx(gr, gi, db)
      else
! Gauntcp calculant Gaunt pour le complexe conjugue, on appele avec
! (-1)**m Y(l1,-m1) = Y(l1,m1)*
        gr = 0._db
        gi = ( gauntcp(l1,-m1,l2,m2,l3,m3)  
     &         - (-1)**m1 * gauntcp(l1,m1,l2,m2,l3,m3) ) / sqrt(2._db)
        gaunttd = cmplx(gr, gi, db)
      end if 

      return
      end

!***********************************************************************

! Transformation harmo comp vers Harmo reel
! La transformation inverse est le conjugue de la transpose

      subroutine Cal_Trans(nlmam_u,Trans)

      use declarations
      implicit none

      integer:: is, l1, l2, lm1, lm2, m1, m2, nlmam_u
 
      complex(kind=db):: r2_r, r2_i
      complex(kind=db),dimension(nlmam_u,nlmam_u):: Trans

      real(kind=db):: r2

      Trans(:,:) = (0._db, 0._db)

      r2 = 1 / sqrt(2._db) 
      r2_r = cmplx( r2,    0._db, db)
      r2_i = cmplx( 0._db, r2,    db)

      lm1 = 0
      boucle_l1: do l1 = 0,100
        do m1 = -l1,l1
          lm1 = lm1 + 1
          if( lm1 > nlmam_u ) exit boucle_l1
          is = (-1)**m1 

          lm2 = 0
          boucle_l2: do l2 = 0,100
            do m2 = -l2,l2
              lm2 = lm2 + 1
              if( l1 /= l2 ) cycle
              if( lm2 > nlmam_u ) exit boucle_l2
                     
              if( m1 == m2 ) then

                if( m1 == 0 ) then
                  Trans(lm1,lm2) = (1._db,0._db)
                elseif( m1 > 0 ) then
                  Trans(lm1,lm2) = r2_r
                else
                  Trans(lm1,lm2) = is * r2_i
                endif

              elseif( m1 == - m2 ) then

                if( m1 > 0 ) then
                  Trans(lm1,lm2) = is * r2_r
                else
                  Trans(lm1,lm2) = - r2_i
                endif

              endif

            end do
          end do boucle_l2

        end do
      end do boucle_l1
   
      return
      end

!***********************************************************************

      subroutine fxcorr(alfpot,fxc,icheck,magnetic,nr,nspin,
     &                  r,rhoato_abs,rsato_abs)

      use declarations
      implicit none

      integer, intent(in):: icheck, nr, nspin

      logical, intent(in):: magnetic

      real(kind=db), intent(in):: alfpot
      real(kind=db),dimension(nr),intent(in):: r, rsato_abs
      real(kind=db),dimension(nr,nspin),intent(in):: rhoato_abs
      real(kind=db),dimension(nr,2,2),intent(out):: fxc

      integer ir, isp, isp1, isp2

      real(kind=db):: f_vonbarth, fprime_vonbarth
      real(kind=db):: a, b, c, c_p, c_f, d, e, f, f1, f2, f3, f4, f5,  
     &                fac, r_p, r_f, rsa, qtr, tr, x, xx

      if( icheck > 2 ) write(3,100)

      fxc(:,:,:) = 0._db

      if( alfpot > eps4 ) then

! Xalpha potential
        fac = - 2 * pi * alfpot / 3
        do isp = 1,2
          fxc(:,isp,isp) = fac * rsato_abs(:)**2
        end do

      else if( alfpot < eps4 ) then

! Pour les valeurs de c_p et r_p on garde les valeurs non polarise
! choisie aussi par Moruzzi Janak et William (1978). Pour r_f et c_f
! on prend aussi leurs parametres plutot que les originaux de Von Barth
! qui sont : c_p = 0.0504, r_p = 30., c_f = 0.0254, r_f = 75.
        c_p =  0.045_db
        r_p = 21._db
        c_f = 0.0254_db
        r_f = 75._db

        tr  =  1._db / 3._db
        qtr = 4._db / 3._db

        f1 = ( 36._db / pi**2 )**tr
        f2 = 4 * pi / 9 
        f3 = (2._db)**(-tr)
        f4 = qtr / ( 1 - f3 )  
        f5 = 1 / ( 1 - f3 )  

        do isp1 = 1,2
          do isp2 = 1,2
            do ir = 1,nr
              if( Magnetic ) then
                x = rhoato_abs(ir,isp2) / sum( rhoato_abs(ir,1:nspin) )
              else
                x = 0.5_db
              end if
              rsa = rsato_abs(ir)

              if( isp1 == isp2 ) then
                xx = 1 - x
                a = f1 * ( x**(-2*tr) * xx + x**tr ) / rsa
              else
                xx = - x
                a = 0._db
              endif                 

              b = c_p * r_p / ( rsa + r_p)          

              c = f4 * rsa * ( x**qtr + (1-x)**qtr - x**tr)
     &          * ( fprime_vonbarth(rsa/r_f) * (c_f/r_f)
     &            - fprime_vonbarth(rsa/r_p) * (c_p/r_p) )

              d = f5 * ( x**qtr + (1-x)**qtr - f3 ) 
     &          * ( c_f*r_f / ( r_f + rsa ) -  c_p*r_p / ( r_p + rsa ) )

              e = f4 * ( x**(-2*tr) + 4*x**tr - 4*(1-x)**tr ) 
     &               * xx * ( c_p * f_vonbarth(rsa/r_p) 
     &                      - c_f * f_vonbarth(rsa/r_f) )

              f = 4 * f5 * xx * ( (1-x)**tr - x**tr )
     &          * ( c_p * log(1 + r_p/rsa) - c_f * log(1 + r_f/rsa) )

              fxc(ir,isp1,isp2) = - ( a + b + c + d + e + f )
     &                          * f2 * rsa**3

            end do
          end do
        end do
      end if

      if( icheck > 2 ) then
        if( magnetic ) then
          write(3,110)
        else
          write(3,120)
        end if
        do ir = 1,nr
          write(3,130) r(ir)*bohr, rsato_abs(ir)*bohr, 
     &            ( rhoato_abs(ir,isp)*r(ir)**2, isp = 1,nspin ),
     &            fxc(ir,:,:)
        end do
      end if

      return
  100 format(/' ---- Fxcorr -------',100('-'))
  110 format(/'     Radius       Rsato   rhoato_up*r**2 ',
     &  'rhoato_dn*r**2  fxc_uu      fxc_ud       fxc_du       fxc_dd')
  120 format(/'     Radius       Rsato     rhoato*r**2     fxc')
  130 format(1p,8e13.5)
      end

!***********************************************************************

! Calcul de Integrale( Y(l1,m1)* Y(l2,m2)* Y(l3,m3) Y(l4,m4) dOmega )

! Utilise : Y(L2)* Y(L3) = Somme_L  Gaunt(L3,L2,L) Y(L)
! avec L = (l,m)

! G4(L1,L2,L3,L4) = Somme_L  G3(L2,L3,L) x G3(L4,L1,L)

      function Gaunt4Y(l1,m1,l2,m2,l3,m3,l4,m4)

      use declarations
      implicit none

      integer, intent(in):: l1, m1, l2, m2, l3, m3, l4, m4
      integer:: l, lmin, lmax, m

      real(kind=db):: Gaunt4Y, Gauntcp

      Gaunt4Y = 0._db

      lmin = max( abs(l2-l3), abs(l1-l4) )
      lmax = min( l2+l3, l1+l4 )
 
      do l = lmin,lmax,2
        do m = -l,l
          Gaunt4Y = Gaunt4Y + Gauntcp(l2,m2,l3,m3,l,m)
     &                      * Gauntcp(l4,m4,l1,m1,l,m)
        end do
      end do

      return 
      end

!***********************************************************************

      function fprime_vonbarth(x)

      use declarations
      implicit none
      real(kind=db):: fprime_vonbarth, x

      fprime_vonbarth = 3*x**2 * log(1 + 1/x) - 1/x -3*x + 1.5_db 

      return
      end
