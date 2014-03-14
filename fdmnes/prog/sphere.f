! FDMNES subroutines
! Resolution de l'equation de Schrodinger dans la partie spherique
! des atomes. Calcul de ces fonctions sur les points en bordures.

!***********************************************************************

      subroutine Sphere(Axe_Atom_grn,Base_ortho,
     &          dcosxyz,Ecinetic,Eimag,Energ,Enervide,Full_atom,
     &          Full_potential,Green,Hubb_a,Hubb_d,iaabsi,iapr,iaprotoi,
     &          ibord,icheck,igreq,igroupi,iopsymr,lmax,m_hubb,n_atom_0,
     &          n_atom_ind,n_atom_proto,natome,nbord,nbtm,neqm,ngroup_m,
     &          nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspin,nspino,
     &          numat,phiato,posi,r,Relativiste,Rmtg,
     &          Rmtsd,Spinorbite,Tau_ato,V_hubb,
     &          V_intmax,V0bd,Vrato,xyz,Ylm_comp,Ylmato)
 
      use declarations
      implicit none

      integer:: ia, iaabsi, iang, iapr, ib, icheck, iga, igr, kgr, l,
     &  l_hubbard, lfin, lm, lm0, lmax, lmp, m, m_hubb, mp, n, n_atom_0, 
     &  n_atom_ind, n_atom_proto, natome, nbtm, neqm, ngroup_m, nlm1,  
     &  nlm2, nlmagm, nlmmax, np, nphiato1, nphiato7, npsom, nr, nrmtg,
     &  nrmtsd, nspin, nspino, numat
      integer, dimension(natome):: iaprotoi, igroupi, nbord
      integer, dimension(nopsm):: iopsymr
      integer, dimension(nbtm,natome):: ibord
      integer, dimension(0:n_atom_proto,neqm):: igreq 
      
      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
      complex(kind=db), dimension(nlmagm,nspin,nlmagm,nspin,
     &                                n_atom_0:n_atom_ind):: Tau_ato
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                       V_hubb
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: u

      logical:: Ylm_comp, Base_ortho, Ecomp, Full_atom, Full_potential,
     &  Green, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm,  
     &  Spinorbite

      real(kind=db):: cosang, Eimag, Energ, Enervide, Rmtg, Rmtsd, 
     &                V_intmax
      real(kind=db), dimension(3):: dcosxyz
      real(kind=db), dimension(nspin):: Ecinetic, V0bd
      real(kind=db), dimension(nr,nspin):: Vrato 
      real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn
      real(kind=db), dimension(nbtm,nlmmax,natome):: Ylmato

      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nr,nspin):: g0, gm, gp, V
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(3,natome):: posi
      real(kind=db), dimension(4,npsom):: xyz
      real(kind=db), dimension(nphiato1,nlmagm,nspin,nspino,natome,
     &                                               nphiato7):: phiato

      real(kind=db), dimension(:,:,:,:,:), allocatable:: ur

      konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

      if( icheck > 1 ) write(3,110) iapr, numat, lmax
      if( icheck > 1 ) write(3,120) Energ*rydb, Enervide*rydb
      if( icheck > 2 ) then
        write(3,130) Ecinetic(:)*rydb
        write(3,140) V0bd(:)*rydb
        write(3,150) konde(:)
      endif

      call mod_V(icheck,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,
     &                 V_intmax,V0bd,Vrato)

      call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nr,nspin,nspino,
     &                        numat,r,Relativiste,Spinorbite,V)

      gp(:,:) = 1 / gp(:,:)

      if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10
     &     .or. Ecinetic(nspin) < eps10 ) then
        Ecomp = .true.
      else
        Ecomp = .false.
      endif
      Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

      if( Full_potential ) then
        lfin = 0
      else
        lfin = lmax
      endif

      Renorm = .true.

      do l = 0,lfin

        if( Hubb_a .and. l == l_hubbard( numat ) )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif  

        if( Full_potential ) then
          nlm1 = ( lmax + 1 )**2
          nlm2 = nlm1
        elseif( Hubb_m .and. .not. Hubb_d ) then
          nlm1 = 2*l + 1
          nlm2 = nlm1
        elseif( Spinorbite .or. Hubb_m ) then
          nlm1 = 2*l + 1
          nlm2 = 1
        else
          nlm1 = 1
          nlm2 = 1
        endif

        allocate( u(nr,nlm1,nlm2,nspin,nspino) )
        allocate( ur(nr,nlm1,nlm2,nspin,nspino) )
        allocate( Tau(nlm1,nspin,nlm1,nspin) )

        call Sch_radial(Ecinetic,Ecomp,Eimag,f2,
     &           Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,
     &           l,lmax,m_hubb,nlm1,nlm2,nr,nrmtg,nspin,nspino,numat,r,
     &           Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,u,
     &           ur,V_hubb)

! Recopie
        if( Full_potential ) then
          do lm = 1,nlm1
            do lmp = 1,nlm2
              Tau_ato(lm,:,lmp,:,iapr) = Tau(lm,:,lmp,:)
            end do
          end do
        else
          lm0 = l**2 + l + 1
           do m = -l,l
            if( nlm1 == 1 ) then
              n = 1
            else
              n = l + 1 + m
            endif
            do mp = -l,l
              if( nlm1 == 1 ) then
                if( m /= mp ) cycle
                np = 1
              else
                np = l + 1 + mp
              endif
              Tau_ato(lm0+m,:,lm0+mp,:,iapr) = Tau(n,:,np,:)
            end do
          end do
        endif

        deallocate( Tau )

! Calcul les fonctions radiales phiato.
        if( .not. Green ) then
          do ib = 1,natome
            if( Full_atom ) then
              if( ib > 1 ) exit
              ia = iapr
            else
              ia = ib
              if( iaprotoi(ia) /= iapr ) cycle
            endif

            if( nspin == 2 ) then
              igr = igroupi(ia)
              iga = igroupi(iaabsi)
              kgr = igreq(iaprotoi(ia),1) 
          
              cosang = sum( Axe_Atom_grn(:,kgr) * Axe_Atom_grn(:,igr) )
              cosang = cosang
     &         * abs( sum( Axe_Atom_grn(:,iga) * Axe_Atom_grn(:,igr) ) )
              if( abs(cosang - 1) < eps4 ) then
                iang = 1
              elseif( abs(cosang + 1) < eps4 ) then
                iang = - 1
              else
                iang = 0
              endif
            else
              iang = 1
            endif

            call cal_phiato(Base_ortho,dcosxyz,
     &           Full_potential,ia,iang,ibord,
     &           icheck,iopsymr,l,lmax,nlm1,nlm2,natome,
     &           nbord(ia),nbtm,nlmagm,nlmmax,
     &           nphiato1,nphiato7,npsom,nr,nspin,nspino,phiato,posi,r,
     &           Radial_comp,Spinorbite,u,ur,xyz,Ylm_comp,Ylmato)
          end do
        endif

        deallocate( u, ur )

      end do   ! fin de la boucle sur l

      return
  110 format(/' iapr =',i3,', Z =',i3,', lmax =',i2)
  120 format(' Energ =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd  =',2f10.3)
  150 format(' konde =',2f12.5)
      end

!***********************************************************************

      Subroutine mod_V(icheck,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,
     &                 V_intmax,V0bd,Vrato)

      use declarations
      implicit none

      integer:: icheck, ir, ispin, nr, nrmtg, nrmtsd, nspin

      real(kind=db):: Rmtg, Rmtsd, V_intmax
      real(kind=db), dimension(nr):: r
      real(kind=db), dimension(nspin):: V0bd
      real(kind=db), dimension(nr,nspin):: V, Vrato

      do ispin = 1,nspin
        V(1:nr,ispin) = Vrato(1:nr,ispin)
      end do

      if( V_intmax < 1000._db ) then
        do ispin = 1,nspin
          do ir = 1,nr
            V(ir,ispin) = min( V(ir,ispin), V_intmax )
          end do
        end do 
      endif

      do ir = 1,nr-1
        if( r(ir) > Rmtg + eps10 ) exit
      end do
      nrmtg = ir
      do ir = 1,nr-1
        if( r(ir) > Rmtsd + eps10 ) exit
      end do
      nrmtsd = ir

! n'a pas d'importance en FDM
      if( abs( V0bd(1) ) > eps10 ) then
        do ir = 1,nr-1
          if( r(ir) > Rmtg - eps10 ) exit
        end do
        do ispin = 1,nspin
          V(ir:nr,ispin) = V0bd(ispin)
        end do
      endif

      if( icheck > 2 ) then
        if( nspin == 2 ) then
          write(3,'(/A)') '    Radius         V(up)         V(dn)'
        else
          write(3,'(/A)') '    Radius           V'
        endif
        do ir = 1,nr
          write(3,110) r(ir)*bohr, V(ir,:)*rydb
        end do
      endif

      return
  110 format(1p,3e14.6)
      end

!***********************************************************************

      subroutine coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nr,nspin,nspino,
     &                        numat,r,Relativiste,Spinorbite,V)

      use declarations
      implicit none

      integer:: ir, ispin, nr, nspin, nspino, numat

      logical:: Relativiste, Spinorbite

      real(kind=db):: a2s4, bder, dr, dvr, Enervide, fac, r0, rm, rp,
     &                Vme

      real(kind=db), dimension(nr):: cgrad0, cgradm, cgradp, clapl0,
     &                               claplm, claplp, f2, r 
      real(kind=db), dimension(nr,nspin):: g0, gm, gp, V 
      real(kind=db), dimension(nr,nspino):: gso 

      gso(:,:) = 0._db

      do ir = 1,nr-1
        rp = r(ir+1)
        if( ir == 1 ) then
          rm = r(1)**2 / r(2)
        else
          rm = r(ir-1)
        endif
        r0 = r(ir)
        dr = 0.5 * ( rp - rm )
        claplm(ir) = 1 / ( ( r0 - rm ) * dr )
        claplp(ir) = 1 / ( ( rp - r0 ) * dr )
        clapl0(ir) = - claplm(ir) - claplp(ir)
        if( Spinorbite .or. Relativiste ) then
          cgradm(ir) = ( rp - r0 ) / ( ( rm - r0 ) * ( rp - rm ) )
          cgradp(ir) = ( rm - r0 ) / ( ( rp - r0 ) * ( rm - rp ) )
          cgrad0(ir) = - cgradm(ir) - cgradp(ir)
        endif
        f2(ir) = 1 / r(ir)**2
      end do

! alfa_sf = constante de structure fine.
      a2s4 = 0.25_db * alfa_sf**2

      do ir = 1,nr-1
        do ispin = 1,nspin

          Vme = V(ir,ispin) - Enervide
          g0(ir,ispin) = - clapl0(ir) + Vme
          gm(ir,ispin) = - claplm(ir)
          gp(ir,ispin) = - claplp(ir)
          if( ir == nr ) cycle

          if( Relativiste .or. Spinorbite ) then

            bder = 1 / ( 1 - a2s4 * Vme )
            if( ir == 1 ) then
              dvr = 2 * numat / r(1)**2
            else
              dvr = cgradm(ir) * V(ir-1,ispin)
     &            + cgrad0(ir) * V(ir,ispin)
     &            + cgradp(ir) * V(ir+1,ispin)
            endif
            fac = a2s4 * bder * dvr

            if( Relativiste ) then
              g0(ir,ispin) = g0(ir,ispin) - a2s4 * Vme**2
     &                     - fac * ( cgrad0(ir) - 1 / r(ir) )
              gm(ir,ispin) = gm(ir,ispin) - fac * cgradm(ir)
              gp(ir,ispin) = gp(ir,ispin) - fac * cgradp(ir)
            endif
            if( Spinorbite ) gso(ir,ispin) = fac / r(ir)

          endif

        end do
      end do

      return
      end

!***********************************************************************

! Resolution de l'equation de Schrodinger radiale dans le cas non
! spherique.
! Applele par Sphere, radial, radial_sd
! Dans ur(ir,m,mp,isp,isol) : (mp, isol) definissent les etats de base
!                             (m, isp) sont le moment et le spin reels

      subroutine Sch_radial(Ecinetic,Ecomp,Eimag,f2,
     &           Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,
     &           ll,lmax,m_hubb,nlm1,nlm2,nr,nrmtg,nspin,nspino,numat,r,
     &           Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,u,
     &           ur,V_hubb)

      use declarations
      implicit none

      integer:: icheck, im, ip, ir, isol, isp, isq, l, l_hubbard, l2,
     &  li, ll, lf, lmax, m, m_hubb, nlm1, nlm2, mp, ms, n, np, nr,
     &  nrmtg, nsol, nspin, nspino, numat

      logical:: Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Hubb_nd, 
     &  Radial_comp, Renorm, Relativiste, Spino_simple, Spinorbite

      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(nlm1,nspin,nlm1,nspin):: Tau
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                        V_hubb
      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u 

      real(kind=db):: br, ci, Eimag, fac, p, Rmtg, td 
      real(kind=db), dimension(nspin):: cr, Ecinetic 
      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nr,nspin):: g0, gm, gp 
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur 
      real(kind=db), dimension(:,:), allocatable:: fac_m 

      if( Radial_comp ) then
        u(:,:,:,:,:) = ( 0._db, 0._db )
      else
        ur(:,:,:,:,:) = 0._db
      endif

      Spino_simple = Spinorbite .and. nlm2 == 1

      if( Full_potential ) then
        li = 0
        lf = lmax
      else 
        li = ll
        lf = ll
      endif

      do l = li,lf
 
        l2 = l * ( l + 1 )

        if( Hubb_a .and. l == l_hubbard( numat ) )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif  
        Hubb_nd = Hubb_m .and. .not. Hubb_d

        br = - numat / ( l + 1._db )
        cr(:) = - ( 2 * numat * br + Ecinetic(:) ) / ( 4 * l + 6 )  
        ci = - Eimag / ( 4 * l + 6 )  

!Je fais l'hypothese que le terme de Hubbard ne croise pas les solutions
! 1 et 2 en cas de spinorbite.

! Developpement a l'origine
        do isol = 1,nspino
          do isp = 1,nspin
            do m = -l,l

              if( Full_potential ) then
                n = l**2 + l + 1 + m
              elseif( Hubb_m .or. Spinorbite ) then
                n = l + 1 + m
              else
                n = 1
              endif

              if( (m == l .and. isp == 1)
     &           .or. (m == -l .and. isp == 2 ) ) then
                nsol = 1
              else
                nsol = 2
              endif
              if( Spinorbite .and. nsol == 1 .and. isp /= isol ) cycle

              if( Spinorbite .and. Relativiste ) then
                if( nsol == 1  .or. isol == 1 ) then          
                  p = sqrt( ( l + 1 )**2 - ( alfa_sf * numat )**2 )
                else
                  p = sqrt( l**2 - ( alfa_sf * numat )**2 )
                endif
              elseif( Spinorbite ) then
                if( nsol == 1 ) then
                  p = l + 1._db
                elseif( isol == 1 ) then
                  p = 0.5_db + 0.5_db * sqrt( 1._db + 4*(l**2) + 8*l )
                else
                  if( l == 1 ) then
! En fait, il n'y a pas de solution pour l=1 avec spin-orbite non
! Relativiste !
!                      p = 1._db * l
                    p = l + 1._db
                  else
                    p = 0.5_db + 0.5_db * sqrt( -5._db + 4*(l**2) )
                  endif
                endif
              elseif( Relativiste ) then
                p = sqrt( l**2 + l + 1 - ( alfa_sf * numat )**2 )
              else
                p = l + 1._db
              endif

              if( Spinorbite .and. nsol /= 1 ) then 
                if( isol == 1 .and. isp == 1 ) then
                  fac = sqrt( ( l - m ) / ( 2*l + 1._db ) )
                elseif( isol == 1 .and. isp == 2 ) then 
                  fac = - sqrt( ( l + m ) / ( 2*l + 1._db ) )
                elseif( isol == 2 .and. isp == 1 ) then 
                  fac = sqrt( ( l + m + 1 ) / ( 2*l + 1._db ) )
                else 
                  fac = sqrt( ( l - m + 1 ) / ( 2*l + 1._db ) )
                endif
              else
                fac = 1._db
              endif
              if( nlm2 == 1 ) then
                np = 1
              elseif( .not. Spinorbite .or. isp == isol ) then
                np = n
              elseif( isp < isol ) then
                np = n + 1
              else
                np = n - 1
              endif
              if( Radial_comp ) then
                u(1:2,n,np,isp,isol) = fac * ( r(1:2)**p ) * cmplx( 
     &                   ( 1 + br * r(1:2) + cr(isp) * r(1:2)**2 ),
     &                   ci * r(1:2)**2, db) 
              else
                ur(1:2,n,np,isp,isol) = fac * ( r(1:2)**p )
     &                    * ( 1 + br * r(1:2) + cr(isp) * r(1:2)**2 )
              endif
            end do
          end do
        end do
      end do

      if( Spinorbite .and. .not. Full_potential ) then
        l = ll
        allocate( fac_m(-l:l,2) )
        do m = -li,li
           fac_m(m,1) = sqrt( ( l - m ) * ( l + m + 1._db ) )
           fac_m(m,2) = sqrt( ( l + m ) * ( l - m + 1._db ) )
        end do
      endif

      do ir = 2,nr-1
        im = ir - 1
        ip = ir + 1

        do isol = 1,nspino
          do isp = 1,nspin
            isq = 3 - isp

            do l = li,lf

              if( Hubb_a .and. l == l_hubbard( numat ) )  then
                Hubb_m = .true.
              else
                Hubb_m = .false.
              endif  

              do m = -l,l

               if( Spino_simple .and. ( ( m == l .and. isp == 1 )
     &           .or. ( m == -l .and. isp == 2 ) ) .and. isp /= isol ) 
     &                                                            cycle
                if( Full_potential ) then
                  n = l**2 + l + 1 + m
                elseif( Hubb_m .or. Spinorbite ) then
                  n = l + 1 + m
                else
                  n = 1
                  if( m /= 0 ) cycle
                endif

                td = g0(ir,isp) + l2 * f2(ir)

                if( Spinorbite ) then
                  if( isp == 1 ) then
                    ms = m
                    mp = m + 1   ! m de l'autre spin                 
                  else
                    ms = - m
                    mp = m - 1                 
                  endif
                  if( Full_potential ) then
                    fac = sqrt( ( l - ms ) * ( l + ms + 1._db ) )
                  else
                    fac = fac_m(m,isp)
                  endif
                  td = td + ms * gso(ir,isp)
                else
                  mp = l + 1 ! pour que le if en dessous soit faux
                endif
                if( Full_potential ) then
                  np = l**2 + l + 1 + mp
                elseif( Hubb_nd .or. Spinorbite ) then
                  np = l + 1 + mp
                else
                  np = 1
                endif

! le deuxieme indice "m" est la composante non nulle a l'origine
                if( Radial_comp ) then
                  u(ip,n,:,isp,isol) = td * u(ir,n,:,isp,isol)
     &                              + gm(ir,isp) * u(im,n,:,isp,isol)

                  if( Ecomp ) u(ip,n,:,isp,isol) = u(ip,n,:,isp,isol)
     &                - cmplx( 0._db, Eimag, db ) * u(ir,n,:,isp,isol)

                  if( abs(mp) <= l ) u(ip,n,:,isp,isol)
     &                      = u(ip,n,:,isp,isol)
     &                      + fac * gso(ir,isp) * u(ir,np,:,isq,isol)

                else
                  ur(ip,n,:,isp,isol) = td * ur(ir,n,:,isp,isol)
     &                              + gm(ir,isp) * ur(im,n,:,isp,isol)

                  if( abs(mp) <= l ) ur(ip,n,:,isp,isol)
     &                      = ur(ip,n,:,isp,isol)
     &                      + fac * gso(ir,isp) * ur(ir,np,:,isq,isol)
                endif

                if( Hubb_m ) then
                  do mp = -l,l
                    if( Hubb_d .and. m /= mp ) cycle
                    if( Full_potential ) then
                      np = l**2 + l + 1 + mp
                    else
                      np = l + 1 + mp
                    endif
                    if( Radial_comp ) then
                      u(ip,n,:,isp,isol) = u(ip,n,:,isp,isol)
     &                      + V_hubb(m,mp,isp)
     &                      * u(ir,np,:,isp,isol)
                    else
                      ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol)
     &                      + real( V_hubb(m,mp,isp), db)
     &                      * ur(ir,np,:,isp,isol)
                    endif
                  end do
                endif
                if( Radial_comp ) then
                  u(ip,n,:,isp,isol) = - u(ip,n,:,isp,isol) * gp(ir,isp)
                else
                  ur(ip,n,:,isp,isol) = - ur(ip,n,:,isp,isol)*gp(ir,isp)
                endif

              end do
            end do
          end do
        end do
      end do

      if( Spinorbite .and. .not. Full_potential ) deallocate( fac_m )

      if( icheck > 2 ) call write_ur(Full_potential,li,lf,nlm1,
     &       nlm2,nr,nspin,nspino,numat,r,Radial_comp,Rmtg,u,ur,1)

      if( .not. Renorm ) return

      do n = 1,nlm1
        do np = 1,nlm2
          do isp = 1,nspin
            do isol = 1 ,nspino
              if( Radial_comp ) then
                u(:,n,np,isp,isol) = u(:,n,np,isp,isol) / r(:)
              else
                ur(:,n,np,isp,isol) = ur(:,n,np,isp,isol) / r(:)
              endif
            end do
          end do
        end do
      end do

      call Renormal(Ecinetic,Radial_comp,Full_potential,icheck,konde,ll,
     &              lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,r,Rmtg,
     &              Tau,u,ur)

      if( icheck > 2 ) then
        do ir = 1,nr
          if( Radial_comp ) then
            u(ir,:,:,:,:) = u(ir,:,:,:,:) * r(ir)
          else
            ur(ir,:,:,:,:) = ur(ir,:,:,:,:) * r(ir)
          endif
        end do
        call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspin,
     &                nspino,numat,r,Radial_comp,Rmtg,u,ur,2)
        do ir = 1,nr
          if( Radial_comp ) then
            u(ir,:,:,:,:) = u(ir,:,:,:,:) / r(ir)
          else
            ur(ir,:,:,:,:) = ur(ir,:,:,:,:) / r(ir)
          endif
        end do
      endif

      return
      end

!***********************************************************************

! Ecriture de la fonction radiale

      subroutine write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspin,
     &                    nspino,numat,r,Radial_comp,Rmtg,u,ur,icom)

      use declarations
      implicit none

      integer:: i, icom, ir, isp, isol, k, kmax, l, li, lf,
     &  lp, m, nlm1, nlm2, mp, n, np, nr, nspin, nspino, numat

      Character(len=14):: mot, mot14
      Character(len=14), dimension(196):: Label, Labeli

      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u 

      logical:: Full_potential, Radial_comp 

      real(kind=db):: Rmtg
      real(kind=db), dimension(nr):: r 
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur 

      if( nlm1 == 1 .and. nspin == 1 ) then
        Label(1) = '     ur       ' 
        Labeli(1) = '     ui       '
        kmax = 1 
      else
        k = 0
        do isol = 1,nspino
          np = 0
          do lp = li,lf 
            do mp = -lp,lp
              if( nlm2 == 1 .and. mp /= 0 ) cycle
              do isp = 1,nspin
                k = k + 1
                mot14 = ' ur(          '
                i = 5
                if( nspin == 2 ) then
                  if( isp == 1 ) then
                    mot14(i:i) = 'u'
                  else
                    mot14(i:i) = 'd'
                  endif
                  i = i + 1
                  if( nlm2 > 1 ) then
                    mot14(i:i) = ','
                    i = i + 1
                  endif
                endif
                if( nlm2 > 1 ) then
                  if( Full_potential ) then
                    call ad_number(abs(lp),mot14,14)
                    i = i + 1
                    mot14(i:i) = ','
                    i = i + 1
                  endif
                  if( mp < 0 ) then
                    mot14(i:i) = '-'
                    i = i + 1
                  endif
                  call ad_number(abs(mp),mot14,14)
                  i = i + 1
                endif
                if( nspino == 2 .and. nlm2 == 1 ) then
                  mot14(i:i) = ','
                  i = i + 1
                  if( isol == 1 ) then
                    mot14(i:i) = 'u'
                  else
                    mot14(i:i) = 'd'
                  endif
                  i = i + 1
                endif
                mot14(i:i) = ')'
                label(k) = mot14
                if( Radial_comp ) then
                  mot14(3:3) = 'i'
                  labeli(k) = mot14
                endif
              end do
            end do
          end do
          if( nlm2 > 1 ) exit
        end do
        kmax = k
        do k = 1,kmax
          mot = label(k)
          Call Center_word( mot, 14 )
          label(k) = mot
          if( Radial_comp ) then
            mot = labeli(k)
            Call Center_word( mot, 14 )
            labeli(k) = mot
          endif
       end do 
      endif  

      n = 0
      np = 0
      do l = li,lf

        if( nlm2 == 1 ) then
          do m = -l,l
            if( nlm1 == 1 .and. m /= 0 ) cycle
            if( nlm1 == 1 ) then
              if( icom == 1 ) then
                write(3,110) numat, l
              elseif( icom == 2 ) then
                write(3,120) numat, l
              else
                write(3,130) numat, l
              endif
            else
              if( icom == 1 ) then
                write(3,140) numat, l, m
              elseif( icom == 2 ) then
                write(3,150) numat, l, m
              else
                write(3,160) numat, l, m
              endif
            endif
            n = n + 1
            if( nspino == 2 ) then
              if( Radial_comp ) then
                write(3,170)
              else
                write(3,180)
              endif
            endif
            if( Radial_comp ) then
              write(3,190) ( Label(k), Labeli(k), k = 1,kmax )
              do ir = 1,nr
                write(3,200) r(ir)*bohr, ( ( ( u(ir,n,np,isp,isol),
     &            isp = 1,nspin), np = 1,nlm2), isol = 1,nspino )
                if( r(ir) > Rmtg ) exit
              end do
            else
              write(3,190) ( Label(k), k = 1,kmax )
              do ir = 1,nr
                write(3,200) r(ir)*bohr, ( ( ( ur(ir,n,np,isp,isol),
     &          isp = 1,nspin), np = 1,nlm2), isol = 1,nspino )
                if( r(ir) > Rmtg ) exit
              end do
            endif
          end do
        else
          do mp = -l,l
            np = np + 1
            do isol = 1,nspino
              if( nspino == 1 ) then
                if( icom == 1 ) then
                  write(3,140) numat, l, mp
                elseif( icom == 2 ) then
                  write(3,150) numat, l, mp
                else
                  write(3,160) numat, l, mp
                endif
              else
                if( icom == 1 ) then
                  write(3,210) numat, l, mp, isol
                elseif( icom == 2 ) then
                  write(3,220) numat, l, mp, isol
                else
                  write(3,230) numat, l, mp, isol
                endif
              endif
              if( Radial_comp ) then
                write(3,190) ( Label(k), Labeli(k), k = 1,kmax )
                do ir = 1,nr
                  write(3,200) r(ir)*bohr, ( ( u(ir,n,np,isp,isol),
     &              isp = 1,nspin), n = 1,nlm1 )
                  if( r(ir) > Rmtg ) exit
                end do
              else
                write(3,190) ( Label(k), k = 1,kmax )
                do ir = 1,nr
                  write(3,200) r(ir)*bohr, ( ( ur(ir,n,np,isp,isol),
     &            isp = 1,nspin), n = 1,nlm1 )
                  if( r(ir) > Rmtg ) exit
                end do
              endif
            end do
          end do
        endif
      end do

      return
  110 format(/' Radial wave function time r:  Z =',i3,', l =',i2,/)
  120 format(/' Radial wave function, time r, after normalization:',
     &        '  Z =',i3,', l =',i2,/)
  130 format(/' Radial singular function time r:  Z =',i3,', l =',i2,/)
  140 format(/' Radial wave function time r:  Z =',i3,', l =',i2,
     &        ', m =',i2,/)
  150 format(/' Radial wave function, time r, after normalization:',
     &        '  Z =',i3,', l =',i2,', m =',i2,/)
  160 format(/' Radial singular function time r:  Z =',i3,', l =',i2,
     &        ', m =',i2,/)
  170 format(37x,' Solution 1',46x,' Solution 2')
  180 format(23x,' Solution 1',18x,' Solution 2')
  190 format('     Radius    ',392a14)
  200 format(1p,491e14.6)
  210 format(/' Radial wave function time r:  Z =',i3,', l =',i2,
     &        ', m =',i2,', Solution',i2,/)
  220 format(/' Radial wave function, time r, after normalization:',
     &        '  Z =',i3,', l =',i2,', m =',i2,', Solution',i2,/)
  230 format(/' Radial singular function time r:  Z =',i3,', l =',i2,
     &        ', m =',i2,', Solution',i2,/)
      end

!***********************************************************************

! Normalisation des fonctions radiales pour que les fonctions de base 
! correspondent aux solutions obtenues par continuite au rayon
! muffin-tin avec les solutions du vide.
! Appele par Sch_radial
 
      subroutine Renormal(Ecinetic,Radial_comp,Full_potential,icheck,
     &                konde,ll,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,r,
     &                Rmtg,Tau,u,ur)

      use declarations
      implicit none

      integer icheck, i, inr, ir, isol, isp, j, l, ll, lmax, lp, m,  
     &  nlm1, nlm2, mp, n, np, nr, nrmtg, nspin, nspino

      complex(kind=db):: ai, bi, ci, fnormc, z
      complex(kind=db), dimension(0:lmax):: bess, neum
      complex(kind=db), dimension(nspin):: konde, s1, s2
      complex(kind=db), dimension(nlm1,nspin,nlm1,nspin)::  Tau
      complex(kind=db), dimension(nlm1,nlm2,nspin,nspino):: Ampl, u1, 
     &                                             u2, Wronske, Wronsks
      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u 
      complex(kind=db), dimension(2,0:lmax,nspin):: bs, nm

      logical Radial_comp, Full_potential

      real(kind=db):: a, b, c, d01, d02, d12, dh, fnorm,
     &     konder, Rmtg, Wronskout, zr
      real(kind=db), dimension(nspin)::  Ecinetic
      real(kind=db), dimension(nr):: r
      real(kind=db), dimension(0:lmax):: bessr, neumr
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur 

      inr = min( nrmtg, nr )
      d12 = r(inr-1) - r(inr-2)
      d01 = r(inr) - r(inr-1)
      d02 = r(inr) - r(inr-2)

      Wronske(:,:,:,:) = (0._db, 0._db)
      Wronsks(:,:,:,:) = (0._db, 0._db)
      u1(:,:,:,:) = (0._db, 0._db)
      u2(:,:,:,:) = (0._db, 0._db)

! Bessel, Neuman et leur derivees au rayon Rmtg
      do isp = 1,nspin

        if( Radial_comp ) then
          konder = real( abs(konde(isp)),db )
        else
          konder = sqrt( Ecinetic(isp) )
        endif
        fnorm = sqrt( konder / pi )
        fnormc = sqrt( konde(isp) / pi )

        j = 0
        dh = 0.0001_db
        do i = -1,1,2
          j = j + 1
          if( Radial_comp ) then
            z = konde(isp) * ( Rmtg + i * 0.5_db * dh )
            call cbessneu(fnormc,z,lmax,lmax,bess,neum)
            bs(j,0:lmax,isp) = bess(0:lmax)
            nm(j,0:lmax,isp) = neum(0:lmax)
          else
            zr = konder * ( Rmtg + i * 0.5_db * dh )
            call cbessneur(fnorm,zr,lmax,lmax,bessr,neumr)
            bs(j,0:lmax,isp) = cmplx( bessr(0:lmax), 0._db,db )
            nm(j,0:lmax,isp) = cmplx( neumr(0:lmax), 0._db,db )
          endif
        end do
        do l = 0,lmax
          bs(2,l,isp) = ( bs(2,l,isp) - bs(1,l,isp) ) / dh
          nm(2,l,isp) = ( nm(2,l,isp) - nm(1,l,isp) ) / dh
        end do

        if( Radial_comp ) then
          z = konde(isp) * Rmtg
          call cbessneu(fnormc,z,lmax,lmax,bess,neum)
          bs(1,0:lmax,isp) = bess(0:lmax)
          nm(1,0:lmax,isp) = neum(0:lmax)
        else
          zr = konder * Rmtg
          call cbessneur(fnorm,zr,lmax,lmax,bessr,neumr)
          bs(1,0:lmax,isp) = cmplx( bessr(0:lmax), 0._db,db )
          nm(1,0:lmax,isp) = cmplx( neumr(0:lmax), 0._db,db )
        endif

      end do

      do isp = 1,nspin 
        do isol = 1,nspino
          do m = 1,nlm1
            do mp = 1,nlm2

              if( Radial_comp ) then
                ai = ( u(inr-2,m,mp,isp,isol) * d01
     &               - u(inr-1,m,mp,isp,isol) * d02
     &               + u(inr,m,mp,isp,isol) * d12 )
     &             / ( d12 * d01 * d02 )
                bi = (u(inr-1,m,mp,isp,isol) - u(inr-2,m,mp,isp,isol))
     &             / d12  - ai * ( r(inr-1) + r(inr-2) )
                ci = u(inr,m,mp,isp,isol) - ai * r(inr)**2 - bi *r(inr)

                u1(m,mp,isp,isol) = ai * Rmtg**2 + bi * Rmtg + ci
                u2(m,mp,isp,isol) = 2 * ai * Rmtg + bi
              else
                a = ( ur(inr-2,m,mp,isp,isol) * d01
     &              - ur(inr-1,m,mp,isp,isol) * d02
     &              + ur(inr,m,mp,isp,isol) * d12 )
     &            / ( d12 * d01 * d02 )
                b = ( ur(inr-1,m,mp,isp,isol) - ur(inr-2,m,mp,isp,isol))
     &            / d12 - a * ( r(inr-1) + r(inr-2) )
                c = ur(inr,m,mp,isp,isol) - a * r(inr)**2 - b * r(inr)

                u1(m,mp,isp,isol) = cmplx( a * Rmtg**2 + b * Rmtg + c,
     &                                     0._db, db )
                u2(m,mp,isp,isol) = cmplx( 2 * a * Rmtg + b, 0._db, db )
              endif

            end do
          end do
        end do
      end do

! Wronskout = 1 / ( pi * rmtg**2 ) a cause de la normalisation en
! rac(k/pi) des fonctions de bessel et hankel
      Wronskout = 1 / ( pi * Rmtg**2 )

      np = 0
      do l = 0,lmax
        if( .not. Full_potential .and. l /= ll ) cycle
        s1(:) = - img * bs(1,l,:) + nm(1,l,:)
        s2(:) = - img * bs(2,l,:) + nm(2,l,:)
        do mp = -l,l
          if( nlm2 == 1 .and. mp > -l ) exit
          np = np + 1     
          do isp = 1,nspin
            Wronsks(:,np,isp,:) = u1(:,np,isp,:) * s2(isp)
     &                          - u2(:,np,isp,:) * s1(isp)
            Wronske(:,np,isp,:) = u1(:,np,isp,:) * bs(2,l,isp)
     &                          - u2(:,np,isp,:) * bs(1,l,isp)
          end do
        end do
      end do

      if( icheck > 2 ) then        
        n = 0
        do l = 0,lmax
          if( .not. Full_potential .and. l /= ll ) cycle
          write(3,110) l
          write(3,120) - img * bs(1,l,:) + nm(1,l,:)
          write(3,130) - img * bs(2,l,:) + nm(2,l,:)
          write(3,140) bs(1,l,:)
          write(3,150) bs(2,l,:)
          write(3,160) Wronskout
          if( nspino == 2 ) then
            write(3,162)
          elseif( nspin == 2 ) then
            write(3,164)
          else
            write(3,166)
          endif
          do m = -l,l
            if( nlm1 == 1 .and. m /= 0 ) cycle            
            n = n + 1
            np = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. lp /= ll ) cycle
              do mp = -lp,lp
                if( nlm2 == 1 .and. mp /= 0 ) cycle            
                np = np + 1
                write(3,170) l, m, lp, mp, ' Wronsks   =',
     &                         ( Wronsks(n,np,:,isol), isol = 1,nspino )
                write(3,170) l, m, lp, mp, ' u1        =',
     &                         ( u1(n,np,:,isol), isol = 1,nspino )
                write(3,170) l, m, lp, mp, ' u2        =',
     &                         ( u2(n,np,:,isol), isol = 1,nspino )
              end do
            end do
          end do
        end do
      endif

      call cal_ampl(Ampl,Full_potential,icheck,ll,lmax,nlm1,nlm2,
     &              nspin,nspino,Tau,Wronske,Wronsks,Wronskout)

      do ir = 1,nr

        if( Radial_comp ) then
          u1(:,:,:,:) =  u(ir,:,:,:,:)
          u(ir,:,:,:,:) = 0._db
        else
          u1(:,:,:,:) =  cmplx( ur(ir,:,:,:,:), 0._db, db )
          ur(ir,:,:,:,:) = 0._db
        endif

        if( nlm2 == 1 .and. nspino == 1 ) then

          do isp = 1,nspin
            do n = 1,nlm1
              z = Ampl(n,1,isp,1) * u1(n,1,isp,1)
              if( Radial_comp ) then
                u(ir,n,1,isp,1) = z
              else
                ur(ir,n,1,isp,1) = real( z, db )
              endif
            end do
          end do

        elseif( nlm2 == 1 .and. nspino == 2 ) then

          do isp = 1,nspin
            n = 0
            do m = -ll,ll
              n = n + 1
              do isol = 1,nspino ! le spin d'attaque devient l'indice de solution

                mp = m + isol - isp
                if( mp > ll .or. mp < -ll ) cycle
                np = n + isol - isp

                z = sum( Ampl(np,1,isol,:) * u1(n,1,isp,:) )

                if( Radial_comp ) then
                  u(ir,n,1,isp,isol) = z
                else
                  ur(ir,n,1,isp,isol) = real( z, db )
                endif

              end do
            end do
          end do

        else

          do isp = 1,nspin
            do n = 1,nlm1
              do isol = 1,nspino ! le spin d'attaque devient l'indice de solution
                do np = 1,nlm2

                  if( nspino == 1 ) then
!                    z = sum( Ampl(n,:,isp,isol) * u1(:,np,isp,isol) )
                    z = sum( Ampl(np,:,isp,isol) * u1(n,:,isp,isol) )
                  else
                    z = ( 0._db, 0._db )
                    do i = 1,nspino    
                      z = z + sum( Ampl(np,:,isol,i)
     &                                 * u1(n,:,isp,i) )
                    end do
                  endif

                  if( Radial_comp ) then
                    u(ir,n,np,isp,isol) = z
                  else
                    ur(ir,n,np,isp,isol) = real( z, db )
                  endif

                end do
              end do
            end do
          end do

        endif

      end do

      return
  110 format(/' Renormalisation parameters, l =',i2) 
  120 format(/' -i.hankel =',1p,8e11.3) 
  130 format('     deriv =',1p,8e11.3) 
  140 format('    bessel =',1p,8e11.3) 
  150 format('     deriv =',1p,8e11.3) 
  160 format(' Wronskout =',1p,8e11.3) 
  162 format('  l  m lp mp',21x,'up up',17x,'dn up',17x,'up dn',17x,
     &        'dn dn') 
  164 format('  l  m lp mp',  23x,'up',19x,'dn') 
  166 format('  l  m lp mp') 
  170 format(4i3,a12,1p,16e11.3) 
      end

!***********************************************************************

! Calcul des amplitudes des fonctions radiales

      subroutine cal_ampl(Ampl,Full_potential,icheck,ll,lmax,nlm1,nlm2,
     &                nspin,nspino,Tau,Wronske,Wronsks,Wronskout)

      use declarations
      implicit none

      integer:: i, icheck, isol, isp, ispp, j, l, ll, lmax, lp,
     &  m, nlm1, nlm2, mp, mq, n, nd, ndim, np, nspin, nu, nspino

      complex(kind=db):: Det
      complex(kind=db), dimension(nlm1,nspin,nlm1,nspin):: Tau
      complex(kind=db), dimension(nlm1,nlm2,nspin,nspino):: Amp, Ampl,
     &                                                 Wronske, Wronsks
      complex(kind=db), dimension(:,:), allocatable:: Mat_A, Mat_e,
     &                                                Mat_s, Mat_T

      logical:: Full_potential, Stop_job

      real(kind=db):: Wronskout

      Stop_job = .true.

      Amp(:,:,:,:) = (0._db, 0._db)
      Ampl(:,:,:,:) = (0._db, 0._db)
      Tau(:,:,:,:) = (0._db, 0._db)

      if( nlm2 == 1 ) then
        if( nspino == 1 ) then
          do isp = 1,nspin
            do n = 1,nlm1
              Ampl(n,1,isp,1) = Wronskout / Wronske(n,1,isp,1)
              Tau(n,isp,n,isp) = - Wronske(n,1,isp,1)
     &                         / Wronsks(n,1,isp,1)
            end do
          end do
        else
! Pour Ampl, isp est le spin d'attaque (donc l'indice de solution apres
! la renormalisation)
! Le m de Amp est m + 1/2 - is = m + isol - 1
          do m = -ll-1,ll
            if( m == -ll-1 ) then
              isp = 2
              isol = 2
              mq = m + isp - 1
              n = ll + 1 + mq
              Ampl(n,1,isp,isol) = Wronskout / Wronske(n,1,isp,isol) 
              Tau(n,isp,n,isp) = - Wronske(n,1,isp,isol)
     &                         / Wronsks(n,1,isp,isol)
            elseif( m == ll ) then
              isp = 1
              isol = 1
              mq = m + isp - 1
              n = ll + 1 + mq
              Ampl(n,1,isp,isol) = Wronskout / Wronske(n,1,isp,isol) 
              Tau(n,isp,n,isp) = - Wronske(n,1,isp,isol)
     &                         / Wronsks(n,1,isp,isol)
            else
              nu = ll + 1 + m 
              nd = ll + 1 + m + 1
              Det = Wronsks(nu,1,1,1) * Wronsks(nd,1,2,2)
     &            - Wronsks(nd,1,2,1) * Wronsks(nu,1,1,2) 
              Amp(nu,1,1,1) = Wronsks(nd,1,2,2) / Det 
              Amp(nu,1,1,2) = - Wronsks(nd,1,2,1) / Det 
              Amp(nd,1,2,1) = - Wronsks(nu,1,1,2) / Det 
              Amp(nd,1,2,2) = Wronsks(nu,1,1,1) / Det 

              do isp = 1,nspin
                n = ll + 1 + m + isp - 1
                do ispp = 1,nspin
                  np = ll + 1 + m + ispp - 1
                  Tau(n,isp,np,ispp) = (0._db,0._db)
                  do isol = 1,nspino
                    Tau(n,isp,np,ispp) = Tau(n,isp,np,ispp)
     &              - Amp(n,1,isp,isol) * Wronske(np,1,ispp,isol)
                  end do
                end do
              end do

              Det = Wronske(nu,1,1,1) * Wronske(nd,1,2,2)
     &            - Wronske(nd,1,2,1) * Wronske(nu,1,1,2) 
              Ampl(nu,1,1,1) = Wronskout * Wronske(nd,1,2,2) / Det 
              Ampl(nu,1,1,2) = - Wronskout * Wronske(nd,1,2,1) / Det 
              Ampl(nd,1,2,1) = - Wronskout * Wronske(nu,1,1,2) / Det 
              Ampl(nd,1,2,2) = Wronskout * Wronske(nu,1,1,1) / Det 
            endif 
          end do
        endif

      else

        ndim = nspino * nlm2
        allocate( Mat_e(ndim,ndim) )
        allocate( Mat_s(ndim,ndim) )
        allocate( Mat_A(ndim,ndim) )
        allocate( Mat_T(ndim,ndim) )

        do isp = 1,nspin ! spin d'attaque
          if( nspino == 1 .or. isp == 1 ) then
            Mat_e(:,:) = ( 0._db, 0._db )
            Mat_s(:,:) = ( 0._db, 0._db )
            Mat_A(:,:) = ( 0._db, 0._db )
            Mat_T(:,:) = ( 0._db, 0._db )
            i = 0
          endif

          n = 0
          do l = 0,lmax
            if( .not. Full_potential .and. l /= ll ) cycle
            do m = -l,l
              n = n + 1
              i = i + 1

              j = 0
              do ispp = 1,nspin  ! solution
                if( nspino == 1 .and. ispp /= isp ) cycle
                isol = min(nspino,ispp)
                np = 0
                do lp = 0,lmax  
                  if( .not. Full_potential .and. lp /= ll ) cycle
                  do mp = -lp,lp
                    np = np + 1
                    j = j + 1
! colonne harmonique reelle, ligne harmonique attaque 
                    Mat_e(j,i) = Wronske(n,np,isp,isol)
                    Mat_s(j,i) = Wronsks(n,np,isp,isol)
                  end do
                end do
              end do
            end do
          end do

          if( nspino == 2 .and. isp == 1 ) cycle

          if( icheck > 3 ) then
            if(nspino == 1 ) then
              write(3,110) 'incoming', isp
            else
              write(3,120) 'incoming'
            endif
            do i = 1,ndim
              write(3,130) i, Mat_e(i,:)
            end do
            if(nspino == 1 ) then
              write(3,110) 'outgoing', isp
            else
              write(3,120) 'outgoing'
            endif
            do i = 1,ndim
              write(3,130) i, Mat_s(i,:)
            end do
          endif

          call invcomp(ndim,Mat_s,ndim,ndim,0,Stop_job)

          Mat_A = - Wronskout * Mat_s
          Mat_T = - matmul( Mat_s, Mat_e )  

          Mat_s = Mat_T
          call invcomp(ndim,Mat_s,ndim,ndim,0,Stop_job)

          Mat_A = matmul( Mat_s, Mat_A )

          if( nspino == 1 ) then
            Ampl(:,:,isp,1) = Mat_A(:,:)
            Tau(:,isp,:,isp) = Mat_T(:,:)
          endif

        end do

        if( nspino == 2 ) then 

          i = 0
          do isp = 1,nspin
            n = 0
            do l = 0,lmax
              if( .not. Full_potential .and. l /= ll ) cycle
              do m = -l,l
                n = n + 1
                i = i + 1 

                j = 0
                do isol = 1,nspino
                  np = 0
                  do lp = 0,lmax  
                    if( .not. Full_potential .and. lp /= ll ) cycle
                    do mp = -lp,lp
                      np = np + 1
                      j = j + 1
                      Ampl(n,np,isp,isol) = Mat_A(i,j)
                      Tau(n,isp,np,isol) = Mat_T(i,j)
                    end do
                  end do
                end do
              end do
            end do
          end do

        endif

        deallocate( Mat_A, Mat_e, Mat_s, Mat_T )

      endif

      if( icheck > 2 ) then        
        if( ll == 0 .or. icheck > 1 ) then
          if( nlm1 == 1 .and. nspin == 1 ) then
            write(3,140)
          elseif( nlm1 == 1 ) then
            write(3,150)
          elseif( nlm2 == 1 .and. nspino == 2 ) then
            write(3,155)
          elseif( nspin == 1 ) then
            write(3,160)
          else
            write(3,170)
          endif
        endif
        do isp = 1,nspin
          n = 0
          do l = 0,lmax
            if( .not. Full_potential .and. l /= ll ) cycle
            do m = -l,l
              if( nlm1 == 1 .and. m /= 0 ) cycle            
              n = n + 1
              if( nlm1 == 1 .and. nspin == 1 ) then
                write(3,220) l, Ampl(n,1,isp,min(isp,nspino))
              elseif( nlm1 == 1 ) then
                write(3,230) l, isp, Ampl(n,1,isp,min(isp,nspino))
              elseif( nlm2 == 1 .and. nspin == 1 ) then
                write(3,230) l, m, Ampl(n,1,isp,:)
              elseif( nlm2 == 1 ) then
                write(3,240) l, m, isp, Ampl(n,1,isp,:)
              elseif( nspin == 1 ) then
                write(3,230) l, m,
     &                         ( Ampl(n,:,isp,ispp), ispp = 1,nspino )
              else
                write(3,240) l, m, isp,
     &                         ( Ampl(n,:,isp,ispp), ispp = 1,nspino )
              endif
            end do
            if( nlm2 == 1 .and. nspino == 1 ) exit
          end do
        end do
      endif

      if( icheck > 1 ) then
        if( ll == 0 .or. icheck > 1 ) then
          if( nlm1 == 1 .and. nspin == 1 ) then
            write(3,180)
          elseif( nlm1 == 1 ) then
            write(3,190)
          elseif( nspin == 1 ) then
            write(3,200)
          else
            write(3,210)
          endif
        endif
        do isp = 1,nspin
          n = 0
          do l = 0,lmax
            if( .not. Full_potential .and. l /= ll ) cycle
            do m = -l,l
              if( nlm1 == 1 .and. m /= 0 ) cycle            
              n = n + 1
              if( nlm1 == 1 .and. nspin == 1 ) then
                write(3,220) l, Tau(n,isp,n,isp)
              elseif( nlm1 == 1 ) then
                write(3,230) l, isp, Tau(n,isp,n,isp)
              elseif( nlm2 == 1 .and. nspin == 1 ) then
                write(3,230) l, m, Tau(n,isp,n,isp)
              elseif( nlm2 == 1 .and. nspino == 1 ) then
                write(3,240) l, m, isp, Tau(n,isp,n,isp)
              else
                write(3,240) l, m, isp,
     &                         ( Tau(n,isp,:,ispp), ispp = 1,nspin )
              endif
            end do
            if( nlm2 == 1 .and. nspino == 1 ) exit
          end do
        end do
      endif

      return
  110 format(/' Wronskian ',a8,' matrix, isp =',i2)
  120 format(/' Wronskian ',a8,' matrix')
  130 format(i3,1p,100(2e11.3,1x))
  140 format(/' Amplitude',/'  l',11x,'Ampl')
  150 format(/' Amplitude',/'  l isp',9x,'Ampl')
  155 format(/' Amplitude',/'  l  m isp',6x,'Ampl(Sol 1)',12x,
     &                                   'Ampl(Sol 2)')
  160 format(/' Amplitude',/'  l  m',10x,'Ampl')
  170 format(/' Amplitude',/'  l  m isp',10x,'Ampl')
  180 format(/' Atomic scattering amplitude',/'  l',11x,'Tau') 
  190 format(/' Atomic scattering amplitude',/'  l isp',10x,'Tau')
  200 format(/' Atomic scattering amplitude',/'  l  m',11x,'Tau') 
  210 format(/' Atomic scattering amplitude',/'  l  m isp',10x,'Tau')
  220 format(i3,1p,16(1x,2e11.3)) 
  230 format(2i3,1p,16(1x,2e11.3)) 
  240 format(3i3,1p,16(1x,2e11.3)) 
      end

!**********************************************************************

! Calcul des integrales radiales et de la solution singuliere
! Appele par tenseur_car

      subroutine radial(Ecinetic,Eimag,Energ,Enervide,Eseuil,
     &       Final_tddft,Full_potential,Green_plus,Hubb_a,Hubb_d,icheck,
     &       initlv,ip_max,ip0,lmax,m_hubb,nbseuil,nlma,nlma2,ninit1,
     &       ninitlv,nr,nrm,nspin,nspino,numat,psii,r,Relativiste,Rmtg,
     &       Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,
     &       V_intmax,V0bd,Vrato,Ylm_comp)
 
      use declarations
      implicit none

      integer:: i, icheck, initl, initlv, ip, ip_max, ip0, ip1,
     &  ip2, iseuil, isp, l, l_hubbard, lfin, lm, 
     &  lmp, lmax, lp, m, m_hubb, nlm1, nlm2, mp, nbseuil, nlma, nlma2, 
     &  ninit1, ninitlv, nr, nrm, nrmtsd, nrmtg, nspin, nspino, numat

      character(len=2):: mot1, mot2      
      character(len=104):: mot      

      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(nlma,nspin,ip0:ip_max,ip0:ip_max,
     &                                                ninitlv):: Singul
      complex(kind=db), dimension(nlma,nlma2,nspin,nspino,ip0:ip_max,
     &                                                   ninitlv):: rof
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                       V_hubb
      complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: u

      logical Ecomp, Final_tddft, Full_potential, Green_plus, Hubb_a,
     &  Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, Solsing,
     &  Spinorbite, Ylm_comp

      real(kind=db):: Eimag, Energ, Enervide, Ephoton,
     &                Rmtg, Rmtsd, V_intmax
      real(kind=db), dimension(nspin):: Ecinetic, V0bd
      real(kind=db), dimension(nr,nspin):: Vrato 

      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nbseuil):: Eseuil, Vecond
      real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gpi, gp, V
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nrm,nbseuil):: psii

      real(kind=db), dimension(:,:,:,:,:), allocatable:: ur

      konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

      if( nbseuil == 2 ) then
        if( Final_tddft .and. nbseuil /= ninitlv ) then
          ninit1 = ( ninitlv - 2 ) / 2
          if( initlv <= ninit1 ) then
            iseuil = 1
          else
            iseuil = 2
          endif
        else ! non utilise en DFT
          iseuil = initlv
        endif
      else
        iseuil = 1
      endif

      if( icheck > 1 ) then
        write(3,110) initlv, iseuil
        write(3,120) Energ*rydb
        write(3,130) Ecinetic(:)*rydb
        write(3,140) V0bd(:)*rydb
        write(3,150) konde(:)
      endif

      do i = 1,nbseuil
        if( Final_tddft ) then
          if( i /= nbseuil ) cycle
          Ephoton = Energ + Eseuil(nbseuil)
        else
          Ephoton = Energ + Eseuil(i)
        endif
! Terme multiplicatif pour les transitions quadrupolaires
! En S.I. vecond = k = E*alfa_sf*4*pi*epsilon0 / (e*e)
! En ua et rydb : k = 0.5 * alfa_sf * E
        Vecond(i) = 0.5 * alfa_sf * Ephoton
      end do

      call mod_V(icheck,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,
     &                 V_intmax,V0bd,Vrato)

      call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nr,nspin,nspino,
     &                        numat,r,Relativiste,Spinorbite,V)

      gpi(:,:) = 1 / gp(:,:)
      if( Solsing ) gmi(:,:) = 1 / gm(:,:)

      if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10
     &     .or. Ecinetic(nspin) < eps10 ) then
        Ecomp = .true.
      else
        Ecomp = .false.
      endif
      Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

      if( Full_potential ) then
        lfin = 0
      else
        lfin = lmax
      endif
      Renorm = .true.

      do l = 0,lfin

        if( Hubb_a .and. l == l_hubbard( numat ) )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif  

        if( Full_potential ) then
          nlm1 = ( lmax + 1 )**2
          nlm2 = nlm1
        elseif( Hubb_m .and. .not. Hubb_d ) then
          nlm1 = 2*l + 1
          nlm2 = nlm1
        elseif( Spinorbite .or. Hubb_m ) then
          nlm1 = 2*l + 1
          nlm2 = 1
        else
          nlm1 = 1
          nlm2 = 1
        endif

        allocate( u(nr,nlm1,nlm2,nspin,nspino) )
        allocate( ur(nr,nlm1,nlm2,nspin,nspino) )
        allocate( Tau(nlm1,nspin,nlm1,nspin) )

        call Sch_radial(Ecinetic,Ecomp,Eimag,f2,
     &          Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,
     &           l,lmax,m_hubb,nlm1,nlm2,nr,nrmtg,nspin,nspino,numat,r,
     &           Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,u,
     &           ur,V_hubb)

! Integrale radiale pour la regle d'or de Fermi
        call radial_matrix(Final_tddft,Green_plus,initlv,ip_max,ip0,
     &         iseuil,l,nlm1,nlm2,nbseuil,
     &         ninitlv,nlma,nlma2,nr,nrm,nrmtsd,nspin,
     &         nspino,psii,r,Radial_comp,Rmtsd,rof,u,ur,Vecond)

! Calcul de la solution singuliere
        if( Solsing )
     &    call Cal_Solsing(Ecomp,Eimag,f2,Final_tddft,
     &       Full_potential,g0,gmi,gp,Green_plus,gso,Hubb_a,Hubb_d,
     &       icheck,initlv,ip_max,ip0,iseuil,konde,l,lmax,m_hubb,nlm1,
     &       nlm2,nbseuil,ninitlv,nlma,nr,nrm,nrmtsd,nspin,nspino,
     &       numat,psii,r,Radial_comp,Rmtsd,Singul,Spinorbite,Tau,u,ur,
     &       V_hubb,Vecond)

        deallocate( Tau, u, ur )

      end do   ! fin de la boucle sur l

      if( icheck > 1 ) then
        mot = ' '
        if( Ecomp .and. nspino == 2 ) then
          mot(11:17) = 'up sol1'
          mot(37:43) = 'up sol2'
          mot(63:69) = 'dn sol1'
          mot(89:95) = 'dn sol2'
        elseif( nspino == 2 ) then
          mot(5:50) = 'up sol1      up sol2      dn sol1      dn sol2'
        elseif( nspin == 2 ) then
          mot = '      up           dn'
        else
          mot = ' '
        endif

        do i = 1,nbseuil
          if( Final_tddft ) then
            if( i /= nbseuil ) cycle
            initl = initlv
          else
            initl = i
            if( nbseuil > 1 ) write(3,155) i
          endif
        
          do ip = ip0,ip_max
            select case(ip)
              case(0)
                write(3,160) '  monopole'
              case(1)
                write(3,160) '    dipole'
              case(2)
                write(3,160) 'quadrupole'
              case(3)
                write(3,160) '  octupole'
            end select
            if( Full_potential .or. Hubb_a ) then
              write(3,170) '  l  m lp mp', mot
            else
              write(3,175) '  l  m', mot
            endif
            lm = 0
            do l = 0,lmax
              do m = -l,l
                lm = lm + 1
                if( .not. ( Hubb_a .or. Full_potential .or. Spinorbite )
     &                                             .and. m /= 0 ) cycle
                lmp = 0
                do lp = 0,lmax
                  if( .not. Full_potential .and. l /= lp ) cycle
                  do mp = -lp,lp
                    if( .not. ( Hubb_a .or. Full_potential
     &                                 .or. m == mp ) ) cycle
                    lmp = lp**2 + lp + 1 + mp
                    lmp = min(nlma2,lmp)
                    if( Full_potential .or. Hubb_a ) then
                      if( Ecomp) then
                        write(3,180) l, m, lp, mp,
     &                      (rof(lm,lmp,isp,:,ip,initl),isp = 1,nspin)
                      else
                        write(3,180) l, m, lp, mp,
     &                  (real(rof(lm,lmp,isp,:,ip,initl)),isp = 1,nspin)
                      endif
                    else
                      if( Ecomp) then
                        write(3,185) l, m,
     &                      (rof(lm,lmp,isp,:,ip,initl),isp = 1,nspin)
                      else
                        write(3,185) l, m,
     &                  (real(rof(lm,lmp,isp,:,ip,initl)),isp = 1,nspin)
                      endif
                    endif
                  end do
                end do
              end do
            end do
          end do

          if( .not. Solsing ) cycle

          do ip1 = ip0,ip_max
            select case(ip1)
              case(0)
                mot1 = 'M1'
              case(1)
                mot1 = 'E1'
              case(2)
                mot1 = 'E2'
              case(3)
                mot1 = 'E3'
            end select

            do ip2 = ip0,ip_max
              if( .not. Full_potential .and. mod(ip1+ip2,2) == 1 ) cycle

              select case(ip2)
                case(0)
                  mot2 = 'M1'
                case(1)
                  mot2 = 'E1'
                case(2)
                  mot2 = 'E2'
                case(3)
                  mot2 = 'E3'
              end select

              write(3,190) mot1, mot2

              if( nspin == 2 ) then
                write(3,210) '  l  m','   up  ','   dn  '
              else
                write(3,210) '  l  m','  Sing '
              endif

              lm = 0
              do l = 0,lmax
                do m = -l,l
                  lm = lm + 1
                  if( .not. ( Hubb_m .or. Full_potential
     &                          .or. Spinorbite ) .and. m /= 0 ) cycle
                  write(3,230) l, m, Singul(lm,:,ip1,ip2,initl)
                end do
              end do
            end do
          end do
        end do

      endif

      return
  110 format(/' initl =',i2,', iseuil =',i2)
  120 format(' Energ    =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd     =',2f10.3)
  150 format(' konde    =',2f12.5)
  155 format(/' iseuil =',i2)
  160 format(/' Radial matrix, ',a10,1x,'term:')
  170 format(a12,a104)
  175 format(a6,a104)
  180 format(4i3,1p,8e13.5)
  185 format(2i3,1p,8e13.5)
  190 format(/' Radial integral of singular solution, ',a2,'-',a2,
     &        ' term:')
  210 format(a6,9x,a7,3(19x,a7))
  230 format(2i3,1p,4e13.5)
      end

!**********************************************************************

! Calcul de l'integrale radiale de la regle d'or de Fermi.
! psii, u et ur sont les fonctions d'onde fois r.

      subroutine radial_matrix(Final_tddft,Green_plus,initlv,ip_max,ip0,
     &         iseuil,l,nlm1,nlm2,nbseuil,
     &         ninitlv,nlma,nlma2,nr,nrm,nrmtsd,nspin,
     &         nspino,psii,r,Radial_comp,Rmtsd,rof,u,ur,Vecond)

      use declarations
      implicit none

      integer initlv, ip, ip_max, ip0, ipp, ir, iseuil, is, isol,
     &  isp, iss, l, lm, lm1, lm2, nlm1, nlm2, n1, n2, nbseuil,
     &  ninitlv, nlma, nlma2, nr, nrm, nrmtsd, ns1, ns2,
     &  nspin, nspino

      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u
      complex(kind=db), dimension(nlma,nlma2,nspin,nspino,ip0:ip_max,
     &                                                   ninitlv):: rof

      logical:: Final_tddft, Green_plus, Radial_comp

      real(kind=db):: f_integr3, fac, radlr, radli, rmtsd
      real(kind=db), dimension(nbseuil):: Vecond
      real(kind=db), dimension(nr):: r, fct
      real(kind=db), dimension(nrm,nbseuil):: psii
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur

      if( Final_tddft ) then
        ns1 = initlv
        ns2 = initlv
      else
        ns1 = 1
        ns2 = nbseuil
      endif

      lm1 = l**2
      lm2 = (l + 1)**2

      do ip = ip0,ip_max

        ipp = ip + 1
  
        do is = ns1,ns2

          if( Final_tddft ) then
            iss = nbseuil   ! c'est le 2eme seuil qui sert d'origine
          else
            iss = is
            iseuil = is 
          endif

          select case(ip)
            case(0)
              fac = 0.5_db * alfa_sf 
            case(1)
              fac = 1._db
            case(2)
              fac = 0.5_db * Vecond(iss)
              if( .not. Green_plus ) fac = - fac
            case(3)
            fac = - ( 1._db / 6 ) * Vecond(iss)**2 
          end select

          do isp = 1,nspin  
            do n1 = 1,nlm1
              do isol = 1,nspino
                do n2 = 1,nlm2

                  if( Radial_comp ) then
                    do ir = 1,nrmtsd
                      fct(ir) = psii(ir,iseuil)
     &                        * real( u(ir,n1,n2,isp,isol), db )
     &                        * r(ir)**ipp
                    end do
                    radlr = fac * f_integr3(r,fct,1,nr,Rmtsd)
                    do ir = 1,nrmtsd
                      fct(ir) = psii(ir,iseuil)
     &                        * aimag( u(ir,n1,n2,isp,isol) )
     &                        * r(ir)**ipp
                    end do
                    radli = fac * f_integr3(r,fct,1,nr,Rmtsd)
                  else
                    do ir = 1,nrmtsd
                      fct(ir) = psii(ir,iseuil) * ur(ir,n1,n2,isp,isol)
     &                        * r(ir)**ipp
                    end do
                    radlr = fac * f_integr3(r,fct,1,nr,Rmtsd)
                    radli = 0._db
                  endif

                  if( nlm1 /= 1 .and. nlm2 /= 1 ) then
                    rof(lm1+n1,lm1+n2,isp,isol,ip,is)
     &                                        = cmplx( radlr, radli,db)
                  elseif( nlm1 /= 1 .and. nlma2 /= 1 ) then
                    rof(lm1+n1,lm1+n1,isp,isol,ip,is)
     &                                        = cmplx( radlr, radli,db)
                  elseif( nlm1 /= 1 ) then
                    rof(lm1+n1,n2,isp,isol,ip,is)
     &                                        = cmplx( radlr, radli,db)
                  elseif( nlm1 == 1 .and. nlma2 /= 1 ) then
                    do lm = lm1+1,lm2 
                      rof(lm,lm,isp,isol,ip,is) = cmplx(radlr, radli,db)
                    end do
                  else
                    do lm = lm1+1,lm2 
                      rof(lm,n2,isp,isol,ip,is) = cmplx(radlr, radli,db)
                    end do
                  endif

                end do
              end do
            end do
          end do
        end do
      end do

      return
      end

!**********************************************************************

! Calcul de l'integrale radiale de la regle d'or de Fermi.
! u et ur sont les fonctions d'onde fois r.

      subroutine radial_matrix_optic(Green_plus,ip_max,ip0,
     &         ne,nlm1g,nlm2g,nr,nrmtsd,nspin,
     &         nspino,r,Radial_comp,Rmtsd,roff_ii,roff_ir,roff_ri,
     &         roff_rr,ui,ur,Vecond)

      use declarations
      implicit none

      integer:: ip, ip_max, ip0, ipp, ir, isol, isoli, isolf,
     &  isp, ispi, ispf, nlm1g, nlm2g, n1i, n1f, n2i, n2f, ne,
     &  nr, nrmtsd, nspin, nspino

      logical:: Green_plus, Radial_comp

      real(kind=db):: f_integr3, fac, Rmtsd, Vecond
      real(kind=db), dimension(nr):: r, rp, fct
      real(kind=db), dimension(nr,nlm1g,nlm2g,nspin,nspino,ne):: ui, ur
      real(kind=db), dimension(nlm1g,nlm1g,nlm2g,nlm2g,nspin**2,
     &       nspino**2,ip0:ip_max):: roff_ii, roff_ir, roff_ri, roff_rr


      do ip = ip0,ip_max

        ipp = ip + 1

        rp(1:nrmtsd) = r(1:nrmtsd)**ipp
  
        select case(ip)
          case(0)
            fac = 0.5_db * alfa_sf 
          case(1)
            fac = 1._db
          case(2)
            fac = 0.5_db * Vecond
            if( .not. Green_plus ) fac = - fac
          case(3)
          fac = - ( 1._db / 6 ) * Vecond**2 
        end select

        do ispi = 1,nspin  
          do n1i = 1,nlm1g
            do isoli = 1,nspino
              do n2i = 1,nlm2g

                do ispf = 1,nspin
                  isp = ispf + ( ispi - 1 ) * nspin  
                  do n1f = 1,nlm1g
                    do isolf = 1,nspino
                      isol = isolf + ( isoli - 1 ) * nspino  
                      do n2f = 1,nlm2g

                        do ir = 1,nrmtsd
                          fct(ir) = ur(ir,n1i,n2i,ispi,isoli,1)
     &                            * ur(ir,n1f,n2f,ispf,isolf,2)
     &                            * rp(ir)
                        end do

                        roff_rr(n1i,n1f,n2i,n2f,isp,isol,ip)
     &                       = fac * f_integr3(r,fct,1,nr,Rmtsd)

                        if( Radial_comp ) then
                          do ir = 1,nrmtsd
                            fct(ir) = ur(ir,n1i,n2i,ispi,isoli,1)
     &                              * ui(ir,n1f,n2f,ispf,isolf,2)
     &                              * rp(ir)
                          end do

                          roff_ri(n1i,n1f,n2i,n2f,isp,isol,ip)
     &                       = fac * f_integr3(r,fct,1,nr,Rmtsd)

                          do ir = 1,nrmtsd
                            fct(ir) = ui(ir,n1i,n2i,ispi,isoli,1)
     &                              * ur(ir,n1f,n2f,ispf,isolf,2)
     &                              * rp(ir)
                          end do

                          roff_ir(n1i,n1f,n2i,n2f,isp,isol,ip)
     &                       = fac * f_integr3(r,fct,1,nr,Rmtsd)

                          do ir = 1,nrmtsd
                            fct(ir) = ui(ir,n1i,n2i,ispi,isoli,1)
     &                              * ui(ir,n1f,n2f,ispf,isolf,2)
     &                              * rp(ir)
                          end do

                          roff_ii(n1i,n1f,n2i,n2f,isp,isol,ip)
     &                       = fac * f_integr3(r,fct,1,nr,Rmtsd)

                        endif

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      return
      end

!***********************************************************************

! Calcul de la solution singuliere pour l'absorption
! Appele par Radial

      Subroutine Cal_Solsing(Ecomp,Eimag,f2,Final_tddft,
     &       Full_potential,g0,gm,gp,Green_plus,gso,Hubb_a,Hubb_d,
     &       icheck,initlv,ip_max,ip0,iseuil,konde,ll,lmax,m_hubb,nlm1,
     &       nlm2,nbseuil,ninitlv,nlma,nr,nrm,nrmtsd,nspin,nspino,
     &       numat,psii,r,Radial_comp,Rmtsd,Singul,Spinorbite,Tau,u,ur,
     &       V_hubb,Vecond)

      use declarations
      implicit none

      integer:: icheck, initlv, ip_max, ip0, ip1, ip2,
     &  is, ise, iseuil, isol, isp, iss,  l, ll, lmax, m_hubb, m, 
     &  ms, mv, n0, n, nbseuil, ninitlv, nlm1, nlm2, np, nlma,
     &  nr, nrm, nrmtsd, ns1, ns2, nspin, nspino, numat, nv 

      complex(kind=db):: integr_sing, Sing
      complex(kind=db), dimension(nspin):: konde 
      complex(kind=db), dimension(nlma,nspin,ip0:ip_max,ip0:ip_max,
     &                                            ninitlv):: Singul
      complex(kind=db), dimension(nrmtsd):: f_reg, f_irg 
      complex(kind=db), dimension(nlm1,nspin,nlm1,nspin):: Tau
      complex(kind=db), dimension(nrmtsd+1,nlm1,nlm2,nspin,nspino):: us 
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                         V_hubb
      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u

      logical:: Ecomp, Final_tddft, Full_potential, Green_plus, Hubb_a,
     &  Hubb_d, Radial_comp, Spinorbite 

      real(kind=db):: Eimag, fac1, fac2, Rmtsd
      real(kind=db), dimension(nbseuil):: Vecond
      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nrmtsd):: rr, phi1, phi2 
      real(kind=db), dimension(nrm,nbseuil):: psii
      real(kind=db), dimension(nr,nspin):: g0, gm, gp
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur

      rr(1:nrmtsd) = r(1:nrmtsd)

      if( Final_tddft ) then
        ns1 = initlv
        ns2 = initlv
      else
        ns1 = 1
        ns2 = nbseuil
      endif
 
! gm et gp inverse dans le sousprogramme
      call Sch_radial_solsing(Ecomp,Eimag,f2,
     &         Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,
     &         ll,lmax,m_hubb,nlm1,nlm2,nr,nrmtsd,nspin,nspino,numat,r,
     &         Radial_comp,Rmtsd,Spinorbite,Tau,us,V_hubb)

      do is = ns1,ns2

        if( Final_tddft ) then
          iss = nbseuil   ! c'est le 2eme seuil qui sert d'origine
          ise = iseuil
        else
          iss = is
          ise = is 
        endif

        do ip1 = ip0,ip_max  ! boucle sur dipole, quadrupole, Octupole

          select case(ip1)
            case(0)
              fac1 = - 0.5_db * alfa_sf 
            case(1)
              fac1 = 1._db
            case(2)
              fac1 = 0.5_db * Vecond(iss)
              if( .not. Green_plus ) fac2 = - fac2
            case(3)
              fac1 = - ( 1._db / 6 ) * Vecond(iss)**2 
          end select

          if( ip1 == 0 ) then
            phi1(1:nrmtsd) = psii(1:nrmtsd,ise)
          else
            phi1(1:nrmtsd) = psii(1:nrmtsd,ise) * r(1:nrmtsd)**ip1
          endif

          do ip2 = ip0,ip_max  ! boucle sur dipole, quadrupole, Octupole

            select case(ip2)
              case(0)
                fac2 = - 0.5_db * alfa_sf 
              case(1)
                fac2 = 1._db
              case(2)
                fac2 = 0.5_db * Vecond(iss)
                if( .not. Green_plus ) fac2 = - fac2
              case(3)
                fac2 = - ( 1._db / 6 ) * Vecond(iss)**2 
            end select

            if( ip2 == 0 ) then
              phi2(1:nrmtsd) = psii(1:nrmtsd,ise)
            else
              phi2(1:nrmtsd) = psii(1:nrmtsd,ise) * r(1:nrmtsd)**ip2
            endif

            n = 0
            do l = 0,lmax
              if( .not. Full_potential .and. l /= ll ) cycle
              n0 = l**2 + l + 1
              do m = -l,l
                if( nlm1 == 1 .and. m /= 0 ) cycle 
                n = n + 1

                do isp = 1,nspin

                  do np = 1,nlm2

                    do isol = 1,nspino
                      if( Spinorbite .and. nlm2 == 1 ) then
                        ms = m + isol - isp
                        if( ms > l .or. ms < -l ) cycle
                      endif                         

                      if( Radial_comp ) then
                        f_reg(1:nrmtsd) = r(1:nrmtsd)
     &                                  * u(1:nrmtsd,n,np,isp,isol)
                      else
                        f_reg(1:nrmtsd) = r(1:nrmtsd)
     &                      * cmplx( ur(1:nrmtsd,n,np,isp,isol),
     &                               0._db, db )
                      endif

                      f_irg(1:nrmtsd) = us(1:nrmtsd,n,np,isp,isol)

                      Sing = fac1 * fac2 * integr_sing(nrmtsd,phi1,phi2,
     &                                      f_reg,f_irg,Rmtsd,rr,icheck)

                      if( nlm1 == 1 ) then
                        do mv = -l,l
                          nv = n0 + mv
                          Singul(nv,isp,ip1,ip2,is)
     &                        = Singul(nv,isp,ip1,ip2,is) + Sing
                        end do
                      else
                        nv = n0 + m
                        Singul(nv,isp,ip1,ip2,is)
     &                        = Singul(nv,isp,ip1,ip2,is) + Sing
                      endif

                    end do
                  end do

                end do    ! fin boucle isp
              end do   ! fin boucle m
            end do  ! fin boucle l
          
          end do ! fin boucle dipole, quadrupole
        end do ! fin boucle dipole, quadrupole

      end do ! fin boucle seuil

      return
      end

!***********************************************************************

! Resolution de l'equation de Schrodinger radiale pour la solution irreguliere

! Recopie de Sch_radial: gm et gp sont inverses dans l'appel

      subroutine Sch_radial_solsing(Ecomp,Eimag,f2,
     &         Full_potential,g0,gp,gm,gso,Hubb_a,Hubb_d,icheck,konde,
     &         ll,lmax,m_hubb,nlm1,nlm2,nr,nrmtsd,nspin,nspino,numat,r,
     &         Radial_comp,Rmtsd,Spinorbite,Tau,us,V_hubb)

      use declarations
      implicit none

      integer:: icheck, im, ip, ir, isol, isp, isq, l,
     &  l_hubbard, l2, li, ll, lf, lmax, lp, m, m_hubb, mp, ms, n,
     &  nlm1, nlm2, np, nr, nrmtsd, ns, nspin, nspino, numat

      logical:: Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Hubb_nd, 
     &  Radial_comp, Spinorbite

      complex(kind=db):: fnormc, z
      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(0:lmax):: bess, neum
      complex(kind=db), dimension(nlm1,nspin,nlm1,nspin):: Tau
      complex(kind=db), dimension(nrmtsd:nrmtsd+1,0:lmax,nspin):: Hankel
      complex(kind=db), dimension(nrmtsd+1,nlm1,nlm2,nspin,nspino):: us 
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                        V_hubb
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: u 

      real(kind=db):: Eimag, fac, fnorm, konder, Rmtsd, td, zr 
      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(0:lmax):: bessr, neumr
      real(kind=db), dimension(nr,nspin):: g0, gm, gp 
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(:), allocatable:: rr 
      real(kind=db), dimension(:,:,:,:,:), allocatable:: ur 

      us(:,:,:,:,:) = ( 0._db, 0._db )

! Hankel est en fait -i.hankel+
      do ir = nrmtsd,nrmtsd+1
        do isp = 1,nspin
          if( Ecomp ) then
            fnormc = sqrt( konde(isp) / pi )
            z = konde(isp) * r(ir)
            call cbessneu(fnormc,z,lmax,lmax,bess,neum)
            Hankel(ir,0:lmax,isp) = neum(0:lmax) - img * bess(0:lmax)
          else
            konder = Real( konde(isp), db )
            fnorm = sqrt( konder / pi )
            zr = konder * r(ir)
            call cbessneur(fnorm,zr,lmax,lmax,bessr,neumr)
            Hankel(ir,0:lmax,isp) = neumr(0:lmax) - img * bessr(0:lmax)
          endif
        end do
      end do

      if( Full_potential ) then
        li = 0
        lf = lmax
      else 
        li = ll
        lf = ll
      endif

      n = 0
      do l = li,lf

        do m = -l,l
          if( nlm1 == 1 .and. m /= 0 ) cycle
          n = n + 1

          do isp = 1,nspin  ! Spin reel

            np = 0
            do lp = li,lf
              if( .not. Full_potential .and. l /= lp ) cycle
              do mp = -lp,lp
                if( nlm2 == 1 .and. mp /= m ) cycle
                np = np  + 1

                do isq = 1,nspin  ! spin d'attaque ( indice solution )
                  isol = min( isq, nspino )  
                  if( nspino == 1 .and. isp /= isq ) cycle
                  if( nlm2 == 1 ) then
                    ms = m + isq - isp
                    if( ms > l .or. ms < -l ) cycle
                    ns = n + isq - isp
                  else
                    ns = np 
                  endif

                  do ir = nrmtsd,nrmtsd+1
                    us(ir,n,np,isp,isol) = Tau(ns,isq,n,isp) * r(ir)
     &                                   * Hankel(ir,l,isp)
                  end do
                end do
              end do
            end do
          end do

        end do
      end do

      do ir = nrmtsd,2,-1
        ip = ir - 1      ! inverse par rapport a Sch_radial
        im = ir + 1

        do isp = 1,nspin
          isq = 3 - isp

          do l = li,lf

            l2 = l * ( l + 1 )

            if( Hubb_a .and. l == l_hubbard(numat) ) then
              Hubb_m = .true.
            else
              Hubb_m = .false.
            endif
            Hubb_nd = Hubb_m .and. .not. Hubb_d

            do m = -l,l

              if( Full_potential ) then
                n = l**2 + l + 1 + m
              elseif( Hubb_m .or. Spinorbite ) then
                n = l + 1 + m
              else
                n = 1
                if( m /= 0 ) cycle
              endif

              td = g0(ir,isp) + l2 * f2(ir)

              if( Spinorbite ) then
                if( isp == 1 ) then
                  ms = m
                  mp = m + 1   ! m de l'autre spin                 
                else
                  ms = - m
                  mp = m - 1                 
                endif
                fac = sqrt( ( l - ms ) * ( l + ms + 1._db ) )
                td = td + ms * gso(ir,isp)
              else
                mp = l + 1 ! pour que le if en dessous soit faux
              endif
              if( Full_potential ) then
                np = l**2 + l + 1 + mp
              elseif( Hubb_nd .or. Spinorbite ) then
                np = l + 1 + mp
              else
                np = 1
              endif

! le deuxieme indice "m" est la composante non nule a l'origine
              us(ip,n,:,isp,:) = td * us(ir,n,:,isp,:)
     &                         + gm(ir,isp) * us(im,n,:,isp,:)
              if( Ecomp ) us(ip,n,:,isp,:) = us(ip,n,:,isp,:)
     &                                 - img * Eimag * us(ir,n,:,isp,:)

              if( abs(mp) <= l ) us(ip,n,:,isp,:) = us(ip,n,:,isp,:)
     &                           + fac * gso(ir,isp) * us(ir,np,:,isq,:)

              if( Hubb_m ) then
                do mp = -l,l
                  if( Hubb_d .and. m /= mp ) cycle
                  if( Full_potential ) then
                    np = l**2 + l + 1 + mp
                  else
                    np = l + 1 + mp
                  endif
                  us(ip,n,:,isp,:) = us(ip,n,:,isp,:)
     &                           + V_hubb(m,mp,isp) * us(ir,np,:,isp,:)
                end do
              endif

              us(ip,n,:,isp,:) = - us(ip,n,:,isp,:) * gp(ir,isp)

            end do
          end do
        end do
      end do

      if( icheck > 2 ) then
        n = nrmtsd + 1
        allocate( u(n,nlm1,nlm2,nspin,nspino) )
        allocate( ur(n,nlm1,nlm2,nspin,nspino) )
        allocate( rr(n) )
        u(1:n,:,:,:,:) = us(1:n,:,:,:,:)
        ur(1:n,:,:,:,:) = Real( us(1:n,:,:,:,:), db )
        rr(1:n) = r(1:n)
        call write_ur(Full_potential,li,lf,nlm1,nlm2,n,nspin,
     &                nspino,numat,rr,Radial_comp,Rmtsd,u,ur,3)
        deallocate( rr, u, ur )
      endif

      return
      end

!**********************************************************************

! Calcul des integrales radiales dans le cas optique
! Appele par tenseur_car

      subroutine radial_optic(Ecinetic_t,Eimag,Energ,
     &       Full_potential,Green_plus,Hubb_a,Hubb_d,icheck,
     &       ip_max,ip0,lmax,m_depend,m_hubb,ne,nlm1g,nlm2g,
     &       No_diag,nr,nspin,nspino,numat,r,Relativiste,Rmtg,Rmtsd,
     &       roff_ii,roff_ir,roff_ri,roff_rr,Solsing,Spinorbite,V_hubb,
     &       V_intmax,V0bd_t,Vrato_t,Ylm_comp)

      use declarations
      implicit none

      integer:: icheck, ie, ip, ip_max, ip0, isp, l, lf, l_hubbard, 
     &  lfin, lm, lmf, lm1, lmax, lmp, lmpf, lp, lpf, m, 
     &  m_hubb, mf, nlm1, nlm1g, nlm2, nlm2g, mp, mpf, ne, 
     &  nr, nrmtsd, nrmtg, nspin, nspino, numat

      character(len=104):: mot      

      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                       V_hubb
      complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: u, us
      complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: us_t

      logical Ecomp, Full_potential, Green_plus, Hubb_a,
     &  Hubb_d, Hubb_m, m_depend, No_diag, Radial_comp, Relativiste,
     &  Renorm, Solsing, Spinorbite, Ylm_comp

      real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax,
     &                Vecond
      real(kind=db), dimension(nspin):: Ecinetic, V0bd
      real(kind=db), dimension(nr,nspin):: Vrato 
      real(kind=db), dimension(nspin,ne):: Ecinetic_t, V0bd_t
      real(kind=db), dimension(nr,nspin,ne):: Vrato_t 

      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gpi, gp, V
      real(kind=db), dimension(nr,nspin,ne):: g0_t, gm_t, gmi_t, gp_t,
     &                                        gpi_t
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nr,nspino,ne):: gso_t 
      real(kind=db), dimension(nlm1g,nlm1g,nlm2g,nlm2g,nspin**2,
     &        nspino**2,ip0:ip_max):: roff_ii, roff_ir, roff_ri, roff_rr

      real(kind=db), dimension(:,:,:,:,:), allocatable:: ur
      real(kind=db), dimension(:,:,:,:,:,:), allocatable:: ui_t, ur_t

      roff_ii(:,:,:,:,:,:,:) = 0._db
      roff_ir(:,:,:,:,:,:,:) = 0._db
      roff_ri(:,:,:,:,:,:,:) = 0._db
      roff_rr(:,:,:,:,:,:,:) = 0._db

      if( icheck > 1 ) then
        write(3,110) Energ*rydb
        write(3,130) Ecinetic_t(:,:)*rydb
        write(3,140) V0bd_t(:,:)*rydb
      endif

! Terme multiplicatif pour les transitions quadrupolaires
! En S.I. vecond = k = E*alfa_sf*4*pi*epsilon0 / (e*e)
! En ua et rydb : k = 0.5 * alfa_sf * E
      Vecond = 0.5 * alfa_sf * Energ

      do ie = 1,ne

        Vrato(:,:) = Vrato_t(:,:,ie)
        V0bd(:) = V0bd_t(:,ie)

        call mod_V(icheck,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,
     &                 V_intmax,V0bd,Vrato)

        Enervide = Ecinetic_t(1,ie) + V0bd_t(1,ie) 

        call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nr,nspin,nspino,
     &                        numat,r,Relativiste,Spinorbite,V)

        g0_t(:,:,ie) = g0(:,:)
        gm_t(:,:,ie) = gm(:,:)
        gso_t(:,:,ie) = gso(:,:)

        gpi(:,:) = 1 / gp(:,:)
        gpi_t(:,:,ie) = gpi(:,:)

        if( Solsing ) then
          gp_t(:,:,ie) = gp(:,:)
          gmi_t(:,:,ie) = 1 / gm(:,:)
        endif

      end do

      Ecomp = .false.
      if( abs(Eimag) > eps10 ) then
        Ecomp = .true.
      else
        do isp = 1,nspin
          do ie = 1,ne
            if( Ecinetic_t(isp,ie) > eps10 ) cycle
            Ecomp = .true.
            exit
          end do
        end do
      endif
      Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

      if( Full_potential ) then
        lfin = 0
      else
        lfin = lmax
      endif
      Renorm = .true.

      allocate( ur_t(nr,nlm1g,nlm2g,nspin,nspino,ne) )
      allocate( ui_t(nr,nlm1g,nlm2g,nspin,nspino,ne) )

      do l = 0,lfin

        if( Hubb_a .and. l == l_hubbard( numat ) )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif  

        if( Full_potential ) then
          nlm1 = ( lmax + 1 )**2
          nlm2 = nlm1
        elseif( Hubb_m .and. .not. Hubb_d ) then
          nlm1 = 2*l + 1
          nlm2 = nlm1
        elseif( Spinorbite .or. Hubb_m ) then
          nlm1 = 2*l + 1
          nlm2 = 1
        else
          nlm1 = 1
          nlm2 = 1
        endif

        allocate( u(nr,nlm1,nlm2,nspin,nspino) )
        allocate( ur(nr,nlm1,nlm2,nspin,nspino) )
        allocate( Tau(nlm1,nspin,nlm1,nspin) )
        allocate( us(nrmtsd+1,nlm1,nlm2,nspin,nspino) )

        if( l == 0 ) allocate( us_t(nrmtsd+1,nlm1g,nlm2g,nspin,nspino,
     &                                                             ne) )

        do ie = 1,ne

          Ecinetic(:) = Ecinetic_t(:,ie)
          konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )
          g0(:,:) = g0_t(:,:,ie) 
          gm(:,:) = gm_t(:,:,ie) 
          gso(:,:) = gso_t(:,:,ie) 
          gpi(:,:) = gpi_t(:,:,ie) 

          call Sch_radial(Ecinetic,Ecomp,Eimag,f2,
     &          Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,
     &          l,lmax,m_hubb,nlm1,nlm2,nr,nrmtg,nspin,nspino,numat,r,
     &          Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,u,
     &          ur,V_hubb)

! Calcul de la solution singuliere
! gm et gp inverse dans le sousprogramme
          if( Solsing ) then
            gp(:,:) = gp_t(:,:,ie) 
            gmi(:,:) = gmi_t(:,:,ie) 
            call Sch_radial_solsing(Ecomp,Eimag,f2,
     &         Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde,
     &         l,lmax,m_hubb,nlm1,nlm2,nr,nrmtsd,nspin,nspino,numat,r,
     &         Radial_comp,Rmtsd,Spinorbite,Tau,us,V_hubb)
          endif

          if( Radial_comp ) then
            if( Full_potential ) then
              ur_t(:,:,:,:,:,ie) = real( u(:,:,:,:,:), db )
              ui_t(:,:,:,:,:,ie) = aimag( u(:,:,:,:,:) )
              us_t(:,:,:,:,:,ie) = u(:,:,:,:,:)
            elseif( Hubb_m .and. .not. Hubb_d ) then
              ur_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) 
     &                           = real( u(:,1:nlm1,1:nlm1,:,:), db )
              ui_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) 
     &                           = aimag( u(:,1:nlm1,1:nlm1,:,:) )
              us_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie)
     &                           = u(:,1:nlm1,1:nlm1,:,:)
            elseif( No_diag ) then
              do lm1 = l**2+1,(l+1)**2
                lm = min( lm1 - l**2, nlm1 )
                ur_t(:,lm1,lm1,:,:,ie) = real( u(:,lm,1,:,:), db )
                ui_t(:,lm1,lm1,:,:,ie) = aimag( u(:,lm,1,:,:) )
                us_t(:,lm1,lm1,:,:,ie) = us(:,lm,1,:,:)
              end do
            elseif( Spinorbite .or. Hubb_m ) then
               do lm1 = l**2+1,(l+1)**2
                 lm = lm1 - l**2
                 ur_t(:,lm1,1,:,:,ie) = real( u(:,lm,1,:,:), db )
                 ui_t(:,lm1,1,:,:,ie) = aimag( u(:,lm,1,:,:) )
                 us_t(:,lm1,1,:,:,ie) = us(:,lm,1,:,:)
               end do
            elseif( m_depend ) then
              do lm1 = l**2+1,(l+1)**2
                ur_t(:,lm1,1,:,:,ie) = real( u(:,1,1,:,:), db )
                ui_t(:,lm1,1,:,:,ie) = aimag( u(:,1,1,:,:) )
                us_t(:,lm1,1,:,:,ie) = us(:,1,1,:,:)
              end do
            else
              ur_t(:,l,1,:,:,ie) = real( u(:,1,1,:,:), db )
              ui_t(:,l,1,:,:,ie) = aimag( u(:,1,1,:,:) )
              us_t(:,l,1,:,:,ie) =  u(:,1,1,:,:)
            endif
          else
            if( Full_potential ) then
              ur_t(:,:,:,:,:,ie) = ur(:,:,:,:,:)
            elseif( Hubb_m .and. .not. Hubb_d ) then
              ur_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) 
     &                           = ur(:,1:nlm1,1:nlm1,:,:)
            elseif( No_diag ) then
              do lm1 = l**2+1,(l+1)**2
                lm = min( lm1 - l**2, nlm1 )
                ur_t(:,lm1,lm1,:,:,ie) = ur(:,lm,1,:,:)
              end do
            elseif( Spinorbite .or. Hubb_m ) then
               do lm1 = l**2+1,(l+1)**2
                 lm = lm1 - l**2
                 ur_t(:,lm1,1,:,:,ie) = ur(:,lm,1,:,:)
               end do
            elseif( m_depend ) then
              do lm1 = l**2+1,(l+1)**2
                ur_t(:,lm1,1,:,:,ie) = ur(:,1,1,:,:)
              end do
            else
              ur_t(:,l+1,1,:,:,ie) = ur(:,1,1,:,:)
            endif
          endif

        end do

        deallocate( Tau, u, ur, us )

      end do   ! fin de la boucle sur l

! Integrale radiale pour la regle d'or de Fermi

      call radial_matrix_optic(Green_plus,ip_max,ip0,
     &         ne,nlm1g,nlm2g,nr,nrmtsd,nspin,
     &         nspino,r,Radial_comp,Rmtsd,roff_ii,roff_ir,roff_ri,
     &         roff_rr,ui_t,ur_t,Vecond)

      deallocate( ur_t, ui_t, us_t )

      if( icheck > 1 ) then
        mot = ' '
        if( Ecomp .and. nspino == 2 ) then
          mot(11:17) = 'up sol1'
          mot(37:43) = 'up sol2'
          mot(63:69) = 'dn sol1'
          mot(89:95) = 'dn sol2'
        elseif( nspino == 2 ) then
          mot(5:50) = 'up sol1      up sol2      dn sol1      dn sol2'
        elseif( nspin == 2 ) then
          mot = '      up           dn'
        else
          mot = ' '
        endif

        do ip = ip0,ip_max
          select case(ip)
            case(0)
              write(3,160) '  monopole'
            case(1)
              write(3,160) '    dipole'
            case(2)
              write(3,160) 'quadrupole'
            case(3)
              write(3,160) '  octupole'
          end select
          if( Full_potential .or. Hubb_a ) then
            write(3,170) '  l  m lp mp', mot
          else
            write(3,175) '  l  m', mot
          endif

          do l = 0,lmax
            do m = -l,l
              if( m_depend ) then
                lm = l**2 + l + 1 + m
              else
                if( m /= 0 ) cycle
                lm = l + 1
              endif
              do lp = 0,lmax
                if( .not. Full_potential .and. l /= lp ) cycle
                do mp = -lp,lp
                  if( .not. ( Hubb_a .or. Full_potential
     &                               .or. m == mp ) ) cycle
                  lmp = min( lp**2 + lp + 1 + mp, nlm2g )

                  do lf = 0,lmax
                    do mf = -lf,lf
                      if( m_depend ) then 
                        lmf = lf**2 + lf + 1 + mf
                      else
                        if( mf /= 0 ) cycle
                        lmf = lf + 1
                      endif
                      lmpf = 0
                      do lpf = 0,lmax
                        if( .not. Full_potential .and. lf /= lpf ) cycle
                        do mpf = -lpf,lpf
                          if( .not. ( Hubb_a .or. Full_potential
     &                                       .or. mf == mpf ) ) cycle
                          lmpf = lpf**2 + lpf + 1 + mpf
                          lmpf = min(nlm2g,lmpf)

                          if( Full_potential .or. Hubb_a ) then
                            if( Ecomp) then
                              write(3,180) l, m, lp, mp, lf, mf, lpf,
     &                                     mpf,
     &                            (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip),
     &                             roff_ri(lm,lmf,lmp,lmpf,isp,:,ip),
     &                             roff_ir(lm,lmf,lmp,lmpf,isp,:,ip),
     &                             roff_ii(lm,lmf,lmp,lmpf,isp,:,ip),
     &                            isp = 1,nspin**2)
                            else
                              write(3,185) l, m, lp, mp, lf, mf, lpf,
     &                                     mpf,
     &                          (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip),
     &                           isp = 1,nspin**2)
                            endif
                          else
                            if( Ecomp) then
                              write(3,180) l, m, lf, mf,
     &                            (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip),
     &                             roff_ri(lm,lmf,lmp,lmpf,isp,:,ip),
     &                             roff_ir(lm,lmf,lmp,lmpf,isp,:,ip),
     &                             roff_ii(lm,lmf,lmp,lmpf,isp,:,ip),
     &                            isp = 1,nspin**2)
                            else
                              write(3,185) l, m, lf, mf,
     &                          (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip),
     &                            isp = 1,nspin**2)
                            endif
                          endif
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
  110 format(/' Photon energy =',f10.3,' eV')
  130 format(' Ecinetic =',4f10.3)
  140 format(' V0bd     =',4f10.3)
  160 format(/' Radial matrix, ',a10,1x,'term:')
  170 format(a12,a104)
  175 format(a6,a104)
  185 format(4i3,1p,16e13.5)
  180 format(4i3,1p,4e13.5)
      end

!***********************************************************************

! Calcul des fonctions de base sur les points du maillage

      subroutine cal_phiato(Base_ortho,dcosxyz,
     &           Full_potential,ia,iang,ibord,
     &           icheck,iopsymr,ll,lmax,nlm1,nlm2,natome,
     &           nbord,nbtm,nlmagm,nlmmax,nphiato1,
     &           nphiato7,npsom,nr,nspin,nspino,phiato,posi,r,
     &           Radial_comp,Spinorbite,u,ur,xyz,Ylm_comp,Ylmato)

      use declarations
      implicit none

      integer:: i, ia, iang, ib, icheck, isol, isp, isym, jsol, jsp, l, 
     &  ll, lm, lm0, lmax, lmv, lp, m, nlm1, nlm2, mp, mpp, mv, n,
     &  natome, nbord, nbtm, nlmagm, nlmmax, np,  
     &  nphiato1, nphiato7, npsom, nr, nspin, nspino
      integer, dimension(nopsm):: iopsymr
      integer, dimension(nbtm,natome):: ibord

      complex(kind=db):: phic, phid, uc0, ucm, ucp, yc
      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u

      logical:: Ylm_comp, Base_ortho, Full_potential,
     &          Radial_comp, Spinorbite

      real(kind=db):: f_interp2, phi, phii, r0, rm, rp, rrel, u0, um, up
      real(kind=db), dimension(3):: dcosxyz, ps, x, w
      real(kind=db), dimension(nr):: r
      real(kind=db), dimension(3,natome):: posi
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur 
      real(kind=db), dimension(4,npsom):: xyz
      real(kind=db), dimension(nbtm,nlmmax,natome):: Ylmato
      real(kind=db), dimension(nphiato1,nlmagm,nspin,nspino,natome,
     &                                             nphiato7):: phiato

      ps(1:3) = posi(1:3,ia)

      do ib = 1,nbord

        i = ibord(ib,ia)
        x(1:3) = xyz(1:3,i)

        call posrel(Base_ortho,dcosxyz,iopsymr,x,ps,w,rrel,isym)

        do i = 2,nr
          if( r(i) > rrel ) exit
        end do
        i = min(i,nr)

        rm = r(i-2)
        r0 = r(i-1)
        rp = r(i)

        n = 0
        do l = 0,lmax
          if( .not. Full_potential .and. l /= ll ) cycle 
          lm0 = l**2 + l + 1
          do m = -l,l
            if( nlm1 == 1 .and. m /= 0 ) cycle
            n = n + 1
            do isp = 1,nspin
              if( iang == 1 ) then
                jsp = isp
              elseif( iang == - 1 ) then
                jsp = nspin + 1 - isp
              endif

              np = 0
              do lp = 0,lmax
                if( .not. Full_potential .and. lp /= l ) cycle
                do mp = -lp,lp
                  if( nlm2 == 1 .and. mp /= m ) cycle
                  np = np + 1
                  do isol = 1,nspino
                    if( iang == 1 ) then
                      jsol = isol
                    elseif( iang == - 1 ) then
                      jsol = nspino + 1 - isol
                    endif

                    if( nlm2 > 1 ) then
                      lmv = lp**2 + lp + 1 + mp
                    elseif( Spinorbite ) then
                      mv = m + isol - isp  ! m  faux
                      if( mv > l .or. mv < -l ) cycle
                      lmv = lm0 + m
                    else
                      lmv = lm0 + m
                    endif
                    mv = m

                    if( Radial_comp ) then
                      if( iang == 0 ) then
                        if( nspino == 2 .and. isp /= isol ) cycle
                        ucm = 0.5_db * ( u(i-2,n,np,1,1)
     &                        + u(i-2,n,np,min(2,nspin),min(2,nspino)) )
                        uc0 = 0.5_db * ( u(i-1,n,np,1,1)
     &                        + u(i-1,n,np,min(2,nspin),min(2,nspino)) )
                        ucp = 0.5_db * ( u(i,n,np,1,1)
     &                        + u(i,n,np,min(2,nspin),min(2,nspino)) )
                      else
                        ucm = u(i-2,n,np,jsp,jsol)
                        uc0 = u(i-1,n,np,jsp,jsol)
                        ucp = u(i,n,np,jsp,jsol)
                      endif
                    else
                      if( iang == 0 ) then
                        if( nspino == 2 .and. isp /= isol ) cycle
                        um = 0.5_db * ( ur(i-2,n,np,1,1)
     &                       + ur(i-2,n,np,min(2,nspin),min(2,nspino)) )
                        u0 = 0.5_db * ( ur(i-1,n,np,1,1)
     &                       + ur(i-1,n,np,min(2,nspin),min(2,nspino)) )
                        up = 0.5_db * ( ur(i,n,np,1,1)
     &                       + ur(i,n,np,min(2,nspin),min(2,nspino)) )
                      else
                        um = ur(i-2,n,np,jsp,jsol)
                        u0 = ur(i-1,n,np,jsp,jsol)
                        up = ur(i,n,np,jsp,jsol)
                      endif
                    endif


                    if( Radial_comp ) then
                      um = real( ucm, db )
                      u0 = real( uc0, db )
                      up = real( ucp, db )
                      phi = f_interp2(rrel,rm,r0,rp,um,u0,up)
                      um = aimag( ucm )
                      u0 = aimag( uc0 )
                      up = aimag( ucp )
                      phii = f_interp2(rrel,rm,r0,rp,um,u0,up)
                    else
                      phi = f_interp2(rrel,rm,r0,rp,um,u0,up)
                      phii = 0._db
                    endif
                    phic = cmplx( phi, phii, db )

                    do mpp = -l,l
                      if( nlm1 /= 1 .and. mpp /= m ) cycle
                      lm = lm0 + mpp
                      if( nlm1 == 1 ) then
                        lmv = lm
                        mv = mpp
                      endif 
                      if( mpp == 0 .or. .not. Ylm_comp ) then
                        phiato(ib,lmv,isp,isol,ia,1)
     &                      = phiato(ib,lmv,isp,isol,ia,1)
     &                      + phi * Ylmato(ib,lm,ia)
                        if( Radial_comp ) phiato(ib,lmv,isp,isol,ia,2)
     &                            = phiato(ib,lmv,isp,isol,ia,2)
     &                            + phii * Ylmato(ib,lm,ia)
                      else
                        phid = phic
     &                   * yc(mv,Ylmato(ib,lm,ia),Ylmato(ib,lm-2*mv,ia))
                        phiato(ib,lmv,isp,isol,ia,1)
     &                  = phiato(ib,lmv,isp,isol,ia,1) + Real( phid, db)
                        phiato(ib,lmv,isp,isol,ia,2)
     &                  = phiato(ib,lmv,isp,isol,ia,2) +  aimag( phid )
                      endif

                    end do
                  end do
                end do
              end do

            end do
          end do
        end do
      end do   

      if( icheck > 2 .and. ( ll == lmax .or. Full_potential ) ) then
        write(3,110) ia
        do l = 0,lmax
          lm0 = l**2 + l + 1
          do mp = -l,l
            if( nphiato7 == 1 ) then
              if( Spinorbite ) then
                write(3,120) l, mp
              else
                write(3,130) l, mp
              endif
            else
              if( Spinorbite ) then
                write(3,140) l, mp
              else
                write(3,150) l, mp
              endif
            endif
            lm = lm0 + mp
            do ib = 1,nbord
              write(3,160) ibord(ib,ia),
     &          ( ( phiato(ib,lm,isp,isol,ia,:), isp = 1,nspin ),
     &                                           isol = 1,nspino )
            end do
          end do
        end do
      endif

      return
  110 format(/' ia =',i3)
  120 format(/' ibord phiato(u,1)  phiato(d,1)  phiato(u,2)  ',
     &        'phiato(d,2)   l =',i3,', m =',i3)
  130 format(/' ibord  phiato(u)    phiato(d)   l =',i3,', m =',i3)
  140 format(/' ibord',6x,' phiato(u,1)',13x,'  phiato(d,1)',13x,
     &        '  phiato(u,2)  ',13x,'phiato(d,2)   l =',i3,', m =',i3)
  150 format(/' ibord',6x,'  phiato(u)',13x,'    phiato(d)   l =',i3,
     &                     ', m =',i3)
  160 format(i5,1p,8e13.4)
      end

!***********************************************************************

! Calcul de l'incrementation de la densite d'etat

      subroutine cal_dens(Ylm_comp,Cal_xanes,drho_self,Ecinetic,
     &          Eimag,Energ,Enervide,Full_atom,Full_potential,
     &          Hubb,Hubb_diag,iaabsi,iaprabs,iaprotoi,icheck,itypei,
     &          itypepr,lla2_state,lla2_state2,lmaxat,m_hubb,mpinodes,
     &          mpirank,n_atom_0,n_atom_0_self,n_atom_ind,
     &          n_atom_ind_self,n_atom_proto,natome,nlmagm,nrato,nrm,
     &          nrm_self,nspin,nspino,ntype,numat,rato,Relativiste,Rmtg,
     &          Rmtsd,Solsing,Solsing_only,Spinorbite,State_all,
     &          Statedens,Taull,V_hubb,V_hubb_abs,V_intmax,V0bd,Vrato)

      use declarations  
      implicit none
 
      integer ia, iaabsi, iapr, iaprabs, icheck, ipr, ir, isp, 
     &  isp1, isp2, it, lla2_state, lla2_state2, lm, lmax, 
     &  m_hubb, mpinodes, mpirank, n_atom_0, 
     &  n_atom_0_self, n_atom_ind, n_atom_ind_self, 
     &  n_atom_proto, natome, nlma, nlma2, nlmagm,
     &  nr, nrs, nrm, nrm_self, nspin, nspino, ntype, Z

      complex(kind=db), dimension(nlmagm,nspin,nlmagm,nspin,natome)::
     &                                                          taull
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                      V_hubb_abs
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin,
     &                        n_atom_0_self:n_atom_ind_self):: V_hubb
      complex(kind=db), dimension(:,:,:), allocatable:: V_hubb_t
      complex(kind=db), dimension(:,:,:,:), allocatable:: Taulla

      integer, dimension(0:ntype):: nrato, numat
      integer, dimension(0:n_atom_proto):: itypepr, lmaxat
      integer, dimension(natome):: iaprotoi, itypei

      logical:: Absorbeur, Cal_xanes, Full_atom, Full_potential, Hubb_a,
     &  Hubb_d, Relativiste, Self, Solsing, Solsing_only,
     &  Spinorbite, State_all, Ylm_comp
      logical, dimension(0:ntype):: Hubb
      logical, dimension(n_atom_0:n_atom_ind):: iapr_done 
      logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag 
 
      real(kind=db):: Eimag, Energ, Enervide, V_intmax 
      real(kind=db), dimension(nspin):: Ecinetic, V0bd 
      real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd 
      real(kind=db), dimension(lla2_state,nspin,lla2_state,nspin,
     &                   n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens
      real(kind=db), dimension(0:nrm,0:ntype):: rato
      real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind)::Vrato 
      real(kind=db), dimension(nrm_self,lla2_state2,nspin,
     &          n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
      real(kind=db), dimension(:,:), allocatable:: Vrato_t
      real(kind=db), dimension(:), allocatable:: r
      real(kind=db), dimension(:,:,:), allocatable:: drho
      real(kind=db), dimension(:,:,:,:), allocatable:: State, State_i

      Self = .not. Cal_xanes
     
      if( Self ) drho_self(:,:,:,:,mpirank) = 0._db
 
      Statedens(:,:,:,:,:,mpirank) = 0._db
      iapr_done(:) = .false.

      boucle_ia: do ia = 1,natome  

        Absorbeur = ( ia == iaabsi .and. Cal_xanes )
        if( Cal_xanes .and. .not. ( State_all .or. Absorbeur ) ) cycle
  
        ipr = iaprotoi(ia)
        it = itypepr(ipr) 
        Hubb_a = Hubb(it)

        Z = numat( it )
        if( Z < 19 ) then
          lmax = min(2,lmaxat(ipr))
        elseif( Z > 18 .and. Z < 55 ) then
          lmax = min(3,lmaxat(ipr))
        else
          lmax = min(4,lmaxat(ipr))
        endif
        nlma = (lmax + 1)**2
        if( Full_potential ) then
          nlma2 = nlma
        else
          nlma2 = 1
        endif

        if( Full_atom ) then
          iapr = ia
        else
          iapr = ipr
          if( iapr_done(iapr) ) cycle 
        endif
        it = itypei(ia)
               
        do ir = 1,nrato(it)
          if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
        end do
        nr = ir + 1

        if( Self ) then
          nrs = nr
        else
          nrs = 0
        endif
        allocate( r(nr) )
        allocate( drho(nrs,nlma2,nspin) )
        allocate( State(nlma,nspin,nlma,nspin))
        allocate( State_i(nlma,nspin,nlma,nspin))
        allocate( Vrato_t(nr,nspin) )
        allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin) )
        allocate( Taulla(nlma,nspin,nlma,nspin) )
        State(:,:,:,:) = 0._db
        State_i(:,:,:,:) = 0._db
        drho(:,:,:) = 0._db
        r(1:nr) = rato(1:nr,it)
        Vrato_t(1:nr,:) = Vrato(1:nr,:,iapr)
        if( Hubb(it) .and. iapr <= n_atom_ind_self ) then
          Hubb_a = .true.
          if( Absorbeur ) then
            V_Hubb_t(:,:,:) = V_Hubb_abs(:,:,:)
            Hubb_d = Hubb_diag(iaprabs)
          else
            V_Hubb_t(:,:,:) = V_Hubb(:,:,:,iapr)
            Hubb_d = Hubb_diag(iapr)
          endif
        else
          Hubb_a = .false.
          Hubb_d = .true.
        end if

        do isp1 = 1,nspin
          do lm = 1,nlma
            do isp2 = 1,nspin
              Taulla(1:nlma,isp1,lm,isp2)
     &                             = Taull(1:nlma,isp1,lm,isp2,ia) 
            end do
          end do
        end do

        call radial_sd(drho,Ecinetic,Eimag,Energ,Enervide,
     &       Full_potential,Hubb_a,Hubb_d,ia,icheck,lmax,m_hubb,nlma,
     &       nlma2,nr,nrs,nspin,nspino,Z,r,Relativiste,Rmtg(ipr),
     &       Rmtsd(ipr),Self,Solsing,Solsing_only,Spinorbite,State,
     &       State_i,Taulla,V_hubb_t,V_intmax,V0bd,Vrato_t,Ylm_comp)

        deallocate( r, Taulla, V_hubb_t, Vrato_t )


        do isp1 = 1,nspin
          do lm = 1,nlma
            do isp2 = 1,nspin
              Statedens(1:nlma,isp1,lm,isp2,iapr,mpirank)
     &                          = State(1:nlma,isp1,lm,isp2)
            end do
          end do
        end do

        if( Self ) then
          do isp = 1,nspin
            do lm = 1,nlma2
              drho_self(1:nr,lm,isp,iapr,mpirank) = drho(1:nr,lm,isp)
            end do
          end do
        endif
 
        deallocate( State, State_i, drho )
        
        iapr_done(iapr) = .true.

      end do boucle_ia
         
      return
      end

!**********************************************************************

! Calcul des integrales radiales et de la solution singuliere pour la
! densite d'etat.
! Appele par cal_data

      subroutine radial_sd(drho,Ecinetic,Eimag,Energ,Enervide,
     &       Full_potential,Hubb_a,Hubb_d,ia,icheck,lmax,m_hubb,nlma,
     &       nlma2,nr,nrs,nspin,nspino,numat,r,Relativiste,Rmtg,
     &       Rmtsd,Self,Solsing,Solsing_only,Spinorbite,State,State_i,
     &       Taull,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)
 
      use declarations
      implicit none

      integer:: ia, icheck, l, l_hubbard, lfin,
     &  lmax, m_hubb, nlm1, nlm2, nlma, nlma2, nr, nrs, 
     &  nrmtg, nrmtsd, nspin, nspino, numat
      
      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(nlma,nspin,nlma,nspin):: Taull
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                          V_hubb
      complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: u

      logical Diagonale, Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, 
     &  Radial_comp, Relativiste, Renorm, Self, Solsing, Solsing_only,
     &  Spinorbite, Ylm_comp

      real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
      real(kind=db), dimension(nspin):: Ecinetic, V0bd
      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gp, gpi, V,
     &                                     Vrato
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nrs,nlma2,nspin):: drho
      real(kind=db), dimension(nlma,nspin,nlma,nspin):: State, State_i
      real(kind=db), dimension(:,:,:,:,:), allocatable:: ur

      konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

      if( icheck > 1 ) then
        write(3,110)
        write(3,120) ia, numat, Rmtg*bohr, Rmtsd*bohr, Energ*rydb
        write(3,130) Ecinetic(:)*rydb
        write(3,140) V0bd(:)*rydb
        write(3,150) konde(:)
      endif

      call mod_V(icheck,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,
     &                 V_intmax,V0bd,Vrato)

      call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nr,nspin,nspino,
     &                        numat,r,Relativiste,Spinorbite,V)

      gpi(:,:) = 1 / gp(:,:)
      if( Solsing ) gmi(:,:) = 1 / gm(:,:)

      if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10
     &     .or. Ecinetic(nspin) < eps10 ) then
        Ecomp = .true.
      else
        Ecomp = .false.
      endif
      Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

      if( Full_potential ) then
        lfin = 0
      else
        lfin = lmax
      endif

      Renorm = .true.

!      if( Ylm_comp .and. Hubb_a ) 
!     &    call Trans_Tau(icheck,0,lmax,nspin,Spinorbite,Taull,.true.)

      do l = 0,lfin

        if( Hubb_a .and. l == l_hubbard( numat ) )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif  

        if( Full_potential ) then
          nlm1 = ( lmax + 1 )**2
          nlm2 = nlm1
        elseif( Hubb_m .and. .not. Hubb_d ) then
          nlm1 = 2*l + 1
          nlm2 = nlm1
        elseif( Spinorbite .or. Hubb_m ) then
          nlm1 = 2*l + 1
          nlm2 = 1
        else
          nlm1 = 1
          nlm2 = 1
        endif
        Diagonale = .not. ( Full_potential .or. Hubb_m )

        allocate( u(nr,nlm1,nlm2,nspin,nspino) )
        allocate( ur(nr,nlm1,nlm2,nspin,nspino) )
        allocate( Tau(nlm1,nspin,nlm1,nspin) )

        call Sch_radial(Ecinetic,Ecomp,Eimag,f2,
     &          Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,
     &           l,lmax,m_hubb,nlm1,nlm2,nr,nrmtg,nspin,nspino,numat,r,
     &           Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,u,
     &           ur,V_hubb)

! Integrale radiale pour la densite d'etat
        if( .not. Solsing_only )
     &    call Radial_matrix_sd(Diagonale,drho,Full_potential,icheck,
     &         l,lmax,nlm1,nlm2,nlma,nlma2,nr,nrs,nrmtsd,nspin,nspino,
     &         r,Radial_comp,Rmtsd,State,State_i,Self,Taull,u,ur,
     &         Ylm_comp)

! Calcul de la solution singuliere
        if( Solsing ) call Cal_Solsing_sd(drho,Ecomp,Eimag,f2,
     &        Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde,l,
     &        lmax,m_hubb,nlm1,nlm2,nlma,nlma2,nr,nrmtsd,nrs,nspin,
     &        nspino,numat,r,Radial_comp,Rmtsd,Self,
     &        Spinorbite,State,Tau,u,ur,V_hubb,Ylm_comp)

        deallocate( Tau )
        deallocate( u, ur )

      end do   ! fin de la boucle sur l

      return
  110 format(/' ---- Radial_sd ----',100('-'))
  120 format(/'     Atom =',i3,/'        Z =',i3,/
     &'     Rmtg =',f10.5,' A,  Rmtsd =',f10.5,' A',/
     &'    Energ =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format('    V0bd  =',2f10.3)
  150 format('    konde =',2f12.5)
      end

!**********************************************************************

! Calcul de l'integrale radiale pour la densite d'etat
! On projete sur la base reelle

      subroutine Radial_matrix_sd(Diagonale,drho,Full_potential,icheck,
     &         l,lmax,nlm1,nlm2,nlma,nlma2,nr,nrs,nrmtsd,nspin,nspino,
     &         r,Radial_comp,Rmtsd,State,State_i,Self,Taull,u,ur,
     &         Ylm_comp)

      use declarations
      implicit none

      integer:: icheck, ir, is1, is2, isg1, isg2, iso, iso1, iso2, isp,
     &  isp1, isp2, l, l1, l2, l3, lm, lm0, lm1, lm2, lma1, lma2, lmax,
     &  lp1, lp2, m1, m2, m3, mp1, mp2, mr1, mr2, mv,
     &  n1, n2, nlm1, nlm2, nlma, nlma2, np1, np2, nr, nr1, nr2, nrmtsd,   
     &  nrs, nspin, nspino

      complex(kind=db):: c_harm, c_harm1, c_harm2, rof_sd, Sta_e, Sta_s
      complex(kind=db), dimension(nspin):: rof_sd_l
      complex(kind=db), dimension(nrmtsd):: fc
      complex(kind=db), dimension(nrmtsd,nspin):: fc_l
      complex(kind=db), dimension(nlma,nspin,nlma,nspin):: Taull
      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u

      logical Diagonale, Full_potential, Radial_comp, Self, Ylm_comp

      real(kind=db):: f_integr3, g, gauntc, gauntcp, rac2, radlr, radli,
     &  Rmtsd, Sta_i, Sta_r
      real(kind=db), dimension(nr):: r 
      real(kind=db), dimension(nrmtsd):: r2, rr, fcr, fci, ri2
      real(kind=db), dimension(nrs,nlma2,nspin):: drho
      real(kind=db), dimension(nlma,nspin,nlma,nspin):: State, State_i
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur

      rac2 = 1 / sqrt(2._db )

      if( Full_potential ) then
        lm0 = 0
      else
        lm0 = l**2
      endif

      rr(1:nrmtsd) = r(1:nrmtsd)
      r2(1:nrmtsd) = r(1:nrmtsd)**2
      ri2(1:nrmtsd) = 1 / ( quatre_pi * r2(1:nrmtsd) )

      if( nlm1 == 1 ) then
        do isp = 1,nspin
          iso = min( isp, nspino )
          if( Radial_comp ) then
            fcr(1:nrmtsd) = r2(1:nrmtsd) * real( 
     &                u(1:nrmtsd,1,1,isp,iso)
     &              * u(1:nrmtsd,1,1,isp,iso), db )
            radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
            fci(1:nrmtsd) = r2(1:nrmtsd) * aimag( 
     &                u(1:nrmtsd,1,1,isp,iso)
     &              * u(1:nrmtsd,1,1,isp,iso) )
            radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)
          else
            fcr(1:nrmtsd) = r2(1:nrmtsd) 
     &                    * ur(1:nrmtsd,1,1,isp,iso)
     &                    * ur(1:nrmtsd,1,1,isp,iso)
            radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
            fci(:) = 0._db
            radli = 0._db
          endif
          rof_sd_l(isp) = cmplx( radlr, radli, db )
          if( Self ) fc_l(1:nrmtsd,isp) = ri2(1:nrmtsd)
     &                     * cmplx( fcr(1:nrmtsd), fci(1:nrmtsd), db )  

        end do
      endif

      if( icheck > 2 ) then
        if( nlm2 > 1 ) then
          write(3,110)
        elseif( nspino ==2 ) then
          write(3,120)
        else
          write(3,130)
        endif
      endif

      do isp1 = 1,nspin
        nr1 = 0  
        do l1 = 0,lmax
          if( .not. Full_potential .and. l /= l1 ) cycle
          do mr1 = -l1,l1
            nr1 = nr1 + 1
            lm1 = lm0 + nr1

            do isp2 = 1,nspin  
              if( nspino == 1 .and. isp1 /= isp2 ) cycle
              nr2 = 0  
              do l2 = 0,lmax
                if( .not. Full_potential .and. l /= l2 ) cycle
                do mr2 = -l2,l2
                  nr2 = nr2 + 1

                  if( Diagonale .and. mr1 /= mr2 ) cycle

                  lm2 = lm0 + nr2

                  do isg1 = -1,1,2
                    if( Ylm_comp ) then
                      if( mr1 < 0 ) then
                        if( isg1 == 1 ) then
                          c_harm1 = cmplx( 0._db, (-1)**mr1 * rac2, db )
                        else
                          c_harm1 = cmplx( 0._db, - rac2, db )
                        endif
                      elseif( mr1 == 0 ) then
                        if( isg1 == -1 ) cycle
                        c_harm1 = ( 1._db, 0._db )
                      else
                        if( isg1 == 1 ) then
                          c_harm1 = cmplx( rac2, 0._db, db )
                        else
                          c_harm1 = cmplx( (-1)**mr1 * rac2, 0._db, db )
                        endif
                      endif
                      m1 = isg1*mr1
                    else
                      if( isg1 == -1 ) cycle
                      m1 = mr1
                      c_harm1 = ( 1._db, 0._db )
                    endif
                    if( isg1 == 1 ) then
                      n1 = nr1
                    else
                      n1 = nr1 - 2 * mr1
                    endif           

                    do isg2 = -1,1,2
                      if( Ylm_comp ) then
                        if( mr2 < 0 ) then
                          if( isg2 == 1 ) then
                            c_harm2 = cmplx( 0._db, (-1)**mr2 * rac2,db)
                          else
                            c_harm2 = cmplx( 0._db, - rac2, db )
                          endif
                        elseif( mr2 == 0 ) then
                          if( isg2 == -1 ) cycle
                          c_harm2 = ( 1._db, 0._db )
                        else
                          if( isg2 == 1 ) then
                            c_harm2 = cmplx( rac2, 0._db, db )
                          else
                            c_harm2 = cmplx( (-1)**mr2 * rac2, 0._db,db)
                          endif
                        endif
                        m2 = isg2*mr2           
                      else
                        if( isg2 == -1 ) cycle
                        m2 = mr2
                        c_harm2 = ( 1._db, 0._db )
                      endif
                      if( isg2 == 1 ) then
                        n2 = nr2
                      else
                        n2 = nr2 - 2 * mr2
                      endif           

                      c_harm = c_harm1 * conjg( c_harm2 )

                      np1 = 0
                      do lp1 = 0,lmax
                        if( .not. Full_potential .and. lp1 /= l1 ) cycle

                        do mp1 = -lp1, lp1
                          if( nlm2 == 1 .and. mp1 /= m1 ) cycle
                          if( nlm2 == 1 ) then
                            np1 = 1
                          else
                            np1 = np1 + 1
                          endif

                          do iso1 = 1,nspino
!       if( iso1 /= isp1 ) cycle

                            if( nspino == 2 .and. nlm2 == 1 ) then
                              mv = mp1 + iso1 - isp1
                              if( mv > lp1 .or. mv < - lp1 ) cycle
                            else
                              mv = mp1
                            endif

                            if( nlma2 /= 1 ) then
                              is1 = iso1
                              lma1 = lp1**2 + lp1 + 1 + mp1
                            elseif( nspino == 2 ) then
                              is1 = iso1
                              lma1 = l1**2 + l1 + 1 + mv
                            else
                              is1 = isp1
                              lma1 = l1**2 + l1 + 1 + mv
                            endif

                            np2 = 0
                            do lp2 = 0,lmax
                              if( .not. Full_potential .and. lp2 /= l2 )
     &                                                             cycle

                              do mp2 = -lp2,lp2
                                if( nlm2 == 1 .and. mp2 /= m2 ) cycle
                                if( nlm2 == 1 ) then
                                  np2 = 1
                                else
                                  np2 = np2 + 1
                                endif

                                do iso2 = 1,nspino
!       if( iso2 /= isp2 ) cycle
                                  if( nspino == 2 .and. nlm2 == 1 ) then
                                    mv = mp2 + iso2 - isp2
                                    if( mv > lp2 .or. mv < - lp2 ) cycle
                                  else
                                    mv = mp2
                                  endif

                                  if( nlma2 /= 1 ) then
                                    is2 = iso2
                                    lma2 = lp2**2 + lp2 + 1 + mp2
                                  elseif( nspino == 2 ) then
                                    is2 = iso2
                                    lma2 = l2**2 + l2 + 1 + mv
                                  else
                                    is2 = isp2
                                    lma2 = l2**2 + l2 + 1 + mv
                                  endif

                                  if( nlm1 == 1 ) then

                                    rof_sd = rof_sd_l(isp1)
                                    if( Self .and. isp1 == isp2 .and.
     &                                lm1 == lm2 ) fc(:) = fc_l(:,isp1)

                                  else

                                    if( Radial_comp ) then
                                      fcr(1:nrmtsd) = r2(1:nrmtsd)
     &                                  * real( 
     &                                  u(1:nrmtsd,n1,np1,isp1,iso1)
     &                                * u(1:nrmtsd,n2,np2,isp2,iso2),db)
                                      radlr
     &                                = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                                      fci(1:nrmtsd) = r2(1:nrmtsd)
     &                                  * aimag( 
     &                                  u(1:nrmtsd,n1,np1,isp1,iso1)
     &                                * u(1:nrmtsd,n2,np2,isp2,iso2) )
                                      radli
     &                                = f_integr3(rr,fci,1,nrmtsd,Rmtsd)
                                    else
                                      fcr(1:nrmtsd) = r2(1:nrmtsd) 
     &                                  * ur(1:nrmtsd,n1,np1,isp1,iso1)
     &                                  * ur(1:nrmtsd,n2,np2,isp2,iso2)
                                      radlr = f_integr3(rr,fcr,1,nrmtsd,
     &                                                  Rmtsd)
                                      fci(:) = 0._db
                                      radli = 0._db
                                    endif

                                    rof_sd = cmplx( radlr, radli, db)

                                    if( Self .and. isp1 == isp2
     &                                 .and. lm1 == lm2 ) fc(1:nrmtsd)
     &                               = cmplx( fcr(1:nrmtsd),
     &                                 fci(1:nrmtsd),db) * ri2(1:nrmtsd) 

                                  endif

                                  Sta_e = c_harm * rof_sd
     &                                  * Taull(lma1,is1,lma2,is2)
                                  Sta_s = conjg( c_harm ) * rof_sd
     &                                  * Taull(lma2,is2,lma1,is1)
 
                                  Sta_r = 0.5_db * real( Sta_e-Sta_s,db)
                                  Sta_i = 0.5_db * aimag( Sta_e + Sta_s)

                                  State(lm1,isp1,lm2,isp2)
     &                              = State(lm1,isp1,lm2,isp2) - Sta_i 
                                  State_i(lm1,isp1,lm2,isp2)
     &                              = State_i(lm1,isp1,lm2,isp2) + Sta_r 

                                  if( Self .and. isp1 == isp2
     &                                .and. lm1 == lm2 ) then
                               
                                    lm = 0
                                    boucle_l3: do l3 = 0,lmax
                                      do m3 = -l3,l3
                                        lm = lm + 1
                                        if( lm > nlma2 ) exit boucle_l3 
  
                                        if( Ylm_comp ) then
                                          g = gauntcp(l1,m1,l2,m2,l3,m3)
                                        else
                                          g = gauntc(l1,m1,l2,m2,l3,m3)
                                        endif

                                        if( abs( g ) < eps10 ) cycle

! plus tard, faire la multiplication par g
                                        if( lm == 1 ) g = 1._db
                                        do ir = 1,nrmtsd 
                                          drho(ir,lm,isp1)
     &                                      = drho(ir,lm,isp1)
     &                                      - g * aimag( c_harm * fc(ir)
     &                                      * Taull(lma1,is1,lma2,is2) )
                                        end do
                                      end do

                                    end do boucle_l3

                                  endif
 
                                  if( icheck > 2 .and.
     &                              abs( Taull(lma1,is1,lma2,is2) )
     &                                            > 1.e-15_db ) then
                                    if( nlm2 > 1 ) then
                                      write(3,140) l1, mr1, m1, isp1,
     &                                  l2, mr2, m2, isp2, lp1, mp1,
     &                                  iso1, l2, mp2, iso2, rof_sd,
     &                                  Taull(lma1,is1,lma2,is2),
     &                                  State(lm1,isp1,lm2,isp2),
     &                                  State_i(lm1,isp1,lm2,isp2)
                                    elseif( nspino == 2 ) then
                                      write(3,150) l1, mr1, m1, isp1,  
     &                                  l2, mr2, m2, isp2, iso1, iso2,
     &                                  rof_sd,Taull(lma1,is1,lma2,is2),
     &                                  State(lm1,isp1,lm2,isp2),
     &                                  State_i(lm1,isp1,lm2,isp2)
                                    else
                                      write(3,160) l1, mr1, m1, isp1,  
     &                                  l2, mr2, m2, isp2, rof_sd,
     &                                  Taull(lma1,is1,lma2,is2),
     &                                  State(lm1,isp1,lm2,isp2),
     &                                  State_i(lm1,isp1,lm2,isp2)
                                    endif
                                  endif 

                                end do  ! boucle iso2
                              end do  ! boucle mp2
                            end do  ! boucle lp2
                          end do  ! boucle iso1
                        end do  ! boucle mp1
                      end do  ! boucle lp1

                    end do  ! boucle isg2
                  end do  ! boucle isg1
                end do  ! boucle mr2  
              end do  ! boucle l2
            end do  ! boucle isp2
          end do  ! boucle mr1
        end do  ! boucle l1
      end do  ! boucle isp1 

      if( icheck > 1 ) then
        write(3,170)
        do l1 = 0,lmax
          if( .not. Full_potential .and. l /= l1 ) cycle
          do mr1 = -l1,l1
            lm1 = lm0 + 1
            do l2 = 0,lmax
              if( .not. Full_potential .and. l /= l2 ) cycle
              do mr2 = -l2,l2
                lm2 = lm0 + 1
                if( Diagonale .and. mr1 /= mr2 ) cycle
                do isp1 = 1,nspin
                  do isp2 = 1,nspin
                    if( nspino == 1 .and. isp1 /= isp2 ) cycle
                    write(3,180) l1, mr1, isp1, l2, mr2, isp2, 
     &                 State(lm1,isp1,lm2,isp2),
     &                 State_i(lm1,isp1,lm2,isp2)
                  end do
                end do
              end do
            end do
          end do
        end do                  
      endif

      if( icheck > 2 .and. Self ) then
        write(3,190)
        do ir = 1,nrmtsd
          write(3,200) r(ir)*bohr, drho(ir,1,:)
        end do
      endif

      return
  110 format(/' Radial integral for the density of state:',/
     &  ' l1 mr1 m1 isp1 l2 mr2 m2 isp2 lp1 mp1 iso1 lp2 mp2 iso2',12x,
     &  'rof_sd',23x,'Taull',17x,'State',7x,'State_i')
  120 format(/' Radial integral for the density of state:',/
     &  ' l1 mr1 m1 isp1 l2 mr2 m2 isp2 iso1 iso2',12x,'rof_sd',23x,
     &  'Taull',17x,'State',7x,'State_i')
  130 format(/' Radial integral for the density of state:',/
     &  ' l1 mr1 m1 isp1 l2 mr2 m2 isp2',12x,'rof_sd',23x,'Taull',17x,
     &  'State',7x,'State_i')
  140 format(4(2i3,2i4,1x),1p,3(2x,2e13.5))
  150 format(2(2i3,2i4,1x),2(i4,1x),1p,3(2x,2e13.5))
  160 format(2(2i3,2i4,1x),1p,3(2x,2e13.5))
  170 format(/' Radial integral for the density of state before singul:'
     &  ,/' l1 m1 isp1 l2 m2 isp2',7x,'State',7x,'State_i')
  180 format(2(2i3,i4,1x),1p,3(2x,2e13.5))
  190 format(/' drho before singul',
     &  /'   Radius         up           dn')
  200 format(1p,3e13.5)
      end

!***********************************************************************

! Calcul de la solution singuliere pour la densite d'etat

      Subroutine Cal_Solsing_sd(drho,Ecomp,Eimag,f2,
     &        Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,ll,
     &        lmax,m_hubb,nlm1,nlm2,nlma,nlma2,nr,nrmtsd,nrs,nspin,
     &        nspino,numat,r,Radial_comp,Rmtsd,Self,
     &        Spinorbite,State,Tau,u,ur,V_hubb,Ylm_comp)

      use declarations
      implicit none

      integer:: icheck, ir, isg, isol, isp, l, ll, lm, lmax, lp, m,  
     &  m_hubb, mp, mr, ms, mv, nlm1, nlm2, n, n0, nlma, nlma2, np, nr, 
     &  nr1, nrmtsd, nrs, nspin, nspino, numat

      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(nlm1,nspin,nlm1,nspin):: Tau
      complex(kind=db), dimension(nrmtsd+1,nlm1,nlm2,nspin,nspino):: us 
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                          V_hubb
      complex(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: u 

      logical:: Ecomp, Full_potential, Hubb_a, Hubb_d, Radial_comp,
     &          Self, Spinorbite, Ylm_comp 

      real(kind=db):: c_harm, Eimag, f_integr3, Rmtsd, Sing
      real(kind=db), dimension(nr):: f2, r
      real(kind=db), dimension(nrmtsd):: fct, ri2, rr
      real(kind=db), dimension(nr,nspin):: g0, gm, gp
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nr,nlm1,nlm2,nspin,nspino):: ur 
      real(kind=db), dimension(nrs,nlma2,nspin):: drho
      real(kind=db), dimension(nlma,nspin,nlma,nspin):: State

      rr(1:nrmtsd) = r(1:nrmtsd) 

! gm et gp inverse dans le sousprogramme
      call Sch_radial_solsing(Ecomp,Eimag,f2,
     &         Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,
     &         ll,lmax,m_hubb,nlm1,nlm2,nr,nrmtsd,nspin,nspino,numat,r,
     &         Radial_comp,Rmtsd,Spinorbite,Tau,us,V_hubb)

      if( icheck == 2 ) write(3,110)

      ri2(1:nrmtsd) = 1 / ( quatre_pi * r(1:nrmtsd)**2 )

! Ici on prend r = r', car on calcule la trace de G(r,r')
! On ne calcule que la partie imaginaire.
      nr1 = 0
      do l = 0,lmax
        if( .not. Full_potential .and. l /= ll ) cycle
        n0 = l**2 + l + 1
        do mr = -l,l
          if( nlm1 == 1 .and. mr /= 0 ) cycle 
          nr1 = nr1 + 1

          do isp = 1,nspin

            fct(:) = 0._db

            do isg = -1,1,2
              if( Ylm_comp .and. mr /= 0 ) then
                c_harm = 0.5_db
                m = isg * mr
                if( isg == 1 ) then
                  n = nr1
                else
                  n = nr1 - 2 * mr
                endif
              else
                if( isg == 1 ) cycle
                c_harm = 1._db
                m = mr
                n = nr1
              endif

            np = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. lp /= l ) cycle
              do mp = -lp,lp
                if( nlm2 == 1 .and. mp /= m ) cycle 
                np = np + 1

                do isol = 1,nspino
                  if( Spinorbite .and. nlm2 == 1 ) then
                    ms = m + isol - isp
                    if( ms > l .or. ms < -l ) cycle
                  endif                         

                  do ir = 1,nrmtsd
                    if( Radial_comp ) then
                      fct(ir) = fct(ir) + c_harm * r(ir) * aimag(
     &                              u(ir,n,np,isp,isol)
     &                            * us(ir,n,np,isp,isol) )
                    else
                      fct(ir) = fct(ir) + c_harm * r(ir)
     &                                 * ur(ir,n,np,isp,isol)
     &                                 * aimag( us(ir,n,np,isp,isol) )
                    endif  
                  end do

                end do
              end do
            end do       
          end do 

            if( icheck > 2 .and. Self ) then
              write(3,115) l, m, isp
              do ir = 1,nrmtsd
                write(3,180) r(ir)*bohr, fct(ir)
              end do
            endif

            Sing = f_integr3(rr,fct,1,nrmtsd,Rmtsd)
            if( Self ) fct(1:nrmtsd) = fct(1:nrmtsd) * ri2(1:nrmtsd)

            if( nlm1 == 1 ) then
              do mv = -l,l
                lm = n0 + mv
                State(lm,isp,lm,isp) = State(lm,isp,lm,isp) + Sing
                if( Self ) drho(1:nrmtsd,1,isp)
     &                         = drho(1:nrmtsd,1,isp) + fct(1:nrmtsd)
                if( icheck > 2 ) write(3,110)
                if( icheck > 1 ) write(3,120) l, m, isp, Sing,
     &                                        State(lm,isp,lm,isp)
              end do
            else
              lm = n0 + m 
              State(lm,isp,lm,isp) = State(lm,isp,lm,isp) + Sing
              if( Self ) drho(1:nrmtsd,1,isp) = drho(1:nrmtsd,1,isp)
     &                                        + fct(1:nrmtsd)
              if( icheck > 2 ) write(3,110)
              if( icheck > 1 ) write(3,120) l, m, isp, Sing,
     &                                        State(lm,isp,lm,isp)
            endif

          end do
        end do
      end do       

      if( icheck > 2 ) then
        write(3,170)
        do ir = 1,nrmtsd
          write(3,180) r(ir)*bohr, drho(ir,1,:)
        end do
      endif

      return
  110 format(/' Radial integral of singular solution for the density',
     &' of state:',/'  l  m isp      Sing         State')
  115 format(/' Singular radial function for the density of state,',
     &' l, m, isp =',3i3)
  120 format(3i3,1p,1x,e13.5,2x,2e13.5)
  170 format(/' drho after singul',
     &  /'   Radius         up           dn')
  180 format(1p,3e13.5)
      end

!***********************************************************************

      subroutine Radial_wave(Ecinetic,Eimag,Energ,Enervide,
     &          Full_potential,Hubb_a,Hubb_d,icheck,initl,
     &          lmax,m_hubb,ninitlu,nlmam,nlmam2,nr,nspin,nspino,
     &          numat,r,Radial_comp,Rmtg,Rmtsd,Spinorbite,V_hubb,
     &          V_intmax,V0bd,Vrato,zet)
 
      use declarations
      implicit none

      integer:: icheck, initl, l, l_hubbard, lfin, lm0, lmax, lmp, m, 
     &  m_hubb, mp, n, ninitlu, nlm1, nlm2, nlmam, nlmam2, np, nr,
     &  nrmtg, nrmtsd, nspin, nspino, numat
      
      complex(kind=db), dimension(nspin):: konde
      complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspin)::
     &                                                       V_hubb
      complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
      complex(kind=db), dimension(:,:,:,:,:), allocatable:: u

      logical:: Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m,  
     &  Radial_comp, Relativiste, Renorm, Spinorbite

      real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
      real(kind=db), dimension(nspin):: Ecinetic, V0bd
      real(kind=db), dimension(nr,nspin):: Vrato 

      real(kind=db), dimension(nr):: f2, r 
      real(kind=db), dimension(nr,nspin):: g0 , gm, gp, V
      real(kind=db), dimension(nr,nspino):: gso 
      real(kind=db), dimension(nr,nlmam,nlmam2,nspin,nspino,ninitlu)::
     &                                                             zet

      real(kind=db), dimension(:,:,:,:,:), allocatable:: ur

      konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

      if( icheck > 1 ) then
        write(3,110) 
        write(3,120) Energ*rydb, Enervide*rydb
        write(3,130) Ecinetic(:)*rydb
        write(3,140) V0bd(:)*rydb
        write(3,150) konde(:)
      endif

      call mod_V(icheck,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,
     &                 V_intmax,V0bd,Vrato)

      call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nr,nspin,nspino,
     &                        numat,r,Relativiste,Spinorbite,V)

      gp(:,:) = 1 / gp(:,:)

      if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10
     &     .or. Ecinetic(nspin) < eps10 ) then
        Ecomp = .true.
      else
        Ecomp = .false.
      endif

      if( Full_potential ) then
        lfin = 0
      else
        lfin = lmax
      endif

      Renorm = .true.

      do l = 0,lfin

        if( Hubb_a .and. l == l_hubbard( numat ) )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif  

        if( Full_potential ) then
          nlm1 = ( lmax + 1 )**2
          nlm2 = nlm1
        elseif( Hubb_m .and. .not. Hubb_d ) then
          nlm1 = 2*l + 1
          nlm2 = nlm1
        elseif( Spinorbite .or. Hubb_m ) then
          nlm1 = 2*l + 1
          nlm2 = 1
        else
          nlm1 = 1
          nlm2 = 1
        endif

        allocate( u(nr,nlm1,nlm2,nspin,nspino) )
        allocate( ur(nr,nlm1,nlm2,nspin,nspino) )
        allocate( Tau(nlm1,nspin,nlm1,nspin) )

        call Sch_radial(Ecinetic,Ecomp,Eimag,f2,
     &           Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,
     &           l,lmax,m_hubb,nlm1,nlm2,nr,nrmtg,nspin,nspino,numat,r,
     &           Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,u,
     &           ur,V_hubb)

! Recopie
        if( Full_potential ) then
          if( Radial_comp ) then
            zet(1:nr,1:nlm1,1:nlm2,:,:,initl)
     &          = real( u(1:nr,1:nlm1,1:nlm2,:,:), db )
          else
            zet(1:nr,1:nlm1,1:nlm2,:,:,initl)
     &                               = ur(1:nr,1:nlm1,1:nlm2,:,:)
          endif
        else
          lm0 = l**2 + l + 1
          do m = -l,l
            if( nlm1 == 1 ) then
              n = 1
            else
              n = l + 1 + m
            endif
            do mp = -l,l
              if( nlm2 == 1 .and. mp /= m ) cycle
              np = min(l + m + 1, nlm2 )
              lmp = min(lm0 + mp, nlmam2 ) 
              if( Radial_comp ) then
                zet(1:nr,lm0+m,lmp,:,:,initl) = real( u(1:nr,n,np,:,:),
     &                                                             db )
              else
                zet(1:nr,lm0+m,lmp,:,:,initl) = ur(1:nr,n,np,:,:)
              endif
            end do
          end do
        endif

        deallocate( Tau )

        deallocate( u, ur )

      end do   ! fin de la boucle sur l

      return
  110 format(/' ---- Radial_wave --',100('-'))
  120 format(/' Energ =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd  =',2f10.3)
  150 format(' konde =',2f12.5)
      end

!***********************************************************************

! Calcul de Integrale_sur_r_et_r' de 
! conj(phi(r) * f_reg( min(r,r') ) * f_irg( max(r,r') ) * phi(r') * dr * dr' ) 
! = Integrale_sur_r ( phi(r) * f_irg(r) * Integrale_sur_r'_de_0_a_r ( f_reg(r') * phi(r') * dr' )
!                   + phi(r) * f_reg(r) * Integrale_sur_r'_de_r_a_Rmax ( f_irg(r') * phi(r') * dr' ) * dr 
!
! f_reg = r * solution reguliere
! f_irg = r * solution irreguliere
! phi = r * fonction initiale * r^p    ou   r * solution reguliere

      function integr_sing(n,phi1,phi2,f_reg,f_irg,Rmtsd,r,icheck) 

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db):: integr_sing

      complex(kind=db), dimension(n):: f_irg, f_reg, fct, phi_irg,
     &                                phi_reg, s_phi_irg, s_phi_reg 
      integer icheck
      real(kind=db), dimension(n):: fct_r, phi1, phi2,  r


      phi_reg(:) = phi1(:) * f_reg(:)
      call ffintegr2(s_phi_reg,phi_reg,r,n,1,Rmtsd)

      phi_irg(:) = phi2(:) * f_irg(:)
      call ffintegr2(s_phi_irg,phi_irg,r,n,-1,Rmtsd)

      fct(:) = phi_irg(:) * s_phi_reg(:) + phi_reg(:) * s_phi_irg(:)

      fct_r(:) = real( fct(:),db )
      fr = f_integr3(r,fct_r,1,n,Rmtsd)

      fct_r(:) = aimag( fct(:) )
      fi = f_integr3(r,fct_r,1,n,Rmtsd)

      integr_sing = - cmplx(fr,fi,db)

      if( icheck > 3 ) then
        write(3,110)
        do ir = 1,n
          write(3,120) r(ir)*bohr, f_reg(ir), phi_reg(ir),s_phi_reg(ir),  
     &                 phi_irg(ir), s_phi_irg(ir), fct(ir)
          if( r(ir) > Rmtsd ) exit
        end do
        write(3,130) integr_sing 
      endif

      return
  110 format(/5x,'Radius',15x,'f_reg',23x,'phi_reg',19x,'s_phi_reg',19x,
     &          'phi_irg',21x,'s_phi_irg',22x,'fct')
  120 format(1p,13e14.6)
  130 format(/' Integr_sing =',1p,2e14.6)
      end

!***********************************************************************

! Calcul l'integrale de 0 a r (is=1) ou r a rmtsd (is=-1) de fct
! Cas complexe

      subroutine ffintegr2(fint,fct,r,n,is,Rmtsd)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      complex(kind=db):: a, a_tiers, b, b_demi, c, f0, fm, fp
      complex(kind=db), dimension(2):: dintegr
      complex(kind=db), dimension(n):: fct, fint

      real(kind=db), dimension(n):: r

      tiers = 1._db / 3._db

      if( is == 1 ) then
        i1 = 1
        i2 = n - 1
      else
        i1 = n - 1
        i2 = 1
        fint(n) = (0._db, 0._db)
      endif

      do i = i1,i2,is
        if( i == 1 ) then
          rm = r(i)
          r0 = r(i+1)
          rp = r(i+2)
          fm = fct(i)
          f0 = fct(i+1)
          fp = fct(i+2)
          xm = 0._db
          x0 = rm
          xp = 0.5 * ( rm + r0 )
        else
          rm = r(i-1)
          r0 = r(i)
          rp = r(i+1)
          fm = fct(i-1)
          f0 = fct(i)
          fp = fct(i+1)
          xm = 0.5 * ( rm + r0 )
          x0 = r0
          xp = 0.5 * ( r0 + rp )
          if( r(i) > Rmtsd ) then
            rm = r(i-2); r0 = r(i-1); rp = r(i)
            fm = fct(i-2); f0 = fct(i-1); fp = fct(i)
          endif 
       endif

        if( is == 1 .and. r0 > Rmtsd ) then
          if( r(i-1) > Rmtsd ) then
            fint(i) = fint(i-1)
            cycle
          elseif( rm > Rmtsd ) then
            fint(i) = fint(i-1) + dintegr(2)
            cycle
          else
            x0 = rmtsd
          endif
        endif
        if( is == - 1 .and. r0 > Rmtsd ) then
          fint(i) = (0._db, 0._db)
          cycle
        endif
        if( xp > Rmtsd ) xp = Rmtsd

        a = ( fm * ( rp - r0 ) - f0 * ( rp - rm ) + fp * ( r0 - rm ) )
     &    / ( ( r0 - rm ) * ( rp - r0 ) * ( rp - rm ) )
        b = ( f0 - fm ) / ( r0 - rm ) - a * ( r0 + rm )
        c = f0 - a * r0**2 - b * r0

        a_tiers = a * tiers
        b_demi = b * 0.5_db

        if( is == 1 ) then
          dintegr(1) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 )
     &               + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
        else
          dintegr(1) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 )
     &               + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
        endif

        if( i == i1 ) then
          fint(i) = dintegr(1)
        else
          fint(i) = fint(i-is) + sum( dintegr(:) )
        endif

        if( is == 1 ) then
          dintegr(2) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 )
     &               + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
        else
          dintegr(2) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 )
     &               + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
        endif

      end do

      if( is == 1 ) fint(n) = fint(n-1)

      return
      end

!***********************************************************************

! Calcul l'integrale de 0 a r (is=1) ou r a rmtsd (is=-1) de fct
! Cas complexe

      subroutine ffintegr2_r(fint,fct,r,n,is,rmtsd)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db):: a, a_tiers, b, b_demi, c, f0, fm, fp
      real(kind=db), dimension(2):: dintegr
      real(kind=db), dimension(n):: fct, fint

      real(kind=db), dimension(n):: r

      tiers = 1._db / 3._db

      if( is == 1 ) then
        i1 = 1
        i2 = n - 1
      else
        i1 = n - 1
        i2 = 1
        fint(n) = 0._db
      endif

      do i = i1,i2,is
        if( i == 1 ) then
          rm = r(i)
          r0 = r(i+1)
          rp = r(i+2)
          fm = fct(i)
          f0 = fct(i+1)
          fp = fct(i+2)
          xm = 0._db
          x0 = rm
          xp = 0.5 * ( rm + r0 )
        else
          rm = r(i-1)
          r0 = r(i)
          rp = r(i+1)
          fm = fct(i-1)
          f0 = fct(i)
          fp = fct(i+1)
          xm = 0.5 * ( rm + r0 )
          x0 = r0
          xp = 0.5 * ( r0 + rp )
        endif

        if( is == 1 .and. r0 > rmtsd ) then
          if( r(i-1) > rmtsd ) then
            fint(i) = fint(i-1)
            cycle
          elseif( rm > rmtsd ) then
            fint(i) = fint(i-1) + dintegr(2)
            cycle
          else
            x0 = rmtsd
          endif
        endif
        if( is == - 1 .and. r0 > rmtsd ) then
          fint(i) = 0._db
          cycle
        endif
        if( xp > rmtsd ) xp = rmtsd

        a = ( fm * ( rp - r0 ) - f0 * ( rp - rm ) + fp * ( r0 - rm ) )
     &    / ( ( r0 - rm ) * ( rp - r0 ) * ( rp - rm ) )
        b = ( f0 - fm ) / ( r0 - rm ) - a * ( r0 + rm )
        c = f0 - a * r0**2 - b * r0

        a_tiers = a * tiers
        b_demi = b * 0.5_db

        if( is == 1 ) then
          dintegr(1) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 )
     &               + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
        else
          dintegr(1) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 )
     &               + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
        endif

        if( i == i1 ) then
          fint(i) = dintegr(1)
        else
          fint(i) = fint(i-is) + sum( dintegr(:) )
        endif

        if( is == 1 ) then
          dintegr(2) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 )
     &               + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
        else
          dintegr(2) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 )
     &               + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
        endif

      end do

      if( is == 1 ) fint(n) = fint(n-1)

      return
      end

!***********************************************************************

! Transformation harmo comp vers Harmo reel si Sens = .true.
!                                 sinon inverse 

      subroutine Trans_Tau(icheck,lmin,lmax,nspin,Spinorbite,Taull,Sens) 

      use declarations
      implicit none

      integer:: icheck, isp1, isp2, l, lm0, lmin, lmax, m, nspin

      complex(kind=db), dimension(:,:), allocatable:: Tau, Trans
      complex(kind=db), dimension((lmax+1)**2-lmin**2,nspin,
     &                            (lmax+1)**2-lmin**2,nspin):: Taull

      logical:: Sens, Spinorbite

      do isp1 = 1,nspin
        do isp2 = 1,nspin
          if( .not. Spinorbite .and. isp1 /= isp2 ) cycle
          do l = lmin,lmax

            if( l == 0 ) cycle

            if( lmin == lmax ) then
              lm0 = l + 1
            else
              lm0 = l**2 + l + 1
            endif
            allocate( Tau(-l:l,-l:l) )   
            allocate( Trans(-l:l,-l:l) )   
            do m = -l,l
              Tau(-l:l,m) = Taull(lm0-l:lm0+l,isp1,lm0+m,isp2)
            end do
        
            call cal_trans_lh(l,Trans)

            if( Sens ) then
              Tau = matmul( Trans,
     &                        matmul( Tau, Transpose( conjg(Trans) ) ) )
            else
              Tau = matmul( Transpose( conjg(Trans) ),
     &                        matmul( Tau, Trans ) )
            endif

            do m = -l,l
              Taull(lm0-l:lm0+l,isp1,lm0+m,isp2) = Tau(-l:l,m)
            end do

            if( icheck > 1 ) then
              if( nspin == 1 ) then
                if( Sens ) then
                  write(3,110) l
                else
                  write(3,120) l
                endif
              elseif( Spinorbite ) then
                if( Sens ) then
                  write(3,130) l, isp1, isp2
                else
                  write(3,140) l, isp1, isp2
                endif
              else
                if( Sens ) then
                  write(3,150) l, isp1
                else
                  write(3,160) l, isp1
                endif
              endif
              do m = -l,l
                write(3,170) m, Tau(m,:)
              end do
            endif 

            deallocate( Tau, Trans )

          end do
        end do
      end do

  110 format(/' Tau in real basis, l =',i2)
  120 format(/' Tau in complex basis, l =',i2)
  130 format(/' Tau in real basis, l =',i2,',isp1, isp2 =',2i2)
  140 format(/' Tau in complex basis, l =',i2,',isp1, isp2 =',2i2)
  150 format(/' Tau in real basis, l =',i2,',isp =',i2)
  160 format(/' Tau in complex basis, l =',i2,',isp =',i2)
  170 format(i3,1p,7(1x,2e13.5))
      return
      end

!***********************************************************************

! Transformation harmo comp vers Harmo reel pour un lh donne.
! La transformation inverse est le conjugue de la transpose

      subroutine Cal_Trans_lh(lh,Trans)

      use declarations
      implicit none

      integer:: is, lh, m1, m2
 
      complex(kind=db):: r2_r, r2_i
      complex(kind=db),dimension(-lh:lh,-lh:lh):: Trans

      real(kind=db):: r2

      Trans(:,:) = (0._db, 0._db)

      r2 = 1 / sqrt(2._db) 
      r2_r = cmplx( r2,    0._db, db)
      r2_i = cmplx( 0._db, r2,    db)

      do m1 = -lh,lh
        is = (-1)**m1 
        do m2 = -lh,lh
                     
          if( m1 == m2 ) then

            if( m1 == 0 ) then
              Trans(m1,m2) = (1._db,0._db)
            elseif( m1 > 0 ) then
              Trans(m1,m2) = r2_r
            else
              Trans(m1,m2) = is * r2_i
            endif

          elseif( m1 == - m2 ) then

            if( m1 > 0 ) then
              Trans(m1,m2) = is * r2_r
            else
              Trans(m1,m2) = - r2_i
            endif

          endif

        end do
      end do 
   
      return
      end
