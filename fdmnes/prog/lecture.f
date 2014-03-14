! FDMNES subroutines

! Sousprogramme de lecture des fichiers d'entree permettant d'etablir
! les dimensions des differents tableaux

      subroutine lectdim(Absauto,Atom_occ_hubb,Atom_nonsph,Bormann,
     &  Extract,
     &  Flapw,Full_self_abs,Hubbard,itape4,Magnetic,Memory_save,
     &  mpinodes,mpirank,n_file_dafs_exp,n_multi_run_e,nb_atom_conf_m,
     &  ncolm,neimagent,
     &  nenerg,ngamme,ngroup,ngroup_neq,nhybm,nklapw,nlatm,nlmlapwm,
     &  nmatsym,norbdil,npldafs,nple,nplrm,nspin,nspino,ntype,
     &  ntype_conf,
     &  Pdb,Readfast,Self_abs,Space_file,Taux,Temperature,Xan_atom)

      use declarations
      implicit none
      include 'mpif.h'

      integer:: i, ie, igamme, igr, igrdat, io, ipr, ipl, ispin, istat,
     &  it, itape4, j, jgr, jpl, kgr, l, ligne, lin_gam, mpinodes,
     &  mpirank, n, n_dic, n_file_dafs_exp, n_multi_run_e, na,
     &  nb_atom_conf_m, ncolm, neimagent,
     &  nenerg, ngamme, ngc, ngroup, ngroup_neq, nhybm, nklapw, nl,
     &  nlatm, nlmlapwm, nmatsym, nn, nnombre, norbdil, norbv, mpierr,
     &  npldafs, nple, nplm, nplrm, nspin, nspino, ntype, ntype_conf

      character(len=2) Chemical_Symbol_c, Symbol
      character(len=3) mot3
      character(len=6) mot6
      character(len=9) grdat
      character(len=10) Space_Group
      character(len=13) Spgr
      character(len=132) Fichier, Fichier_pdb, identmot, mots, motsb,
     &            nomstruct, nomvcoul, Space_file

      integer, dimension(:), allocatable :: igra, neq, numat

      logical Absauto, Atom_conf, Atom_nonsph, Atom_occ_hubb, Axe_loc,
     &   Bormann, Cylindre,    
     &   Extract, Flapw, Full_self_abs, Hubbard, Magnetic, Matper,   
     &   Memory_save, Pdb, Quadrupole, Readfast, Recup_potlapw, 
     &   Self_abs, Spherique, Taux, Temperature, Xan_atom

      real(kind=db):: de, def, E, phi, popatc, r, Theta
      real(kind=db), dimension(3) :: p
      real(kind=db), dimension(3,3) :: Mat
      real(kind=db), dimension(:), allocatable :: Egamme
      real(kind=db), dimension(:,:), allocatable :: posn, posout
    
      common/Axe_loc/ Axe_loc

      Absauto = .true.
      Atom_conf = .false.
      Atom_nonsph = .false.
      Atom_occ_hubb = .false.
      Axe_loc = .false.
      Cylindre = .false.
      Extract = .false.
      Flapw = .false.
      Full_self_abs = .false.
      Hubbard = .false.
      Matper = .false.
      Memory_save = .false.
      n_dic = 0
      n_file_dafs_exp = 0
      n_multi_run_e = 1
      nb_atom_conf_m = 0
      neimagent = 0
      nenerg = 131
      ngamme = 3
      nhybm = 1
      nklapw = 1
      nlmlapwm = 1
      nlatm = 0
      nmatsym = 1
      norbdil = 0
      nple = 0
      nspin = 1
      nspino = 1
      ntype = 0
      ntype_conf = 0
      Pdb = .false.
      Quadrupole = .false.
      Readfast = .false.
      recup_potlapw = .false. 
      Self_abs = .false.
      Space_Group = ' '
      Spherique = .false.
      Taux = .false.
      Temperature = .false.
      Xan_atom = .false.

      if( Bormann ) then
        npldafs = 36
      else
        npldafs = 0
      endif

      if( mpirank /= 0 ) goto 2000

      Rewind(itape4)

      do igrdat = 1,100000

        read(itape4,'(A)',end=1000) mots
        grdat = identmot(mots,9)
        if( grdat(1:1) == '!' ) cycle

        select case(grdat)

          case('absorbeur')
            absauto = .false.
            n = nnombre(itape4,132)
            if( n > 0 ) n_multi_run_e = 0
            do i = 1,100000
              n = nnombre(itape4,132)
              if( n < 1 ) exit
              n_multi_run_e = n_multi_run_e + n
              read(itape4,*)
            end do

          case('end')
            exit

          case('extract') 
            Extract = .true.

          case('hubbard')
            Hubbard = .true.

          case('spgroup')
            n = nnombre(itape4,132)
            read(itape4,'(A)') mots
            if( mots(1:1) == ' ' ) mots = adjustl(mots)
            Space_group = mots(1:10)

          case('temperatu')
            Temperature = .true.

          case('xan_atom') 
            Xan_atom = .true.

          case('range','rangel')
            ngamme = nnombre(itape4,132)
            if( mod(ngamme,2) == 0 ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,110)
                write(ipr,120) ngamme
              end do
              stop
            endif
            if( grdat == 'rangel' ) then
              lin_gam = 1
              if( ngamme /= 1 ) ngamme = 3
            else
              lin_gam = 0
            endif
            allocate( egamme(ngamme) )
            read(itape4,*,err=9999) egamme(1:ngamme)
            do igamme = 2,ngamme,2
              if( egamme(igamme) > eps6 ) cycle
              call write_error
              do ipr = 6,9,3
                write(ipr,110)
                write(ipr,130)
              end do
              stop
            end do

            E = egamme(1)
            if( ngamme == 1 ) then
              nenerg = 1
            elseif( egamme(3) <= egamme(1) ) then
              nenerg = 1
            elseif( lin_gam == 1 ) then
              def = 10 / rydb
              do ie = 2,10000000
                r = 1 + e / def
                r = max( r, 0.25_db )
                de = sqrt( r ) * egamme(2)
                e = e + de
                if( e > egamme(ngamme) + eps10 ) then
                  nenerg = ie - 1
                  exit
                endif
              end do
            else
              ngc = 2
              do ie = 2,1000000
                e = e + egamme(ngc)                     
                if( e > egamme(ngamme) + eps10 ) then
                  nenerg = ie - 1                     
                  exit
                elseif( e >= egamme(ngc+1) - eps10 ) then
                  if( ngc+1 == ngamme ) then 
                    nenerg = ie
                    exit
                  endif
                  if( egamme(ngc+3) <= egamme(ngc+1) ) then  
                    nenerg = ie
                    exit
                  endif
                  ngc = ngc + 2
                endif
              end do
            endif
            deallocate( egamme )

          case('eimag')
            do ie = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              read(itape4,*)
            end do
            neimagent = ie - 1          

          case('quadrupol')
            Quadrupole = .true.

          case('e1e2')
            Quadrupole = .true.

          case('e2e2')
            Quadrupole = .true.

          case('spinorbit')
            nspin = 2
            nspino = 2

          case('magnetism')
            nspin = 2

          case('dilatorb')
            do i = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              norbdil = norbdil + 1
              read(itape4,*)
            end do

          case('polarized')
            do jpl = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              read(itape4,*,err=9999) p(:)
              if( sum( p(:) )**2 < eps10 ) n_dic = n_dic + 1 
            end do
            nple = jpl - 1

          case('dafs')
            do ipl = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              select case(n)
                case(3)
                  read(itape4,*)
                  nn = nnombre(itape4,132)
                  select case(nn)
                    case(2,3,4,5)
                      read(itape4,*)
                    case(6)
                      read(itape4,*)
                      read(itape4,*)
                    case default
                      call write_error
                      do ipr = 6,9,3
                        write(ipr,140) ipl
                      end do
                      stop
                  end select
                case(5,6,7,8)
                  read(itape4,*)
                case default
                  call write_error
                  do ipr = 6,9,3
                     write(ipr,140) ipl
                  end do
                  stop
              end select
            end do
            npldafs = ipl - 1

          case('dafs_exp')
            do i = 1,3
              read(itape4,*)
            end do
            npldafs = 0
            do i = 1,1000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              read(itape4,'(/A)') Fichier
              Fichier = Adjustl(Fichier)
              l = len_trim(Fichier)
              if( l > 4 ) then
                if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.txt' 
              endif
              open(99, file=Fichier, status='old',iostat=istat)
              if( istat /= 0 ) call write_open_error(Fichier,istat,1)
              n = nnombre(99,100000)
              Close(99)
              npldafs = npldafs + n / 3
            end do
            n_file_dafs_exp = i - 1

          case('readfast')
            Readfast = .true.

          case('self_abs')
            Self_abs = .true.

          case('full_self')
            Full_self_abs = .true.

          case('atom_conf')

            Atom_conf = .true.
            do it = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              ntype_conf = ntype_conf + 1
              read(itape4,*) n
              nb_atom_conf_m = max( n, nb_atom_conf_m )
              backspace(itape4) 
              read(itape4,*) ( nl, i = 1,n+2 )
              nlatm = max( nlatm, nl )
            end do
            ntype = ntype_conf

          case('atom')

            do it = 1,100000
              n = nnombre(itape4,132)
              if( n == 3 ) cycle

              if( n > 0 ) then

                ntype = ntype + 1
                if( n == 1 ) then
                  read(itape4,*)
                else
                  read(itape4,*) n, nl
                  nlatm = max( nlatm, nl )
                endif

              else

                read(itape4,'(A)',end=1000) motsb
                open(8, file = motsb, status='old', iostat=istat)
                if( istat /= 0 ) then
                  backspace(itape4)
                  exit
                endif
                if( it /= 1 ) ntype = ntype + 1
                read(8,*)
                do i = 1,100000
                  read(8,'(A)') mot3
                  if(mot3 == '---') exit
                end do
                read(8,*,err=9999) numat, popatc, nl
                nlatm = max( nlatm, nl )
                Close(8)
                if( nl > 0 ) read(itape4,*)
              endif
            end do

! Pour le cas des atomes charge ou il faut ajouter une orbitale:
            nlatm = nlatm + 1

! Description de l'agregat :
          case('crystal','molecule','crystal_t','molecule_')
            if( grdat == 'crystal_t' .or. grdat == 'molecule_' )
     &                                 Taux = .true.

            Matper = grdat(1:5) == 'cryst'

            n = nnombre(itape4,132)
            if( n == 0 ) then
              call write_error
              read(itape4,'(A)') motsb
              do ipr = 6,9,3
                write(ipr,150) grdat
                write(ipr,'(A)') motsb
                write(ipr,160)
              end do
              stop
            endif
            read(itape4,*)

            if( Readfast .or. Taux ) then
              do igr = 1,100000
                read(itape4,*,err=999,end=999) i, p(:)
              end do
  999         continue
            else
              do igr = 1,100000
                n = nnombre(itape4,132)
                if( n == 0 ) exit
                if( n == 2 .or. n == 3 ) then
                  Axe_loc = .true.
                  read(itape4,*)
                  n = nnombre(itape4,132)
                endif
                norbv = 0
                select case(n)
                  case(4)
                    read(itape4,*)
                  case(5)
                    read(itape4,*,err=9999) i, p(:), norbv
                    if( norbv < 0 ) then
                      nhybm = max( nhybm, - norbv - 1 )
                    else
                      nhybm = max( nhybm, norbv )
                    endif
                  case default
                    call write_error
                    read(itape4,'(A)') motsb
                    do ipr = 6,9,3
                      write(ipr,150) grdat
                      write(ipr,'(A)') motsb
                      write(ipr,160)
                    end do
                    stop
                end select
                if( norbv == 0 ) cycle
                do io = 1,abs(norbv)
                  read(itape4,*)
                end do
                if( norbv /= -1 ) Atom_nonsph = .true.
                if( norbv < 0 ) then 
                  Atom_occ_hubb = .true.
                  Hubbard = .true.
                endif
              end do
            endif
            ngroup = igr - 1

! Description de l'agregat :
          case('crystal_p')

            Pdb = .true.
            Taux = .true.
            Matper = .true.

            Fichier_pdb = ' '
            read(itape4,'(A)') Fichier_pdb
            Fichier_pdb = Adjustl(Fichier_pdb)
            l = len_trim(Fichier_pdb)
            if( l > 4 ) then
              if( Fichier_pdb(l-3:l-3) /= '.' )
     &                             Fichier_pdb(l+1:l+4) = '.pdb'
            endif 
            open(8, file = Fichier_pdb, status='old', iostat=istat) 
            if( istat /= 0 ) call write_open_error(Fichier_pdb,istat,1)

            igr = 0

            do ligne = 1,100000
            
              read(8,'(a6)',end=1030,err=1030) mot6

              select case(mot6)

                case('END   ')

                  exit

                case('CRYST1' )

                  backspace(8) 
                  read(8,'(54x,a13)') Spgr

                  l = len_trim(Spgr)
                  j = 0
                  do i = 1,l
                    if( Spgr(i:i) == ' ' ) cycle
                    j = j + 1
                    Space_group(j:j) = Spgr(i:i)
                  end do 

                case('ATOM  ','HETATM') 

                  read(8,*)

                  igr = igr + 1

                case default

                  cycle

                end select

            end do

 1030       continue
       
            Close(8)
            ngroup = igr

          case('flapw','flapw_s','flapw_r','flapw_s_p','flapw_psi',
     &         'flapw_n','flapw_n_p','flapw_s_n')
            Flapw = .true.
            if( grdat(6:7) == '_s' ) then
              n = nnombre(itape4,132)
              read(itape4,*)
            elseif( grdat(6:7) == '_r' ) then
              recup_potlapw = .true.
              n = nnombre(itape4,132)
              read(itape4,*)
            endif
            n = nnombre(itape4,132)
            read(itape4,'(A)') nomstruct
            if( .not. recup_potlapw ) then
              n = nnombre(itape4,132)
              read(itape4,'(A)') nomvcoul
              if( nomvcoul(1:1) == ' ' ) nomvcoul = adjustl( nomvcoul )
              n = nnombre(itape4,132)
              read(itape4,*)
              if( grdat(6:7) == '_n' .and. nspin == 2 ) read(itape4,*) 
              do ispin = 1,2*nspin-1
                n = nnombre(itape4,132)
                read(itape4,*)
              end do
            endif
            n = nnombre(itape4,132)
            if( grdat /= 'flapw_s_p' .and. ( grdat /= 'flapw_psi' .and.
     &         grdat /= 'flapw_n_p' ) ) read(itape4,*)
            call lectpot_dim(ngroup,nklapw,nlmlapwm,nmatsym,
     &                       nomstruct,nomvcoul,ntype,recup_potlapw)

          case('memory_sa')
            Memory_save = .true.

        end select

      end do

 1000 continue

      if( Extract ) xan_atom = .false.
      if( Flapw .or. Extract ) Hubbard = .false. 
      if( npldafs == 0 ) Self_abs = .false.
      if( npldafs == 0 ) Full_self_abs = .false.
      if( Self_abs ) n_dic = n_dic + 2 * npldafs
      if( Full_self_abs ) n_dic = n_dic + 4 * npldafs
      if( nspin == 2 ) then
        Magnetic = .true.
      else
        Magnetic = .false.
      endif
      if( .not. Matper ) Space_Group = ' ' 
      if( Pdb ) Temperature = .true.

      ngroup_neq = ngroup

      if( Space_Group /= ' ' .or. ntype == 0 .or. Atom_conf ) then

        allocate( neq(ngroup_neq) )
        allocate( posn(3,ngroup_neq) )
        allocate( posout(3,ngroup) )
        allocate( numat(ngroup) )

        if( Pdb ) then

          open(8, file = Fichier_pdb, status='old', iostat=istat) 

          igr = 0

          do ligne = 1,100000
            
            read(8,'(a6)') mot6

            select case(mot6)

              case('END   ')

                exit

              case('SCALE1' )

                backspace(8)
                do i = 1,3
                  read(8,'(10x,3f10.6)') Mat(i,:)
                end do

              case('ATOM  ','HETATM') 

                backspace(8)
                read(8,'(30x,3f8.3,22x,a2)') p(:), Symbol
                read(8,*)

                igr = igr + 1
                p = Matmul( Mat, p )
                posn(:,igr) = p(:)    

                Symbol = adjustl(Symbol)              
                do i = 1,103
                  if( Chemical_Symbol_c(i) /= Symbol ) cycle
                  numat(igr) = i
                  exit
                end do
                if( i == 104 ) then
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,170) igr, Symbol
                  end do
                  stop
                endif
                if( igr == ngroup_neq ) exit
  
              case default

                cycle

            end select

          end do

          Close(8)
          
        else

          Rewind(itape4)

          do igrdat = 1,100000
            read(itape4,'(A)') mots
            grdat = identmot(mots,9)
            if( grdat(1:7) /= 'crystal' .and. grdat(1:7) /= 'molecul' )
     &                                                             cycle 

            n = nnombre(itape4,132)
            if( n == 2 ) then
              cylindre = .true.
            elseif( n == 1 ) then
              spherique = .true.
            endif
            read(itape4,*)
            do igr = 1,ngroup_neq
              if( Readfast .or. Taux .or. .not. ( Atom_nonsph
     &                       .or. Atom_occ_hubb .or. Axe_loc) ) then
                read(itape4,*) numat(igr), posn(:,igr)
              else
                n = nnombre(itape4,132)
                if( n == 0 ) exit
                if( n == 2 .or. n == 3 ) then
                  read(itape4,*)
                  n = nnombre(itape4,132)
                endif
                norbv = 0
                select case(n)
                  case(4)
                    read(itape4,*) numat(igr), posn(:,igr) 
                  case default
                    read(itape4,*,err=9999) numat(igr), posn(:,igr),
     &                                      norbv
                end select
                if( norbv == 0 ) cycle
                do io = 1,abs(norbv)
                  read(itape4,*)
                end do
              endif
            end do

            exit
          end do

        endif

        if( ntype == 0 ) then

          ntype = 1
          boucle_1: do igr = 2,ngroup_neq
            do jgr = 1,igr-1
              if( numat(igr) == numat(jgr) ) cycle boucle_1
            end do
            ntype = ntype + 1
          end do boucle_1

        elseif( Atom_conf ) then

          allocate( igra(ngroup_neq) )

          Rewind(itape4)

          do igrdat = 1,100000
            read(itape4,'(A)') mots
            grdat = identmot(mots,9)
            if( grdat /= 'atom_conf' ) cycle 

            na = 0
            do it = 1,ntype
              read(itape4,*) n, igra(na+1:na+n)
              na = na + n
            end do
            exit
          end do
          
          boucle_2: do igr = 1,ngroup_neq
            do kgr = 1,na
              if( igr == igra(kgr) ) cycle boucle_2
            end do
            boucle_3: do jgr = 1,igr-1
              do kgr = 1,na
                if( jgr == igra(kgr) ) cycle boucle_3
              end do
              if( numat(igr) == numat(jgr) ) cycle boucle_2
            end do boucle_3
            ntype = ntype + 1
          end do boucle_2

          deallocate( igra )

        endif

        if( Space_Group /= ' ' ) then

          do igr = 1,ngroup_neq
            if( cylindre ) then
              r = posn(1,igr)
              Theta = pi * posn(2,igr) / 180 
              posn(1,igr) = r * cos( theta )
              posn(2,igr) = r * sin( theta )
            elseif( spherique ) then
              r = posn(1,igr)
              theta = pi * posn(2,igr) / 180 
              phi = pi * posn(3,igr) / 180 
              posn(1,igr) = r * sin( theta ) * cos( phi)
              posn(2,igr) = r * sin( theta ) * sin( phi)
              posn(3,igr) = r * cos( theta )
            endif
          end do

          call spgroup(0,neq,ngroup,ngroup_neq,posn,posout,
     &                       Space_file,space_group)

        endif

        deallocate( posn )
        deallocate( posout )
        deallocate( neq )
        deallocate( numat )

      endif

 2000 continue

      if( mpinodes > 1 ) then
        call MPI_Bcast(Absauto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Atom_nonsph,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Atom_occ_hubb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Axe_loc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Memory_save,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Extract,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Flapw,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Full_self_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Hubbard,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Magnetic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(n_dic,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(n_multi_run_e,1,MPI_INTEGER,0,
     &                                            MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nb_atom_conf_m,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(neimagent,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nenerg,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(ngamme,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(ngroup,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(ngroup_neq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nhybm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nklapw,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nlatm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nlmlapwm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nmatsym,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(norbdil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(npldafs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nple,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nspin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(nspino,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(ntype,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(ntype_conf,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Pdb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Quadrupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Recup_potlapw,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Self_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Taux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Temperature,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Xan_atom,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 

        call MPI_BARRIER(MPI_COMM_WORLD,mpierr) 
      endif

      if( nple == 0 ) then
        if( Quadrupole ) then
          nplm = 6
        else
          nplm = 3
        endif
      else
        nplm = nple
      endif

      nplrm = nplm + n_dic
      ncolm = nplm + 2*npldafs + 2*n_dic + 1
      if( xan_atom ) ncolm = ncolm + 1
      if( self_abs ) ncolm = ncolm  + 2*npldafs
      if( Full_self_abs ) ncolm = ncolm  + 4*npldafs

      return

 9999 call write_err_form(itape4,grdat)

      return
  110 format(//'  Error in the indata file :')
  120 format(//' After keyword "Range", the number of value given for',
     &' the energies and steps',/' is even (',i2,'), it must be odd !',/
     &' If it is not the case, check the presence of',
     &' extra characters or tabulations.',/' They are forbidden !'//)
  130 format(//' Energy step is zero or negatif after keyword "Range",',
     &' forbidden !'//)
  140 format(//' After the keyword "Dafs", for the reflection number',
     &         i3,',',/' the polarization or the indexes are not',
     &        ' well set.',/' Check the format !'//)
  150 format(///' Just after the keyword "',A,'"',
     &  ', the following line is red :',/)
  160 format(/' It must not be there or it contains unwanted',
     &  ' characters !'//) 
  170 format(/' Error in the Pdb file :',//
     &        ' The chemical symbol of atom number',i6,' is: ',
     &        a2,', what is not known!'//)
      end

!*********************************************************************

      subroutine write_err_form(irec,keyword)

      character(len=9):: keyword
      character(len=132):: mot

      logical:: pb_line

      pb_line = .false.
      call write_error
      backspace(irec)
      read(irec,'(A)',err=1000) mot
      goto 1010
 1000 pb_line = .true.
 1010 continue

      do ipr = 6,9,3
        write(ipr,110) keyword
        if( pb_line ) then
          write(ipr,'(//5x,A//)')
     &     ' Check if the line is terminated by a cariage return !'
        else
          write(ipr,'(//A/)') ' The line is:'
          write(ipr,'(A)') mot
          write(ipr,120)
        endif
      end do
      stop

      return
  110 format(//' Format error when reading in the indata file under',
     &  ' the keyword:',//,5x,A)
  120 format(//' Check :',/
     &         '  - How many numbers must be in the line ?',/
     &         '  - Are there spaces between the numbers ?',/
     &         '  - Are there unwanted characters, extra  points... ?',/
     &         '  - Tabulations are forbidden !'//)
      end

!***********************************************************************

! Routine de lecture des potentiels et densites electroniques venant de
! WIEN

      subroutine lectpot_dim(ngroup,nklapw,nlmlapwm,nmatsym,
     &                       nomstruct,nomvcoul,ntype,recup_potlapw)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      character(len=1) Trans
      character(len=132) nomstruct, nomvcoul

      integer, dimension(:), allocatable :: nrato_lapw

      logical recup_potlapw

      open(8, file = nomstruct, status='old', iostat=istat) 
      if( istat /= 0 ) call write_open_error(nomstruct,istat,1)

      read(8,*)
      read(8,'(a1,27x,i2)') Trans, ntype
      select case(Trans)
        case('F')
          ntrans = 3
        case('B','C')
          ntrans = 1
        case default
          ntrans = 0
      end select
      read(8,*)
      read(8,*)

      allocate( nrato_lapw(ntype) )

      index = 0
      do ia = 1,ntype
        index = index + ntrans + 1
        read(8,*) 
        read(8,'(15x,i2)') mult
        do mu = 1,mult-1
          index = index + ntrans + 1
          read(8,*) 
        end do
        read(8,'(15x,i5)') nrato_lapw(ia)
        read(8,*) 
        read(8,*) 
        read(8,*) 
      end do

      read(8,'(i4)') nmatsym

      close(8)

      ngroup = index

      if( .not. recup_potlapw ) then

        open(8, file = nomvcoul, status='old', iostat=istat) 
        if( istat /= 0 ) call write_open_error(nomvcoul,istat,1)

        do ia = 1,ntype

          do i = 1,4
            read(8,*)
          end do
          read(8,'(16x,i2)') ll
          nlmlapwm = max( abs(ll), nlmlapwm )

          do l = 1,ll
            do i = 1,4
              read(8,*)
            end do
            do j = 1,nrato_lapw(ia),4
              read(8,*)
            end do
          end do

          do i = 1,3
            read(8,*)
          end do

        end do ! fin de la boucle sur les atomes

        do i = 1,5
          read(8,*)
        end do
        read(8,'(14x,i5)') nklapw

        close(8)

      endif

      deallocate ( nrato_lapw )

      return
      end

!*********************************************************************

!   Sousprogramme de lecture.
!   Toutes les entrees sauf les densites electroniques et les fonctions
! d'onde sont en Angstroem et en eV. Elles sont converties, pour tout le
! programme, en unites atomiques et Rydberg dans ce sous-programme.

      subroutine lectur(Absauto,adimp,alfpot,Allsite,Ang_borm,Angle_or,
     &  Angpoldafs,Angxyz,Atom_occ_hubb,Atom_nonsph,Atom_nsph_e,
     &  Atomic_scr,Axe_atom_gr,axyz,Base_spin,basereel,Bormann,       
     &  BSE,Charge_free,Clementi,com,comt,Core_resolved,
     &  coupelapw,Cubmat,Current,D_max_pot,Dafs,Dafs_bio,Density,dipmag,
     &  dpos,dyn_eg,dyn_g,Eclie,Ecrantage,Eeient,Egamme,Eimagent,
     &  Delta_En_conv,Delta_Epsii,E1E1,E1E2e,E1E3,E1M1,E1M2,E2E2,E3E3,
     &  Eneg_i,Eneg_n_i,Energphot,Extract,Fit_cal,Flapw,Flapw_new,
     &  Force_ecr,Full_atom_e,Full_potential,Full_self_abs,Gamma_hole,
     &  Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_int,Green_s,
     &  Green_self,hkl_borm,hkl_dafs,Hubb,Hubbard,Hybrid,iabsm,iabsorig,
     &  icheck,icom,indice_par,iscratch,isigpi,itdil,its_lapw,iord,
     &  itape4,itype,jseuil,Kern_fac,Kgroup,korigimp,l_selec_max,
     &  lamstdens,ldil,lecrantage,lin_gam,lmaxfree,lmaxso0,lmaxat0,
     &  lmoins1,lplus1,lseuil,lvval,m_hubb_e,M1M1,M1M2,M2M2,
     &  Magnetic,Mat_or,Matper, mix_repr,mpinodes,mpirank,Muffintin,
     &  multrmax,n_atom_proto,n_file_dafs_exp,n_multi_run_e,
     &  nb_atom_conf_m,nbseuil,nchemin,necrantage,neimagent,nenerg,
     &  ngamh,ngamme,ngroup,ngroup_hubb,ngroup_lapw,ngroup_m,
     &  ngroup_neq,ngroup_nonsph,ngroup_par,ngroup_pdb,ngroup_taux,
     &  ngroup_temp,nhybm,nlat,nlatm,nnlm,No_solsing,nom_fich_Extract,
     &  nomfich,nomfich_tddft_data,nomfichbav,nomclm,nomfile_atom,
     &  nompsii,nomr2v,nomr2vdn,nomstruct, nomvcoul,Noncentre,
     &  nonexc,norbdil,norbv,normaltau,normrmt,npar,nparm,nphi_dafs,
     &  nphim,npldafs,nple,nposextract,nrato,nrato_dirac,nrato_lapw,nrm,
     &  nself,nseuil,nspin,nsymextract,ntype,ntype_conf,numat,numat_abs,
     &  nvval,occ_hubb_e,Octupole,One_run,Optic,overad,overlap,p_self0,
     &  param,pdpolar,PointGroup,polar,poldafsem,poldafssm, 
     &  pop_nonsph,popatc,popats,popatv,popval,posn,quadmag,Quadrupole,
     &  r0_lapw,Radius_current,rchimp,Readfast,Recup_tddft_data,
     &  Relativiste,r_self,rlapw,rmt,rmtimp,Rot_Atom_gr,rotloc_lapw,
     &  roverad,RPALF,rpotmax,rydberg,rsorte_s,Save_tddft_data,
     &  scf_elecabs,SCF_mag_free,Self_abs,
     &  self_cons,self_nonexc,Solsing_s,solsing_only,Space_file,
     &  Spinorbite,state_all,state_all_out,Struct,supermuf,symauto,Taux,
     &  Taux_oc,tddft,Tddft_mix,Tddft_so,Temp,Temp_coef,Temperature,
     &  Test_dist_min,Trace_format_wien,Typepar,V_hubbard,V_intmax,
     &  Vec_orig,vecdafsem,vecdafssm,veconde,v0bdcFimp,Ylm_comp_inp,
     &  Z_nospinorbite)

      use declarations
      implicit none
      include 'mpif.h'

      integer:: i, ia, ie, igr, igrdat, io, iord, ip, ipar, ipl,
     &  ipl0, ipr, ipr0, iscratch, isp, ispin, istat, 
     &  istop, isymeq, it, itape4, j, jgr, jpl, jseuil, jt, k, kgr,
     &  l, l_hubbard, l_level_val, l_selec_max, l1, l2, lamstdens,
     &  lecrantage, ligne,
     &  lin_gam, lmaxat0, lmaxso0, long, lseuil, ltrace, m, 
     &  m_hubb_e, mpierr, mpinodes, mpirank, multi_run, multrmax,
     &  n, n_atom_proto, n_file_dafs_exp, n_multi_run_e, n_orb_rel, 
     &  n1, n2, natomsym, nb_atom_conf_m, nbseuil,
     &  nchemin, necrantage, neimagent, nenerg, ngamme, ngamh, 
     &  ngroup, ngroup_hubb, ngroup_lapw, ngroup_m, ngroup_neq,
     &  ngroup_nonsph,
     &  ngroup_par, ngroup_pdb, ngroup_taux, ngroup_temp, nhybm, nlatm,
     &  nn, nnlm, nnombre, non_relat, norb, norbdil, normrmt, 
     &  nparm, nphim, nphimt, npldafs, nple, nrato_dirac, nrm, nscan,
     &  nself, 
     &  nseuil, nspin, ntype, ntype_conf, numat_abs, Z_nospinorbite

      character(len=1):: Let
      character(len=2):: Chemical_Symbol_c, Symbol
      character(len=3):: mot3, seuil
      character(len=5):: struct
      character(len=6):: mot6
      character(len=8):: PointGroup
      character(len=9):: grdat
      character(len=10):: Space_Group
      character(len=11):: motpol
      character(len=13):: Chemical_Name, mot13, Spgr 
      character(len=35), dimension(0:ntype):: com
      character(len=50):: com_date, com_time, Revision 
      character(len=132), dimension(0:ntype):: nomfile_atom
      character(len=132), dimension(n_file_dafs_exp):: File_dafs_exp
      character(len=132) comt, fdmnes_error, Fichier,
     &     identmot, file_sauve_potlapw, mots, motsb, nomfich,
     &     nomfich_tddft_data, nom_fich_Extract, nomfichbav, nomstruct,
     &     nomvcoul, 
     &     nomr2v, nomr2vdn, nomclm(2*nspin-1), nompsii, Space_file
      character(len=9), dimension(ngroup_par,nparm):: typepar

      integer, dimension(2):: mix_repr
      integer, dimension(3):: hkl_borm
      integer, dimension(30):: icheck
      integer, dimension(n_multi_run_e):: iabsm, iabsorig, nposextract,
     &                                    nsymextract  
      integer, dimension(ngroup_par) :: npar
      integer, dimension(ngroup_lapw):: its_lapw
      integer, dimension(norbdil):: itdil, ldil
      integer, dimension(0:ntype):: icom, nlat, nrato, nrato_lapw,
     &                                 numat
      integer, dimension(ngroup):: itype
      integer, dimension(ngroup_pdb):: Kgroup 

      integer, dimension(0:ngroup_nonsph):: norbv
      integer, dimension(npldafs):: nphi_dafs
      integer, dimension(3,npldafs):: hkl_dafs
      integer, dimension(npldafs,2):: isigpi
      integer, dimension(ntype_conf):: nb_atom_conf
      integer, dimension(nb_atom_conf_m,ntype_conf):: igra
      integer, dimension(ngroup_par,nparm):: indice_par
      integer, dimension(0:ntype,nlatm):: lvval, nvval
      integer, dimension(:), allocatable :: itZ, neq

      complex(kind=db), dimension(3,npldafs):: poldafsem, poldafssm
      complex(kind=db), dimension(nhybm,16,ngroup_nonsph):: Hybrid

      logical Absauto, Allsite, Atom, Atom_conf, Atom_nonsph,
     &  Atom_occ_hubb, Atomic_scr, Axe_loc, Base_spin, Basereel, 
     &  Bormann, BSE, Cartesian_tensor, Charge_free,   
     &  Centre_auto, Centre_auto_abs, Clementi, Core_resolved,  
     &  Core_resolved_e, Coupelapw, Cylindre, Current, Dafs, Dafs_bio,
     &  Density, Dipmag, dyn_eg, dyn_g,   
     &  E1E1, E1E2e, E1E3, E1M1, E1M2, E2E2, E3E3, 
     &  eneg_i, eneg_n_i, Energphot,      
     &  exc_imp, Extract, Fermi_auto, Fit_cal, Flapw, Flapw_new,  
     &  Force_ecr, Full_atom_e, Full_potential, Full_self_abs,
     &  Gamma_hole_imp, Gamma_tddft, Green_s, Green_self, Green_int, 
     &  Hedin, Hubbard, korigimp, lmaxfree, lmoins1, lplus1, M1M1,  
     &  M1M2, M2M2, magnetic, matper, muffintin, noncentre, nonexc,  
     &  normaltau, no_core_resolved, no_dipquad, no_e1e3, no_e2e2,  
     &  no_e3e3, No_solsing, Octupole, old_reference, One_run,
     &  Optic, Overad, Pdb, Perdew, 
     &  PointGroup_Auto, polarise, quadmag, Quadrupole, r_self_imp,
     &  Readfast, recup_potlapw, Recup_tddft_data, Relativiste,
     &  rydberg, sauve_potlapw, Save_tddft_data, scf_elecabs,     
     &  SCF_mag_free, Self_abs, self_cons, self_exc_imp, self_nonexc,  
     &  self_nonexc_imp, single_prec, solsing_only, solsing_s, 
     &  spherical_signal, spherical_tensor, spherique, 
     &  Spinorbite, State_all, State_all_out, Supermuf, Symauto,   
     &  symmol, Taux, Tddft, Tddft_mix, Tddft_so, Temperature,
     &  Trace_format_wien, rpalf, Ylm_comp_inp  

      logical, dimension(ngroup):: Atom_nsph_e
      logical, dimension(0:ntype):: Hubb

      real(kind=db):: Adimp, Alfpot, Ang_borm, D_max_pot,
     &  Delta_En_conv, Delta_Epsii, Eclie, f_no_res_mag, f_no_res_mom,
     &  g1, g2, Gamma_max, Overlap, p_self0, phi, pop_nsph, pp, q, r,
     &  rad, Radius_current, Rmtt, rn, Roverad, Rpotmax, rrydb,
     &  rsorte_s, Step_azim, t, tc, Temp, Test_dist_min, Theta,
     &  V_intmax, vv

      real(kind=db), dimension(3):: Ang, Ang_spin, angxyz, Axe,
     &  Axe_spin, axyz, Centre, dpos, p, Vec_orig
      real(kind=db), dimension(10):: Gamma_hole
      real(kind=db), dimension(3,3):: Cubmat, Cubmati, Mat, Mat_or, Rot,
     &                                Rot_gen
      real(kind=db), dimension(norbdil):: cdil
      real(kind=db), dimension(0:ntype):: popatc, r0_lapw, rchimp,
     &  rlapw, rmt, rmtimp, V_hubbard
      real(kind=db):: Kern_fac, r_self
      real(kind=db), dimension(neimagent):: eeient, eimagent
      real(kind=db), dimension(ngamme):: egamme
      real(kind=db), dimension(nspin):: ecrantage, V0bdcFimp
      real(kind=db), dimension(3,nple):: polar, veconde
      real(kind=db), dimension(nple,2):: pdpolar
      real(kind=db), dimension(3,0:ntype) :: Ang_base_loc
      real(kind=db), dimension(n_file_dafs_exp):: Angle_dafs_exp
      real(kind=db), dimension(ngroup_taux):: Taux_oc 
      real(kind=db), dimension(ngroup_temp):: Temp_coef 
      real(kind=db), dimension(3,ngroup):: Ang_base_loc_gr, posn
      real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
      real(kind=db), dimension(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,
     &                                  nspin,ngroup_taux):: occ_hubb_e
      real(kind=db), dimension(3,3,ngroup_m):: Rot_atom_gr
      real(kind=db), dimension(3,3,ngroup_lapw):: rotloc_lapw
      real(kind=db), dimension(npldafs):: Angle_or 
      real(kind=db), dimension(3,npldafs):: angpoldafs, poldafse, 
     &  poldafsei, poldafss, poldafssi, vecdafsem, vecdafssm
      real(kind=db), dimension(ngroup_par,nparm):: param 
      real(kind=db), dimension(nhybm,ngroup_nonsph):: pop_nonsph
      real(kind=db), dimension(ngroup,nlatm,nspin):: popats
      real(kind=db), dimension(0:ntype,nlatm):: popatv
      real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
      real(kind=db), dimension(:), allocatable:: Hyb, x
      real(kind=db), dimension(:,:), allocatable:: pos
      real(kind=db), dimension(:,:,:), allocatable:: Hybrid_i, Hybrid_r

      integer, dimension(3):: ldipimp 
      integer, dimension(3,3):: lquaimp
      real(kind=db), dimension(3):: ang_rotsup, vectrace, ptrace 
 
      common/ang_rotsup/ ang_rotsup
      common/Axe_loc/ Axe_loc
      common/cartesian/ cartesian_tensor 
      common/com_out/ com_date, com_time, fdmnes_error, Revision
      common/f_no_res/ f_no_res_mag, f_no_res_mom
      common/file/ file_sauve_potlapw
      common/ldipimp/ ldipimp, lquaimp
      common/old_reference/ old_reference
      common/PointGroup_Auto/ PointGroup_Auto
      common/polarise/ polarise
      common/recup/ recup_potlapw, sauve_potlapw
      common/rrydb/ rrydb
      common/single_prec/ single_prec
      common/spheric/ spherical_tensor 
      common/spherical_signal/ spherical_signal 
      common/symmol/ symmol
      common/trac1/ vectrace, ptrace
      common/trac2/ ltrace

! Parametres par defaut
      adimp = 0.25_db
      alfpot = 0._db
      Allsite = .false.
      Ang_base_loc(1,:) = -10000._db; Ang_base_loc(2:3,:) = 0._db
      Ang_base_loc_gr(1,:) = -10000._db; Ang_base_loc_gr(2:3,:) = 0._db
      ang_rotsup = 0._db
      Angpoldafs(:,:) = 0._db
      Atom = .false.
      Atom_conf = .false.
      Atomic_scr = .false.
      Ang_spin(:) = 0._db
      Axe_spin(1) = 0._db; Axe_spin(2) = 0._db; Axe_spin(3) = 1._db
      Base_spin = .false.
      Basereel = .true.
      BSE = .false.
      cdil(:) = 0._db
      Cartesian_tensor = .false.
      Centre(:) = 0._db
      Centre_auto = .false.
      Centre_auto_abs = .false.
      Clementi = .false.
      com(:) = ' Dirac'
      Current = .false.
      E1E1 = .true.
      E1E2e = .false.
      E1E3 = .false.
      E1M1 = .false.
      E1M2 = .true.
      E2E2 = .false.
      E3E3 = .false.
      M1M1 = .false.
      M1M2 = .false.
      M2M2 = .false.
      Coupelapw = .false.
      Cylindre = .false.
      D_max_pot = 2.5_db
      Dafs = .false.
      Dafs_bio = .false.
      dipmag = .false.
      dyn_eg = .false.
      dyn_g = .false.
      dpos(:) = 0._db
      Ecrantage(:) = 1._db / nspin
      Eclie = 1._db
      Delta_En_conv = 1._db
      Delta_Epsii = 1000000._db * Rydb
      Density = .false.
      Eneg_i = .false.
      Eneg_n_i = .false.
      Energphot = .false.
      Exc_imp = .false.
      f_no_res_mag = 1._db
      f_no_res_mom = -100._db
      Fermi_auto = .true.
      Flapw_new = .false.
      Force_ecr = .false.
      Full_atom_e = .false.
      Full_potential = .false.
      Gamma_tddft = .false.
      Green_int = .false.
      Green_s = .false.
      Green_self = .true.
      Hedin = .false.
      iabsm(1) = 1
      Charge_free = .false.
      icom(:) = 1
      iord = 4
      isigpi(:,:) = 0
      istop = 0
      Hubb(:) = .false.
      Kern_fac = 2._db
      korigimp = .false.
      lamstdens = -1
      ldipimp(:) = -1
      lquaimp(:,:) = -1
      lin_gam = -1
      lmaxso0 = -5
      lmaxat0 = -1
      lmaxfree = .false.
      lmoins1 = .false.
      lplus1 = .false.
      ltrace = 0
      matper = .true.
      muffintin = .false.
      multrmax = 1
      nchemin = - 1
      nlat(:) = 0
      no_dipquad = .false.
      no_e1e3 = .false.
      no_e3e3 = .false.
      no_e2e2 = .false.
      No_solsing = .false.
      nonexc = .false.
      noncentre = .false.
      non_relat = 0
      norbv(:) = 0
      normaltau = .false.
      normrmt = 1
      nphim = 180
      nrato(:) = 0
      nrato_dirac = 600
      nrm = 0
      nself = 0
      nsymextract(:) = 1
      do i = 1,n_multi_run_e
        nposextract(i) = i
      end do
      numat_abs = 0
      if( Atom_occ_hubb ) occ_hubb_e(:,:,:,:) = 0._db
      Octupole = .false.
      old_reference = .true.
      One_run = .false.
      Optic = .false.
      overad = .false.
      overlap = 0.1_db
      p_self0 = 0.1_db
      Pdb = .false.
      pdpolar(:,:) = 0._db
      perdew = .false.
      polar(:,:) = 0._db
      polarise = .false.
      popatc(:) = 0._db
      popats(:,:,:) = 0._db
      quadmag = .false.
      Quadrupole = .false.
      Radius_current = 2.5_db
      recup_potlapw = .false.
      Recup_tddft_data = .false.
      relativiste = .false.
      rchimp(:) = 0._db
      r_self_imp = .false.
      rmt(:) = 0._db
      rmtimp(:) = 0._db
      roverad = 0._db
      rpotmax = 0._db
      rsorte_s = 3._db / bohr
      rydberg = .false.
      rrydb = 1._db
      sauve_potlapw = .false.
      Save_tddft_data = .false.
      scf_elecabs = .false.
      SCF_mag_free = .false.
      self_cons = .false.
      self_exc_imp = .false.
      self_nonexc = .true. 
      self_nonexc_imp = .false.
      seuil = 'K1'
      single_prec = .false.
      solsing_s = .false.
      solsing_only = .false.
      Space_Group = ' '
      spherical_tensor = .false.
      spherical_signal = .false.
      spherique = .false.
      Spinorbite = .false.
      core_resolved_e = .false.
      no_core_resolved = .false.
      state_all = .false.
      state_all_out = .false.
      supermuf = .false.
      PointGroup = ' '
      PointGroup_Auto = .true.
      RPALF = .false.
      symauto = .true.
      symmol = .false.
      if( Taux ) Taux_oc(:) = 1._db
      if( Temperature ) Temp_coef(:) = 0._db
      Tddft = .false.
      Tddft_mix = .false.
      Tddft_so = .true.
      Temp = 0._db
      Test_dist_min = 0.7_db * bohr ! distance minimum entre 2 atomes
      trace_format_wien = .false.
      veconde(:,:) = 0._db
      v0bdcFimp(:) = 0._db
      V_hubbard(:) = 0._db
      v_intmax = 1000000 * rydb
      Vec_orig(1:2) = 0._db; Vec_orig(3) = 1._db 
      Ylm_comp_inp = .false.
      Z_nospinorbite = 0

      if( mpirank > 0 ) goto 1320

      Rewind(itape4)

      write(6,*)

      do igrdat = 1,100000

        read(itape4,'(A)',end=1310) mots
        grdat = identmot(mots,9)
        if( grdat(1:1) /= ' ' ) write(6,'(3x,A)') grdat

        select case(grdat)

          case('full_pote') 

            Full_potential = .true.

          case('optic') 

            Optic = .true.

          case('clementi') 
            clementi = .true.

          case('current')
            Current = .true.
            n = nnombre(itape4,132)
            if( n > 0 ) read(itape4,*,err=9999) Radius_current

          case('scf_mag_f') 
            SCF_mag_free = .true.

          case('xan_atom') 

          case('spherical') 
            spherical_tensor = .true.

          case('sphere_al') 
            spherical_tensor = .true.
            spherical_signal = .true.

          case('cartesian') 
            cartesian_tensor = .true.

          case('extract')
            n = nnombre(itape4,132)
            read(itape4,'(A)') nom_fich_Extract
            if( nom_fich_Extract(1:1) == ' ' )
     &          nom_fich_Extract = adjustl( nom_fich_Extract )
            l = len_trim(nom_fich_Extract)
            if( nom_fich_Extract(l-3:l) == '_bav' )
     &         nom_fich_Extract(l+1:l+4) = '.txt'   

          case('extractsy')
            n = nnombre(itape4,132)
            n = min(n,n_multi_run_e)
            read(itape4,*) nsymextract(1:n)

          case('extractpo')
            n = nnombre(itape4,132)
            n = min(n,n_multi_run_e)
            read(itape4,*) nposextract(1:n)

          case('trace')
            n = nnombre(itape4,132)
            coupelapw = .true.
            read(itape4,*,err=9999) ltrace, vectrace(1:3), ptrace(1:3)
            if( grdat == 'trace_for' .or. grdat == 'trace_wie' )
     &        trace_format_wien = .true.

          case('range','rangel')
            n = nnombre(itape4,132)
            if( grdat == 'rangel' ) then
              lin_gam = 1               
            else
              lin_gam = 0
            endif
            read(itape4,*,err=9999) egamme(1:ngamme)

          case('energphot')
            energphot = .true.

          case('density')
            Density = .true.

          case('state_all')
            Density = .true.
            state_all = .true.
            state_all_out = .true.

          case('supermuf')
            supermuf = .true.

          case('old_refer')
            old_reference = .true.

          case('new_refer')
            old_reference = .false.

          case('etatlie')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) eclie

          case('eneg')
            Eneg_i = .true.
            Eclie = 0._db

          case('not_eneg')
            Eneg_n_i = .true.

          case('rydberg')
            rydberg = .true.
            n = nnombre(itape4,132)
            if( n > 0 ) read(itape4,*,err=9999) rrydb

          case('nchemin')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) nchemin
            normaltau = .true.

          case('lmoins1')
            lmoins1 = .true.

          case('lplus1')
            lplus1 = .true.

          case('base_reel')
            basereel = .true.

          case('base_comp')
            basereel = .false.

          case('base_spin')
            base_spin = .true.

          case('spinorbit')
            Spinorbite = .true.

          case('core_reso')
            core_resolved_e = .true.

          case('no_core_r')
            no_core_resolved = .true.

          case('ang_spin')
            n = nnombre(itape4,132)
            n = min(3,n)
            read(itape4,*,err=9999) Ang_spin(1:n)
            Ang_spin(1:n) = Ang_spin(1:n) * pi / 180

          case('axe_spin')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) Axe_spin(1:3)

          case('rot_sup')
            n = nnombre(itape4,132)
            n = min(3,n)
            read(itape4,*,err=9999) ang_rotsup(1:n)

          case('relativis')
            relativiste = .true.

          case('non_relat')
            non_relat = 1

          case('no_res_ma')
            read(itape4,*,err=9999) f_no_res_mag

          case('no_res_mo')
            read(itape4,*,err=9999) f_no_res_mom

          case('z_nospino')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) Z_nospinorbite

          case('debye')
            read(itape4,*,err=9999) temp
              
            if(temp > 10000-eps10) then
              call write_error
              do ipr = 6,9,3
                write(ipr,'(/A/)')
     &                 ' The temperature must be less than 10000K!'
              end do
              stop  
            elseif(temp <= - eps10) then
              call write_error
              do ipr = 6,9,3
                write(ipr,'(//A//)')
     &                      ' Negative temperatures are not allowed!' 
              end do
              stop
            end if

          case('radius')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) rsorte_s

          case('nrato')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) nrato_dirac

          case('multrmax')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) multrmax

          case('rpotmax')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) rpotmax

          case('over_rad')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) roverad
            overad = .true.

          case('screening')
            n = nnombre(itape4,132)
            if( n < 3 ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,'(/A)')
     &            ' Wrong format under the keyword screening'
              end do
              stop
            endif
            n = n - 2
            read(itape4,*,err=9999) necrantage, lecrantage,
     &                     ecrantage(1:min(nspin,n))
            if( n == 1 .and. magnetic ) then
              ecrantage(nspin) = ecrantage(1) / nspin
              ecrantage(1) = ecrantage(nspin)
            endif
            Force_ecr = .true.
            Charge_free = .true.

          case('tddft')
            tddft = .true.
            
          case('rpalf')
            tddft = .true.
            rpalf = .true.

          case('dyn_eg')   ! noyau TDDFT fxc diag sur seuils
            tddft = .true.
            dyn_eg = .true.

          case('dyn_g')  ! noyau TDDFT fxc diag sur etats initiaux
            tddft = .true.
            dyn_g = .true.

          case('bse')
            BSE = .true.
            tddft = .true.

          case('kern_fac')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) Kern_fac

          case('gamma_tdd')
            Gamma_tddft = .true.

          case('tddft_dat')
            Recup_tddft_data = .true.
            read(itape4,'(A)',err=9999) nomfich_tddft_data

          case('save_tddf')
            Save_tddft_data = .true.

          case('tddft_sca')
            tddft = .true.
            tddft_so = .false.

          case('mix_repr')
            tddft_mix = .true.
            state_all = .true.  ! to ensure that all repr are calculated
            tddft = .true.
            rpalf = .true.
            read(itape4,*) mix_repr(1:2) 

          case('atomic_sc')
            Atomic_scr = .true.

          case('edge')
            n = nnombre(itape4,132)
            read(itape4,'(A)') motsb
            seuil = identmot(motsb,3)
            select case( seuil(1:1) )
              case('k')
                seuil = 'K1'
              case('l')
                seuil(1:1) = 'L'
              case('m')
                seuil(1:1) = 'M'
              case('n')
                seuil(1:1) = 'N'
              case('o')
                seuil(1:1) = 'O'
              case('p')
                seuil(1:1) = 'P'
            end select

          case('quadrupol')
            Quadrupole = .true.

          case('octupole','dipole_oc','dip_oct')
            Octupole = .true.

          case('e1e2')
            E1E2e = .true.

          case('e1e3')
            E1E3 = .true.

          case('e2e2')
            E2E2 = .true.

          case('e1m1')
            E1M1 = .true.

          case('e1m2')
            E1M2 = .true.

          case('e3e3')
            E3E3 = .true.

          case('m1m1')
            M1M1 = .true.

          case('m2m2')
            M2M2 = .true.

          case('dipmag')
            dipmag = .true.
            E1M1 = .true.
            M1M1 = .true.

          case('quadmag')
            quadmag = .true.
            E1M2 = .true.
            M2M2 = .true.
            M1M2 = .true.

          case('no_e1e1')
            E1E1 = .false.

          case('no_e2e2')
            no_e2e2 = .true.

          case('no_e1e3')
            no_e1e3 = .true.

          case('no_e1e2')
            no_dipquad = .true.

          case('no_e3e3')
            no_e3e3 = .true.

          case('ldipimp')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) ldipimp(1:3)

          case('lquaimp')
            do i = 1,3
              n = nnombre(itape4,132)
              read(itape4,*,err=9999) lquaimp(i,1:3)
            end do

          case('normaltau')
            normaltau = .true.

          case('noncentre')
            noncentre = .true.
            state_all = .true.

          case('center')
            noncentre = .true.
            state_all = .true.
            n = nnombre(itape4,132)
            select case(n)
              case(0)
                Centre_auto = .true.
              case(1,2)
                Centre_auto = .true.
                read(itape4,*)
              case default
                read(itape4,*,err=9999) Centre(1:3)
            end select 
 
          case('center_ab')
            noncentre = .true.
            State_all = .true.
            n = nnombre(itape4,132)
            select case(n)
              case(0)
                Centre_auto = .true.
                Centre_auto_abs = .true.
              case(1,2)
                Centre_auto = .true.
                Centre_auto_abs = .true.
                read(itape4,*)
              case default
                read(itape4,*,err=9999) Centre(1:3)
            end select 
            Centre_auto = .true.

          case('polarized')
            polarise = .true.
            do ipl = 1,nple
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              select case(n)
                case(3)
                  read(itape4,*,err=9999) polar(:,ipl)
                case(4)
                  read(itape4,*,err=9999) polar(:,ipl), pdpolar(ipl,1)
                case(6)
                  read(itape4,*,err=9999) polar(:,ipl), veconde(:,ipl)
                case(7)
                  read(itape4,*,err=9999) polar(:,ipl), veconde(:,ipl),
     &                           pdpolar(ipl,1)
                  pdpolar(ipl,2) = pdpolar(ipl,1) 
                case(8)
                  read(itape4,*,err=9999) polar(:,ipl), veconde(:,ipl),
     &                           pdpolar(ipl,:)
                case default

                  call write_error
                  do ipr = 6,9,3
                    write(ipr,100)
                    write(ipr,110)
                  end do
                  stop
              end select
            end do

          case('allsite')
            allsite = .true.

          case('step_azim')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) Step_azim
            if( abs( Step_azim ) < eps10 ) then
              nphim = 1
            else
              nphim = nint( 360._db / abs( Step_azim ) )
            endif

          case('symsite')
            symauto = .false.
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) n_atom_proto
            igr = 0
            do ipr = 1,n_atom_proto
              n = nnombre(itape4,132)
              read(itape4,*,err=9999) natomsym 
              write(iscratch,*) natomsym
              do i = 1,natomsym
                n = nnombre(itape4,132)
                igr = igr + 1
                if( n == 1 ) then
                  read(itape4,*,err=9999) isymeq
                  p(:) = 0._db
                else
                  read(itape4,*,err=9999) isymeq, p(:)
                endif
                write(iscratch,*) igr, p(:), isymeq
              end do
            end do

          case('dafs')
            Dafs = .true.
            do ipl = 1,npldafs
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              select case(n)
                case(3)
                  read(itape4,*) hkl_dafs(:,ipl)
                  nn = nnombre(itape4,132)
                  select case(nn)
                    case(6)
                      read(itape4,*,err=9999) poldafse(1:3,ipl),
     &                                        vecdafsem(1:3,ipl)
                      read(itape4,*,err=9999) poldafss(1:3,ipl),
     &                                        vecdafssm(1:3,ipl)
                      angpoldafs(3,ipl) = 10000._db
                    case(9)
                      read(itape4,*,err=9999) (poldafse(i,ipl),
     &                    poldafsei(i,ipl), i = 1,3), vecdafsem(1:3,ipl)
                      read(itape4,*,err=9999) (poldafss(1:3,ipl),
     &                    poldafssi(i,ipl), i = 1,3), vecdafssm(1:3,ipl)
                      angpoldafs(3,ipl) = 10000._db
                    case default
                      call write_error
                      do ipr = 6,9,3
                        write(ipr,120) ipl
                      end do
                      stop
                  end select
                case(5)
                  read(itape4,*,err=9999) hkl_dafs(:,ipl), 
     &                                    isigpi(ipl,1:2)
                  angpoldafs(3,ipl) = - 10000._db
                  do i = 1,2
                    if( isigpi(ipl,i) == 1 ) then
                      angpoldafs(i,ipl) = 0._db
                    elseif( isigpi(ipl,i) == 2 ) then
                      angpoldafs(i,ipl) = 90._db
                    endif
                  end do
                case(6)
                  read(itape4,*,err=9999) hkl_dafs(:,ipl),
     &                                isigpi(ipl,1:2), angpoldafs(3,ipl)
                  do i = 1,2
                    if( isigpi(ipl,i) == 1 ) then
                      angpoldafs(i,ipl) = 0._db
                    elseif( isigpi(ipl,i) == 2 ) then
                      angpoldafs(i,ipl) = 90._db
                    endif
                  end do
                case(7)
                  read(itape4,*,err=9999) hkl_dafs(:,ipl), 
     &                       (isigpi(ipl,i), angpoldafs(i,ipl), i = 1,2)
                  angpoldafs(3,ipl) = - 10000._db
                case(8)
                  read(itape4,*,err=9999) hkl_dafs(:,ipl),
     &               (isigpi(ipl,i), angpoldafs(i,ipl), i = 1,2),
     &                     angpoldafs(3,ipl)
                case default
                  call write_error
                  do ipr = 6,9,3
                     write(ipr,100)
                     write(ipr,120) ipl
                  end do
                  stop
              end select
            end do
            do ipl = 1,npldafs
              do i = 1,2
                if( isigpi(ipl,i) == 3 .or. isigpi(ipl,i) == 4 .or.
     &              isigpi(ipl,i) == 10 ) cycle
                if( abs( angpoldafs(i,ipl) ) < eps10 ) then
                  isigpi(ipl,i) = 1
                elseif( abs( angpoldafs(i,ipl) - 90 ) < eps10 ) then
                  isigpi(ipl,i) = 2
                else
                  isigpi(ipl,i) = 5
                endif   
              end do
              do i = 1,2
                if( isigpi(ipl,i) == 10 ) angpoldafs(i,ipl) = -10000._db
              end do
              nscan = 0
              do i = 1,3
                if( angpoldafs(i,ipl) < - 9999._db ) nscan = nscan + 1
              end do
              if( nscan > 1 ) then
                call write_error
                do ipr = 6,9,3
                   write(ipr,100)
                   write(ipr,122) ipl
                end do
                stop
              endif
            end do

          case('dafs_exp')
            Dafs_bio = .true.
            Dafs = .true.
            do i = 1,3
              read(itape4,*,err=9999) Mat_or(i,:)
            end do
            ipl0 = 0
            do i = 1,n_file_dafs_exp
              read(itape4,*,err=9999) Angle_dafs_exp(i)
              n = nnombre(itape4,132)
              read(itape4,'(A)') Fichier
              Fichier = Adjustl(Fichier)
              l = len_trim(Fichier)
              if( l > 4 ) then
                if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.txt'
              endif 
              open(99, file=Fichier, status='old',iostat=istat)
              File_dafs_exp(i) = Fichier
              n = nnombre(99,100000)
              n = n / 3
              read(99,*,err=1010) (hkl_dafs(:,ipl),ipl=ipl0+1,ipl0+n)
              goto 1020
 1010         call write_err_form(99,grdat)
 1020         Close(99)
              Angle_or(ipl0+1:ipl0+n) = Angle_dafs_exp(i)
              ipl0 = ipl0 + n
            end do

          case('green')
            Green_s = .true.

          case('green_int')
            Green_s = .true.
            Green_int = .true.

          case('zero_azim')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) Vec_orig(:)

          case('solsing')
            solsing_only = .true.
            solsing_s = .true.

          case('no_solsin')
            No_solsing = .true.

          case('norman')
            normrmt = 2

          case('rchimp')
            do it = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              if( it <= ntype ) then
                read(itape4,*,err=9999) rchimp(it)
              else
                read(itape4,*)
              endif
            end do

          case('raydem')
            normrmt = 3

          case('rmtg')
            normrmt = 4
            n1 = 1
            do i = 1,ntype
              n = nnombre(itape4,132)
              n2 = min( n1+n-1, ntype )
              if( n2 >= n1 ) then
                read(itape4,*,err=9999) rmtimp(n1:n2)
                n1 = n2 + 1
              else
                exit
              endif
            end do
            if( n2 < ntype ) normrmt = - n2
            overlap = 0._db

          case('rmtv0')
            normrmt = 5
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) v0bdcFimp(1:min(n,nspin))
            if( nspin > n ) v0bdcFimp(nspin) = v0bdcFimp(1)
            korigimp = .true.
            overlap = 0._db

          case('overlap')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) overlap

          case('muffintin')
            muffintin = .true.

          case('iord')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) iord

          case('adimp','interpoin','inter_poi')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) adimp

          case('rmt')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) rmtt
            rmt(:) = rmtt

          case('lmax')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) lmaxat0

          case('lmaxfree')
            lmaxfree = .true.

          case('lmaxstden')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) lamstdens

          case('lmaxso')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) lmaxso0

          case('chlib')
            Charge_free = .true.

          case('hedin')
            hedin = .true.

          case('perdew')
            perdew = .true.

          case('xalpha')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) alfpot

          case('flapw','flapw_s','flapw_r','flapw_s_p','flapw_psi',
     &         'flapw_n','flapw_n_p','flapw_s_n')
            if( grdat(6:7) == '_n' ) flapw_new = .true.
            if( grdat(6:7) == '_s' ) then
              sauve_potlapw = .true.
              n = nnombre(itape4,132)
              read(itape4,'(A)') file_sauve_potlapw
            elseif( grdat(6:7) == '_r' ) then
              recup_potlapw = .true.
              n = nnombre(itape4,132)
              read(itape4,'(A)') file_sauve_potlapw
            endif
            com(:) = ' Come from FLAPW'
            icom(:) = 3
            n = nnombre(itape4,132)
            read(itape4,'(A)') nomstruct
            if( .not. recup_potlapw ) then
              n = nnombre(itape4,132)
              read(itape4,'(A)') nomvcoul
              if( nomvcoul(1:1) == ' ' ) nomvcoul = adjustl( nomvcoul )
              n = nnombre(itape4,132)
              read(itape4,'(A)') nomr2v
              if( nomr2v(1:1) == ' ' ) nomr2v = adjustl( nomr2v )
              if( flapw_new .and. nspin == 2 )  then
                read(itape4,'(A)') nomr2vdn
                if( nomr2vdn(1:1) == ' ' ) nomr2vdn = adjustl( nomr2vdn)
              endif
              do ispin = 1,2*nspin-1
                n = nnombre(itape4,132)
                read(itape4,'(A)') nomclm(ispin)
                if( nomclm(ispin)(1:1) == ' ' )
     &                   nomclm(ispin) = adjustl( nomclm(ispin) )
              end do
            endif
            n = nnombre(itape4,132)
            if( grdat /= 'flapw_s_p' .and. grdat /= 'flapw_psi' .and.
     &         grdat /= 'flapw_n_p' ) then
              read(itape4,'(A)') nompsii
              if( nompsii(1:1) == ' ' ) nompsii = adjustl( nompsii )
            else
              nompsii = 'dirac'
            endif

          case('delta_en_')
            read(itape4,*,err=9999) Delta_En_conv
         
          case('scf') 
            self_cons = .true.
            n = nnombre(itape4,132)   
            if( n == 1 ) then
              read(itape4,*,err=9999) nself 
            elseif( n == 2 ) then
              read(itape4,*,err=9999) nself, p_self0
            elseif( n > 2 ) then
              read(itape4,*,err=9999) nself, p_self0, Delta_En_conv
            endif 

          case('scf_abs') 
            scf_elecabs = .true.

          case('p_self') 
            read(itape4,*,err=9999) p_self0

          case('n_self') 
            read(itape4,*,err=9999) nself

          case('r_self')
            read(itape4,*,err=9999) r_self
            r_self_imp = .true.

          case('scf_exc')
            self_exc_imp = .true.
            
          case('scf_non_e')
            self_nonexc_imp = .true.
            
          case('nonexc')
            nonexc = .true.

          case('excited')
            exc_imp = .true.
            
          case('no_fermi') 
            fermi_auto = .false.

          case('hubbard')
            n1 = 1
            do i = 1,100
              n = nnombre(itape4,132)
              n2 = min( n1+n-1, ntype)
              if( n2 >= n1 ) then
                read(itape4,*,err=9999) V_hubbard(n1:n2)
                do j = n1,n2   
                 if( abs(V_hubbard(j)) > eps10 ) Hubb(j) = .true.
                end do
                n1 = n2 + 1
              else
                exit
              endif
            end do

          case('full_atom') 
            Full_atom_e = .true.

          case('absorbeur')
            k = 0
            do i = 1,100000
              n = nnombre(itape4,132)
              if( n < 1 ) exit
              read(itape4,*,err=9999) iabsm(k+1:k+n)
              k = k + n
            end do

          case('atom_conf')

            Atom_conf = .true.

            do it = 1,ntype_conf
              read(itape4,*) nb_atom_conf(it),
     &                       igra(1:nb_atom_conf(it),it), nlat(it)
              backspace(itape4)
              n = nnombre(itape4,132)
              if( n == nb_atom_conf(it) + 2 + 3*nlat(it) ) then
                read(itape4,*,err=9999) nb_atom_conf(it),
     &                       igra(1:nb_atom_conf(it),it), nlat(it),
     &                        ( nvval(it,l), lvval(it,l), 
     &                         popval(it,l,1), l = 1,nlat(it) )
                if( nspin == 2 ) then         
                  do l = 1,nlat(it)          
                    popval(it,l,1) = 0.5_db * popval(it,l,1)
                    popval(it,l,2) = popval(it,l,1)
                  end do
                endif
              elseif( n == nb_atom_conf(it) + 2 + 4*nlat(it) ) then
                if( nspin == 1 ) then
                  allocate( x(nlat(it)) )
                  read(itape4,*,err=9999) nb_atom_conf(it),
     &                       igra(1:nb_atom_conf(it),it), nlat(it),
     &                    ( nvval(it,l),
     &                      lvval(it,l), popval(it,l,1), x(l),
     &                      l = 1,nlat(it) )
                  do l = 1,nlat(it)
                    popval(it,l,1) = popval(it,l,1) + x(l)
                  end do
                  deallocate( x )
                else
                  read(itape4,*,err=9999) numat(it), nlat(it),
     &               ( nvval(it,l),
     &               lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
                endif
              else
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,100)
                    write(ipr,125)  it
                  end do
                  stop
              endif  
            end do

! Lecture des densites electroniques
          case('atom')
            atom = .true.
            it = 0
            do jt = 1,100000
              n = nnombre(itape4,132)

              if( n == 3 ) then
                read(itape4,*,err=9999) Ang_base_loc(:,it) 
                n = nnombre(itape4,132)
              endif

              if( n > 0 ) then

                nrm = max( nrm, nrato_dirac )
                if( jt <= ntype ) it = it + 1
                if( n == 1 ) then
                  read(itape4,*,err=9999) numat(it)
                  nlat(it) = 0
                else
                  read(itape4,*,err=9999) numat(it), nlat(it) 
                  if( nlat(it) > 0 ) then
                    backspace(itape4)
                    n = nnombre(itape4,132)
                    if( n == 2 + 3*nlat(it) ) then
                      read(itape4,*,err=9999) numat(it), nlat(it),
     &                              ( nvval(it,l), lvval(it,l), 
     &                               popval(it,l,1), l = 1,nlat(it) )
                      if( nspin == 2 ) then         
                        do l = 1,nlat(it)          
                          popval(it,l,1) = 0.5_db * popval(it,l,1)
                          popval(it,l,2) = popval(it,l,1)
                        end do
                      endif
                    elseif( n == 2 + 4*nlat(it) ) then
                      if( nspin == 1 ) then
                        allocate( x(nlat(it)) )
                        read(itape4,*,err=9999) numat(it),nlat(it),
     &                          ( nvval(it,l),
     &                            lvval(it,l), popval(it,l,1), x(l),
     &                            l = 1,nlat(it) )
                        do l = 1,nlat(it)
                          popval(it,l,1) = popval(it,l,1) + x(l)
                        end do
                        deallocate( x )
                      else
                        read(itape4,*,err=9999) numat(it), nlat(it),
     &                     ( nvval(it,l),
     &                     lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
                      endif
                    else
                      call write_error
                      do ipr = 6,9,3
                        write(ipr,100)
                        write(ipr,130)  it
                      end do
                      stop
                    endif
                  endif
                endif

              else

                read(itape4,'(A)',end=1310) motsb
                open(8, file = motsb, status='old', iostat=istat)
                if( istat /= 0 ) then
                  backspace(itape4)
                  exit
                endif
                if( jt <= ntype .and. jt /= 1 ) it = it + 1
                nomfile_atom(it) = motsb
                read(8,'(A)') com(it)
                icom(it) = 4
                do i = 1,100000
                  read(8,'(A)') mot3
                  if(mot3 == '---') exit
                end do
                read(8,*,err=9999) numat(it), popatc(it), nlat(it)
                backspace(8)
                read(8,*,err=9999) numat(it), popatc(it), nlat(it),
     &         (nvval(it,l), lvval(it,l), popatv(it,l), l = 1,nlat(it))
                read(8,*,err=9999) nrato(it)
                nrm = max( nrm, nrato(it) )
                Close(8)
                if( nlat(it) > 0 ) then
                  n = nnombre(itape4,132)
                  if( n == nlat(it) ) then
                    read(itape4,*,err=9999) popval(it,1:nlat(it),1)
                  else
                    read(itape4,*,err=9999) 
     &                             ( popval(it,l,:), l = 1,nlat(it) )
                  endif
                endif
              endif
            end do

          case('dilatorb')
            do i = 1,norbdil
              n = nnombre(itape4,132)
              read(itape4,*,err=9999) itdil(i), ldil(i), cdil(i)
            end do

          case('v0imp')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) v0bdcFimp(1:min(n,nspin))
            if( nspin > n ) v0bdcFimp(nspin) = v0bdcFimp(1)
            korigimp = .true.

          case('vmax')
            n = nnombre(itape4,132)
            read(itape4,*,err=9999) v_intmax

          case('eimag')
            do ie = 1,neimagent
              n = nnombre(itape4,132)
              if( n == 1 ) then
                read(itape4,*,err=9999) eimagent(ie)
                eeient(ie) = 0._db
              else
                read(itape4,*,err=9999) eeient(ie), eimagent(ie)
              endif
            end do

          case('pointgrou')
            n = nnombre(itape4,132)
            read(itape4,'(A)') motsb
            PointGroup = identmot(motsb,8)
            select case( PointGroup(1:1) )
              Case('_')
                PointGroup(1:1) = '-'
              Case('c')
                PointGroup(1:1) = 'C'
              Case('d')
                PointGroup(1:1) = 'D'
              Case('s')
                PointGroup(1:1) = 'S'
              Case('o')
                PointGroup(1:1) = 'O'
              Case('t')
                PointGroup(1:1) = 'T'
            end select
            PointGroup_Auto = .false.
 
          case('symmol')
            symmol = .true.

          case('spgroup')
            n = nnombre(itape4,132)
            read(itape4,'(A)') motsb
            if( motsb(1:1) == ' ' ) motsb = adjustl(motsb)
            Space_group = motsb(1:10)

! Description de l'agregat :
          case('crystal','molecule','crystal_t','molecule_')
            if( grdat(1:8) == 'molecule') then
              matper = .false.
            else
              matper = .true.
            endif
            n = nnombre(itape4,132)
            if( n == 3 ) then
              read(itape4,*,err=9999) axyz(1:3)
              angxyz(1:3) = 90._db
            elseif( n == 2 ) then
              cylindre = .true.
              read(itape4,*,err=9999) axyz(1:3:2)
              axyz(2) = axyz(1)
              angxyz(1:3) = 90._db
            elseif( n == 1 ) then
              spherique = .true.
              read(itape4,*,err=9999) axyz(1)
              axyz(2) = axyz(1)
              axyz(3) = axyz(1)
              angxyz(1:3) = 90._db
            else
              read(itape4,*,err=9999) axyz(1:3), angxyz(1:3)
            endif
            do igr = 1,ngroup_neq
              if( Temperature .and. Taux ) then
                n = nnombre(itape4,132)
                if( n > 5 ) then
                  read(itape4,*) itype(igr), posn(:,igr), Taux_oc(igr),
     &                           Temp_coef(igr)
                else
                  read(itape4,*) itype(igr), posn(:,igr)
                endif
              elseif( Temperature ) then
                n = nnombre(itape4,132)
                if( n > 4 ) then
                  read(itape4,*) itype(igr), posn(:,igr), Temp_coef(igr)
                else
                  read(itape4,*) itype(igr), posn(:,igr)
                endif
              elseif( Taux ) then
                n = nnombre(itape4,132)
                if( n > 4 ) then
                  read(itape4,*) itype(igr), posn(:,igr), Taux_oc(igr)
                else
                  read(itape4,*) itype(igr), posn(:,igr)
                endif
              elseif( Readfast .or. .not.( Atom_nonsph
     &                          .or. Atom_occ_hubb .or. Axe_loc) ) then
                read(itape4,*) itype(igr), posn(:,igr)
              else
                n = nnombre(itape4,132)
                if( n == 0 ) exit
                if( n == 2 .or. n == 3 ) then
                  read(itape4,*) Ang_base_loc_gr(1:n,igr)
                  Ang_base_loc_gr(1:n,igr) = Ang_base_loc_gr(1:n,igr)
     &                                     * pi / 180
                  n = nnombre(itape4,132)
                endif
                select case(n)
                  case(4)
                    read(itape4,*) itype(igr), posn(:,igr)
                    norb = 0
                  case(5)
                    read(itape4,*) itype(igr), posn(:,igr), norb
                end select
                if( norb == 0 ) cycle
                if( norb /= -1 ) then          
                  if( norb < 0 ) then
                    norbv(igr) = - norb - 1
                  else 
                    norbv(igr) = norb
                  endif
                  do io = 1,norbv(igr)
                    hybrid(io,:,igr) = 0._db
                    n = nnombre(itape4,132)
                    allocate( Hyb(n-1) )
                    read(itape4,*) Hyb(1:n-1), pop_nonsph(io,igr)
                    select case(n)
                      Case(5,6,8,17)
                        Hybrid(io,1:n-1,igr)
     &                      = cmplx( Hyb(1:n-1), 0._db, db )
                      Case(9,11,15,33)
                        do i = 1,n-2,2
                          Hybrid(io,(i+1)/2,igr)
     &                      = cmplx( Hyb(i), Hyb(i+1), db )
                        end do
                    end select
                    deallocate ( Hyb )
                  end do
                endif 
                if( norb < 0 ) then
                  n = nnombre(itape4,132)
                  l = (n/nspin - 1)/2
                  read(itape4,*) 
     &             ((occ_hubb_e(m,m,isp,igr), m = -l,l ), isp = 1,nspin)
                endif 
              endif
            end do

          case('crystal_p')

            Pdb = .true.
            Fichier = ' '
            read(itape4,'(A)') Fichier
            Fichier = Adjustl(Fichier)
            l = len_trim(Fichier)
            if( l > 4 ) then
              if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.pdb'
            endif 
            open(8, file = Fichier, status='old', iostat=istat) 
            if( istat /= 0 ) call write_open_error(Fichier,istat,1)

            Matper = .true.

            igr = 0

            do ligne = 1,100000
            
              read(8,'(a6)',end=1030,err=1030) mot6

              select case(mot6)

                case('END   ')

                  exit

                case('SCALE1' )

                  backspace(8)
                  do i = 1,3
                    read(8,'(10x,3f10.6)') Mat(i,:)
                  end do

                case('CRYST1' )

                  backspace(8) 
                  read(8,'(6x,3f9.3,3f7.2,a13)') axyz(1:3), angxyz(1:3),
     &                                           Spgr

                  l = len_trim(Spgr)
                  j = 0
                  do i = 1,l
                    if( Spgr(i:i) == ' ' ) cycle
                    j = j + 1
                    Space_group(j:j) = Spgr(i:i)
                  end do 

                case('ATOM  ','HETATM') 

                  backspace(8)
                  read(8,'(16x,a1,13x,3f8.3,2f6.2,10x,a2)') Let ,p(:),
     &                                           t, tc,  Symbol
                  read(8,*)

                  igr = igr + 1
                  p = Matmul( Mat, p )
                  posn(:,igr) = p(:)
                  if( Taux ) Taux_oc(igr) = t    
                  if( Temperature ) Temp_coef(igr) = tc    

                  select case(Let)
                    case(' ')
                      Kgroup(igr) = 0
                    case('A')
                      Kgroup(igr) = 1
                    case('B')
                      Kgroup(igr) = 2
                    case('C')
                      Kgroup(igr) = 3
                    case('D')
                      Kgroup(igr) = 4
                    case default
                      Kgroup(igr) = 5
                  end select

                  Symbol = adjustl( Symbol )              
                  do i = 1,103
                    if( Chemical_Symbol_c(i) /= Symbol ) cycle
                    itype(igr) = i
                    exit
                  end do
                  if( i == 104 ) then
                    call write_error
                    do ipr = 6,9,3
                      write(ipr,140) igr, Symbol
                    end do
                    stop
                  endif
 
                case default

                  cycle

                end select

            end do

 1030       continue

            Close(8)

          case('dpos')
            n = nnombre(itape4,132)
            read(itape4,*) dpos(1:3)

          case('single_pr')
            single_prec = .true.

          case('test_dist')
            n = nnombre(itape4,132)
            read(itape4,*) Test_dist_min

          case('ylm_comp')
            Ylm_comp_inp = .true.

          case('delta_eps')
            n = nnombre(itape4,132)
            read(itape4,*) Delta_Epsii

          case('one_run')
            One_run = .true.
            State_all = .true.

          case('z_absorbe')
            n = nnombre(itape4,132)
            if( Absauto ) then
              read(itape4,*) numat_abs
            else 
              read(itape4,*)
            endif

          case('d_max_pot')
            n = nnombre(itape4,132)
            read(itape4,*) D_max_pot

! Parametres deja lus dans lectdim
          case('full_self','magnetism','memory_sa','self_abs',
     &         'readfast','temperatu')

          case default

            if( igrdat == 1 ) then
              comt = mots
            elseif( grdat(1:1) /= ' ' ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,100)
                write(ipr,150) mots
              end do
              stop
            endif

        end select

      end do
 1310 continue

! Fin de la lecture.

      l = len_trim(nomfich)
      write(6,'(/a9,A)') ' Filout: ',nomfich(1:l)
      if( l > 4 ) then
        if( nomfich(l-3:l) == '.txt' .or. nomfich(l-3:l) == '.dat')
     &    nomfich(l-3:l) = '    '
      endif
      nomfichbav = nomfich
      long = len_trim(nomfich)
      nomfichbav(long+1:long+8) = '_bav.txt'

      do i = 1,28
        if( icheck(i) > 1 ) exit
      end do
      if( Extract .and. i > 28 ) icheck(1:28) = 0
  
      i = sum( icheck(1:28) )
      if( i > 0 ) then
        open(3, file = nomfichbav, status='unknown',iostat=istat)
        if( istat /= 0 ) call write_open_error(nomfichbav,istat,1)
      endif

      if( icheck(1) > 0 ) then
        write(3,'(A/A/A)') Revision, com_date, com_time
        if( comt /= ' ') write(3,'(/A)') comt
        ipr0 = 3
      else
        ipr0 = 6
      endif

! Modification en cas de fit.
      if( fit_cal ) then
        do igr = 2,ngroup_par
          istop = 0
          do ipar = 1,npar(igr)
            if( typepar(igr,ipar) /= 'dposx' .and.
     &          typepar(igr,ipar) /= 'dposy' .and.
     &          typepar(igr,ipar) /= 'dposz' .and.
     &          typepar(igr,ipar) /= 'posx' .and.
     &          typepar(igr,ipar) /= 'posy' .and.
     &          typepar(igr,ipar) /= 'posz' .and.
     &          typepar(igr,ipar) /= 'theta' .and.
     &          typepar(igr,ipar) /= 'phi'   ) cycle
            if( indice_par(igr,ipar) > ngroup ) then
              call write_error
              do ipr = ipr0,9,3
                write(ipr,100)
                write(ipr,160) typepar(igr,ipar), indice_par(igr,ipar),
     &                         ngroup
              end do
              istop = 1
             endif
          end do
        end do
        if( istop == 1 ) stop

        do i = 2,ngroup_par
          do ip = 1,npar(i)
            select case( typepar(i,ip) )
              case('dposx')
                posn(1,indice_par(i,ip)) = posn(1,indice_par(i,ip)) 
     &                                   + param(i,ip)
              case('dposy')
                posn(2,indice_par(i,ip)) = posn(2,indice_par(i,ip))
     &                                   + param(i,ip)
              case('dposz')
                posn(3,indice_par(i,ip)) = posn(3,indice_par(i,ip))
     &                                   + param(i,ip)
              case('posx')
                posn(1,indice_par(i,ip)) = param(i,ip)
              case('posy','theta')
                posn(2,indice_par(i,ip)) = param(i,ip)
              case('posz','phi')
                posn(3,indice_par(i,ip)) = param(i,ip)
              case('abc')
                axyz(1:3) = axyz(1:3) * (1 + 0.01*param(i,ip))
              case('a')
                axyz(1) = axyz(1) * (1 + 0.01*param(i,ip))
              case('b')
                axyz(2) = axyz(2) * (1 + 0.01*param(i,ip))
              case('c')
                axyz(3) = axyz(3) * (1 + 0.01*param(i,ip))
              case('angx','anga')
                angxyz(1) = param(i,ip)
              case('angy','angb')
                angxyz(2) = param(i,ip)
              case('angz','angc')
                angxyz(3) = param(i,ip)
              case('poporb')
                io = 0
                boucle_it: do it = 1,ntype
                  do l = 1,nlat(it)
                    do ispin = 1,nspin
                      io = io + 1
                      if( io /= indice_par(i,ip) ) cycle 
                      popval(it,l,ispin) = param(i,ip)
                      exit boucle_it
                    end do
                  end do
                end do boucle_it
                if( it > ntype ) then
                  call write_error
                  do ipr = ipr0,9,3
                    write(ipr,100)
                    write(ipr,170) typepar(i,ip), indice_par(i,ip)
                  end do
                  stop 
                endif
            end select
          end do
        end do
      endif

      do igr = 1,ngroup
        if( cylindre ) then
          r = posn(1,igr)
          theta = pi * posn(2,igr) / 180 
          posn(1,igr) = r * cos( theta )
          posn(2,igr) = r * sin( theta )
        elseif( spherique ) then
          r = posn(1,igr)
          theta = pi * posn(2,igr) / 180 
          phi = pi * posn(3,igr) / 180 
          posn(1,igr) = r * sin( theta ) * cos( phi)
          posn(2,igr) = r * sin( theta ) * sin( phi)
          posn(3,igr) = r * cos( theta )
        endif
      end do

      if( Atom_conf ) then

        do it = 1,ntype_conf
          numat(it) = itype(igra(1,it))
        end do 

        allocate( itZ(103) )

        jt = ntype_conf

        boucle_2: do igr = 1,ngroup_neq
          do it = 1,ntype_conf
            do kgr = 1,nb_atom_conf(it)
              if( igr == igra(kgr,it) ) cycle boucle_2
            end do
          end do
          boucle_3: do jgr = 1,igr-1
            do it = 1,ntype_conf
              do kgr = 1,nb_atom_conf(it)
                if( jgr == igra(kgr,it) ) cycle boucle_3
              end do
            end do
            if( itype(igr) == itype(jgr) ) cycle boucle_2
          end do boucle_3
          jt = jt + 1
          numat(jt) = itype(igr)
          itZ( numat(jt) ) = jt
        end do boucle_2

        boucle_igr: do igr = 1,ngroup_neq
          do it = 1,ntype_conf
            do kgr = 1,nb_atom_conf(it)
              if( igr /= igra(kgr,it) ) cycle
              itype(igr) = it
              cycle boucle_igr
            end do
          end do
          itype(igr) = itZ( itype(igr) )
        end do boucle_igr

        deallocate( itZ )

        Atom = .true.

      endif

      iabsorig(:) = iabsm(:)
      if( Space_group /= ' ' .and. matper ) then
        allocate( neq(ngroup_neq) )
        allocate( pos(3,ngroup_neq) )
        do igr = 1,ngroup_neq
          pos(:,igr) = posn(:,igr) 
        end do
        call spgroup(1,neq,ngroup,ngroup_neq,pos,posn,Space_file,
     &               space_group)
        ia = ngroup + 1
        do igr = ngroup_neq,1,-1
          do i = 1,neq(igr)
            ia = ia - 1
            itype(ia) = itype(igr)
            if( Taux ) Taux_oc(ia) = Taux_oc(igr)
            if( Temperature ) Temp_coef(ia) = Temp_coef(igr)
            if( ngroup_pdb > 0 ) Kgroup(ia) = Kgroup(igr)

            Ang_base_loc_gr(:,ia) = Ang_base_loc_gr(:,igr )

            if( ia > ngroup_nonsph ) cycle

            norbv(ia) = norbv(igr)
            if( norbv(ia) /= 0 ) then
              hybrid(:,:,ia) = hybrid(:,:,igr)
              pop_nonsph(:,ia) = pop_nonsph(:,igr)
            endif
            if( Atom_occ_hubb ) occ_hubb_e(:,:,:,ia)
     &                             = occ_hubb_e(:,:,:,igr)

          end do
          do multi_run = 1,n_multi_run_e
            if( iabsorig(multi_run) == igr ) iabsm(multi_run) = ia
          end do
        end do
        deallocate( pos )
        deallocate( neq )
      endif

      if( flapw ) then
        nrm = nrato_dirac
        call lect_struct_lapw(angxyz,axyz,icheck(1),its_lapw,itype,
     &      ngroup,ngroup_lapw,nomstruct,nrato_lapw,nrm,ntype,
     &      numat,posn,r0_lapw,rlapw,rotloc_lapw)
        if( normrmt == 1 ) normrmt = 2
      elseif( nrm == 0 ) then
        nrm = nrato_dirac
      endif

      if( Pdb ) One_run = .false.

      if( Atom ) then
        istop = 0
        do igr = 1,ngroup
          if( itype(igr) <= ntype ) cycle
          if( istop == 0 ) then
            call write_error
            do ipr = ipr0,9,3
              write(ipr,100)
            end do
          endif
          do ipr = ipr0,9,3
            write(ipr,172) igr, itype(igr), ntype
          end do
          istop = 1
        end do
        if( istop == 1 ) then 
          do ipr = ipr0,9,3
            write(ipr,174)
          end do
          stop
        endif
      endif

      if( numat_abs /= 0 ) then
        do igr = 1,ngroup
          if( Atom .or. Flapw ) then
            if( numat( abs( itype( igr ) ) ) == numat_abs ) exit
          else
            if( itype( igr )  == numat_abs ) exit
          endif
        end do
        if( igr > ngroup ) then
          call write_error
          do ipr = ipr0,9,3
            write(ipr,175) numat_abs
          end do
          stop
        endif        
      else
        if( Atom .or. Flapw ) then
          numat_abs = numat( abs( itype( iabsm(1) ) ) )
        else
          numat_abs = itype( iabsm(1) )
        endif
      endif

      if( clementi ) then
        nompsii = 'clementi'
        nrm = nrato_dirac
        Fermi_auto = .false.
      endif

      if( Atom_occ_hubb ) Fermi_auto = .false. 
      if( Flapw .or. Extract .or. .not. Fermi_auto ) then
        Self_cons = .false. 
        Fermi_auto = .false.
        nself = 0
      end if
      if( .not. r_self_imp ) r_self = rsorte_s
      if( Self_cons ) then
        if( nself == 0 ) nself = 30
      elseif( ( Fermi_auto .or. ( Hubbard .and. .not. Atom_occ_hubb) )
     &          .and. nself == 0 ) then
        Self_cons = .true.
        p_self0 = 0._db 
        nself = 1
        if( .not. r_self_imp ) r_self = min( rsorte_s, 3.5_db )
      endif

      if( Quadrupole ) then
        if( .not. no_dipquad ) E1E2e = .true.
        if( .not. no_e2e2 ) E2E2 = .true.
      endif
      if( Octupole .and. .not. no_e1e3 ) E1E3 = .true.
      if( Octupole .and. .not. no_e3e3 ) E3E3 = .true.
      if( Dipmag ) then
        E1M1 = .true.
        M1M1 = .true.
      endif
      if( E1E2e .or. E2E2 ) Quadrupole = .true.
      if( E1M1 .or. M1M1 ) Dipmag = .true.
      if( E1M2 .or. M2M2 ) Quadmag = .true.
      if( E1E3 .or. E3E3 ) Octupole = .true.

      if( Recup_tddft_data ) Save_tddft_data = .false.

      if( Extract ) then
        Density = .false.
        state_all = .false.
        state_all_out = .false.

        open(1, file = nom_fich_Extract, status='old', iostat=istat) 
        if( istat /= 0 )
     &    call write_open_error(nom_fich_Extract,istat,1)
        do i = 1,100000
          read(1,'(A)') mots
          if( mots(2:10) == 'Threshold' ) then
            l = len_trim(mots)
            if( mots(l:l) == 'e' ) then
              seuil = mots(l-7:l-5)
              seuil = adjustl( seuil ) 
            else 
              seuil = mots(14:16)
            endif
            exit
          endif
          if( Seuil == 'Opt' ) Optic = .true.
        end do
        if( Quadrupole ) then
          do i = 1,100000
            read(1,'(A)') mots
            if( mots(2:6) == 'Dipol'  ) then
              read(1,'(A)') mots
              if( mots(2:6) /= 'Quadr' ) Quadrupole = .false.
              read(1,'(A)') mots
              if( mots(2:6) /= 'Octup' ) octupole = .false.
              exit
            endif
          end do
        endif
        Rewind(1)
        Core_resolved = .false.
        do i = 1,100
          read(1,'(A)') mots
          if( mots(2:5) /= 'Core'  ) cycle
          Core_resolved = .true.
          exit
        end do
        Close(1)
      endif

      if( Spinorbite .and. non_relat == 0 ) relativiste = .true.

      if( .not. flapw .and. abs(alfpot) < eps6 .and. .not. perdew )
     &  hedin = .true.
      if( hedin ) then
        alfpot = 0._db
      elseif( flapw ) then
        alfpot = 0.33
      elseif( perdew ) then
        alfpot = -1._db
      endif

      if( Optic ) Seuil = 'Opt'

      if( seuil /= 'K1' .and. seuil /= 'L1' .and. seuil /= 'L2' .and.
     &    seuil /= 'L3' .and. seuil /= 'M1' .and. seuil /= 'M2' .and.
     &    seuil /= 'M3' .and. seuil /= 'M4' .and. seuil /= 'M5' .and.
     &    seuil /= 'N1' .and. seuil /= 'N2' .and. seuil /= 'N3' .and.
     &    seuil /= 'N4' .and. seuil /= 'N5' .and. seuil /= 'N6' .and.
     &    seuil /= 'N7' .and. seuil /= 'O1' .and. seuil /= 'O2' .and.
     &    seuil /= 'O3' .and. seuil /= 'O4' .and. seuil /= 'O5' .and.
     &    seuil /= 'P1' .and. seuil /= 'P2' .and. seuil /= 'P3' .and.
     &    seuil /= 'L23' .and. seuil /= 'M23' .and. seuil /= 'M45' .and.
     &    seuil /= 'N23' .and. seuil /= 'N45' .and. seuil /= 'N67' .and.
     &    seuil /= 'O23' .and. seuil /= 'O45' .and. seuil /= 'P23' .and.
     &    seuil /= 'Opt' ) then
        call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,180) seuil
        end do
        stop
      endif

      if( Optic ) then
        nseuil = 0; lseuil = 0; jseuil = 0; nbseuil = 1
      else
        select case( seuil(1:1) )
          case('K')
            nseuil = 1
          case('L')
            nseuil = 2
          case('M')
            nseuil = 3
          case('N')
            nseuil = 4
          case('O')
            nseuil = 5
          case('P')
            nseuil = 6
        end select
        select case( seuil(2:2) )
          case('1')
            lseuil = 0; jseuil = 1
          case('2')
            lseuil = 1; jseuil = 2
          case('3')
            lseuil = 1; jseuil = 3
          case('4')
            lseuil = 2; jseuil = 4
          case('5')
            lseuil = 2; jseuil = 5
          case('6')
            lseuil = 3; jseuil = 6
          case('7')
            lseuil = 3; jseuil = 7
        end select
        nbseuil = len_trim( seuil ) - 1
      endif

      if( Flapw .or. lseuil > 0 .or. Optic ) nonexc = .true.
      if( exc_imp .and. .not. flapw ) nonexc = .false. 
      if( nonexc ) symmol = .true.
      if( self_nonexc_imp ) self_nonexc = .true.
      if( self_exc_imp .and. .not. Optic ) self_nonexc = .false.
      if( nonexc ) self_nonexc = .true.
      l = l_level_val(numat_abs)
      if( .not. ( nonexc .and. self_nonexc ) 
     &      .and. ( l == 2 .or. l == 3 ) ) scf_elecabs = .true.

      if( .not. flapw ) then

        if( .not. atom ) then
! Dans ce cas itype est pour l'instant le numero atomique
          jt = 0
          boucle_1: do igr = 1,ngroup
            n = itype(igr)
            do it = 1,jt
              if( numat(it) /= n ) cycle
              itype(igr) = it
              cycle boucle_1
            end do
            jt = jt + 1
            itype(igr) = jt
            numat(jt) = n
          end do boucle_1
          nlat(1:jt) = 0
        endif

        do igr = 1,ngroup
          it = abs( itype(igr) )
          do l = 1,nlat(it)
            do ispin = 1,nspin
              popats(igr,l,ispin) = popval(it,l,ispin)
            end do
          end do
        end do

      endif

      do it = 1,ntype
        if( numat(it) == 0 ) then
          com(it) = ' Empty sphere'
          icom(it) = 5
        endif
      end do

      if( Extract ) then
        if( no_core_resolved ) Core_resolved = .false.
      elseif( ( nspin == 1 .or. no_core_resolved )
     &                     .and. .not. Core_resolved_e ) then
        Core_resolved = .false.
      else
        Core_resolved = .true. 
      endif
! S'il y a plusieurs seuils, le Gamma doit etre pris dans le calcul de
! Chi_0.  
      if( nbseuil > 1 .or. ngamh > 1 ) then
        if( Gamma_max > eps10 ) Gamma_tddft = .true.
        if( Gamma_hole_imp ) then
          g1 = Gamma_hole(1) 
          g2 = Gamma_hole(1) 
          do i = 2,ngamh
            g1 = min( g1, Gamma_hole(i) )
            g2 = max( g2, Gamma_hole(i) )
          end do
          if( abs( g1 - g2 ) > eps10 ) Gamma_tddft = .true.
        else
          if( lseuil == 1 .and. nseuil == 2 ) Gamma_tddft = .true.          
        endif
      endif 

! Verification des entrees :

      if( .not. ( E1E1 .or. E1E2e .or. E1E3 .or. E2E2 .or. E1M1
     &           .or. E1M2 .or. M1M1 .or. M1M2 .or. M2M2 ) ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,190)
        end do
        istop = 1
      endif

      if( Bormann ) then
        Dafs = .true.
        Angpoldafs(3,:) = Ang_borm
        do ipl = 1,npldafs/4
         hkl_dafs(:,ipl) = hkl_borm(:)
        end do            
        do ipl = npldafs/4 + 1,3*npldafs/4
          hkl_dafs(:,ipl) = 0
        end do            
        do ipl = 3*npldafs/4 + 1,npldafs
          hkl_dafs(:,ipl) = - hkl_borm(:)
        end do            
      endif

      if( ngroup == 0 ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,200)
        end do
        istop = 1
      endif
      if( ntype == 0 .and. ngroup /= 0 ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,210)
        end do
        istop = 1
      endif
      if( iord /= 2 .and. iord /= 4 ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,220) iord
        end do
        istop = 1
      endif
      do igr = 1,ngroup
        if( abs(itype(igr)) <= ntype ) cycle
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,230) igr, itype(igr), ntype
        end do
        istop = 1
      end do
      do ipl = 1,nple
        pp = sum( polar(1:3,ipl)**2 )
        q = sum( veconde(1:3,ipl)**2 )
        if( pp < eps10 .and. q < eps10 ) then
          if( istop == 0 ) call write_error
          do ipr = ipr0,9,3
            write(ipr,240)
          end do
          istop = 1
        endif
      end do
      if( normrmt < 0 ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,250) -normrmt, ntype
        end do
        istop = 1
      endif
      if( Full_self_abs ) then
        do ipl = 1,npldafs
          jpl = mod(ipl-1,4) + 1
          select case(jpl)
            case(1)
              if( isigpi(ipl,1) == 1 .and. isigpi(ipl,2) == 1 ) cycle 
            case(2)
              if( isigpi(ipl,1) == 1 .and. isigpi(ipl,2) == 2 ) cycle 
            case(3)
              if( isigpi(ipl,1) == 2 .and. isigpi(ipl,2) == 1 ) cycle 
            case(4)
              if( isigpi(ipl,1) == 2 .and. isigpi(ipl,2) == 2 ) cycle
          end select
          if( istop == 0 ) call write_error
          do ipr = ipr0,9,3
            write(ipr,100)
            write(ipr,260) ipl, hkl_dafs(:,ipl), isigpi(ipl,:)
          end do
          istop = 1
          exit
        end do
      endif

      if( istop == 1 ) stop

! Normalisation des vecteurs hybridations :
      if( Atom_nonsph ) then
        do igr = 1,ngroup
          if( norbv(igr) == 0 ) cycle
          do io = 1,norbv(igr)
            rn = sqrt( sum( abs(hybrid(io,:,igr))**2 ) )
            hybrid(io,:,igr) = hybrid(io,:,igr) / rn
          end do
        end do
      endif

      if( Tddft ) eimagent(:) = 0._db
      if( Spinorbite ) then
        tddft_so = .false.
        Ylm_comp_inp = .true.
      endif

      if( lin_gam == - 1 ) then            
        lin_gam = 1 
        egamme(1) = -5.0_db; egamme(2) =  0.5_db;
        egamme(3) = egamme(1) + ( nenerg - 1 ) * egamme(2)
      endif
      if( .not. green_s .and. neimagent > 0 ) eimagent(:) = 0._db

! Ecriture des entrees :
      if( icheck(1) > 0 ) then

        if( Extract ) then
          write(3,300) nom_fich_Extract
          write(6,300) nom_fich_Extract
          open(1, file = nom_fich_Extract, status='old',iostat=istat) 
          if( istat /= 0 )
     &      call write_open_error(nom_fich_Extract,istat,1)
          do i = 1,100000
            read(1,'(A)') mots
            write(3,'(A)') mots
            if( mots(2:6) == 'Dipol'  ) exit
          end do
        else
          mot13 = Chemical_Name(numat_abs) 
          l1 = len_trim(mot13)
          l2 = len_trim(seuil)
          mots = ' '
          mots =' Threshold:'
          mots(13:12+l1) = mot13(1:l1)
          mots(14+l1:13+l1+l2) = seuil(1:l2)
          mots(15+l1+l2:18+l1+l2) = 'edge'
          write(3,'(/A)') mots
          write(6,'(A)') mots(12:18+l1+l2)
          write(3,320) rsorte_s
          if( .not. green_s .and. overad ) write(3,325) roverad
          write(3,330) icheck(:)
          if( lin_gam == 1 ) write(3,340)
          write(3,350) egamme(1:ngamme)
          write(3,360)
        endif
        if( Quadrupole ) write(3,365)
        if( Octupole ) write(3,370)
        if( Dipmag ) write(3,371)
        if( .not. E1E1 ) write(3,372)
        if( no_e2e2 ) write(3,373)
        if( no_e1e3 ) write(3,374)
        if( no_e3e3 ) write(3,375)
        if( no_dipquad ) write(3,376)
        if( tddft ) then
          if( rpalf ) then
            write(3,377)
          else
            write(3,378)
          endif
          if( Gamma_tddft ) then
            write(3,379)
          else
            write(3,381)
          endif
        endif
        if( Core_resolved ) write(3,382)

        if( Extract ) then
          do i = 1,100000
            read(1,'(A)') mots
            if( mots(2:6) /= 'Relat' .and. mots(6:10) /= 'relat' ) cycle
            write(3,'(A)') mots
            exit          
          end do
          do i = 1,100000
            read(1,'(A)') mots
            if(  mots(2:7) == 'ngroup' .or. mots(2:6) == 'XANES'
     &         .or. mots(2:5) == 'DAFS' ) then
              backspace(1)
              exit
            endif
            write(3,'(A)') mots
          end do
        else
          if( relativiste ) then
            write(3,384)
          else
            write(3,385)
          endif
          if( Spinorbite ) then
            write(3,390)
          elseif( Magnetic ) then
            write(3,400)
          else
            write(3,410)
          endif
          if( Z_nospinorbite /= 0 ) write(3,416) Z_nospinorbite
          if( lmoins1 ) write(3,419)
          if( lplus1 ) write(3,420)
          if( basereel ) then
            write(3,421)
          else
            write(3,422)
          endif
          if( green_s ) then
            write(3,425)
            if( Green_int ) write(3,426)
            if( supermuf ) write(3,427)
            if( nchemin /= - 1 ) write(3,430) nchemin
            if( normaltau ) write(3,435)
            select case(normrmt)
              case(1)
                write(3,436)
              case(2)
                write(3,437)
              case(3)
                write(3,438)
              case(4)
                write(3,439) rmtimp(1:ntype)
              case(5)
                write(3,440)
            end select
            if( abs(overlap) > eps6 ) write(3,445) overlap
            if( lmaxfree ) then
              write(3,452)
            else
              write(3,454)
            endif
          else
            write(3,460) iord, adimp
            write(3,470) lmaxso0
            if( muffintin ) write(3,492)
            if( rydberg ) write(3,494) rrydb
            if( noncentre ) write(3,496)
            if( sum( abs(Centre(:)) ) > epspos ) write(3,497) Centre(:)
            if( .not. Eneg_i ) write(3,498) eclie
            if( v_intmax < 100000._db ) write(3,499) v_intmax
          endif
        endif
        if( Atom_nonsph ) write(3,500)
        if( single_prec ) then
          write(3,501)
        else
          write(3,502)
        endif
        if( temp > eps10 ) write(3,503) temp
        if( Temperature ) write(3,504) 

        if( nple > 0 ) then
          if( sum( abs( pdpolar(:,:) ) ) > eps10 ) then
            write(3,505)
            do ipl = 1,nple
              write(3,509) polar(1:3,ipl), veconde(1:3,ipl),
     &                     pdpolar(ipl,:)
            end do
          else
            write(3,510)
            do ipl = 1,nple
              write(3,509) polar(1:3,ipl), veconde(1:3,ipl)
            end do
          endif
        endif

        if( Dafs_bio ) then

          write(3,'(1x,A)')
     &      'DAFS reflections selected using the experimental files:'
          do i = 1,n_file_dafs_exp
            write(3,'(3x,A)') File_dafs_exp(i)
            write(3,'(5x,a12,f10.3)') 'with angle =', Angle_dafs_exp(i)
          end do

          write(3,'(A)') ' Experimental orientation matrix:'
          do i = 1,3
            write(3,'(3x,3f12.6)') Mat_or(i,:)
          end do

          Write(3,'(a23,i6)') ' Number of reflexions =', npldafs

        elseif( Bormann ) then

          write(3,511) hkl_dafs(:,1), Angpoldafs(3,1)

        else

          do ipl = 1,npldafs
            motpol = ' '
            l = 0
            do i = 1,2
              if( isigpi(ipl,i) == 1 ) then
                motpol(l+1:l+5) = 'sigma'
              elseif( isigpi(ipl,i) == 2 ) then
                motpol(l+1:l+2) = 'pi'
              elseif( isigpi(ipl,i) == 3 ) then
                motpol(l+1:l+5) = 'right'
              elseif( isigpi(ipl,i) == 4 ) then
                motpol(l+1:l+4) = 'left'
              elseif( isigpi(ipl,i) == 5 .or. isigpi(ipl,i) == 10 ) then
                motpol(l+1:l+5) = 'recti'
              endif
              l = l + 6
            end do
 
            if( angpoldafs(3,ipl) < -9999._db ) then
              if( ipl == 1 ) write(3,512)
              write(3,514) hkl_dafs(:,ipl), motpol,
     &                     angpoldafs(1:2,ipl)
            elseif( isigpi(ipl,1) == 10 ) then
              if( ipl == 1 ) write(3,512)
              write(3,517) hkl_dafs(:,ipl), motpol,
     &                     angpoldafs(2:3,ipl)
            elseif( isigpi(ipl,2) == 10 ) then
              if( ipl == 1 ) write(3,512)
              write(3,518) hkl_dafs(:,ipl), motpol,
     &                     angpoldafs(1:3:2,ipl)
            elseif( angpoldafs(3,ipl) < 9999._db ) then
              if( ipl == 1 ) write(3,512)
              write(3,519) hkl_dafs(:,ipl), motpol,
     &                     angpoldafs(1:3,ipl)
            else
              if( ipl == 1 ) write(3,520)
              write(3,530) poldafse(1:3,ipl), vecdafsem(1:3,ipl)
              write(3,540) poldafss(1:3,ipl), vecdafssm(1:3,ipl)
            endif
          end do

        endif

        if( Self_abs ) write(3,522)
        if( Full_self_abs ) write(3,523)

        if( Extract ) then

          Close(1)

        else

          if( matper ) then
            write(3,'(/A)') ' Crystal'
          else
            write(3,'(/A)') ' Molecule'
          endif
          write(3,542) ngroup, ntype

          if( nlatm > 0 ) then
            write(3,'(A)') '   Typ  Z   n  l  Popul.'
            do it = 1,ntype
              if( nspin == 1 ) then
                write(3,543) it, numat(it), ( nvval(it,l), lvval(it,l),
     &                               popval(it,l,:), l = 1,nlat(it) ) 
              else
                write(3,544) it, numat(it), ( nvval(it,l), lvval(it,l),
     &                               popval(it,l,:), l = 1,nlat(it) ) 
              endif
            end do
          endif

          if( abs( dpos(1) ) > epspos .or. abs( dpos(2) ) > epspos .or.
     &        abs( dpos(3) ) > epspos ) write(3,545) dpos(1:3)
          if( norbdil > 0 ) write(3,548)
          do io = 1,norbdil
            write(3,549) itdil(io), ldil(io), cdil(io)
          end do
          if( nonexc ) write(3,550)
          if( .not. PointGroup_Auto ) write(3,552) PointGroup
          if( symmol ) write(3,553)
          if( Space_group /= ' ' .and. matper ) write(3,554) Space_group

          write(3,555) axyz(1:3)
          write(3,556) angxyz(1:3)
          if( Pdb .and. Taux .and. Temperature ) then
            write(3,557)
          elseif( Pdb .and. Taux ) then
            write(3,558)
          elseif( Pdb .and. Temperature ) then
            write(3,559)
          elseif( Pdb ) then
            write(3,560)
          elseif( Taux .and. Temperature ) then
            write(3,561)
          elseif( Taux ) then
            write(3,562)
          elseif( Temperature ) then
            write(3,563)
          elseif( Atom_nonsph ) then
            write(3,564)
          else
            write(3,565)
          endif
          do igr = 1,ngroup
            it = abs(itype(igr))
            if( Pdb .and. Taux .and. Temperature ) then
              write(3,570) numat(it), it, posn(:,igr), Kgroup(igr),
     &                     Taux_oc(igr), Temp_coef(igr) 
            elseif( Pdb .and. Taux ) then
              write(3,570) numat(it), it, posn(:,igr), Kgroup(igr),
     &                     Taux_oc(igr) 
            elseif( Pdb .and. Temperature ) then
              write(3,570) numat(it), it, posn(:,igr), Kgroup(igr),
     &                     Temp_coef(igr) 
            elseif( Pdb ) then
              write(3,570) numat(it), it, posn(:,igr), Kgroup(igr) 
            elseif( Taux .and. Temperature ) then
              write(3,572) numat(it), it, posn(:,igr), Taux_oc(igr),
     &              Temp_coef(igr) 
            elseif( Taux ) then
              write(3,572) numat(it), it, posn(:,igr), Taux_oc(igr) 
            elseif( Temperature ) then
              write(3,572) numat(it), it, posn(:,igr), Temp_coef(igr) 
            elseif( .not. Atom_nonsph ) then
              write(3,570) numat(it), it, posn(:,igr)
            else
              write(3,580) numat(it), it, posn(:,igr), norbv(igr),
     &                       (popats(igr,l,1:nspin), l = 1,nlat(it))
              if( norbv(igr) == 0 ) cycle
              write(3,600) (hybrid(io,:,igr), pop_nonsph(io,igr), io,
     &                      io = 1,norbv(igr))
            endif
            l = l_hubbard( numat(it) )
            if( Atom_occ_hubb .and. Hubb(it) )
     &          write(3,610) ( ( occ_hubb_e(m,m,isp,igr), m = -l,l),
     &                               isp = 1,nspin )  
          end do
          write(3,*)

! About potental
          if( flapw ) then
            if( hedin .or. perdew ) then
              write(3,630)
            else
              write(3,640)
            endif
          else
            if( hedin ) then
              write(3,650)
            elseif( perdew ) then
              write(3,655)
            else
              write(3,660) alfpot
            endif
          endif
          if( Full_potential ) write(3,670)
          if( neimagent == 1 ) then
            write(3,680) eimagent(1)
          elseif( neimagent > 1 ) then
            write(3,685)
            write(3,690) (eeient(ie), eimagent(ie), ie = 2,neimagent)
          endif
          if( multrmax /= 1 ) write(3,702) multrmax
          if( abs( rpotmax ) > eps4 ) write(3,703) rpotmax
          write(3,704) D_max_pot
        endif

        if( Hubbard ) then
          write(3,'(/A)') ' Hubbard calculation '
          write(3,'(/A)') ' Type  Z  Hubbard parameter (eV)'
          do it = 1,ntype
            if( Hubb(it) ) write(3,'(2i4,f12.3)') it, numat(it), 
     &                                                   V_hubbard(it)
          end do
        endif

        if( self_cons .or. fermi_auto ) then
          if( Fermi_auto .and. abs(p_self0) < eps10 .and. nself == 1 )
     &       then
            write(3,'(/A)')' One cycle for the Fermi level calculation'
          else
            write(3,'(/A)') ' Self consistent calculation'
            write(3,708) nself, p_self0, Delta_En_conv
          endif
          write(3,709) r_self
          if( scf_elecabs ) write(3,'(A,A)')
     &       '   Cuting energy criterium using the number of electron',
     &       ' in the absorbing atom'
          if( self_nonexc ) then 
            write(3,'(A)') '   Non excited absorbing atom in this part'
          else
            write(3,'(A)') '   Excited absorbing atom in this part'
          endif
        end if

        if( Current ) write(3,710) Radius_current

      endif

 1320 continue  ! Point d'arrivee mpirank /= 0

      if( mpinodes > 1 ) then
        call MPI_Bcast(adimp,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(alfpot,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(allsite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        if( npldafs > 0 ) call MPI_Bcast(angpoldafs,3*npldafs, 
     &                              MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(angxyz,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Ang_base_loc,3*(ntype+1),MPI_REAL8,0,
     &                                          MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Ang_base_loc_gr,3*ngroup,MPI_REAL8,0,
     &                                          MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(ang_rotsup,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Ang_spin,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Axe_spin,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(axyz,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(basereel,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(base_spin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(BSE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(cartesian_tensor,1,MPI_LOGICAL,0,
     &                                           MPI_COMM_WORLD,mpierr)
        if( norbdil > 0 ) call MPI_Bcast(cdil,norbdil,MPI_REAL8,0,
     &                                           MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Centre,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Centre_auto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(Centre_auto_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(clementi,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Current,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(D_max_pot,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(dyn_eg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(dyn_g,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Bormann,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Charge_free,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(Dafs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Dafs_bio,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Delta_Epsii,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Delta_En_conv,1,MPI_REAL8,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(Density,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Dipmag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(dpos,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E1E1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E1E2e,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E1E3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E1M1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E1M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E2E2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(E3E3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(M1M1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(M1M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(M2M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Quadmag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Eclie,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Ecrantage,nspin,MPI_REAL8,0, 
     &                 MPI_COMM_WORLD,mpierr)
        if( neimagent > 0 ) call MPI_Bcast(Eeient,neimagent,MPI_REAL8,0, 
     &                                     MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(egamme,ngamme,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        if( neimagent > 0 ) call MPI_Bcast(eimagent,neimagent,MPI_REAL8, 
     &                                     0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Eneg_i,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Eneg_n_i,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Energphot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(f_no_res_mag,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(f_no_res_mom,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(force_ecr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Full_atom_e,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(Full_potential,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(Gamma_tddft,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(Green_int,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Green_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Green_self,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Hedin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        if( npldafs > 0 ) call MPI_Bcast(hkl_dafs,3*npldafs, 
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if( ngroup_nonsph > 0 ) then
          allocate( Hybrid_r(nhybm,16,ngroup_nonsph) )
          allocate( Hybrid_i(nhybm,16,ngroup_nonsph) )
          if( mpirank == 0 ) then
            hybrid_r(:,:,:) = real( hybrid(:,:,:),db )
            hybrid_i(:,:,:) = aimag( hybrid(:,:,:) )
          endif
          call MPI_Bcast(hybrid_r, 
     &         nhybm*16*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(hybrid_i, 
     &         nhybm*16*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

          call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
          if( mpirank /= 0 ) Hybrid(:,:,:)
     &                = cmplx( Hybrid_r(:,:,:), Hybrid_i(:,:,:),db )
          deallocate( Hybrid_i, Hybrid_r )
        endif
        if( Atom_occ_hubb ) call MPI_Bcast(occ_hubb_e,nspin*ngroup_hubb 
     &            *(2*m_hubb_e+1)**2,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(iabsm,n_multi_run_e,MPI_INTEGER,0,
     &                                MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(iabsorig,n_multi_run_e,MPI_INTEGER,0,
     &                                MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(iord,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if( ngroup_pdb > 0 ) call MPI_Bcast(Kgroup,ngroup_pdb,
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if( npldafs > 0 ) call MPI_Bcast(isigpi,npldafs*2, 
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if( norbdil > 0 ) call MPI_Bcast(itdil,norbdil, 
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(itype,ngroup,MPI_INTEGER,0,
     &                                            MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(jseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Kern_fac,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(korigimp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lamstdens,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        if( norbdil > 0 ) call MPI_Bcast(ldil,norbdil,MPI_INTEGER,0,
     &                                           MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(ldipimp,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lecrantage,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lin_gam,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lmaxat0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(lmaxfree,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lmaxso0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(lmoins1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lplus1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lquaimp,9,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(lseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        if( nlatm > 0 ) call MPI_Bcast(lvval,(ntype+1)*nlatm, 
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(matper,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(mix_repr,2,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(muffintin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(multrmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(n_atom_proto,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(nbseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nchemin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(necrantage,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(No_solsing,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(noncentre,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nonexc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        if( ngroup_nonsph > 0 ) call MPI_Bcast(norbv,ngroup_nonsph+1,
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(normaltau,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(normrmt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nphim,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nself,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nsymextract,n_multi_run_e,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nposextract,n_multi_run_e,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,mpierr)
        if( nlatm > 0 ) call MPI_Bcast(nvval,(ntype+1)*nlatm, 
     &                              MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(numat_abs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(octupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(old_reference,1,MPI_LOGICAL,0, 
     &                                            MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(One_run,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Optic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(overad,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(overlap,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(p_self0,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        if( nple > 0 ) then
          call MPI_Bcast(pdpolar,nple*2,MPI_REAL8,0,     
     &                                            MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(polar,3*nple,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        endif
        call MPI_Bcast(polarise,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        if( npldafs > 0 ) then
          call MPI_Bcast(poldafse,3*npldafs,MPI_REAL8,0,
     &                                            MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(poldafsei,3*npldafs,MPI_REAL8,0,
     &                                            MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(poldafss,3*npldafs,MPI_REAL8,0,
     &                                            MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(poldafssi,3*npldafs,MPI_REAL8,0,
     &                                            MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(vecdafsem,3*npldafs,MPI_REAL8,0,
     &                                            MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(vecdafssm,3*npldafs,MPI_REAL8,0,
     &                                            MPI_COMM_WORLD,mpierr)
        endif
        if( ngroup_nonsph > 0 ) call MPI_Bcast(pop_nonsph, 
     &            nhybm*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        if( nlatm > 0 ) then 
          call MPI_Bcast(popats,ngroup*nlatm*nspin,MPI_REAL8, 
     &                                        0,MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(popval,(ntype+1)*nlatm*nspin,MPI_REAL8, 
     &                                        0,MPI_COMM_WORLD,mpierr)
        endif
        call MPI_Bcast(posn,3*ngroup,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Quadrupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(r_self,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Radius_current,1,MPI_REAL8,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(Recup_tddft_data,1,MPI_LOGICAL,0, 
     &                                            MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Relativiste,1,MPI_LOGICAL,0, 
     &                                            MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(roverad,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(rpotmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(rrydb,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(rsorte_s,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Rydberg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Save_tddft_data,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(scf_elecabs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(self_cons,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(self_nonexc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                                                   mpierr)
        call MPI_Bcast(single_prec,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        do i = 1,3
          if( mpirank == 0 ) j = iachar( seuil(i:i) )
          call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
          call MPI_BARRIER(MPI_COMM_WORLD,mpierr) 
          if( mpirank /= 0 ) seuil(i:i) = achar( j ) 
        end do  
        call MPI_Bcast(solsing_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(solsing_only,1, 
     &                              MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(spherical_tensor,1, 
     &                              MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(spherical_signal,1, 
     &                              MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(core_resolved,1,MPI_LOGICAL,0, 
     &                                            MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Spinorbite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(state_all,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(state_all_out,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(supermuf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
          do i = 1,8
            if( mpirank == 0 ) j = iachar( PointGroup(i:i) )
            call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
            call MPI_BARRIER(MPI_COMM_WORLD,mpierr) 
            if( mpirank /= 0 ) PointGroup(i:i) = achar( j ) 
          end do  
        call MPI_Bcast(PointGroup_Auto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(symauto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(symmol,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Tddft,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Tddft_mix,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Atomic_scr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Tddft_so,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Test_dist_min,1,MPI_REAL8,0,MPI_COMM_WORLD,
     &                 mpierr) 
        call MPI_Bcast(rpalf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(temp,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(Taux_oc,ngroup_taux,MPI_REAL8,0,MPI_COMM_WORLD,
     &                                                          mpierr) 
        call MPI_Bcast(Vec_orig,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(v_intmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(V_hubbard,ntype+1,MPI_REAL8,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(V0bdcFimp,nspin,MPI_REAL8,0, 
     &                                            MPI_COMM_WORLD,mpierr)
        if( nple > 0 ) call MPI_Bcast(veconde,3*nple,MPI_REAL8,0,     
     &                                            MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Ylm_comp_inp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &                                            mpierr)
        call MPI_Bcast(Z_nospinorbite,1,MPI_INTEGER,0, 
     &                                            MPI_COMM_WORLD,mpierr)

        n = ntype + 1
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(icom,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nlat,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(Hubb,n,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nrato,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(nrato_dirac,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &                 mpierr)
        call MPI_Bcast(nrm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(numat,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_Bcast(popatc,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(rchimp,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr) 
        call MPI_Bcast(rmt,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)   
        call MPI_Bcast(rmtimp,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr) 
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr) 
        if( nlatm > 0 ) call MPI_Bcast(popatv,n*nlatm,MPI_REAL8,0,
     &                                         MPI_COMM_WORLD,mpierr) 

        if( flapw ) then
          call MPI_Bcast(its_lapw,ngroup_lapw,MPI_INTEGER,0,
     &                                          MPI_COMM_WORLD,mpierr)   
          call MPI_Bcast(nrato_lapw,n,MPI_INTEGER,0,
     &                                          MPI_COMM_WORLD,mpierr)   
          call MPI_Bcast(r0_lapw,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)   
          call MPI_Bcast(rlapw,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)   
          call MPI_Bcast(rotloc_lapw,9*ngroup_lapw,MPI_REAL8,0,
     &                                          MPI_COMM_WORLD,mpierr)   
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr) 

        if( mpirank /= 0 ) then
          if( Extract ) return
          icheck(:) = 0
        endif

      endif

      if( mpirank == 0 ) then
        istop = 0
        do k = 1,3
          if( abs(axyz(k)) < epspos ) then
            if( istop == 0 ) call write_error
            do ipr = ipr0,9,3
              write(ipr,100)
              write(ipr,740) k
            end do
            istop = 1
          endif
        end do
        if( istop == 1 ) stop
      endif
      
! Conversion en unites atomiques.
      rad = pi / 180

      adimp = adimp / bohr

      where( abs(angpoldafs) < 9999._db ) angpoldafs = angpoldafs * rad
      axyz(1:3) = axyz(1:3) / bohr

      D_max_pot = D_max_pot / bohr
      Delta_En_conv = Delta_En_conv / Rydb
      Delta_Epsii = Delta_Epsii / Rydb
      dpos(:) = dpos(:) / bohr
      eclie = eclie / Rydb
      eeient(:) = eeient(:) / Rydb
      egamme(:) = egamme(:) / Rydb
      eimagent(:) = eimagent(:) / Rydb
      if( Hubbard ) V_hubbard(:) = V_hubbard(:) / Rydb
      r_self = r_self / bohr
      rchimp(:) = rchimp(:) / bohr
      rmtimp(:) = rmtimp(:) / bohr
      roverad = roverad / bohr
      rpotmax = rpotmax / bohr
      rrydb = rrydb / bohr
      rsorte_s = rsorte_s / bohr
      Test_dist_min = Test_dist_min / bohr
      v_intmax = v_intmax / Rydb
      v0bdcFimp(:) = v0bdcFimp(:) / Rydb

! Nombre d'orbitales atomiques maximum
      nnlm = 0
      do it = 1,ntype
        n = n_orb_rel( numat(it) )
        if( numat(it) == numat_abs )then
          n = n + 2 * max( 1, nlat(it) )
        elseif( nlat(it) > 0 ) then
          n = n + 2 * nlat(it)
        endif 
        nnlm = max( nnlm, n )
      end do

      call cal_cubmat(angxyz,cubmat,struct)

      if( magnetic .or. Atom_nonsph ) then
        call invermat(cubmat,cubmati)

        if( abs( Axe_spin(1) ) < eps10 .and. abs( Axe_spin(2) ) < eps10
     &       .and. abs( Axe_spin(3) - 1 ) < eps10 ) then
          Axe_spin(1:2) = 0._db; Axe_spin(3) = 1._db
          Axe_spin = matmul(cubmat, Axe_spin)
          call mat_euler( Ang_spin, Rot_gen )
          Axe_spin = matmul( Rot_gen, Axe_spin )
        else
          Axe_spin(:) = Axe_spin(:) * axyz(:) * bohr
          Axe_spin = matmul( cubmat, Axe_spin )
          vv = sqrt( sum( Axe_spin(:)**2 ) )
          Axe_spin = Axe_spin / vv
          if( abs( Axe_spin(3) ) < 1 - eps10 ) then
            Ang_spin(1) = datan2( Axe_spin(2), Axe_spin(1) )
            Ang_spin(2) = acos( Axe_spin(3) )
            Ang_spin(3) = 0._db
          else
            Ang_spin(:) = 0._db
          endif
          call mat_euler( Ang_spin, Rot_gen )
        endif
        Axe_spin = matmul( cubmati, Axe_spin )
        Axe_spin(:) = Axe_spin(:) / axyz(:)

        do igr = 1,ngroup
          if( Ang_base_loc_gr(1,igr) > - 1000._db ) cycle
          if( Ang_base_loc(1, abs(itype(igr)) ) < -1000._db ) cycle
          Ang_base_loc_gr(:,igr) = Ang_base_loc(:, abs(itype(igr)) )
        end do
        do igr = 1,ngroup
          if( Ang_base_loc_gr(1,igr) < - 1000._db )
     &             Ang_base_loc_gr(:,igr) = Ang_spin(:)
          Ang(:) = Ang_base_loc_gr(:,igr)
          call mat_euler( Ang, Rot )
          Rot_Atom_gr(:,:,igr) = Rot(:,:)
          if( struct == 'trigo' ) then
            Axe(1:3) = 1._db / sqrt(3._db)
          else
            Axe(1:2) = 0._db; Axe(3) = 1._db
          endif
          Axe = matmul(cubmat, Axe)
          Axe = matmul(Rot, Axe)
          Axe = matmul( cubmati, Axe )
          Axe_atom_gr(:,igr) = Axe(:) / axyz(:)
          if( itype(igr) < 0 ) Axe_atom_gr(:,igr) = - Axe_atom_gr(:,igr)

          p(:) = Axe_atom_gr(:,igr) * axyz(:)
          p = matmul( cubmat, p )
          vv = sum( p(:)**2 )
          if( vv > eps10 ) p(:) = p(:) / sqrt( vv ) 
          p = matmul( cubmati, p )
! Axe_atom_gr est maintenant normalise a 1.
          Axe_atom_gr(:,igr) = p(:) / axyz(:)
        end do
        Radius_current = Radius_current / bohr

        if( icheck(1) > 0 ) then
          write(3,750) Axe_spin(1:3) * axyz(1:3)
          write(3,760) Ang_spin(1:3) * 180 / pi
          do igr = 1,ngroup
            write(3,770) igr 
            do i = 1,3
              write(3,780) Rot_atom_gr(i,:,igr), 
     &                                    Axe_atom_gr(i,igr) * axyz(i)
            end do
          end do
        endif
      endif

      Atom_nsph_e(:) = .false.
      if( Atom_nonsph ) then
        do igr = 1,ngroup
          if( norbv(igr) == 0 ) cycle
          pop_nsph = sum( pop_nonsph(1:norbv(igr),igr) )
          if( pop_nsph > eps10 ) Atom_nsph_e(igr) = .true.
        end do
      endif

      if( octupole) then
        l_selec_max = lseuil + 3
      elseif( Quadrupole) then
        l_selec_max = lseuil + 2
      else
        l_selec_max = lseuil + 1
      endif

      if( Dafs_bio ) then
        nphim = 4
        nphi_dafs(:) = nphim
      else
        nphimt = 1
        nphi_dafs(:) = 1
        do ipl = 1,npldafs
          if( angpoldafs(1,ipl) > -9999._db .and. 
     &        angpoldafs(2,ipl) > -9999._db .and.
     &        angpoldafs(3,ipl) > -9999._db ) cycle
          nphimt = nphim
          nphi_dafs(ipl) = nphim
        end do
        nphim = nphimt
      endif

      if( Dafs ) then
        poldafsem(:,:) = cmplx( poldafse(:,:), poldafsei(:,:), db ) 
        poldafssm(:,:) = cmplx( poldafss(:,:), poldafssi(:,:), db ) 
      endif

      if( Centre_auto ) call Auto_center(axyz,Centre,
     &                   Centre_auto_abs,cubmat,icheck(1),itype,
     &                   ngroup,ntype,numat,numat_abs,posn,struct)

      do igr = 1,ngroup
        posn(:,igr) = posn(:,igr) - Centre(:)
      end do

      D_max_pot = D_max_pot / adimp

      return

 9999 call write_err_form(itape4,grdat)

      return
  100 format(//'  Error in the indata file :')
  110 format(//' Wrong number of elements under the card Polarized !'/)
  120 format(//' After the keyword Dafs (or RXS),'/,
     &         ' for the reflexion number',i3,/
     &         ', the polarization is not defined !'//)
  122 format(//' After the keyword Dafs (or RXS),'/,
     &         ' for the reflexion number',i3,/
     &         ' there are at least 2 directions for an angular scan.',/
     &         ' Just one is allowed !'//)
  125 format(/' After the keyword Atom_conf, for the atom number',i2,
     &        ','/,' check how is written the corresponding electronic',
     &        ' configuration !'//)
  130 format(/' After the keyword Atom, for the atom number',i2,','/
     &        ' check how is written the corresponding electronic',
     &        ' configuration !'//)
  140 format(/' Error in the Pdb file :',//
     &        ' The chemical symbol of atom number',i6,' is: ',
     &        a2,', what is not known!'//)
  150 format(//' The following line is not understood :',/A,//
     &        ' If it is a keyword, check the spelling.'/,
     &        ' If the line is not supposed to be a keyword but',
     &        ' contains numbers, check:'/
     &        5x,' - How many numbers must be in the line ?'/,
     &        5x,' - Are there spaces between the numbers ?'/,
     &        5x,' - Tabulations are forbidden !'//)
  160 format(///' Error under the keyword Par_',a6,/
     &         ' The wanted atom is the number',i3,' !',/
     &         ' There are only',i3,' atoms in the job !'//)
  170 format(//'  A parameter index for the fit is not possible !'/,
     &         '  Check your indata file under the keyword ',a9,/
     &         '  The index is',i4//)
  172 format(/' Atom number',i4,' is indexed with',i3,' what is more',
     & ' than',/' the number of atom type =',i2,' declared under',
     & ' keyword Atom !')
  174 format(//' Under keyword Crystal or Molecule, when the index is',/
     &         ' supposed to be the atomic number, the simultaneous',/
     &         ' use of the keyword Atom is forbidden !'//)
  175 format(//' The atomic number given under keyword Z_absorber =',
     &           i4,' is not in the list of atoms !'//)
  180 format(/' Edge = ',a3,' not programmed !'//)
  190 format(//'  No transitions are allowed in your calculation !'/,
     &         ' Check the keywords Quadrupole, Octupole, No_dipole,',
     &         ' No_dipquad, No_E2E2, Dipmag, Quadmag, ...')
  200 format(//'  There is no atom in your calculation !'/,
     &         ' Some necessary keywords as molecule or crystal could',
     &         ' be missing'//)
  210 format(//'  There is no chemical species specified in your',
     &         ' calculation !'//)
  220 format(//' iord =',i2, ' must be equal to 2 or 4 !'//)
  230 format(/' itype(igr=',i2,') =',i2,' > ntype =',i2,', forbidden !')
  240 format(//' Bad value in the indata file !'/,
     &     ' Under the keyword polarise, polarisation and wave vector',
     &     ' (if specified) cannot be both zero !')
  250 format(//' Bad number of values in the indata file !'/,
     &     ' Under the keyword Rmtg, the number of radius values must',
     &     ' be equal to the number of atom type !',/
     &     ' In the indata file number of radius =',i3,/
     &     '                 number of atom type =',i3,/)
  260 format(//' Bad polarization definition for DAFS when using',
     &         ' Full_self_abs option !'//,
     &     ' For reflection number',i3,', (h,k,l) = (',3i3,')',/
     &     ' polarization indexes are',i2,' and',i2,' !'//,
     &     ' For each reflection (h,k,l) one must have in order',/,
     &     2x,' sigma-sigma, sigma-pi, pi-sigma and pi-pi,',/,
     &     2x,' that is respectively indexes 1 1, 1 2, 2 1 and 2 2',/)
  300 format(/' Tensors Extracted from the file :'/,A,/) 
  320 format(' Radius =',f6.2)
  325 format('    Roverad =',f6.2)
  330 format(' icheck =',30i2)
  340 format(' Linear range :')
  350 format(' Range =',9f8.3,5(/9x,9f8.3))
  360 format(' Dipole component')
  365 format(' Quadrupole component')
  370 format(' Octupole component')
  371 format(' Magnetic dipole component')
  372 format(' ... but E1E2 component neglected in the output')
  373 format(' ... but E2E2 component neglected in the output')
  374 format(' ... but E1E3 component neglected in the output')
  375 format(' ... but E3E3 component neglected in the output')
  376 format(' ... but E1E2 component neglected in the output')
  377 format(' TDDFT calculation using the RPA LF approximation')
  378 format(' TDDFT calculation')
  379 format('    Broadening in the Chi_0 calculation')
  381 format('    No broadening in the Chi_0 calculation')
  382 format(' Core resolved in outputs')
  384 format(' Relativistic calculation')
  385 format(' Non-relativistic calculation')
  390 format(' Magnetic calculation with spin-orbit interaction')
  400 format(' Magnetic calculation without spin-orbit interaction')
  410 format(' Non-magnetic calculation')
  416 format(' Spin-orbit not taken into account for the atomic number',
     &         i3)
  419 format(' Approximation l-1')
  420 format(' Approximation l+1')
  421 format(' Real bases')
  422 format(' Complex bases')
  425 format(' Multiple scattering calculation (Green)')
  426 format('    Full Green mode')
  427 format('    Continuous potential (Supermuf)')
  430 format('    Path expansion, n =',i3)
  435 format('    Tau normalization')
  436 format('    Optimized type muffin-tin radius')
  437 format('    Norman type muffin-tin radius')
  438 format('    Half interatomic distance type muffin-tin radius')
  439 format('    Imposed type muffin-tin radius,',
     &       ' Rmtimp =',10f6.3,/9x,10f6.3)
  440 format('    Potential imposed type muffin-tin radius')
  445 format('    Overlap of the muffin-tin radius  =',f6.2)
  452 format('    No limitation on the maximum value of l')
  454 format('    Limitation on the maximum value of l')
  460 format(' Finite difference method calculation',/
     &       '   iord =',i2,', adimp =',f6.2)
  470 format('   lmaxso0 =',i3)
  492 format('   Muffin-tin potential')
  494 format('   Rrydb =',f7.3,' A')
  496 format('   Non centered absorbing atom')
  497 format('     Center =',3f7.3)
  498 format('   Eclie =',f7.3,' eV')
  499 format('   V_intmax =',f7.3,' eV')
  500 format(/' Calculation with non spherical orbitals')
  501 format(/' Calculation in single precision')
  502 format(/' Calculation in double precision')
  503 format(/' Temperature =',f6.1,' K')
  504 format(/' Temperature coefficients taken into account for',
     &        ' diffraction')
  505 format(/' XANES :    Polarization             Wave vector ',
     &       '    Weight_dip Weight_quad')
  509 format(6x,3f7.3,3x,3f7.3,3x,f7.3,4x,f7.3)
  510 format(/' XANES :    Polarization             Wave vector')
  511 format(/' Bormann',/
     &  '  (h, k, l) = (',i3,',',i3,',',i3,')   Azimuth =',f7.2)
  512 format(/' DAFS : (h, k, l)  Polarization   Angle_i   Angle_o',
     &'  Azimuth')
  514 format(7x,3i3,3x,a11,2f10.3,'      scan')
  517 format(7x,3i3,3x,a11,'      scan',2f10.3)
  518 format(7x,3i3,3x,a11,f10.3,'          scan',f10.3)
  519 format(7x,3i3,3x,a11,3f10.3)
  520 format(/' DAFS :     Polarization             Wave vector ' )
  522 format(/' XANES for the RXS polarizations')
  523 format(/' Self absorption taken into account with',
     &        ' birefringence effect')
  530 format(6x,3f7.3,3x,3f7.3,3x,' incoming')
  540 format(6x,3f7.3,3x,3f7.3,3x,' outcoming')
  542 format('   ngroup =',i5,', ntype =',i2)
  543 format(1x,2i4,10(i4,i3,f7.3))
  544 format(1x,2i4,10(i4,i3,2f7.3))
  545 format('    dpos =',3f7.3)
  548 format(' Orbital dilatation :',/'   it    l   cdil')
  549 format(2i5,f7.3)
  550 format('   Non excited absorbing atom')
  552 format('   Point Group = ',a8)
  553 format('   Absorbing atom taken as not excited for the symmetry',
     &       ' calculation')
  554 format('   Space Group = ',a10)
  555 format('   a, b, c =',3f10.5)
  556 format('   alfa, beta, gamma =',3f9.3)
  557 format('    Z  Typ    posx      posy      posz   Kgroup   Occup.',
     &       ' Temp_cf')
  558 format('    Z  Typ    posx      posy      posz   Kgroup   Occup.')
  559 format('    Z  Typ    posx      posy      posz   Kgroup',
     &       '   Temp_cf')
  560 format('    Z  Typ    posx      posy      posz   Kgroup')
  561 format('    Z  Typ    posx      posy      posz     Occup.',
     &       ' Temp_cf')
  562 format('    Z  Typ    posx      posy      posz     Occup.')
  563 format('    Z  Typ    posx      posy      posz    Temp_cf')
  564 format('    Z  Typ    posx      posy      posz   norbv   popats')
  565 format('    Z  Typ    posx      posy      posz      popats')
  570 format(i5,i4,3f10.5,i6,2x,2f8.2)
  572 format(i5,i4,3f10.5,2x,2f8.2)
  580 format(i5,i4,3f10.5,i5,2x,12f8.4)
  600 format(4x,17(1x,2f7.3),1x,'= hybrid, pop_nonsph(',i1,')')
  610 format('    Occ. matrix :', 14f5.2)
  630 format(' FLAPW potential energy dependant')
  640 format(' FLAPW potential not energy dependant')
  650 format(' Hedin and Lundqvist exchange-correlation potential')
  655 format(' Perdew and Zunger exchange-correlation potential')
  660 format(' Xalfa potential , Xalfa =',f8.5)
  670 format(' Full potential inside the atomic spheres')
  680 format(' E_imag =',f9.3,' eV')
  685 format('   Energy   E_imag    (eV)')
  690 format(2f9.3)
  702 format('   multrmax =',i2)
  703 format('   rpotmax =',f7.3)
  704 format('   D_max_pot =',f7.3)
  708 format('   Maximum number of iteration = ',i3,/
     &       '   Weight =',f6.3,/
     &       '   Delta energy for convergence =',f7.3,' eV / atom')
  709 format('   Cluster radius used for this part = ',f7.3,' A')
  710 format(/' Calculation taking into account the current.'/,
     &        ' inside the sphere of radius R =',f6.2,
     &        ' A, around the absorbing atoms.')
  740 format(//' The mesh parameter',i2,' is zero !'//)
  750 format(/' General Z axis =',3f9.5)
  760 format(' Euler angles   =',3f9.3)
  770 format(/6x,' local matrix rotation',7x,'local Z axis   Atom =',i3)
  780 format(3x,3f9.5,5x,f9.5)

      end
     
!***********************************************************************

! Routine de lecture de la structure venant de WIEN

      subroutine lect_struct_lapw(angxyz,axyz,icheck,its_lapw,itype,
     &      ngroup,ngroup_lapw,nomstruct,nrato_lapw,nrm,ntype,
     &      numat,posn,r0_lapw,rlapw,rotloc_lapw)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      character(len=1) Trans
      character(len=2) Plan
      character(len=132) nomstruct

      integer, dimension(ngroup) :: numprot
      integer, dimension(0:ntype) :: nrato_lapw, numat
      integer, dimension(ngroup_lapw) :: its_lapw
      integer, dimension(ngroup) :: itype

      real(kind=db), dimension(3):: angxyz, axyz, dp, p
      real(kind=db), dimension(3,3):: Mas, Mat, ptrans, Rot
      real(kind=db), dimension(3,ngroup):: posn
      real(kind=db), dimension(3,3,ntype):: rotloc
      real(kind=db), dimension(3,3,ngroup_lapw):: rotloc_lapw
      real(kind=db), dimension(0:ntype):: r0_lapw, rlapw

      common/lapwksym/ matsym(3,3,nslapwm)
      common/lapwtau/ taulap(3,nslapwm)
! ll             = nombre de (l,m) par atome
! ntype          = nombre d'atomes inequivalents
! nmatsym        = nombre d'op. de symetrie
!
! Unite 8 : 'case.struct'  (donnees structurales)

      open(8, file = nomstruct, status='old', iostat=istat) 
      if( istat /= 0 ) call write_open_error(nomstruct,istat,1)

! Lecture de xxxx.struct
! Read et formats pris en partie dans la routine main1.f de
! Wien97/SRC_lapw5

      read(8,*)
      read(8,'(a1,a2)') Trans, Plan
      demi = 0.5_db
      ptrans = 0._db
      select case(Trans)
        case('F')
          ntrans = 3
          ptrans(1,1) = demi; ptrans(2,1) = demi; ptrans(3,1) = 0.0
          ptrans(1,2) = demi; ptrans(2,2) = 0.0;  ptrans(3,2) = demi
          ptrans(1,3) = 0.0;  ptrans(2,3) = demi; ptrans(3,3) = demi
        case('B')
          ntrans = 1
          ptrans(1:3,1) = demi
        case('C')
          ntrans = 1
          if( Plan == 'XY' ) then
            ptrans(1:2,1) = demi
          elseif( Plan == 'YZ' ) then
            ptrans(1:3:2,1) = demi
          else
            ptrans(2:3,1) = demi
          endif  
      case default
          ntrans = 0
          ptrans(1:3,1) = 0._db
      end select

      if( Trans == 'H' ) then
        Mat(1,1) = 1._db; Mat(1,2) = -1._db / 2;      Mat(1,3) = 0._db
        Mat(2,1) = 0._db; Mat(2,2) = sqrt(3._db) / 2; Mat(2,3) = 0._db
        Mat(3,1) = 0._db; Mat(3,2) = 0._db;           Mat(3,3) = 1._db
        Mas(1,1) = 1._db; Mas(1,2) = 1 / sqrt(3._db); Mas(1,3) = 0._db
        Mas(2,1) = 0._db; Mas(2,2) = 2 / sqrt(3._db); Mas(2,3) = 0._db
        Mas(3,1) = 0._db; Mas(3,2) = 0._db;           Mas(3,3) = 1._db
      endif

      read(8,*)
      read(8,'(6f10.7)') axyz(1:3), angxyz(1:3)
      axyz(1:3) = axyz(1:3) * bohr

      index = 0
      do jatom = 1,ntype

        index = index + 1
        read(8,'(5x,i3,1x,3(3x,f10.7))') its, posn(:,index)
        it = abs( its )
        itype(index) = it
        its_lapw(index) = its
        numprot(index) = index
        iprot = index
        do itr = 1,ntrans
          index = index + 1
          posn(1:3,index) = posn(1:3,iprot) + ptrans(1:3,itr)
          itype(index) = it
          its_lapw(index) = its
          numprot(index) = iprot
        end do

        read(8,'(15x,i2)') mult
        do mu = 1,mult-1
          index = index + 1
          read(8,'(5x,i3,1x,3(3x,f10.7))') ittt, posn(1:3,index)
          itype(index) = it
          its_lapw(index) = its
          numprot(index) = iprot
          do itr = 1,ntrans
            index = index + 1
            posn(1:3,index) = posn(1:3,index-itr) + ptrans(1:3,itr)
            itype(index) = it
            its_lapw(index) = its
            numprot(index) = iprot
          end do
        end do

! Maillage radial dans les spheres atomiques
        read(8,'(15x,i5,2(5x,f10.5),5x,f5.2)')
     &                    nrato_lapw(it), r0_lapw(it), rlapw(it), rZ
        nrm = max( nrm, nrato_lapw(it) + 200 )

        numat(it) = nint( rZ )

! Matrice de rotation
        read(8,'(20x,3f10.7)') (rotloc(i,1:3,it), i = 1,3)

      end do

      where( posn > 1._db - eps10 ) posn = posn - 1._db

! Operations de symetrie
      read(8,'(i4)') nmatsym
      do is = 1,nmatsym
        read(8,'(3(3i2,f10.5,/))') (matsym(i,1:3,is), taulap(i,is),
     &                              i = 1,3)
      end do
      close(8)

! Recherche de l'operation de symetrie qui renvoit a l'atome
! prototypique
      nt = ntrans + 1

      boucle_exter: do igr = 1,ngroup
        it = itype(igr)
        iprot = numprot(igr)

        if( igr == iprot ) then
          do i = 1,3
            rotloc_lapw(i,1:3,igr) = rotloc(1:3,i,it)
          end do
          cycle
        elseif( ntrans > 0 .and. mod(igr,nt) /= 1 ) then
          do i = 1,3
            rotloc_lapw(i,1:3,igr) = rotloc_lapw(i,1:3,igr-1)
          end do
          cycle
        endif

        do is = 1,nmatsym
          do j = 1,3
            p(j) = sum( matsym(j,:,is) * (posn(:,igr) - taulap(:,is)) )
            if( p(j) < -epspos ) then
              p(j) = p(j) + 1._db
            elseif( p(j) >= 1._db-epspos ) then
              p(j) = p(j) - 1._db
            endif
          end do

          dp(1:3) = abs( p(1:3) - posn(1:3,iprot) )
          if( dp(1) < epspos .and. dp(2) < epspos .and. dp(3) < epspos )
     &        then
            Rot(:,:) = Matsym(:,:,is)
            if( Trans == 'H' ) Rot = Matmul( Mat, Matmul( Rot, Mas ) ) 
            do i = 1,3
              do j = 1,3
                Rotloc_lapw(j,i,igr) =
     &             sum( Rot(i,1:3) * Rotloc(1:3,j,it) )
              end do
            end do
            cycle boucle_exter
          endif
        end do

        call write_error
        do ipr = 3,9,3
          if( icheck == 0 .and. ipr == 3 ) cycle
          write(ipr,110)
        end do
        stop

      end do boucle_exter

      if( icheck > 1 ) then
        do ia = 1,ngroup
          write(3,130) ia
          write(3,'(3f7.3)') (rotloc_lapw(i,1:3,ia), i = 1,3)
        end do
      endif

      return
  110 format(/' Atome symetrique non trouve !')
  130 format(/' rotloc(ia=',i2,')')
      end

!***********************************************************************

! Calcul de la matrice de changement de repere maille - orthogonale

      subroutine cal_cubmat(angxyz,cubmat,struct)

      use declarations
      implicit none
      include 'mpif.h'

      character(len=5) struct

      logical ang(3), ange(3)

      real(kind=db):: a, alfa, b, beta, cosa, cosb, cosg, gamma, rad,
     &                sina, sinb
      real(kind=db), dimension(3):: angxyz(3)
      real(kind=db), dimension(3,3):: cubmat

! Matrice de changement de repere cristallo, cubique
      ang(:) = abs( angxyz(:) - 90._db ) < eps4
      ange(1) = abs( angxyz(2) - angxyz(3) ) < eps4
      ange(2) = abs( angxyz(3) - angxyz(1) ) < eps4
      ange(3) = abs( angxyz(1) - angxyz(2) ) < eps4
      if( ange(1) .and. ange(2) .and. ange(3) ) then
        if( ang(1) ) then
          struct = 'cubic'
        else
          struct = 'trigo'
        endif
      elseif( ( abs( angxyz(3) - 120. ) < eps4 ) .and. ang(1)
     &        .and. ang(2) ) then
        struct = 'hexag'
      else
        struct = 'autre'
      endif

      if( struct /= 'cubic' ) then

        rad = pi / 180
        if( struct == 'trigo' ) then
          alfa = angxyz(1) * rad
          cosa = sqrt( ( 1._db + 2 * cos( alfa ) ) / 3. )
          sina = sqrt( 1._db - cosa**2 )
          cubmat(1,1) = sina;  cubmat(1,2:3) = -0.5_db * sina
          cubmat(2,1) = 0._db;  cubmat(2,2) = sqrt(3._db) * sina / 2
          cubmat(2,3) = - cubmat(2,2)
          cubmat(3,1:3) = cosa
        else
          alfa = angxyz(1) * rad
          beta = angxyz(2) * rad
          gamma = angxyz(3) * rad
          sina = sin( alfa )
          cosa = cos( alfa )
          sinb = sin( beta )
          cosb = cos( beta )
          cosg = cos( gamma )
          a = ( cosg - cosa*cosb ) / sinb
          b = sqrt( sina**2 - a**2 )
          cubmat(1,1) = sinb;  cubmat(1,2) = a;     cubmat(1,3) = 0._db
          cubmat(2,1) = 0._db; cubmat(2,2) = b;     cubmat(2,3) = 0._db
          cubmat(3,1) = cosb;  cubmat(3,2) = cosa;  cubmat(3,3) = 1._db
        endif

      else

         cubmat(:,:) = 0._db
         cubmat(1,1) = 1._db; cubmat(2,2) = 1._db; cubmat(3,3) = 1._db

      endif

      return
      end

!***********************************************************************

! Calcule la matrice de rotation en fonction des angles d'Euler.

      subroutine mat_euler(Ang,Rot)

      use declarations
      implicit real(kind=db) (a-h,o-z)

      real(kind=db), dimension(3):: Ang
      real(kind=db), dimension(3,3):: mat, Rot

      do l = 1,3
  
        Angr = Ang(4-l) 

        cs = cos( Angr )
        ss = sin( Angr )
        i = mod(l,3) + 1
        j = mod(l+1,3) + 1
        k = mod(l-1,3) + 1 
        mat(i,i) = cs;    mat(i,j) = -ss;   mat(i,k) = 0._db
        mat(j,i) = ss;    mat(j,j) = cs;    mat(j,k) = 0._db
        mat(k,i) = 0._db;  mat(k,j) = 0._db;  mat(k,k) = 1._db
        if( l == 1 ) then
          rot = mat
        else
          rot = matmul( mat, rot )
        endif
      end do

      return
      end

!***********************************************************************

! Calcul du centre de l'agregat.

      subroutine Auto_center(axyz,Centre,Centre_auto_abs,
     &                       cubmat,icheck,itype,ngroup,ntype,
     &                       numat,numat_abs,posn,struct)

      use declarations
      implicit none

      integer:: ia1, ia2, ia3, ia4, icheck, igr, jgr, ngr, ngroup,
     &          ntype, numat_abs, Z
      integer, dimension(0:ntype):: numat
      integer, dimension(ngroup):: itype

      character(len=2):: Chemical_Symbol
      character(len=5) struct

      logical:: Centre_auto_abs

      real(kind=db):: dist, dist_max, Radius
      real(kind=db), dimension(3):: axyz, b, Centre, p, q, v
      real(kind=db), dimension(3,3):: Cubmat, Mat, Mati
      real(kind=db), dimension(3,ngroup):: posn, pos

      dist_max = 0._db

      jgr = 0
      do igr = 1,ngroup
        Z = numat( abs( itype(igr) ) )
        if( Centre_auto_abs .and. Z /= numat_abs ) cycle
        jgr = jgr + 1
        if( struct /= 'cubic' ) then
          p(:) = posn(:,igr)
          p = matmul( Cubmat, p )
          pos(:,jgr) = p(:) * axyz(:) 
        else
          pos(:,jgr) = posn(:,igr) * axyz(:) 
        endif
      end do

      ngr = jgr

      do igr = 1,ngr
        Z = numat( abs( itype(igr) ) )
        if( igr == 1 ) Centre(:) = pos(:,igr)
        do jgr = igr+1,ngr
          p(:) = pos(:,igr) - pos(:,jgr)
          dist = sqrt( sum( ( p(:) )**2 ) )
          if( dist < dist_max ) cycle
          dist_max = dist
          Centre(:) = 0.5_db * ( pos(:,igr) + pos(:,jgr) ) 
          ia1 = igr
          ia2 = jgr
        end do
      end do

      Radius = dist_max / 2

      dist_max = Radius
      ia3 = 0
      do igr = 1,ngr
        p(:) = pos(:,igr) - Centre(:)
        dist = sqrt( sum( ( p(:) )**2 ) )
        if( dist < dist_max + eps10 ) cycle
        dist_max = dist
        ia3 = igr
      end do

      if( ia3 /= 0 ) then
! Recherche du cercle circonscrit

! Plan hauteur 1
        v(:) = pos(:,ia2) - pos(:,ia1)
        Mat(1,:) = v(:)
        b(1) = sum( v(:)*Centre(:) )

! Plan hauteur 2
        v(:) = pos(:,ia3) - pos(:,ia1)
        Mat(2,:) = v(:)
        b(2) = 0.5_db * sum( v(:) * ( pos(:,ia3) + pos(:,ia1) ) )

! Plan du triangle
        p(:) = pos(:,ia2) - pos(:,ia1)
        q(:) = pos(:,ia3) - pos(:,ia1)
        call prodvec(v,p,q)
        Mat(3,:) = v(:)
        b(3) = sum( v(:)*Centre(:) )

        call invermat(Mat,Mati)

        Centre = Matmul( Mati, b )

        v(:) = pos(:,ia3) - Centre(:)       
        Radius = sqrt( sum( v(:)**2 ) )

        dist_max = Radius
        ia4 = 0
        do igr = 1,ngr
          p(:) = pos(:,igr) - Centre(:)
          dist = sqrt( sum( ( p(:) )**2 ) )
          if( dist < dist_max + eps10 ) cycle
          dist_max = dist
          ia4 = igr
        end do

        if( ia4 /= 0 ) then
! Recherche de la sphere circonscrite a 4 points

! Plan hauteur 3
          v(:) = pos(:,ia4) - pos(:,ia1)
          Mat(3,:) = v(:)
          b(3) = 0.5_db * sum( v(:) * ( pos(:,ia4) + pos(:,ia1) ) )

          call invermat(Mat,Mati)

          Centre = Matmul( Mati, b )

          v(:) = pos(:,ia4) - Centre(:)       
          Radius = sqrt( sum( v(:)**2 ) )

        endif

      endif

      call invermat(Cubmat,Mati)

      Centre = Matmul( Mati, Centre )
      Centre(:) = Centre(:) / axyz(:)

      if( icheck > 0 ) then
        write(3,110) Centre(:)
        if( Centre_auto_abs ) then
          write(3,120) Chemical_Symbol(numat_abs), Radius * bohr
        else
          write(3,130) Radius * bohr
        endif
      endif

      return
  110 format(/' Center set at  =',3f9.5,'  (cell unit)')
  120 format(' Farest ',a2,' atom at =',f9.5,' A')
  130 format(' Farest atom at =',f7.3,' A')
      end



