!***********************************************************************
! Read_DVR_Module
!**begin prologue     Read_DVR_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Read in FEDVRata
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Read_DVR_Module
!***********************************************************************
!***********************************************************************
                           MODULE Read_DVR_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Legendre_Data
                    IMPLICIT NONE
!
!***********************************************************************
                           Contains
!***********************************************************************
!***********************************************************************
!deck Set_General_Keywords
!***begin prologue     Set_General_Keywords
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input routines for dvr basis sets
!***description        user interface for dvr library
!                      Note that the various Xkey functions enable keyword
!                      input, with default values if the keyword does not appear 
!                      in the input line.  There are Xkey routine for real, integer,
!                      logical and character type variables.
!                      The general process is that the card variable is read and
!                      the string decoded by the Xkey routine which set the variables.
!
!                      The prnkey parameter is;
!                      'sector_points', 'sector_factors', 'sector_polynomials', 
!                      'sector_matrices' ,'global_points', 'global_polynomials', 
!                      'potential','global_matrices', 'hamiltonian', 'eigenvalues',
!                      'eigenvectors', 'all'
!***references
!***routines called    
!***end prologue       Set_General_Keywords
  SUBROUTINE Set_General_Keywords ( sets ) 
  IMPLICIT NONE
  REAL(idp)                   :: fpkey
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: i
  INTEGER                     :: intkey
  INTEGER                     :: sets
  CHARACTER(LEN=80)           :: chrkey
  CHARACTER(LEN=8)            :: itoc
  CHARACTER(LEN=3)            :: pointer
  pointer=itoc(sets)
  len=lenth(pointer)
  IF ( dollar('$general_keywords_'//pointer(1:len),card,cpass,inp) )THEN
       keyword = chrkey(card,'coordinate_system','cartesian',' ')
       spdim = intkey(card,'number_of_spatial_dimensions',1,' ')
       units=chrkey(card,'units','atomic_units',' ')
       mass=fpkey(card,'mass',1.d0,' ')
       typwt=chrkey(card,'weighting_function','one',' ')
       refwt=chrkey(card,'reference_weight','legendre',' ')
       bc=fpkey(card,'boundary_condition',0.d0,' ')
       ALLOCATE(reg_grid(1:spdim))
       DO i = 1, 12
          prn(i) = logkey(card,prnkey(i),.false.,' ')
       END DO
       DO i=1,spdim
          reg_grid(i)%label = chrkey(card,'coordinate_label_'//itoc(i),'x',' ')
       END DO
       form_nabla=logkey(card,'form_nabla',.true.,' ')
       form_hamiltonian=logkey(card,'form_hamiltonian',.true.,' ')
       IF (form_nabla == .true. ) THEN
           read_nabla=logkey(card,'read_nabla',.false.,' ')
           write_nabla=logkey(card,'write_nabla',.true.,' ')
           diagonalize_nabla=logkey(card,'diagonalize_nabla',.false.,' ')
       END IF
       IF (form_hamiltonian == .true. ) THEN
           diagonalize_hamiltonian=logkey(card,'diagonalize_hamiltonian',.false.,' ')
           write_hamiltonian=logkey(card,'write_hamiltonian',.true.,' ')
           read_hamiltonian=logkey(card,'read_hamiltonian',.false.,' ')
       END IF
       poisson=logkey(card,'solve_poisson_equation',.false.,' ')
       atomic=logkey(card,'atomic_calculation',.false.,' ')
       typke='dvr'
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
       proj=.false.
!
!      Z_a is located on the negative z axis, z_b on the positive.
!
       Z_a=fpkey(card,'charge_nucleus_a',1.d0,' ')
       Z_b=fpkey(card,'charge_nucleus_b',1.d0,' ')
       R_ab=fpkey(card,'internuclear_distance',1.d0,' ')
       two_electron=logkey(card,'two_electron_integrals',.false.,' ')
  ELSE
       call lnkerr('no keyword for '//keyword)
  END If
  IF ( form_nabla == .false. .and. form_hamiltonian == .false. ) THEN
       call lnkerr('neither form_nabla or form_hamiltonian task set to true')
  END IF
  Write(iout,1) keyword, spdim, units, mass, typwt, refwt, form_nabla, &
                form_hamiltonian, write_nabla, write_hamiltonian,      &
                read_nabla, read_hamiltonian, diagonalize_nabla,       &
                diagonalize_hamiltonian, poisson, atomic, drctv
1 Format(/,15x,'coordinate system             = ',a16,/15x,            &
               'physical dimensions           = ',i2,/15x,             &
               'units                         = ',a16,/15x,            &
               'mass                          = ',f18.8/15x,           &
               'weighting function            = ',a16,/15x,            &
               'reference weight              = ',a16,/15x,            &
               'compute Nabla                 = ',l1,/15x,             &
               'compute Hamiltonian           = ',l1,/15x,             &
               'write Nabla                   = ',l1,/15x,             &
               'write Hamiltonian             = ',l1,/15x,             &
               'read Nabla                    = ',l1,/15x,             &
               'read Hamiltonian              = ',l1,/15x,             &
               'diagonalize Nabla             = ',l1,/15x,             &
               'diagonalize Hamiltonian       = ',l1,/15x,             &
               'solve poisson equation        = ',l1,/15x,             &
               'is this an atomic calculation = ',l1,/15x,             &        
               'type calculation              = ',a16)
  END SUBROUTINE Set_General_Keywords
!***********************************************************************
!***********************************************************************
!deck Read_Data
!***begin prologue     Read_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Data
  SUBROUTINE Read_Data(grid)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
!
  Call FEDVR_File_Open(grid,'new')
  IF(keyword == 'cartesian') THEN
     Call Read_Cartesian(grid)
  ELSE IF(keyword == 'spherical') THEN
     Call Read_Spherical(grid)
  ELSE IF(keyword == 'cylindrical') THEN
     Call Read_Cylindrical(grid)
  ELSE IF(keyword == 'spheroidal') THEN
     Call Read_Spheroidal(grid)
  ELSE IF(keyword == 'fourier') THEN
     Call Read_Fourier(grid)
  END IF  
  Call IOsys('write integer "number of physical grid points" to '//FEDVR_File,1,physical_points,0,0)
  Call IOsys('write integer "drop left" to '//FEDVR_File,1,drop(1),0,0)
  Call IOsys('write integer "drop right" to '//FEDVR_File,1,drop(2),0,0)
  Call IOsys('write integer "L Max" to '//FEDVR_File,1,l_max,0,' ')
  Call IOsys('write integer "M Max" to '//FEDVR_File,1,m_max,0,' ')
  END SUBROUTINE Read_Data
!***********************************************************************
!***********************************************************************
!deck FEDVR_File_Open
!***begin prologue     FEDVR_File_Open
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       FEDVR_File_Open
  SUBROUTINE FEDVR_File_Open(grid,new_old)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
  CHARACTER(LEN=*)                 :: new_old
  INTEGER, DIMENSION(5)            :: len
!
  len(1)=lenth(keyword)
  len(2)=lenth(grid%label)
!
!        Note a separate file is written for each coordinate
!        to avoid confusion of labels when data is written to IOsys
!
  FEDVR_File = keyword(1:len(1))//'_'//grid%label(1:len(2))
  len(3)=lenth(FEDVR_File)
  write(iout,*) 'Opening Data File = '//FEDVR_File(1:len(3))
  file_loc=File_Directory(4)(1:len_dir(4))//'/'//FEDVR_File(1:len(3))
  len(4)=lenth(file_loc)
  write(iout,*) 'File Location = ',file_loc(1:len(4))
  Call IOsys('open '//FEDVR_File(1:len(3))//' as '//new_old,0,0,0,file_loc(1:len(4)))
  END SUBROUTINE FEDVR_File_Open
!***********************************************************************
!***********************************************************************
!deck Read_Grid_Parameters
!***begin prologue     Read_Grid_Parameters
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Grid_Parameters
  SUBROUTINE Read_Grid_Parameters(grid)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: intkey
  INTEGER                          :: i
  INTEGER                          :: j
  INTEGER                          :: nply
  INTEGER                          :: nfun
  INTEGER                          :: nsubr
  INTEGER                          :: nblock
  INTEGER                          :: begin
  INTEGER                          :: ntest
  INTEGER                          :: n_sub
  LOGICAL                          :: logkey
  LOGICAL                          :: automte
  LOGICAL                          :: dollar
  REAL(idp)                        :: fpkey
  REAL(idp)                        :: step
  REAL(idp)                        :: boundl
  REAL(idp)                        :: boundr
  len=lenth(grid%label)
  automte=logkey(card,'automate',.false.,' ')
  reuse_sector_information=logkey(card,'reuse_sector_information',.false.,' ')
  nfix=intkey(card,'number_of_fixed_points',0,' ')
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  IF(nfix /= 0) THEN
     fix(1)=logkey(card,'left_fixed_point',.false.,' ')
     fix(2)=logkey(card,'right_fixed_point',.false.,' ')
     drop(1)=logkey(card,'drop_left_function',.false.,' ')
     drop(2)=logkey(card,'drop_right_function',.false.,' ')
  END IF
  bcl=1
  bcr=1
  IF(drop(1)) THEN
     bcl=0
  END IF
  IF(drop(2)) THEN
     bcr=0
  END IF
  IF(.not.automte) THEN
     nreg=intkey(card,'number_of_regions',1,' ')
     CALL fparr(card,'region_boundaries',edge,nreg+1,' ')
     CALL intarr(card,'quadrature_order_per_region',n,nreg,' ')
     npt=n+1
     nrq=npt 
     CALL intarr(card,'number_of_reference_quadrature_'//  &
                      'points_per_region',nrq,nreg,' ')
  ELSE
     write(iout,*)
     write(iout,*)
     write(iout,*) '                   Automated Selection of', &
                   ' Steps   '
     nreg=0
     nblock=intkey(card,'number_of_major_blocks',1,' ')       
     DO i=1,nblock 
        IF ( dollar('$'//grid%label(1:len)//'_block_'//itoc(i),card, cpass,inp) ) THEN
             skip=logkey(card,'skip',.false.,' ')
             IF(skip) THEN
                write(iout,2) i
             ELSE 
                nply=intkey(card,'default_order',10,' ')
                nsubr=intkey(card,'number_of_subregions',1,' ')
                boundl=fpkey(card,'left_boundary',0.d0,' ')
                boundr=fpkey(card,'right_boundary',0.d0,' ')
                n_sub=intkey(card,'number_of_reference_quadrature_'//  &
                                  'points_per_sub_region',2*nply,' ')
                write(iout,1) i, nply, nsubr, boundl, boundr
                step = ( boundr - boundl ) / nsubr
                nreg = nreg + 1
                begin = nreg
                edge(nreg)=boundl
                DO j=1,nsubr
                   nrq(nreg) = n_sub
                   nreg = nreg + 1
                   edge(nreg) = edge(nreg-1) + step
                END DO
                nreg = nreg - 1
                n(begin:nreg)=nply
             END IF
        END IF
     END DO     
     npt = n + 1
  END IF
1   FORMAT(/,1x, 'block = ',i3, &
            /,15x,'quadrature order              = ',i4,    &
            /,15x,'number of subregions          = ',i4,    & 
            /,15x,'left hand boundary            = ',e15.8, &
            /,15x,'right hand boundary           = ',e15.8)
2   FORMAT(/,1x,'skipping input block            = ',i4)
  END SUBROUTINE Read_Grid_Parameters
!***********************************************************************
!***********************************************************************
!deck Read_Potential_Data
!***begin prologue     Read_Potential_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for spheroidal dvr basis sets
!***description        sets up the unperturbed dv hamiltonian for a two center
!***                   problem where charge Z_a is at -R/2 and Z_b at +R/2
!***references
!***routines called    
!***end prologue       Read_Potential_Data
  SUBROUTINE Read_Potential_Data (grid)
  IMPLICIT NONE
  TYPE(coordinates)                             :: grid  
  LOGICAL                                       :: dollar
  CHARACTER(LEN=80)                             :: chrkey
  CHARACTER(LEN=3)                              :: itoc
  INTEGER                                       :: i
!
  ALLOCATE(reg_pot(1:nreg))
  IF ( dollar('$potential_'//grid%label,card,cpass,inp) )THEN
!
       write(iout,1)
       i = 1
       reg_pot(1)%type = chrkey(card,'potential_type_region_1','none',' ')
       write(iout,2) i, reg_pot(i)%type
       IF ( nreg > 1 ) THEN
            IF ( reuse == .true. ) THEN
                 reg_pot(2:nreg)%type = reg_pot(1)%type
            ELSE
                 DO i = 2, nreg
                    reg_pot(i)%type = chrkey(card,'potential_type_region_'//itoc(i),'none',' ')        
                    write(iout,2) i, reg_pot(i)%type
                 END DO
            END IF
       END IF
  ELSE
       reg_pot(1:nreg)%type='none'
  END IF
1 FORMAT(/,20x,'Potential Data')
2 FORMAT(2x,'Region = ',i3,2x,'Potential = ',a24)
  END SUBROUTINE Read_Potential_Data
!***********************************************************************
!***********************************************************************
!deck Read_Cartesian
!***begin prologue     Read_Cartesian
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Cartesian
  SUBROUTINE Read_Cartesian(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  len=lenth(grid%label)
  IF ( dollar('$cartesian_'//grid%label(1:len),card,cpass,inp) ) THEN
       atomic=logkey(card,'atomic_calculation',.false.,' ')
       typke='dvr'
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
       proj=.false.
       angmom=intkey(card,'angular_momentum',0,' ')
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       IF (typwt == 'hermite' ) THEN
            Call Read_Hermite(grid)
       ELSE
            Call Read_Grid_Parameters(grid)
       END IF
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label(1:len))
  END If
!
!      Determine nphy
!
!     THere are unused but may be read so zero them
!
  l_max=int_zero
  m_max=int_zero
!
  call ptcal(physical_points,global_points,typwt)
  write(iout,1) typwt, refwt, global_points, physical_points, nfix, fix, drop
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg), grid%reg_pt_wt(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
1 Format(/,15x,'type weight               = ',a16,/15x,                  &
               'reference weight          = ',a16,/15x,                  &
               'number of global points   = ',i5,/15x,                   &
               'number of physcal points  = ',i5,/15x,                   &
               'number of fixed points    = ',i1,/15x,                   &
               'left point fixed          = ',l1,/,15x,                  &
               'right point fixed         = ',l1,/,15x,                  &
               'left point dropped        = ',l1,/15x,                   &
               'right point dropped       = ',l1)
  END SUBROUTINE Read_Cartesian
!***********************************************************************
!***********************************************************************
!deck Read_Spherical
!***begin prologue     Read_Spherical
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Spherical
  SUBROUTINE Read_Spherical(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  len=lenth(grid%label)
  IF ( dollar('$spherical_'//grid%label(1:len),card,cpass,inp) ) THEN
       l_max=intkey(card,'maximum_l_value',0,' ')
       m_max=intkey(card,'maximum_m_value',0,' ')
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
       IF(drctv == 'poisson') THEN
          dentyp=chrkey(card,'type_density','exponential',' ')
       END IF
       proj=.false.
       IF ( typwt == 'laguerre') THEN
            Call Read_Laguerre(grid)
       ELSE IF ( typwt == 'spherical_hermite') THEN
            Call Read_Hermite(grid)
       ELSE IF ( typwt == 'one' .or. typwt == 'legendre') THEN
            Call Read_Grid_Parameters(grid)
       ELSE
            Call lnkerr('error in weighting function '//typwt)
       END IF
!
!      Determine nphy
!
       call ptcal(physical_points,global_points,typwt)
       IF (grid%label(1:len) == 'r') THEN
           write(iout,1) typwt, refwt, global_points, physical_points, nfix, fix, drop, l_max
       ELSE IF (grid%label(1:len) == 'theta') THEN
           write(iout,2) typwt, refwt, global_points, physical_points, nfix, fix, drop, m_max
       END IF
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg), grid%reg_pt_wt(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label(1:len))
  END If
1 Format(/,15x,'type weight              = ',a20,/15x,                  &
               'reference weight         = ',a20,/15x,                  &
               'number of global points  = ',i5,/15x,                   &
               'number of physcal points = ',i5,/15x,                   &
               'number of fixed points   = ',i1,/15x,                   &
               'left point fixed         = ',l1,/,15x,                  &
               'right point fixed        = ',l1,/,15x,                  &
               'left point dropped       = ',l1,/15x,                   &
               'right point dropped      = ',l1,/15x,                   &
               'maximum l value          = ',i4,/15x)
2 Format(/,15x,'type weight              = ',a20,/15x,                  &
               'reference weight         = ',a20,/15x,                  &
               'number of global points  = ',i5,/15x,                   &
               'number of physcal points = ',i5,/15x,                   &
               'number of fixed points   = ',i1,/15x,                   &
               'left point fixed         = ',l1,/,15x,                  &
               'right point fixed        = ',l1,/,15x,                  &
               'left point dropped       = ',l1,/15x,                   &
               'right point dropped      = ',l1,/15x,                   &
               'maximum m value          = ',i4) 
  END SUBROUTINE Read_Spherical
!***********************************************************************
!***********************************************************************
!deck Read_Cylindrical
!***begin prologue     Read_Cylindrical
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Cylindrical
  SUBROUTINE Read_Cylindrical(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  len=lenth(grid%label)
  IF ( dollar('$cylindrical_'//grid%label(1:len),card,cpass,inp) ) THEN
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       l_max=intkey(card,'maximum_l_value',0,' ')
       m_max=intkey(card,'maximum_m_value',0,' ')
       Call Read_Grid_Parameters(grid)
!
!      Determine nphy
!
       call ptcal(physical_points,global_points,typwt)
       write(iout,1) typwt, refwt, global_points, physical_points, nfix, fix, drop, m_max
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg), grid%reg_pt_wt(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label(1:len))
  END If
1 Format(/,15x,'type weight              = ',a20,/15x,                  &
               'reference weight         = ',a20,/15x,                  &
               'number of global points  = ',i5,/15x,                   &
               'number of physcal points = ',i5,/15x,                   &
               'number of fixed points   = ',i1,/15x,                   &
               'left point fixed         = ',l1,/,15x,                  &
               'right point fixed        = ',l1,/,15x,                  &
               'left point dropped       = ',l1,/15x,                   &
               'right point dropped      = ',l1,/15x,                   &
               'maximum m value          = ',i4)  
  END SUBROUTINE Read_Cylindrical
!***********************************************************************
!***********************************************************************
!deck Read_Fourier
!***begin prologue     Read_Fourier
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Fourier
  SUBROUTINE Read_Fourier(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$fourier',card,cpass,inp) ) THEN
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       nfix=intkey(card,'number_of_fixed_points',0,' ')
       fix(1)=.false.
       fix(2)=.false.
       drop(1)=.false.
       drop(2)=.false.
       IF(nfix /= 0) THEN
          fix(1)=logkey(card,'left_fixed_point',.false.,' ')
          fix(2)=logkey(card,'right_fixed_point',.false.,' ')
          drop(1)=logkey(card,'drop_left_function',.false.,' ')
          drop(2)=logkey(card,'drop_right_function',.false.,' ')
       END IF
       bcl=1
       bcr=1
       IF(drop(1)) THEN
          bcl=0
       END IF
       IF(drop(2)) THEN
          bcr=0
       END IF
       nreg=1
       edge(1) = 0.d0
       edge(2) = two * pi
       n(1) = intkey(card,'number_of_points',5,' ')
       npt(1) = n(1) - 2 * ( n(1) / 2 )
       IF(npt(1) == 0 ) THEN
          n(1) = n(1) + 1
       END IF
       npt(1) = n(1)
       physical_points = npt(1)
       global_points = npt(1)
       typwt = 'fourier'
       write(iout,1) typwt, global_points, physical_points, nfix, fix,   &
                     drop, edge(1), edge(2)
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg), grid%reg_pt_wt(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = fourier')
  END IF
  l_max=int_zero
  m_max=int_zero
1 Format(/,15x,'type weight               = ',a16,/15x,                  &
               'number of global points   = ',i5,/15x,                   &
               'number of physcal points  = ',i5,/15x,                   &
               'number of fixed points    = ',i1,/15x,                   &
               'left point fixed          = ',l1,/,15x,                  &
               'right point fixed         = ',l1,/,15x,                  &
               'left point dropped        = ',l1,/15x,                   &
               'right point dropped       = ',l1,/15x,                   &
               'left edge                 = ',e15.8,/15x,                &
               'right edge                = ',e15.8)
  END SUBROUTINE Read_Fourier
!***********************************************************************
!***********************************************************************
!deck Read_Hermite
!***begin prologue     Read_Hermite
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Hermite
  SUBROUTINE Read_Hermite(grid)
  IMPLICIT NONE
  TYPE(coordinates)                     :: grid  
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  LOGICAL                               :: dollar
  CHARACTER(LEN=16)                     :: chrkey
  nreg=1
  nfix=0
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  bcl=1
  bcr=1
  npt(1) = n(1)
  IF (typwt == 'hermite') THEN
      write(iout,1)
      CALL intarr(card,'number_of_points',npt(1),nreg,' ')
  ELSE IF(typwt == 'spherical_hermite') THEN
      write(iout,2)     
      Call Read_Grid_Parameters(grid)
  END IF
1 FORMAT(/,1x,'     Hermite : Region is - infinity to + infinity')
2 FORMAT(/,1x,'     Hermite : Region is zero infinity')
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Read_Hermite
!***********************************************************************
!***********************************************************************
!deck read_laguerre
!***begin prologue     read_laguerre
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_laguerre
  SUBROUTINE Read_Laguerre(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  INTEGER                     :: intkey
  LOGICAL                     :: logkey
  LOGICAL                     :: dollar
  nreg=1
  nfix=intkey(card,'number_of_fixed_points',0,' ')
  CALL intarr(card,'number_of_points',n,nreg,' ')
  CALL fparr(card,'region_boundaries',edge,2,' ')
  npt(1) = n(1)
  IF(nfix == 1) THEN
     fix(1)=logkey(card,'left_fixed_point',.true.,' ')
     drop(1)=logkey(card,'drop_left_function',.false.,' ')
  END IF
  bcl=1
  IF(drop(1)) THEN
     bcl=0
  END IF
!*************************************************************************************
  END SUBROUTINE Read_Laguerre
!***********************************************************************
!***********************************************************************
!deck Read_Theta
!***begin prologue     Read_Theta
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Theta
  SUBROUTINE Read_Theta(grid)
  IMPLICIT NONE
  TYPE(coordinates)                     :: grid  
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  LOGICAL                               :: dollar
  l_max=int_zero
  IF ( dollar('$theta',card,cpass,inp) ) THEN
       nreg=1
       m_max=intkey(card,'maximum_m_value',0,' ')
       write(iout,1) m_val
!
!     This will do what needs to be done for one region
!     If there is more than one region the user needs to call
!     the grid_parameters subroutine and explicitly set things up.
!
       skip=logkey(card,'read_grid_parameters',.false.,' ')
       IF (skip) THEN
           grid%label='theta'
           CALL Read_Grid_Parameters(grid)
       ELSE
           CALL intarr(card,'number_of_points',n,nreg,' ')
           npt(1)=n(1)
           nrq(1)=npt(1) 
           edge(1)=0.d0
           edge(2)=pi
           IF(m_val == 0 ) THEN
              nfix=0
              fix(1)=.false.
              fix(2)=.false.
              drop(1)=.false.
              drop(2)=.false.
              bcl=1
              bcr=1
           ELSE 
              nfix=2
              fix(1)=.true.
              fix(2)=.true.
              drop(1)=.true.
              drop(2)=.true.
              bcl=0
              bcr=0
           END IF 
       END IF
       WRITE(iout,2) npt(1:nreg)
  ELSE
       call lnkerr('no coordinate keyword found for label = theta')
  END IF
  call ptcal(physical_points,global_points,'theta')
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg), grid%reg_pt_wt(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
1 Format(/,1x,'m value = ',i2)
2 FORMAT(/,15x,'    Polynomial order          = ',(/,15x,5(i4,1x)))
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Read_Theta
!***********************************************************************
!***********************************************************************
!deck Read_Spheroidal
!***begin prologue     Read_Spheroidal
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for spheroidal dvr basis sets
!***description        sets up the unperturbed dv hamiltonian for a two center
!***                   problem where charge Z_a is at -R/2 and Z_b at +R/2
!***references
!***routines called    
!***end prologue       Read_Spheroidal
  SUBROUTINE Read_Spheroidal(grid)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
  LOGICAL                          :: dollar
  INTEGER                          :: intkey
  len=lenth(grid%label)
  IF ( dollar('$spheroidal_'//grid%label(1:len),card,cpass,inp) ) THEN
!
       l_max=intkey(card,'maximum_l_value',0,' ')
       m_max=intkey(card,'maximum_m_value',0,' ')
!
       Call Read_Grid_Parameters(grid)
!
!      Determine nphy
!
       call ptcal(physical_points,global_points,'spheroidal')
!      Compute the Lobatto regional polynomials for even and odd types.
!
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg), grid%reg_pt_wt(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label)
  END IF
  write(iout,1) typwt, refwt, global_points, physical_points, nfix, fix, &
                drop, m_max, z_a, z_b, R_ab
1 Format(/,15x,'type weight               = ',a16,/15x,                  &
               'reference weight          = ',a16,/15x,                  &
               'number of global points   = ',i5,/15x,                   &
               'number of physcal points  = ',i5,/15x,                   &
               'number of fixed points    = ',i1,/15x,                   &
               'left point fixed          = ',l1,/,15x,                  &
               'right point fixed         = ',l1,/,15x,                  &
               'left point dropped        = ',l1,/15x,                   &
               'right point dropped       = ',l1,/15x,                   &
               'maximum m value           = ',i4,/15x,                   &
               'charge on nucleus A       = ',f10.5,/15x,                & 
               'charge on nucleus b       = ',f10.5,/15x,                &
               'internuclear distance     = ',e15.8)
  END SUBROUTINE Read_Spheroidal
!***********************************************************************
!***********************************************************************
!deck ptcal
!***begin prologue     ptcal
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            determine the size of the gid with all constraints
!***description
!***references
!***routines called    
!***end prologue       ptcal
  SUBROUTINE ptcal(nphy,nglobal,typwt)
  IMPLICIT NONE
  INTEGER                    :: nphy
  INTEGER                    :: nglobal
  INTEGER                    :: i
  CHARACTER(LEN=*)           :: typwt
  IF(typwt == 'hermite') THEN
     nphy = npt(1)
     nglobal=nphy
  ELSE IF(typwt == 'laguerre') THEN
     nphy = npt(1)
     nglobal=nphy
     IF(bcl == 0) THEN
        nphy=nphy-1
     END IF
  ELSE IF(typwt == 'fourier') THEN
     nphy=npt(1)
     nglobal=nphy
  ELSE
     nphy=0
     DO  i=1,nreg
!
!        number of internal functions
! 
         nphy = nphy + npt(i) - 2
     END DO
! 
!     add one bridge function between each interval.
!
     nphy = nphy + nreg - 1
! 
!     add the extreme left and extreme right points
!     we have not yet dropped any functions at the endpoints.
!
     nphy=nphy + 2
     nglobal=nphy
     IF(bcl == 0) THEN
        nphy=nphy-1
     END IF
     IF(bcr == 0) THEN
        nphy=nphy-1
     END IF
  END IF
END SUBROUTINE ptcal
!***********************************************************************
!***********************************************************************
           END MODULE Read_DVR_Module
!***********************************************************************
!***********************************************************************
