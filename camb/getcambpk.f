!!!! getcambpk.f   Time-stamp: <2017-12-22 08:50:15 hamana>

! based on getpkmattercamb.F90

! gfortran -O3 -fopenmp -ffast-math -fmax-errors=4 -march=native -JRelease -IRelease/  -c getcambpk.f -o Release/getcambpk.o


      module getmatterpkcamb

      contains
      
      subroutine getcambpk(par_s8,par_ts100,par_h0,par_oc,par_ob,par_ov,
     $     par_on,par_as,par_ns,par_tu,par_ok,par_om,
     $     nzdim,nkdim,nkmax,cambpk,par_zdat)

      use IniFile
      use CAMB
      use LambdaGeneral
      use Lensing
      use AMLUtils
      use Transfer
      use constants
      use Bispectrum
      use CAMBmain
      use NonLinear
      use KeepMatterPower

      implicit none

      Type(CAMBparams) P

      real par_s8,par_ts100,par_h0,par_oc,par_ob,par_ov,par_on,
     $     par_as,par_ns,par_tu,par_ok,par_om

      integer nzdim,nkdim,nkmax
      real cambpk(nkmax,0:nzdim)
      real par_zdat(nzdim)
      
      character(LEN=Ini_max_string_len) numstr
      integer i,j


!!! parameters from recfast.f90
      real(dl), parameter :: RECFAST_fudge_default = 1.14_dl !1.14_dl
      real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 + 0.02d0
      real(dl) :: AGauss1 =      -0.14D0 !Amplitude of 1st Gaussian
      real(dl) :: AGauss2 =       0.079D0 ! 0.05D0  !Amplitude of 2nd Gaussian
      real(dl) :: zGauss1 =       7.28D0 !ln(1+z) of 1st Gaussian
      real(dl) :: zGauss2=        6.73D0 !ln(1+z) of 2nd Gaussian
      real(dl) :: wGauss1=        0.18D0 !Width of 1st Gaussian
      real(dl) :: wGauss2=        0.33D0 !Width of 2nd Gaussian
    
      real, dimension(:,:), allocatable :: keepmatterpk

      character(LEN=Ini_max_string_len) outroot, version_check
      real(dl) output_factor, nmassive

!      logical bad

      integer narg
      character*256 opt,arg

!      Ini_fail_on_not_found = .false.

!      outroot ='getcambpk'      !!!Ini_Read_String('output_root')

      highL_unlensed_cl_template=
     $     'HighLExtrapTemplate_lenspotentialCls.dat' !!! Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)

      call CAMB_SetDefParams(P)

      P%WantScalars =.false.    !!! Ini_Read_Logical('get_scalar_cls')
      P%WantVectors =.false.    !!! Ini_Read_Logical('get_vector_cls',.false.)
      P%WantTensors =.false.    !!! Ini_Read_Logical('get_tensor_cls',.false.)

      P%OutputNormalization=outNone
      output_factor = 1.0d0     !!!Ini_Read_Double('CMB_outputscale',1.d0)

      P%WantCls=.false.         !!! P%WantScalars .or. P%WantTensors .or. P%WantVectors

      P%PK_WantTransfer=.true.  !!!Ini_Read_Logical('get_transfer')

      AccuracyBoost  =1.0d0     !!! Ini_Read_Double('accuracy_boost',AccuracyBoost)
      lAccuracyBoost =1.0       !!! Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
      HighAccuracyDefault =.true. !!! Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

      P%NonLinear =1            !!! Ini_Read_Int('do_nonlinear',NonLinear_none)

      P%DoLensing = .false.

    !  Read initial parameters.

!
!!!!!!!!!!!!!!!!!!! read options
!      
c      par_h0=70.0               ! H0
c      par_ob=0.0462             !Omega_b
c      par_oc=0.2538             !Omega_CDM
c      par_ov=0.698694           !Omega_Lambda
c      par_on=0.001306           !Omega_nutrino
c      par_as=2.0e-9             !A_s
c      par_ns=0.96               !n_s
c      par_tu=0.09               !optical depth tau
c
c      narg=iargc()
c      do i=1,narg
c         call getarg(i,opt)
c         call getarg(i+1,arg)
c         select case (opt)
c         case ('-h0')
c            read(arg,*)par_h0
c         case ('-oc')
c            read(arg,*)par_oc
c         case ('-ob')
c            read(arg,*)par_ob
c         case ('-ov')
c            read(arg,*)par_ov
c         case ('-on')
c            read(arg,*)par_on
c         case ('-as')
c            read(arg,*)par_as
c         case ('-ns')
c            read(arg,*)par_ns
c         case ('-tau')
c            read(arg,*)par_tu
c         end select
c      enddo


!!!    call DarkEnergy_ReadParams(DefIni)       
      w_lam =-1.0d0             !!!Ini_Read_Double_File(Ini,'w', -1.d0)
      cs2_lam =1.0d0            !!!Ini_Read_Double_File(Ini,'cs2_lam',1.d0)

      P%h0     = par_h0         !!! Ini_Read_Double('hubble')

      P%omegab = par_ob        !!! Ini_Read_Double('omega_baryon')
      P%omegac = par_oc        !!!  Ini_Read_Double('omega_cdm')	
      P%omegav = par_ov      !!!  Ini_Read_Double('omega_lambda	')
      P%omegan = par_on      !!! Ini_Read_Double('omega_neutrino')
      
!    if (Ini_Read_Logical('use_physical',.false.)) then
!        P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
!        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
!        P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
!        P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
!    else
!        P%omegab = Ini_Read_Double('omega_baryon')
!        P%omegac = Ini_Read_Double('omega_cdm')
!        P%omegav = Ini_Read_Double('omega_lambda')
!        P%omegan = Ini_Read_Double('omega_neutrino')
!    end if

      P%tcmb   =2.7255d0        !!! Ini_Read_Double('temp_cmb',COBE_CMBTemp)
      P%yhe    =0.24d0          !!! Ini_Read_Double('helium_fraction',0.24_dl)
      P%Num_Nu_massless =2.046d0 !!! Ini_Read_Double('massless_neutrinos')
      
      P%Nu_mass_eigenstates =1  !!! Ini_Read_Int('nu_mass_eigenstates',1)
      if (P%Nu_mass_eigenstates > max_nu) 
     $     error stop 'too many mass eigenstates'

      numstr ='1'               !!! Ini_Read_String('massive_neutrinos')
      read(numstr, *) nmassive
      read(numstr,*, end=100, err=100) 
     $     P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
      P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))
      
      if (P%Num_Nu_massive>0) then
         P%share_delta_neff =.true. !!! Ini_Read_Logical('share_delta_neff', .true.)
         numstr =''             !!! Ini_Read_String('nu_mass_degeneracies')
         if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 
     $    'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
         else
            if (numstr=='') error stop 
     $'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) 
     $           P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
         endif
         numstr ='1'            !!! Ini_Read_String('nu_mass_fractions')
         if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) error stop 
     $      'must give nu_mass_fractions for the eigenstates'
            P%Nu_mass_fractions(1)=1
         else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
         end if
      end if

      if (((P%NonLinear==NonLinear_lens .or. 
     $     P%NonLinear==NonLinear_both) .and. P%DoLensing) 
     $     .or. P%PK_WantTransfer) then
         P%Transfer%high_precision=.false. !!!  Ini_Read_Logical('transfer_high_precision',.false.)
      else
         P%transfer%high_precision = .false.
      endif

      DebugParam =0.0d0         !!! Ini_Read_Double('DebugParam',DebugParam)
      ALens =1.0d0              !!! Ini_Read_Double('Alens',Alens)

!!!    call Reionization_ReadParams(P%Reion, DefIni)
      P%Reion%Reionization =.true. !!! Ini_Read_Logical_File(Ini,'reionization')
      if (P%Reion%Reionization) then
         P%Reion%use_optical_depth =.true. !!! Ini_Read_Logical_File(Ini,'re_use_optical_depth')
         if (P%Reion%use_optical_depth) then
!!!!!!!!!!! need parametrization !!!!!!!!!!!!!!!!!!!!
            P%Reion%optical_depth = par_tu !0.09d0 !!! Ini_Read_Double_File(Ini,'re_optical_depth')
         else
            P%Reion%redshift =11.0d0 !!! Ini_Read_Double_File(Ini,'re_redshift')
         end if

         P%Reion%delta_redshift =1.5d0 !!! Ini_Read_Double_File(Ini,'re_delta_redshift', 0.5_dl) !default similar to CMBFAST original
         P%Reion%fraction =-1.0 !!! Ini_Read_Double_File(Ini,'re_ionization_frac',Reionization_DefFraction)
         
         P%Reion%helium_redshift  =3.5d0 !!! Ini_Read_Double_File(Ini,'re_helium_redshift', 3.5_dl)
         P%Reion%helium_delta_redshift  =0.5d0 !!! Ini_Read_Double_File(Ini,'re_helium_delta_redshift', 0.5_dl)
         P%Reion%helium_redshiftstart  =  P%Reion%helium_redshift 
     $        + 3*P%Reion%helium_delta_redshift !!! Ini_Read_Double_File(Ini,'re_helium_redshiftstart', P%Reion%helium_redshift + 3*P%Reion%helium_delta_redshift)
      endif

!!!    call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
      P%InitPower%k_0_scalar =0.05 !!! Ini_Read_Double_File(Ini,'pivot_scalar',P%InitPower%k_0_scalar)
      P%InitPower%k_0_tensor =0.05 !!! Ini_Read_Double_File(Ini,'pivot_tensor',P%InitPower%k_0_tensor)
      P%InitPower%nn =1         !!! Ini_Read_Int_File(Ini,'initial_power_num',1)
      P%InitPower%rat(:) = 1
      do i=1, P%InitPower%nn
!!!!!!!!!!! need parametrization !!!!!!!!!!!!!!!!!!!!
         P%InitPower%an(i) =par_ns ! 0.96d0 !!! Ini_Read_Double_Array_File(Ini,'scalar_spectral_index', i)
         P%InitPower%n_run(i) =0.0d0 !!! Ini_Read_Double_Array_File(Ini,'scalar_nrun',i,0._dl)
         P%InitPower%n_runrun(i) =0.0d0 !!! Ini_Read_Double_Array_File(Ini,'scalar_nrunrun',i,0._dl)
         
!!!!!!!!!!! need parametrization !!!!!!!!!!!!!!!!!!!!
         P%InitPower%ScalarPowerAmp(i) =par_as ! 2.1d-9 !!! Ini_Read_Double_Array_File(Ini,'scalar_amp',i,1._dl)
!Always need this as may want to set tensor amplitude even if scalars not computed
      enddo
      
!!!   call Recombination_ReadParams(P%Recomb, DefIni)
      P%Recomb%RECFAST_fudge_He =0.86d0 !!! Ini_Read_Double_File(Ini,'RECFAST_fudge_He',RECFAST_fudge_He_default)
      P%Recomb%RECFAST_Heswitch =6 !!! Ini_Read_Int_File(Ini, 'RECFAST_Heswitch',RECFAST_Heswitch_default)
      P%Recomb%RECFAST_Hswitch =.true. !!! Ini_Read_Logical_File(Ini, 'RECFAST_Hswitch',RECFAST_Hswitch_default)
      P%Recomb%RECFAST_fudge =1.14 !!! Ini_Read_Double_File(Ini,'RECFAST_fudge',RECFAST_fudge_default)
      if (P%Recomb%RECFAST_Hswitch) then
         P%Recomb%RECFAST_fudge = P%Recomb%RECFAST_fudge - 
     $        (RECFAST_fudge_default - RECFAST_fudge_default2)
      end if

      Ini_fail_on_not_found = .false.
      
    !optional parameters controlling the computation

      P%AccuratePolarization =.true. !!! Ini_Read_Logical('accurate_polarization',.true.)
      P%AccurateReionization =.true. !!! Ini_Read_Logical('accurate_reionization',.false.)
      P%AccurateBB =.false.     !!! Ini_Read_Logical('accurate_BB',.false.)
      P%DerivedParameters =.true. !!! Ini_Read_Logical('derived_parameters',.true.)

      version_check ='Jan17'    !!! Ini_Read_String('version_check')

      if (HighAccuracyDefault) then
         DoTensorNeutrinos = .true.
      else
         DoTensorNeutrinos = .true. !!! Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
      end if
      FeedbackLevel =0          ! 1 for parameter output to display !!! Ini_Read_Int('feedback_level',FeedbackLevel)

      output_file_headers =.true. !!! Ini_Read_Logical('output_file_headers',output_file_headers)

      P%MassiveNuMethod  =1     !!! Ini_Read_Int('massive_nu_approx',Nu_best)

      ThreadNum      = 0        !!!Ini_Read_Int('number_of_threads',ThreadNum)
      use_spline_template =.true. !!! Ini_Read_Logical('use_spline_template',use_spline_template)


!!!!!!!!!!!!!!!! need set-up begin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (P%PK_WantTransfer)  then
         P%WantTransfer  = .true.
         P%transfer%kmax          =5.0d0 !!!  Ini_Read_Double('transfer_kmax')
         P%transfer%k_per_logint  =0 !!!  Ini_Read_Int('transfer_k_per_logint')
         P%transfer%PK_num_redshifts =nzdim !!!  Ini_Read_Int('transfer_num_redshifts')
         
         transfer_interp_matterpower =.false. !!! Ini_Read_Logical('transfer_interp_matterpower', transfer_interp_matterpower)
!which variable to use for defining the matter power spectrum and sigma8
!main choices are 2: CDM, 7: CDM+baryon+neutrino, 8: CDM+baryon, 9: CDM+baryon+neutrino+de perts
         transfer_power_var =7  !!! Ini_read_int('transfer_power_var',transfer_power_var)
         if (P%transfer%PK_num_redshifts > max_transfer_redshifts) 
     $        error stop 'Too many redshifts'
!!!!!!!!!!!!!! need input-scheme for redshifs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do i=1,nzdim
            P%transfer%PK_redshifts(i)=par_zdat(i)
         enddo
!
         P%transfer%PK_num_redshifts=P%transfer%PK_num_redshifts+1
         P%transfer%PK_redshifts(P%transfer%PK_num_redshifts)=0.0d0
!!!!!!!!!!!!!! need input-scheme for redshifs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        do i=1, P%transfer%PK_num_redshifts
!!!            P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
!!!            transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
!!!            MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
!!!            if (TransferFileNames(i) == '') then
!!!                TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
!!!            end if
!!!            if (MatterPowerFilenames(i) == '') then
!!!                MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
!!!            end if
!!!            if (TransferFileNames(i)/= '') &
!!!                TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
!!!            if (MatterPowerFilenames(i) /= '') &
!!!                MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
!!!        end do
!!!!!!!!!!!!!!!! need set-up end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
         P%Transfer%PK_num_redshifts = 1
         P%Transfer%PK_redshifts = 0
      end if

      call Transfer_SortAndIndexRedshifts(P%Transfer)

      P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)
      
      Ini_fail_on_not_found = .false.
      

      call Ini_Close

      if (.not. CAMB_ValidateParams(P)) 
     $     error stop 'Stopped due to parameter error'

      if (global_error_flag==0) call CAMB_GetResults(P)
      if (global_error_flag/=0) then
         write(*,*) 'Error result '//trim(global_error_message)
         error stop
      endif
      
      if (P%PK_WantTransfer) then
!         call Transfer_SaveToFiles(MT,TransferFileNames)
!         call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
	 call Transfer_GetMatterPowerDim(MT)
!         write(*,*)pkmatterdim,pkredshiftsdim
         pkredshiftsdim=pkredshiftsdim         
	 allocate(keepmatterpk(pkmatterdim,0:pkredshiftsdim))
         call Transfer_KeepMatterPower(MT,keepmatterpk)	
!         call Transfer_output_sig8(MTk)
      endif
    
!     do i=1,P%Transfer%PK_num_redshifts
!         write(*,*)MT%sigma_8(i,1)
!     enddo

      nkdim=pkmatterdim
      do i=0,nzdim
         do j=1,nkdim            
            cambpk(j,i)=keepmatterpk(j,i)
         enddo
      enddo

c      write(*,'("# Om_b                = ",f9.6)') par_ob
c      write(*,'("# Om_c h^2            = ",f9.6)') par_oc
c      write(*,'("# Om_nu h^2           = ",f9.6)') par_on
c      write(*,'("# Om_Lambda           = ",f9.6)') par_ov
      par_ok=P%omegak
c      write(*,'("# Om_K                = ",f9.6)') par_ok
      par_om=1-P%omegak-CP%omegav
c      write(*,'("# Om_m (1-Om_K-Om_L)  = ",f9.6)') par_om
c      write(*,'("# Scalar amplitude    = ",1pe13.6)') par_as
c      write(*,'("# Spectral index      = ",f9.6)') par_ns
      par_ts100=100*CosmomcTheta()
c      write(*,'("# 100 theta (CosmoMC) = ",f9.6)') par_ts100
      par_s8=MT%sigma_8(P%Transfer%PK_num_redshifts,1)
c      write(*,'("# sigma8 at z=0       = ",f9.6)') par_s8
      
      call CAMB_cleanup
      
      return
    
 100  stop 'Must give num_massive number 
     $     of integer physical neutrinos for each eigenstate'
      end subroutine getcambpk

      end module getmatterpkcamb
