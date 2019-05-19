    module wlhsc
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use IniObjects
    use MatrixUtils
    use MiscUtils
    implicit none
    private

    !------------------------------------------------------------
    ! SM: This block has all private variables required by the
    ! subroutines in the likelihood function
    !------------------------------------------------------------
    integer npzbin       ! Number of tomographic redshift bins
    parameter (npzbin=4) ! Number of tomographic redshift bins
    real pdftot(npzbin)  ! Array to normalize pdf
    real zmean(npzbin)   ! Array with mean redshift of sources
    
    integer nspdf,nspdfmax    ! Array size for P(z)
    parameter (nspdfmax=1024)
    real delnspdf,delnspdf2   ! Deltaz for P(z), Deltaz/2
    real znspdf0(nspdfmax)    ! Redshifts in the P(z) array
    real pnspdf0(npzbin,nspdfmax) ! P(z) per tomographic bin
    real znspdf(npzbin,nspdfmax)  ! Redshifts in the P(z+deltaz) array
    real pnspdf(npzbin,nspdfmax)  ! Modified P(z+deltaz)
    real*8 db0,db1,db2,db3,db4    ! Variables to read in the file for P(z)s

    ! These are cosmological variables that will be shared with other
    ! routines. This is legacy from Hamana-san's previous code
    ! structure.
    ! Omegam, OmegaLambda, Omega_curvature, sqrt(|Omega_curvature|)
    real cosmo_om,cosmo_ov,cosmo_crv,cosmo_scrv

    ! Parameters for neutrino mass, these are now handled by the
    ! cosmology module
    ! Neutrino mass, Neutrino fraction
    real par_mneu,fneu     

    ! CAMB related parameters, These could in principle be borrowed
    ! from theory
    integer nkdim,nkmax
    ! Sigma8, T/S ratio, h0, Omegacdm, Omegab, OmegaLambda, Omega_nu,
    ! As, ns, tau, Omegak, Omegam
    real par_sig8,par_ts100,par_h0,par_oc,par_ob,par_ov,par_on,par_as,par_ns,par_tu,par_ok,par_om

    ! Parameters for FFTlog
    integer nsbin
    real*8 slogmax,slogmin

    ! Parameters for IA model
    real par_aia,par_etaia,par_z0ia

    ! Parameters for PSF error
    real par_apsf,par_bpsf

    ! Derived parameters, I have not yet put these in. My guess is they
    ! go in an array called derived.
    real par_s8,par_s845

    ! Nuisance parameters to account for P(z) uncertainties
    ! Deltam, Deltaz1, Deltaz2, Deltaz3, Deltaz4
    real par_delm,par_delz1,par_delz2,par_delz3,par_delz4

    ! Correction for constant-e (xi_+ only)
    real par_conste

    ! Variable describing size of two point correlation function
    integer ntpcf,ndvec      
    parameter (ntpcf=10) ! number of Auto/cross TPCF, 11,12,13,14,22,23,24,33,34,44

    ! Allocate these arrays describing the data only once when the HSC
    ! WL likelihood is added to list of likelihoods to calculate
    ! Model vector, data vector, (model-data), (model-data)^TCinv,
    ! data vector uncorrected for systematic issues
    real, allocatable, dimension(:)::modelvec,datavec,delvec,tmpvec,datavec0 !(ndvec)
    integer, allocatable, dimension(:)::datavecflag
    ! I suppose this is actually an inverse covariance matrix
    real, allocatable, dimension(:,:)::covmat !(ndvec,ndvec)      
    real, allocatable, dimension(:,:)::psfvec !(ndvec,3)


    ! Number of theta bins, theta of each bin
    real tbin
    real, allocatable:: thetabin(:)
    ! File describing the data vectors
    character(LEN=:), allocatable :: fdatavec,fdatavecflag,fcovmat,fnspdf,ftbinin,fpsfvec
    integer ntbin,ntp1,ntp2,ntm1,ntm2
    real cor_z_sigmae(4)

    !------------------------------------------------------------
    ! SM: The block above has all private variables required by the
    ! subroutines in the likelihood function
    !------------------------------------------------------------

    !------------------------------------------------------------
    ! SM: This defines the type of our likelihood TCosmoCalcLikelihood
    ! can borrow the CMB parameters, as well as the theory
    ! calculations of sigma8 from the main modules
    !------------------------------------------------------------
    type, extends(TCosmoCalcLikelihood) :: HSCWL_Likelihood
    contains
    procedure :: LogLike => HSCWL_lnlike
    end type HSCWL_Likelihood

    PUBLIC :: HSCWL_Likelihood_Add, HSCWL_lnlike

    contains

    !------------------------------------------------------------
    ! SM: This subroutined allows the addition of the likelihood code
    ! to the list of likelihoods from cosmomc
    !------------------------------------------------------------
    subroutine HSCWL_Likelihood_Add(LikeList, Ini)
    implicit none
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(HSCWL_Likelihood), pointer :: this

    ! Check that use_HSCWL is true in the corresponding ini file
    if (Ini%Read_Logical('use_HSCWL',.false.)) then
        allocate(this)
        call this%loadParamNames(trim(DataDir)//'HSCWL.paramnames')
        call LikeList%Add(this)
        this%LikelihoodType = 'HSCWL'
        this%name='HSCWL'
        this%needs_background_functions = .true.
        this%needs_powerspectra = .false.
        this%needs_sigmaR = .true.

        !------------------------------------------------------------
        ! Initialize any constant variables from ini file
        !------------------------------------------------------------
        ntbin = Ini%Read_Int('ntbin',15)
        write(*,*)'ntbin = ',ntbin
        ntp1 = Ini%Read_Int('ntp1',10)
        write(*,*)'ntp1 = ',ntp1
        ntm1 = Ini%Read_Int('ntm1',13)
        write(*,*)'ntm1 = ',ntm1
        ntp2 = Ini%Read_Int('ntp2',23)
        write(*,*)'ntp2 = ',ntp2
        ntm2 = Ini%Read_Int('ntm2',23)
        write(*,*)'ntm2 = ',ntm2

        ! Read in correction factors for redshift dependent sigma_e
        cor_z_sigmae(1) = Ini%Read_Real('cor_z_sigmae1',1.0)
        write(*,*)'cor_z_sigmae(1) = ',cor_z_sigmae(1)
        cor_z_sigmae(2) = Ini%Read_Real('cor_z_sigmae2',1.0)
        write(*,*)'cor_z_sigmae(2) = ',cor_z_sigmae(2)
        cor_z_sigmae(3) = Ini%Read_Real('cor_z_sigmae3',1.015)
        write(*,*)'cor_z_sigmae(3) = ',cor_z_sigmae(3)
        cor_z_sigmae(4) = Ini%Read_Real('cor_z_sigmae4',1.03)
        write(*,*)'cor_z_sigmae(4) = ',cor_z_sigmae(4)

        ! Read in filename containing the data vector flags
        fdatavecflag = Ini%ReadFileName('file_datavecflag', NotFoundFail = .true.)
        if (fdatavecflag.eq.'') then
           write(*,*)'file_datavecflag not exist'
           write(*,*)'STOP'
        endif
        write(*,*)'file_datavecflag = ',fdatavecflag

        ! Read in filename containing the data vector
        fdatavec = Ini%ReadFileName('file_datavec', NotFoundFail = .true.)
        if (fdatavec.eq.'') then
           write(*,*)'file_datavec not exist'
           write(*,*)'STOP'
        endif
        write(*,*)'file_datavec = ',fdatavec

        ! Read in filename containing the covariance matrix
        fcovmat = Ini%ReadFileName('file_covmat', NotFoundFail = .true.)
        if (fcovmat.eq.'') then
           write(*,*)'file_covmat not exist'
           write(*,*)'STOP'
        endif
        write(*,*)'file_covmat = ',fcovmat

        ! Read in filename containing the P(z)s
        fnspdf = Ini%ReadFileName('file_nspdf', NotFoundFail = .true.)
        if (fnspdf.eq.'') then
           write(*,*)'file_nspdf not exist'
           write(*,*)'STOP'
        endif
        write(*,*)'file_nspdf = ',fnspdf

        ! Read in filename containing the theta bins
        ftbinin = Ini%ReadFileName('file_tbinin', NotFoundFail = .true.)
        if (ftbinin.eq.'') then
           write(*,*)'file_tbinin not exist'
           write(*,*)'STOP'
        endif
        write(*,*)'file_tbinin = ',ftbinin

        ! Read in filename containing the PSF systematics corrections
        fpsfvec = Ini%ReadFileName('file_psfvec', NotFoundFail = .true.)
        if (fpsfvec.eq.'') then
           write(*,*)'file_psfvec not exist'
           write(*,*)'STOP'
        endif
        write(*,*)'file_psfvec = ',fpsfvec

        ! Now read in all the relevant files
        CALL HSCWL_init

    endif
    end subroutine HSCWL_Likelihood_Add

    !------------------------------------------------------------
    ! This function is the main likelihood function and should return
    ! logLike, so chisq/2
    !------------------------------------------------------------
    function HSCWL_lnlike(this,CMB,Theory,DataParams)
    Class(HSCWL_Likelihood) :: this
    Class (CMBParams):: CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) DataParams(:)
    real(mcp)  HSCWL_lnlike

    !********* Now here we should have all the variables required for wltpcf

    ! These variables here are old variables that are not really
    ! required, and can be commented out
    ! real*8 prm_h0,prm_ob,prm_oc,prm_ov,prm_on,prm_om
    ! real*8 prm_as9,prm_ns,prm_tu,prm_sig8,prm_ts100,prm_ok
    ! real*8 prm_aia,prm_etaia,prm_z0ia
    ! real*8 prm_s8,prm_s845
    ! real*8 prm_delm,prm_delz1,prm_delz2,prm_delz3,prm_delz4
    ! real*8 prm_apsf,prm_bpsf
    ! real*8 prm_ce !!! TH20190130
    ! real*8 prm_mneu
    ! The once above may not be really required, but we will see! SM SM

    !real tloglikeli
    real*8 chisq
    integer i, j, k

    !
    ! Map internal variables to cosmological parameters
    ! and to other parameters required for
    ! likelihood calculation
    par_oc = CMB%omc ! Omegacdm
    par_as = 1e-10*(CMB%InitPower(As_index)) ! SMSM: DANGER: check docs for what InitPower stores, I think it is (As/1e10)
    par_h0=CMB%H0 !H0
    par_ob=CMB%omb  ! Omega_b
    par_on = CMB%omnu ! Omega_nu
    fneu=par_on/(par_oc+par_ob+par_on) ! frac_neu

    par_ns=CMB%InitPower(ns_index) !ns
    par_tu=CMB%tau !optical depth tau
    par_ok=CMB%omk ! Omega_curvature


    !Map the nuisance parameters
    par_aia=DataParams(1)
    par_etaia=DataParams(2)
    par_z0ia=DataParams(3)

! nuisance parameters
    par_delm=(1.0+0.01*DataParams(4))**2 ! note par_delm = (1.0 + 0.01 * prm_delm)**2
    
    par_delz1=0.01*DataParams(5) 
    par_delz2=0.01*DataParams(6) 
    par_delz3=0.01*DataParams(7) 
    par_delz4=0.01*DataParams(8) 

!!! TH20180802 correction for constant-e (xi_+ only)
    par_conste=1.0e-4*DataParams(9)
    par_conste=par_conste*par_conste

    par_apsf=0.01*DataParams(10) ! note par_apsf = 0.01 * prm_apsf
    par_bpsf=0.01*DataParams(11) ! note par_bpsf = 0.01 * prm_bpsf

    write(*, *) "Please remove the below 4 lines after debugging is complete"
    write(*, *) "par_oc, par_as, par_h0, par_ob, par_on, par_ns, par_tu, par_ok"
    write(*, *) par_oc, par_as, par_h0, par_ob, par_on, par_ns, par_tu, par_ok
    write(*, *) "par_aia, par_etaia, par_z0ia, par_delm, par_delz1, par_delz2, par_delz3, par_delz4, par_conste, par_apsf, par_bpsf"
    write(*, *) par_aia, par_etaia, par_z0ia, par_delm, par_delz1, par_delz2, par_delz3, par_delz4, par_conste, par_apsf, par_bpsf

    !-----------------------------------------------------------
    ! derived parameters: Do not know if this is really required and
    ! whether these can be output by the chain. My guess is we have to
    ! use the array derived
    !-----------------------------------------------------------
    par_ov=CMB%omv
    par_sig8=Theory%sigma_8
    par_om=CMB%omc+CMB%omb+CMB%omnu
    par_ts100=0.0
    par_s8=Theory%sigma_8*(par_om/0.3)**0.5
    par_s845=Theory%sigma_8*(par_om/0.3)**0.45

    ! correction for PSF error xi => xi - par_apsf**2 * xi_pp
    !                                - 2 * par_apsf * par_bpsf * xi_pq
    !                                - par_bpsf**2 * xi_qq
    do i=1,ndvec
       datavec(i)=datavec0(i)-par_apsf**2*psfvec(i,1)-2.0*par_apsf*par_bpsf*psfvec(i,2)-par_bpsf**2*psfvec(i,3)
    enddo

    !!! Get the Photoz PDF by shifting it by Deltazi
    ! Reassign mean and integral P(z)
    do j=1,npzbin
       pdftot(j)=0.0
       zmean(j)=0.0
    enddo

    ! Shift the redshifts of the pdf
    do i=1,nspdf
       znspdf(1,i)=znspdf0(i)+par_delz1
       znspdf(2,i)=znspdf0(i)+par_delz2
       znspdf(3,i)=znspdf0(i)+par_delz3
       znspdf(4,i)=znspdf0(i)+par_delz4
       !
       ! Assign the redshift non-negative prior
       !
       do k = 1, 4
           if (znspdf(k,i).le.0.0) then
               pnspdf(k,i)=0.0
           else
               pnspdf(k,i)=pnspdf0(k, i)
           endif
       enddo
       do j=1,npzbin
          pdftot(j)=pdftot(j)+pnspdf(j,i)*delnspdf
          zmean(j)=zmean(j)+znspdf(j,i)*pnspdf(j,i)*delnspdf
       enddo
    enddo

    ! Assign mean redshift
    do j=1,npzbin
       zmean(j)=zmean(j)/pdftot(j)
    enddo      

    ! Normalize the PDF
    do i=1,nspdf
       do j=1,npzbin
          pnspdf(j,i)=pnspdf(j,i)/pdftot(j)
       enddo
    enddo

    ! setup parameters for calling getcambpklin, could these go in modelvec, perhapsDataParams(1)
    nkdim=0
    nkmax=1024

    nsbin=2048
    slogmin=-4.0
    slogmax=8.0

    ! This call will write up the datavec, no parameters need to be
    ! given as all parameters are private to this module and can be
    ! accessed from make_modelvec
    call make_modelvec

    ! Some debugging statements to verify the datavectors, and
    ! modelvectors
    !do i=1, ndvec
    !    write(*, *) ">>> Is it this one?:", datavec(i), modelvec(i)
    !enddo

    ! Assign Deltvec = theoryvec-modelvector(1+deltam)
    do i=1,ndvec
       delvec(i)=float(datavecflag(i))*(datavec(i)-modelvec(i)*par_delm) ! note par_delm = (1.0 + 0.01 * prm_delm)**2
    enddo
   
    ! Deltavec^T C^-1 
    do j=1,ndvec
       tmpvec(j)=0.0
       do i=1,ndvec
          tmpvec(j)=tmpvec(j)+delvec(i)*covmat(i,j)
       enddo
    enddo

    ! Deltavec^T C^-1 Deltavec
    chisq=0.0
    do i=1,ndvec
       chisq=chisq+tmpvec(i)*delvec(i)
    enddo
    ! Debug HSCWL likelihood
    !write(*, *)"HSCWL:", chisq
    ! log Likelihood is chisq/2
    HSCWL_lnlike = 0.5*chisq

    end function HSCWL_lnlike

      subroutine HSCWL_init
      integer i, j, k, indvec
      !----------------------------------------------------
      ! This subroutine should initialize all data vectors and
      ! covariance matrices
      !----------------------------------------------------

      ! Allocate vectors
      ndvec=ntpcf*(ntp2-ntp1+1+ntm2-ntm1+1)
      allocate(modelvec(ndvec))
      allocate(datavec(ndvec))
      allocate(datavec0(ndvec))
      allocate(datavecflag(ndvec))
      allocate(delvec(ndvec))
      allocate(tmpvec(ndvec))
      allocate(covmat(ndvec,ndvec))   
      allocate(psfvec(ndvec,3))
      allocate(thetabin(ntbin))

      ! write(*,*)'in',prm_oc,prm_as9


      ! Now read the data vectors
      open(10,file=fdatavec,status='old',form='unformatted')
      read(10)indvec
      if (indvec.ne.ndvec) then
         write(*,*)'Bad data vector input'
         write(*,*)'STOP'
         stop
      endif         
      read(10)datavec0
      close(10)

      open(10,file=fdatavecflag,status='old',form='unformatted')
      read(10)indvec
      if (indvec.ne.ndvec) then
         write(*,*)'Bad data vector flag input'
         write(*,*)'STOP'
         stop
      endif         
      read(10)datavecflag      
      close(10)
!     stop

!cc   re  ad psf vector
!          write(*,'("### input data vector: ")')
!          write(*,'("# ",a)')trim(fdatavec)
      open(10,file=fpsfvec,status='old',form='unformatted')
      read(10)indvec
      if (indvec.ne.ndvec) then
         write(*,*)'Bad psf vector input'
         write(*,*)'STOP'
         stop
      endif         
      read(10)psfvec
      close(10)
!          stop
!         
!          write(*,'("### input covariance matrix: ")')
!          write(*,'("# ",a)')trim(fcovmat)         
      open(10,file=fcovmat,status='old',form='unformatted')
      read(10)indvec
      if (indvec.ne.ndvec) then
         write(*,*)'Bad covariance input'
         write(*,*)'STOP'
         stop
      endif         
      read(10)covmat
      close(10)
       
!cc   read theta-bin
      open(10,file=ftbinin,status='old')
      do i=1,ntbin
         read(10,*)tbin
         thetabin(i)=tbin
      enddo
      close(10)

      open(11,file=fnspdf,status='unknown')
      read(11,*)nspdf
      read(11,*)delnspdf
      ! This fixes bug in previous code
      delnspdf2=0.5*delnspdf
      do i=1,nspdf
         read(11,*)db0,db1,db2,db3,db4
         znspdf0(i)=db0
         pnspdf0(1,i)=db1
         pnspdf0(2,i)=db2
         pnspdf0(3,i)=db3
         pnspdf0(4,i)=db4
      enddo
      close(11)

      ! This completes all the initializations

      end subroutine HSCWL_init

      subroutine make_modelvec

      ! TH20181101 lin->neu !!!      use getmatterpklincamb
      use getmatterpkneucamb

      implicit none
      
      real e
      
      integer i,j,k,jl,is,ks,ix,jx,l,n
      real ss,cdl,al,zl,fzia,as,xx,x,y,dx,z,ak3
      real zlmin,zlmax,pnsfac,zpsfmin,zpsfmax !!! TH-20180702 P(a)da !!!
      
      integer ia,igrflag,kcv
      
      real*8 pi,pk2del,del2pk,twopi,tps,cv

      real hubp,cosmo_ob,cosmo_e,asf2
      
      real qcoef,zmax,amin,dela,dela2,dsinv,dd,dod,esa
      real yp1,ypn,pow,ak,dps2,spfit,phi,twphi,phiam

!cc parameters for FFTlog
      integer fldir,flkropt
      logical flok
      real*8, allocatable,dimension(:)::fla !fla(nsbin)
      real*8 fldlnr,fldlogr,flkr,fllogkc,fllogrc,flmu,flnc,flq,flrk
      real*8, allocatable,dimension(:)::wsave !wsave(2*nsbin+3*(nsbin/2)+19)      

      real*8, allocatable,dimension(:)::s !s(nsbin)

!!! CAMB
      real*8, dimension(:), allocatable :: pkredshifts!(ia) !!! TH20181101 lin->neu !!!
      real, dimension(:,:), allocatable :: cambpk!(nkmax,0:ia) !!! TH20181101 lin->neu !!!
      real, dimension(:), allocatable :: rklin,pklin,dpklin!(nkdim)
      !real plin
      !real omega_m,omega_v

! halofit
      real aexp,om_m,om_v,grow,grow0,amp,amp2
      real*8 xlogr1,xlogr2,rmid,sig,d1,d2,diff
      real rknl,rneff,rncur      
      real*8 xlogr2start,xlogr2max
      
      real pq,ph,plinz,pnl
      integer nonlinearflag

      real am2rad,rad2am,cln2log,am2deg

      real*8, allocatable,dimension(:,:,:)::dpsnl,pskapnl!(npzbin,npzbin,nsbin)
      real*8, allocatable,dimension(:,:,:)::dpsgi,pskapgi!(npzbin,npzbin,nsbin)
      real*8, allocatable,dimension(:,:,:)::dpsii,pskapii!(npzbin,npzbin,nsbin)
      real, allocatable,dimension(:,:,:)::xipnl,ximnl !(npzbin,npzbin,nsbin)
      real, allocatable,dimension(:,:,:)::xipgi,ximgi !(npzbin,npzbin,nsbin)
      real, allocatable,dimension(:,:,:)::xipii,ximii !(npzbin,npzbin,nsbin)

      real, allocatable,dimension(:,:,:)::xipsum,ximsum !(npzbin,npzbin,ntbin)
      
      real deffs(npzbin),deffs0(npzbin),deffs1(npzbin)
      real weightnl(npzbin,npzbin)
      real weightgi(npzbin,npzbin)
      real weightii(npzbin,npzbin)

      real cornlp,cornlm
      real corgip,corgim
      real coriip,coriim

!!! z-bin
      parameter (ia=100)
      real asf(0:ia),cd(0:ia),d(0:ia)
      real panspdf(npzbin,ia) !!! TH-20180702 P(a)da !!!
      real effs(npzbin,ia)
      real effn(npzbin,ia)

      integer nrkmin,nrkmax,nrk
      real delrk
      parameter (nrkmin=-40,nrkmax=60)
      parameter (nrk=nrkmax-nrkmin+1)
      parameter (delrk=0.1)      

      character*256 fname,fnspdf,ftbinin,flog

      !external dw

      integer itertest
      

      include 'omp_lib.h'
      
      pi=4.0*atan(1.0)
      twopi=2.0*pi
      tps=2.0*pi**2
      am2deg=pi/(180.0*60.0)
      rad2am=180.0*60.0/pi
      am2rad=pi/(180.0*60.0)
      cln2log=log(10.0)
      cv=2997.92458 ! c/H0 = 3e5[km/s] / 100h[km/s/Mpc]  = 3000 /h Mpc
                    ! c/H0 unit:  d <=> w
                    ! Mpc/h unit: r <=> k
                    ! d: distance in c/H0 unit => r = d*cv [Mpc/h unit]
                    ! w: wave number in H0/c unit => k = w/cv [h/Mpc unit]
      pk2del=cv**3/(2.0*pi**2)
      del2pk=2.0*pi**2/(cv**3)

      igrflag=-1

      flq=0.0d0 ! bias exponent: q = 0 is unbiased        
      flkr=1.0d0 ! sensible approximate choice of k_c r_c
!        
      flkropt=0 ! tell fhti to change kr to low-ringing value
                ! 0 to use input kr as is;
                ! 1 to change kr to nearest low-ringing kr, quietly;
                ! 2 to change kr to nearest low-ringing kr, verbosely;
                ! 3 for option to change kr interactively.
      fldir=-1  ! 1/-1 (=forward/inverse) transform
 
      allocate(wsave(2*nsbin+3*(nsbin/2)+19))
      allocate(fla(nsbin))      
      allocate(s(nsbin))
      allocate(dpsnl(npzbin,npzbin,nsbin))
      allocate(pskapnl(npzbin,npzbin,nsbin))
      allocate(dpsgi(npzbin,npzbin,nsbin))
      allocate(pskapgi(npzbin,npzbin,nsbin))
      allocate(dpsii(npzbin,npzbin,nsbin))
      allocate(pskapii(npzbin,npzbin,nsbin))
      allocate(xipnl(npzbin,npzbin,nsbin))
      allocate(ximnl(npzbin,npzbin,nsbin))
      allocate(xipgi(npzbin,npzbin,nsbin))
      allocate(ximgi(npzbin,npzbin,nsbin))
      allocate(xipii(npzbin,npzbin,nsbin))
      allocate(ximii(npzbin,npzbin,nsbin))
      allocate(xipsum(npzbin,npzbin,ntbin))
      allocate(ximsum(npzbin,npzbin,ntbin))

!cc FFTlog parameters
      fllogrc=(slogmin+slogmax)/2.0d0 !central point log10(r_c) of periodic interval
      flnc=dble(nsbin+1)/2.0d0 !central index (1/2 integral if n is even)
      fldlogr=(slogmax-slogmin)/nsbin ! log spacing of points
      fldlnr=fldlogr*log(10.0d0)
      
!
!!!!!!!!! Cosmic shear powewr spectrum
!

!cc source redshift
      zmax=nspdf*delnspdf
      amin=1.0/(1.0+zmax)
      dela=(1.0-amin)/ia
      dela2=0.5*dela
      
!cc  complute a-w a-r relation 
      asf(0)=1.0
!$omp parallel
!$omp do
      do j=1,ia
         asf(j)=1.0-j*dela
      enddo
!$omp end do
!$omp end parallel

!!! TH20181101 lin->neu !!! begin
!
!!!!!! call getcambpk to get linear Pk(at z=0)
!
      allocate(pkredshifts(ia))
      allocate(cambpk(nkmax,0:ia))
!!! TH20181101 lin->neu !!!      allocate(cambpk(nkmax,0:1)) 

      do i=1,ia
         j=ia-i+1
	 pkredshifts(i)=1.0/asf(j)-1.0
      enddo

!!! TH20181101 lin->neu !!!      call getcambpklin(par_sig8,par_ts100,par_h0,par_oc,par_ob,par_ov,par_on,par_as,par_ns,par_tu,par_ok,par_om,nkdim,nkmax,cambpk)
      print*, par_sig8,par_ts100,par_h0,par_oc,par_ob,par_ov,par_on,par_as,par_ns,par_tu,par_ok,par_om,nkdim,nkmax,ia!,pkredshifts,cambpk 
      !stop
      call getcambpkneu(par_sig8,par_ts100,par_h0,par_oc,par_ob,par_ov,par_on,par_as,par_ns,par_tu,par_ok,par_om,nkdim,nkmax,ia,pkredshifts,cambpk)

      allocate(rklin(nkdim))
      allocate(pklin(nkdim))
      allocate(dpklin(nkdim))
!$omp parallel
!$omp do
      do i=1,nkdim
         rklin(i)=cambpk(i,0) !!! TH20181101 lin->neu !!!
      enddo
!$omp end do
!$omp end parallel
!!! TH20181101 lin->neu !!! end

!
!!!!! set cosmological parameters for private functions
!
      cosmo_om=par_om
      cosmo_ov=par_ov
      hubp=0.01*par_h0
      cosmo_ob=par_ob
      
      cosmo_crv=cosmo_om+cosmo_ov-1.0
      cosmo_scrv=sqrt(abs(cosmo_crv))
      kcv=0
      if (cosmo_crv.gt.0.001) kcv=1
      if (cosmo_crv.lt.-0.001) kcv=-1
      
      qcoef=3.0/2.0*cosmo_om

!cc  complute a-r relation 
      cd(0)=0.0
      d(0)=0.0
!$omp parallel private(ss)
!$omp do
      do j=1,ia
         call qromb(dw,asf(j),1.0,ss)
         cd(j)=ss
         d(j)=fk(ss,kcv)
      enddo
!$omp end do
!$omp end parallel

!!! TH-20180702 P(a)da begin
!$omp parallel private(zl,zlmin,zlmax,zpsfmin,zpsfmax,pnsfac)
!$omp do
      do j=1,npzbin ! loop for photo-z bin (1-4)
         do is=1,ia ! loop for lens redshift
            zl=1.0/asf(is)-1.0
            zlmin=1.0/(asf(is)+0.5*dela)-1.0
            zlmax=1.0/(asf(is)-0.5*dela)-1.0
            panspdf(j,is)=0.0  
            do i=1,nspdf ! loop for PDF(z)
               if (znspdf(j,i)+delnspdf2.ge.zlmin.and.znspdf(j,i)-delnspdf2.lt.zlmax) then                  
                  zpsfmin=znspdf(j,i)-delnspdf2
                  if (zpsfmin.lt.zlmin) zpsfmin=zlmin
                  zpsfmax=znspdf(j,i)+delnspdf2
                  if (zpsfmax.ge.zlmax) zpsfmax=zlmax
                  pnsfac=(zpsfmax-zpsfmin)/delnspdf
!
                  panspdf(j,is)=panspdf(j,is)+pnsfac*pnspdf(j,i)*delnspdf/dela
               endif
               if (znspdf(j,i)-delnspdf2.gt.zlmax) exit
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel

!      do j=1,npzbin
!         pdftot(j)=0.0
!         do is=1,ia
!            pdftot(j)=pdftot(j)+panspdf(j,is)*dela
!
!        write(*,505)'# ',j,pdftot(j)
!      enddo
! 505  format(a2,i2,f9.5)

!      do is=1,ia         
!         zl=1.0/asf(is)-1.0
!         write(*,*)zl,panspdf(1,is),
!     $        panspdf(2,is),panspdf(3,is),panspdf(4,is)
!      enddo
!      stop

!!! TH-20180702 P(a)da end

!cc efficiency function
      do jl=1,ia-1 ! loop for zl
         cdl=cd(jl)
         al=1.0-jl*dela
         zl=1.0/al-1.0
         fzia=-0.01387683*par_aia*cosmo_om/gr(al,igrflag)*((1.0+zl)/(1.0+par_z0ia))**par_etaia
!cc source
!$omp parallel private(as,xx,dsinv,dd,dod)
!$omp do
         do j=1,npzbin
            deffs(j)=0.0
            deffs0(j)=0.0
            deffs1(j)=0.0
            do is=1,nspdf       ! loop for zs
               if (znspdf(j,is).le.zl) cycle
               as=1.0/(1.0+znspdf(j,is))
               call qromb(dw,as,1.0,xx)
               dsinv=1.0/fk(xx,kcv)
               dd=xx-cdl
               dod=fk(dd,kcv)*dsinv
!     
               deffs1(j)=dod*pnspdf(j,is)
               deffs(j)=deffs(j)+(deffs0(j)+deffs1(j))*delnspdf2
               deffs0(j)=deffs1(j)
            enddo
            effs(j,jl)=qcoef*deffs(j)/al
            effn(j,jl)=fzia*panspdf(j,jl)/d(jl) !!! TH-20180702 P(a)da !!!
         enddo        
!$omp end do
!$omp end parallel
      enddo
      do j=1,npzbin
         effs(j,ia)=0.0
         effn(j,ia)=0.0
      enddo
      
!cc clear arrays
!$omp parallel
!$omp do
      do i=1,nsbin
         s(i)=10.0d0**(fllogrc+(i-flnc)*fldlogr)
         do k=1,npzbin
            do j=k,npzbin    
               pskapnl(j,k,i)=0.0
               dpsnl(j,k,i)=0.0
               pskapgi(j,k,i)=0.0
               dpsgi(j,k,i)=0.0
               pskapii(j,k,i)=0.0
               dpsii(j,k,i)=0.0
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel

!
! compute power-spectrum
!
      grow0=gg(par_om,par_ov)     
      do i=1,ia

!!! TH20181101 lin->neu !!! begin
!
!!!!! spline modeling of the matter power spectrum
!
      j=ia-i+1
!$omp parallel
!$omp do
      do k=1,nkdim
         pklin(k)=pk2del*cambpk(k,j)*rklin(k)**3
      enddo
!$omp end do
!$omp end parallel

      yp1=(pklin(2)-pklin(1))/(rklin(2)-rklin(1))
      ypn=(pklin(nkdim)-pklin(nkdim-1))/(rklin(nkdim)-rklin(nkdim-1))
      call nr_spline(rklin,pklin,nkdim,yp1,ypn,dpklin)
      pow=(log(pklin(nkdim))-log(pklin(nkdim-1)))/(log(rklin(nkdim))-log(rklin(nkdim-1)))
!!! TH20181101 lin->neu !!! end

         z=1.0/asf(i)-1.0
         asf2=asf(i)*asf(i)
         cosmo_e=sqrt(asf(i)*cosmo_om-asf2*cosmo_crv+asf2*asf2*cosmo_ov)
         esa=1.0/cosmo_e
         do k=1,npzbin
            do j=k,npzbin               
               weightnl(j,k)=effs(j,i)*effs(k,i)*esa
               weightgi(j,k)=(effn(j,i)*effs(k,i)+effn(k,i)*effs(j,i))*esa
               weightii(j,k)=effn(j,i)*effn(k,i)*esa
            enddo
         enddo
! expansion factor from redshift

         aexp=1.0/(1.0+z)

! calculate matter density, vacuum density at desired redshift

         om_m=omega_m(aexp,par_om,par_ov) 
         om_v=omega_v(aexp,par_om,par_ov)

! calculate the amplitude of the power spectrum at desired redshift 
! using linear growth factors (gg given by Caroll, Press & Turner 1992, ARAA, 30, 499)

         grow=gg(om_m,om_v) 
!!! TH20181101 lin->neu !!!         amp=aexp*grow/grow0
         amp=1.0

! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and 
! curvature (rncur) of the power spectrum at the desired redshift, using method 
! described in Smith et al (2002).

!      write(*,*) 'computing effective spectral quantities:'

         xlogr2start=3.5d0
         xlogr2max=xlogr2start-0.0001d0
         xlogr1=-2.0d0
         xlogr2=xlogr2max
         itertest=0
         do 
            itertest=itertest+1            
            rmid=(xlogr2+xlogr1)/2.0d0
            rmid=10.0d0**rmid
            call wint(rmid,sig,d1,d2,rklin,pklin,dpklin,nkdim,pow,par_ns,amp)      
            diff=sig-1.0d0
            if (mod(itertest,20).eq.0) write(*,*)itertest,xlogr1,xlogr2,diff
            if (abs(diff).le.0.001d0) then
               rknl=1./rmid
               rneff=-3.0d0-d1
               rncur=-d2
               nonlinearflag=1
               exit
            elseif (diff.gt.0.001d0) then            
               xlogr1=log10(rmid)
            elseif (diff.lt.-0.001d0) then
               xlogr2=log10(rmid)
            endif
!ccc from CAMB halofit_ppf.f90
            if (xlogr2.lt.-1.9999d0) then !is still linear, exit
               nonlinearflag=0
               exit
            else if (xlogr2.gt.xlogr2max) then ! Totally crazy non-linear
!cc TH-20180125
               xlogr2start=xlogr2start+0.5d0
               xlogr2max=xlogr2start-0.0001d0
               xlogr2=xlogr2max
!               write(*,*)'*** Error in halofit ***'
!               write(*,*)'*** STOP ***************'
!               stop
            end if
         end do

!      write(*,20) 'rknl [h/Mpc] =',rknl,'rneff=',rneff, 'rncur=',rncur
! 20   format(a14,f12.6,2x,a6,f12.6,2x,a6,f12.6)

!!! TH20181101 lin->neu !!!         amp2=amp*amp
         amp2=1.0
         asf2=asf(i)*asf(i)
         cosmo_e=sqrt(asf(i)*cosmo_om-asf2*cosmo_crv+asf2*asf2*cosmo_ov)
!$omp parallel private(ak,ak3,plinz,dps2,pnl,pq,ph,y)
!$omp do
         do ks=1,nsbin
            ak=s(ks)/(d(i)*cv)
            ak3=1.0/(ak**3)
!!! linear
            plinz=amp2*plin(rklin,pklin,dpklin,nkdim,pow,par_ns,ak)
            if (nonlinearflag.eq.1) then
!!! non-linear
               call halofit_class(ak,rneff,rncur,rknl,plinz,pnl,pq,ph,om_m,om_v,fneu,par_om) ! halo fitting formula with neutrino effect
!               call halofit(ak,rneff,rncur,rknl,plinz,pnl,pq,ph,om_m,om_v) ! halo fitting formula
               y=del2pk*pnl*ak3
            else
               y=del2pk*plinz*ak3
            endif
            do k=1,npzbin
               do j=k,npzbin    
                  dps2=weightnl(j,k)*y
                  pskapnl(j,k,ks)=pskapnl(j,k,ks)+(dpsnl(j,k,ks)+dps2)*dela2
                  dpsnl(j,k,ks)=dps2
!
                  dps2=weightgi(j,k)*y
                  pskapgi(j,k,ks)=pskapgi(j,k,ks)+(dpsgi(j,k,ks)+dps2)*dela2
                  dpsgi(j,k,ks)=dps2
!
                  dps2=weightii(j,k)*y
                  pskapii(j,k,ks)=pskapii(j,k,ks)+(dpsii(j,k,ks)+dps2)*dela2
                  dpsii(j,k,ks)=dps2
                 
               enddo
            enddo            
         enddo
!$omp end do
!$omp end parallel
      enddo

! Removed lines which do nothing
!!!!$omp parallel
!!!!$omp do
!!!      do ks=1,nsbin
!!!         do k=1,npzbin
!!!            do j=k,npzbin
!!!               pskapnl(j,k,ks)=pskapnl(j,k,ks)
!!!               pskapgi(j,k,ks)=pskapgi(j,k,ks)
!!!               pskapii(j,k,ks)=pskapii(j,k,ks)
!!!            enddo
!!!         enddo
!!!      enddo
!!!!$omp end do
!!!!$omp end parallel

! 
! compute TPCF
!     

!
!!!!!!  for xi_+
!
!--------initialize FFTLog transform - note fhti resets flkr
      flmu=0.0d0 ! order of Bessel function
      call fhti(nsbin,flmu,flq,fldlnr,flkr,flkropt,wsave,flok)
      if (.not.flok) then
         write(*,*)'failed in fhti'
         write(*,*)'STOP'
         stop
      endif      
      fllogkc=log10(flkr)-fllogrc
!      write(*,*)'central point in k-space at log10(k_c) =',fllogkc
!
      flrk=10.0d0**(fllogrc-fllogkc) !rk = r_c/k_c

!!!! shear-shear
      do k=1,npzbin
         do j=k,npzbin 
!$omp parallel
!$omp do
            do ks=1,nsbin
               fla(ks)=s(ks)*pskapnl(j,k,ks)
            enddo
!$omp end do
!$omp end parallel
            call fht(nsbin,fla,fldir,wsave)
!$omp parallel
!$omp do
            do ks=1,nsbin
               xipnl(j,k,ks)=fla(ks)
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo

!!!! G-I
      do k=1,npzbin
         do j=k,npzbin 
!$omp parallel
!$omp do
            do ks=1,nsbin
               fla(ks)=s(ks)*pskapgi(j,k,ks)
            enddo
!$omp end do
!$omp end parallel
            call fht(nsbin,fla,fldir,wsave)
!$omp parallel
!$omp do
            do ks=1,nsbin
               xipgi(j,k,ks)=fla(ks)
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo

!!!! I-I
      do k=1,npzbin
         do j=k,npzbin 
!$omp parallel
!$omp do
            do ks=1,nsbin
               fla(ks)=s(ks)*pskapii(j,k,ks)
            enddo
!$omp end do
!$omp end parallel
            call fht(nsbin,fla,fldir,wsave)
!$omp parallel
!$omp do
            do ks=1,nsbin
               xipii(j,k,ks)=fla(ks)
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo

!
!!!!!!  for xi_-
!
!--------initialize FFTLog transform - note fhti resets flkr
      flmu=4.0d0 ! order of Bessel function
      call fhti(nsbin,flmu,flq,fldlnr,flkr,flkropt,wsave,flok)
      if (.not.flok) then
         write(*,*)'failed in fhti'
         write(*,*)'STOP'
         stop
      endif      
      fllogkc=log10(flkr)-fllogrc
!      write(*,*)'central point in k-space at log10(k_c) =',fllogkc
!
      flrk=10.0d0**(fllogrc-fllogkc) !rk = r_c/k_c

!!! shear-shear
      do k=1,npzbin
         do j=k,npzbin 
!$omp parallel
!$omp do
            do ks=1,nsbin
               fla(ks)=s(ks)*pskapnl(j,k,ks)
            enddo
!$omp end do
!$omp end parallel
            call fht(nsbin,fla,fldir,wsave)
!$omp parallel
!$omp do
            do ks=1,nsbin
               ximnl(j,k,ks)=fla(ks)
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo

!!! G-I
      do k=1,npzbin
         do j=k,npzbin 
!$omp parallel
!$omp do
            do ks=1,nsbin
               fla(ks)=s(ks)*pskapgi(j,k,ks)
            enddo
!$omp end do
!$omp end parallel
            call fht(nsbin,fla,fldir,wsave)
!$omp parallel
!$omp do
            do ks=1,nsbin
               ximgi(j,k,ks)=fla(ks)
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo

!!! I-I
      do k=1,npzbin
         do j=k,npzbin 
!$omp parallel
!$omp do
            do ks=1,nsbin
               fla(ks)=s(ks)*pskapii(j,k,ks)
            enddo
!$omp end do
!$omp end parallel
            call fht(nsbin,fla,fldir,wsave)
!$omp parallel
!$omp do
            do ks=1,nsbin
               ximii(j,k,ks)=fla(ks)
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo


!$omp parallel private(phi,twphi)
!$omp do
      do ks=1,nsbin
         do k=1,npzbin
            do j=k,npzbin
               phi=10.0d0**(fllogkc+(ks-flnc)*fldlogr)
               twphi=1.0/(twopi*phi)
               xipnl(j,k,ks)=xipnl(j,k,ks)*twphi
               ximnl(j,k,ks)=ximnl(j,k,ks)*twphi
               xipgi(j,k,ks)=xipgi(j,k,ks)*twphi
               ximgi(j,k,ks)=ximgi(j,k,ks)*twphi
               xipii(j,k,ks)=xipii(j,k,ks)*twphi
               ximii(j,k,ks)=ximii(j,k,ks)*twphi
            enddo
         enddo
      enddo
!$omp end do
!$omp end parallel

!cc for xi_+
      do n=ntp1,ntp2
         phiam=thetabin(n)
         phi=phiam*am2rad
         x=(log10(phi)-fllogkc)/fldlogr+flnc
         ix=int(x)
         jx=ix+1
         dx=x-1.0*ix
         do k=1,npzbin
            do j=k,npzbin
               cornlp=xipnl(j,k,ix)+(xipnl(j,k,jx)-xipnl(j,k,ix))*dx
               corgip=xipgi(j,k,ix)+(xipgi(j,k,jx)-xipgi(j,k,ix))*dx
               coriip=xipii(j,k,ix)+(xipii(j,k,jx)-xipii(j,k,ix))*dx
               xipsum(k,j,n)=cornlp+corgip+coriip
            enddo
         enddo
      enddo
      
!cc for xi_-
      do n=ntm1,ntm2
         phiam=thetabin(n)
         phi=phiam*am2rad
         x=(log10(phi)-fllogkc)/fldlogr+flnc
         ix=int(x)
         jx=ix+1
         dx=x-1.0*ix
         do k=1,npzbin
            do j=k,npzbin
               cornlm=ximnl(j,k,ix)+(ximnl(j,k,jx)-ximnl(j,k,ix))*dx
               corgim=ximgi(j,k,ix)+(ximgi(j,k,jx)-ximgi(j,k,ix))*dx
               coriim=ximii(j,k,ix)+(ximii(j,k,jx)-ximii(j,k,ix))*dx
               ximsum(k,j,n)=cornlm+corgim+coriim
            enddo
         enddo
      enddo

      n=0
      do k=1,npzbin
         do j=k,npzbin
            do i=ntp1,ntp2
               n=n+1
!!! TH20180802 correction for constant-e (xi_+ only)
               modelvec(n)=xipsum(k,j,i)*cor_z_sigmae(k)*cor_z_sigmae(j)+par_conste
            enddo
            do i=ntm1,ntm2
               n=n+1
               modelvec(n)=ximsum(k,j,i)*cor_z_sigmae(k)*cor_z_sigmae(j)
            enddo
         enddo
      enddo

      deallocate(wsave)
      deallocate(fla)
      deallocate(s)
      deallocate(dpsnl)
      deallocate(pskapnl)
      deallocate(dpsgi)
      deallocate(pskapgi)
      deallocate(dpsii)
      deallocate(pskapii)
      deallocate(xipnl)
      deallocate(ximnl)
      deallocate(xipgi)
      deallocate(ximgi)
      deallocate(xipii)
      deallocate(ximii)
      deallocate(xipsum)
      deallocate(ximsum)
      deallocate(rklin)
      deallocate(pklin)
      deallocate(dpklin)
      deallocate(pkredshifts)
      deallocate(cambpk)
      end
!cc

!!! Pk_lin(k) by spline fitting
      function plin(rklin,pklin,dpklin,nkdim,pow,par_ns,x)
      integer nkdim
      real rklin(nkdim),pklin(nkdim),dpklin(nkdim)
      real pow,par_ns,x,y,plin     
      if (x.lt.rklin(1)) then
         y=pklin(1)*(x/rklin(1))**(par_ns+3.0)
      elseif (x.gt.rklin(nkdim)) then
         y=pklin(nkdim)*(x/rklin(nkdim))**pow
      else
         call nr_splint(rklin,pklin,dpklin,nkdim,x,y)
      endif
      plin=y
      return
      end

!
!!!!! subroutines & functions from halofit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

! halo model nonlinear fitting formula as described in 
! Appendix C of Smith et al. (2002)

      subroutine halofit(rk,rn,rncur,rknl,plin,pnl,pq,ph,om_m,om_v)

      implicit none

      real gam,a,amod,b,c,xmu,xnu,alpha,beta,f1,f2,f3,f4
      real rk,rn,plin,pnl,pq,ph
      real om_m,om_v,h,p_index,gams,sig8,amp
      real rknl,y,rncur,p_cdm
      real f1a,f2a,f3a,f4a,f1b,f2b,f3b,f4b,frac
      real yy
      real rn2,rn3,rn4

!      gam=0.86485+0.2989*rn+0.1631*rncur
!      a=1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+      
!     &     0.1670756*rn*rn*rn*rn-0.620695*rncur
!      b=10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
!      c=10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
!      xmu=10**(-3.54419+0.19086*rn)
!      xnu=10**(0.95897+1.2857*rn)
!      alpha=1.38848+0.3701*rn-0.1452*rn*rn
!      beta=0.8291+0.9854*rn+0.3400*rn**2
!
      rn2=rn*rn
      rn3=rn*rn2
      rn4=rn*rn3
!     
      gam=0.1971-0.0843*rn+0.8460*rncur
      a=1.5222+2.8553*rn+2.3706*rn2+0.9903*rn3+0.2250*rn4-0.6038*rncur
      a=10**a
      b=10**(-0.5642+0.5864*rn+0.5716*rn2-1.5474*rncur)
      c=10**(0.3698+2.0404*rn+0.8161*rn2+0.5869*rncur)
      xmu=0.0
      xnu=10**(5.2105+3.6902*rn)
      alpha=abs(6.0835+1.3373*rn-0.1959*rn2-5.5274*rncur)
      beta=2.0379-0.7354*rn+0.3757*rn2+1.2490*rn3+0.3980*rn4-0.1682*rncur

      if(abs(1-om_m).gt.0.01) then ! omega evolution 
         f1a=om_m**(-0.0732)
         f2a=om_m**(-0.1423)
         f3a=om_m**(0.0725)
         f1b=om_m**(-0.0307)
         f2b=om_m**(-0.0585)
         f3b=om_m**(0.0743)       
         frac=om_v/(1.-om_m) 
         f1=frac*f1b + (1-frac)*f1a
         f2=frac*f2b + (1-frac)*f2a
         f3=frac*f3b + (1-frac)*f3a
      else         
         f1=1.0
         f2=1.
         f3=1.
      endif

      y=(rk/rknl)

      ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
      ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))
      yy=-y/4.0-y**2/8.0
      if (yy.gt.-88.0) then
         pq=plin*(1+plin)**beta/(1+plin*alpha)*exp(-y/4.0-y**2/8.0)
      else
         pq=0.0
      endif
      pnl=pq+ph

      return
      end       

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!! TH20180825 neutrino mass
!!! Notice this is the version in arXiv:1109.4416v3 which is diferrent from the MNRAS version in some points

      subroutine halofit_bird(rk,rn,rncur,rknl,plin,pnl,pq,ph,om_m,om_v,fneu,par_om)

      implicit none

      real gam,a,amod,b,c,xmu,xnu,alpha,beta,f1,f2,f3,f4
      real rk,rn,plin,pnl,pq,ph
      real om_m,om_v,h,p_index,gams,sig8,amp
      real rknl,y,rncur,p_cdm
      real f1a,f2a,f3a,f4a,f1b,f2b,f3b,f4b,frac
      real yy
      real rn2,rn3,rn4
      real fneu,par_om
      real qneu,pldash

      rn2=rn*rn
      rn3=rn*rn2
      rn4=rn*rn3
      
      gam=0.86485+0.2989*rn+0.1631*rncur
      a=1.4861+1.83693*rn+1.67618*rn2+0.7940*rn3+0.1670756*rn4-0.620695*rncur
      a=10**a
      b=10**(0.9463+0.9466*rn+0.3084*rn2-0.940*rncur)
      c=10**(-0.2807+0.6669*rn+0.3214*rn2-0.0793*rncur)
      xmu=10**(-3.54419+0.19086*rn)
      xnu=10**(0.95897+1.2857*rn)
      alpha=1.38848+0.3701*rn-0.1452*rn2
      beta=0.8291+0.9854*rn+0.3400*rn2

!!! Takahashi et al's modified version
!      gam=0.1971-0.0843*rn+0.8460*rncur
!      a=1.5222+2.8553*rn+2.3706*rn2+0.9903*rn3+0.2250*rn4-0.6038*rncur
!      a=10**a
!      b=10**(-0.5642+0.5864*rn+0.5716*rn2-1.5474*rncur)
!      c=10**(0.3698+2.0404*rn+0.8161*rn2+0.5869*rncur)
!      xmu=0.0
!      xnu=10**(5.2105+3.6902*rn)
!      alpha=abs(6.0835+1.3373*rn-0.1959*rn2-5.5274*rncur)
!      beta=2.0379-0.7354*rn+0.3757*rn2+1.2490*rn3+0.3980*rn4-0.1682*rncur

      if(abs(1-om_m).gt.0.01) then ! omega evolution 
         f1a=om_m**(-0.0732)
         f2a=om_m**(-0.1423)
         f3a=om_m**(0.0725)
         f1b=om_m**(-0.0307)
         f2b=om_m**(-0.0585)
         f3b=om_m**(0.0743)       
         frac=om_v/(1.-om_m) 
         f1=frac*f1b + (1-frac)*f1a
         f2=frac*f2b + (1-frac)*f2a
         f3=frac*f3b + (1-frac)*f3a
      else         
         f1=1.0
         f2=1.
         f3=1.
      endif

      y=(rk/rknl)

!!! Bird et al 2012 but the version in arXiv:1109.4416v3
      qneu=fneu*(2.080-12.4*(par_om-0.3))/(1.0+0.0012*y**3)
      gam=gam+0.316-0.0765*rn-0.835*rncur
      beta=beta+fneu*(-6.49+1.44*rn2)
!      write(*,*)rk,qneu

      ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
      ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))
      ph=ph*(1.0+qneu)

      yy=-y/4.0-y**2/8.0

      if (yy.gt.-88.0) then
!         pq=plin*(1+plin)**beta/(1+plin*alpha)*exp(-y/4.0-y**2/8.0)
         pldash=plin*(1.0+26.3*fneu*rk*rk/(1.0+1.5*rk*rk))
         pq=plin*((1.0+pldash)**beta)/(1.0+alpha*pldash)*exp(yy)
      else
         pq=0.0
      endif
      pnl=pq+ph

      return
      end       

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!! TH20180829 neutrino mass
!!! based on Bird et al 2012 but the version implemented in halofit_ppf.f90 in camb
!!! which is different from the published versions.

      subroutine halofit_bird_camb(rk,rn,rncur,rknl,plin,pnl,pq,ph,om_m,om_v,fneu)

      implicit none

      real gam,a,amod,b,c,xmu,xnu,alpha,beta,f1,f2,f3,f4
      real rk,rn,plin,pnl,pq,ph
      real om_m,om_v,h,p_index,gams,sig8,amp
      real rknl,y,rncur,p_cdm
      real f1a,f2a,f3a,f4a,f1b,f2b,f3b,f4b,frac
      real yy
      real rn2,rn3,rn4
      real fneu,par_om
      real plinaa

      rn2=rn*rn
      rn3=rn*rn2
      rn4=rn*rn3
      
!      gam=0.86485+0.2989*rn+0.1631*rncur  !!! TH20180829
      gam=0.86485+0.2989*rn+0.1631*rncur+0.3159-0.0765*rn-0.8350*rncur
      a=1.4861+1.83693*rn+1.67618*rn2+0.7940*rn3+0.1670756*rn4-0.620695*rncur
      a=10**a
      b=10**(0.9463+0.9466*rn+0.3084*rn2-0.940*rncur)
      c=10**(-0.2807+0.6669*rn+0.3214*rn2-0.0793*rncur)
      xmu=10**(-3.54419+0.19086*rn)
      xnu=10**(0.95897+1.2857*rn)
      alpha=1.38848+0.3701*rn-0.1452*rn2
!      beta=0.8291+0.9854*rn+0.3400*rn2 !!! TH20180829
      beta=0.8291+0.9854*rn+0.3400*rn2+fneu*(-6.4868+1.4373*rn2)

      if(abs(1-om_m).gt.0.01) then ! omega evolution 
         f1a=om_m**(-0.0732)
         f2a=om_m**(-0.1423)
         f3a=om_m**(0.0725)
         f1b=om_m**(-0.0307)
         f2b=om_m**(-0.0585)
         f3b=om_m**(0.0743)       
         frac=om_v/(1.-om_m) 
         f1=frac*f1b + (1-frac)*f1a
         f2=frac*f2b + (1-frac)*f2a
         f3=frac*f3b + (1-frac)*f3a
      else         
         f1=1.0
         f2=1.
         f3=1.
      endif

      y=(rk/rknl)

      ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
      ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+fneu*0.977) !!! TH20180829

      yy=-y/4.0-y**2/8.0
      if (yy.gt.-88.0) then         
         plinaa=plin*(1+fneu*47.48*rk**2/(1+1.5*rk**2)) !!! TH20180829
         pq=plin*(1.0+plinaa)**beta/(1.0+plinaa*alpha)*exp(yy)
      else
         pq=0.0
      endif
      pnl=pq+ph

      return
      end       

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!! TH20180829 neutrino mass
!!! based on Bird et al 2012 but the version implemented in nonlinear.c in CLASS
!!! which is different from the published versions.
!!! In this version, effect of neutrino is added on Takahashi et al 2012's halofit.

      subroutine halofit_class(rk,rn,rncur,rknl,plin,pnl,pq,ph,om_m,om_v,fneu,par_om)

      implicit none

      real gam,a,amod,b,c,xmu,xnu,alpha,beta,f1,f2,f3,f4
      real rk,rn,plin,pnl,pq,ph
      real om_m,om_v,h,p_index,gams,sig8,amp
      real rknl,y,rncur,p_cdm
      real f1a,f2a,f3a,f4a,f1b,f2b,f3b,f4b,frac
      real yy
      real rn2,rn3,rn4
      real fneu,par_om,plinaa

      rn2=rn*rn
      rn3=rn*rn2
      rn4=rn*rn3
      
!      gam=0.86485+0.2989*rn+0.1631*rncur
!      a=1.4861+1.83693*rn+1.67618*rn2+0.7940*rn3+0.1670756*rn4-0.620695*rncur
!      a=10**a
!      b=10**(0.9463+0.9466*rn+0.3084*rn2-0.940*rncur)
!      c=10**(-0.2807+0.6669*rn+0.3214*rn2-0.0793*rncur)
!      xmu=10**(-3.54419+0.19086*rn)
!      xnu=10**(0.95897+1.2857*rn)
!      alpha=1.38848+0.3701*rn-0.1452*rn2
!      beta=0.8291+0.9854*rn+0.3400*rn2

!!! Takahashi et al's modified version + neutrino effect implemented in CLASS
      gam=0.1971-0.0843*rn+0.8460*rncur
      a=1.5222+2.8553*rn+2.3706*rn2+0.9903*rn3+0.2250*rn4-0.6038*rncur
      a=10.0**a
      b=10.0**(-0.5642+0.5864*rn+0.5716*rn2-1.5474*rncur)
      c=10.0**(0.3698+2.0404*rn+0.8161*rn2+0.5869*rncur)
      xmu=0.0
      xnu=10**(5.2105+3.6902*rn)
      alpha=abs(6.0835+1.3373*rn-0.1959*rn2-5.5274*rncur)
      beta=2.0379-0.7354*rn+0.3757*rn2+1.2490*rn3+0.3980*rn4-0.1682*rncur+fneu*(1.081+0.395*rn2) !!! TH20180829 from CLASS

      if(abs(1-om_m).gt.0.01) then ! omega evolution 
         f1a=om_m**(-0.0732)
         f2a=om_m**(-0.1423)
         f3a=om_m**(0.0725)
         f1b=om_m**(-0.0307)
         f2b=om_m**(-0.0585)
         f3b=om_m**(0.0743)       
         frac=om_v/(1.0-om_m) 
         f1=frac*f1b + (1.0-frac)*f1a
         f2=frac*f2b + (1.0-frac)*f2a
         f3=frac*f3b + (1.0-frac)*f3a
      else         
         f1=1.0
         f2=1.0
         f3=1.0
      endif

      y=(rk/rknl)

      ph=a*y**(f1*3)/(1.0+b*y**(f2)+(f3*c*y)**(3-gam))
      ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1.0+fneu*(0.977-18.015*(par_om-0.3))) !!! TH20180829 from CLASS

      yy=-y/4.0-y**2/8.0
      if (yy.gt.-88.0) then         
         plinaa=plin*(1.0+fneu*47.48*rk**2/(1+1.5*rk**2)) !!! TH20180829
         pq=plin*(1.0+plinaa)**beta/(1.0+plinaa*alpha)*exp(yy)
      else
         pq=0.0
      endif
      pnl=pq+ph

      return
      end       

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! The subroutine wint, finds the effective spectral quantities
! rknl, rneff & rncur. This it does by calculating the radius of 
! the Gaussian filter at which the variance is unity = rknl.
! rneff is defined as the first derivative of the variance, calculated 
! at the nonlinear wavenumber and similarly the rncur is the second
! derivative at the nonlinear wavenumber. 
      subroutine wint(r,sig,d1,d2,rklin,pklin,dpklin,nkdim,pow,par_ns,amp)
      implicit none
      integer nkdim
      real rklin(nkdim),pklin(nkdim),dpklin(nkdim)
      real pow,par_ns,sum1,sum2,sum3,t,y,x,w1,w2,w3
      real rk, p_cdm
      real*8 r, sig, d1,d2
      real amp
      integer i,nint

      nint=3000
      sum1=0.
      sum2=0.
      sum3=0.
      do i=1,nint
         t=(float(i)-0.5)/float(nint)
         y=-1.+1./t
         rk=y
         d2=amp*amp*plin(rklin,pklin,dpklin,nkdim,pow,par_ns,rk)
         x=y*r
         w1=exp(-x*x)
         w2=2*x*x*w1
         w3=4*x*x*(1-x*x)*w1
         sum1=sum1+w1*d2/y/t/t
         sum2=sum2+w2*d2/y/t/t
         sum3=sum3+w3*d2/y/t/t
      enddo
      sum1=sum1/float(nint)
      sum2=sum2/float(nint)
      sum3=sum3/float(nint)
      sig=sqrt(sum1)
      d1=-sum2/sum1
      d2=-sum2*sum2/sum1/sum1 - sum3/sum1

      return
      end
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! evolution of omega matter with expansion factor

      function omega_m(aa,om_m0,om_v0)
      implicit none
      real omega_m,omega_t,om_m0,om_v0,aa
      omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*aa*aa+om_m0/aa)
      omega_m=omega_t*om_m0/(om_m0+om_v0*aa*aa*aa)
      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! evolution of omega lambda with expansion factor

      function omega_v(aa,om_m0,om_v0)      
      implicit none
      real aa,omega_v,om_m0,om_v0,omega_t
      omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*aa*aa+om_m0/aa)
      omega_v=omega_t*om_v0/(om_v0+om_m0/aa/aa/aa)
      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! growth factor for linear fluctuations 

      function gg(om_m,om_v)        
      implicit none
      real gg,om_m,om_v
      gg=2.5*om_m/(om_m**(4./7.)-om_v+(1.+om_m/2.)*(1.+om_v/70.))
      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!cc linear growth factor
      function gr(a,igrflag)
      implicit none
      real gr,a,et,ot,vt,grnorm,graw,e,a2
      integer igrflag
      save grnorm
      if (igrflag.eq.-1) then
         a2=a*a
         e=a*cosmo_om-a2*cosmo_crv+a2*a2*cosmo_ov
         et=1.0/e
         ot=cosmo_om*et
         vt=cosmo_ov*et
         grnorm=0.4*(ot**(0.571428571)-vt+(1.0+0.5*ot)*(1.0+vt*0.014285714))/ot
         igrflag=1
      endif
      a2=a*a
      e=a*cosmo_om-a2*cosmo_crv+a2*a2*cosmo_ov
      et=1.0/e
      ot=cosmo_om*a*et
      vt=cosmo_ov*et*a**4
      graw=2.5*ot/(ot**0.571428571-vt+(1.0+0.5*ot)*(1.0+vt*0.014285714))
      gr=graw*grnorm
      return
      end

!cc
      function dw(a)
      implicit none
      real dw,a,e,a2
      a2=a*a
      e=a*cosmo_om-a2*cosmo_crv+a2*a2*cosmo_ov
      dw=1.0/sqrt(e)
      return
      end

!cc      
      function fk(w,kcv)
      implicit none
      real fk,w
      integer kcv
      if (kcv.eq.-1) then
         fk=sinh(cosmo_scrv*w)/cosmo_scrv
      elseif (kcv.eq.0) then
         fk=w
      else 
         fk=sin(cosmo_scrv*w)/cosmo_scrv
      endif
      return
      end

!ccccc integration pakege cccccccccccccccccccccccccccccccc
      SUBROUTINE qromb(func,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
!U    USES polint,trapzd
      INTEGER j
      REAL dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      write(*,*)'too many steps in qromb'
      stop
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 41m.
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 41m.
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*)'failure in polint'
             stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 41m.

!ccccc spline fitting ccccccccccccccccccccccccccccccccc
      SUBROUTINE nr_splint(xa,ya,y2a,n,x,y)           
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         write(*,*)'bad xa input in nr_splint'
         write(*,*)'STOP'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 41m.
      SUBROUTINE nr_spline(x,y,n,yp1,ypn,y2)                 
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 41m.

    end module wlhsc


