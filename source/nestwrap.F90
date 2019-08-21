! Author:Will Handley (wh260@cam.ac.uk) Mar 2015

module nestwrap

    use settings
    use MpiUtils, only: IsMainMPI
    use BaseParameters
    use GeneralSetup

    use priors_module
    use settings_module
    use read_write_module, only: prior_info_file
    use random_module,          only: initialise_random
    use feedback_module
    use nested_sampling_module,   only: NestedSampling



    implicit none

    type(program_settings) :: nest_settings  ! The program settings 
    type(prior),allocatable, dimension(:) :: priors

    !retain CurParams throughout, for the code to be run in the only fast params mode
    Type(ParamSet) :: CurParams

    integer :: size_derived

    ! Prior temporary variables
    character(LEN=:), allocatable :: sorted_uniform_parameter_names ! Store of Sorted parameter names from file 
    real(dp) sorted_uniform_max ! The maximum for the sorted uniform parameters
    real(dp) sorted_uniform_min ! The minimum for the sorted uniform parameters

    integer :: num_repeats_slow
    integer :: num_repeats_semislow
    integer :: num_repeats_semifast
    integer :: num_repeats_fast

    real(dp), dimension(4) :: output_info
    real(dp) :: frac_fast,frac_semi_fast,frac_semi_slow,frac_slow


    integer :: num_slow
    integer :: num_semi_slow
    integer :: num_semi_fast
    integer :: num_fast

    contains

    !-------------------------------------------------------------------------

    subroutine Initialise_PolyChord_Settings(Ini)

        implicit none
        Type(TSettingIni) :: Ini

        nest_settings%boost_posterior = Ini%Read_Double('boost_posterior',5d0)

        ! The number of live points
        nest_settings%nlive                =  Ini%Read_Int('nlive',500)    

        ! The number of slow repeats
        nest_settings%num_repeats = Ini%Read_Int('num_slow_repeats', 8)
        frac_slow      = Ini%Read_Double('nest_frac_slow', 0.75d0)
        frac_semi_slow = Ini%Read_Double('nest_frac_semi_slow', 0.20d0)
        frac_semi_fast = Ini%Read_Double('nest_frac_semi_fast', 0d0)
        frac_fast      = Ini%Read_Double('nest_frac_fast',0.05d0)

        ! The precision criterion for when nested sampling stops
        nest_settings%precision_criterion  =  Ini%Read_Double('nest_precision',1d-3)

        ! The root name for any output files
        nest_settings%file_root            =  Ini%ReadFileName('file_root') 
        nest_settings%base_dir             =  Ini%ReadFileName('root_dir')  

        ! The degree of feedback
        nest_settings%feedback             =  Ini%Read_Int('feedback',1)     

        ! No maximum number of dead points
        nest_settings%max_ndead            =  Ini%Read_Int('max_ndead',-1)      

        ! Whether to resume from a previous run
        nest_settings%read_resume          = Ini%Read_Logical('checkpoint',.false.) 
        ! Write a resume file
        nest_settings%write_resume         = .true.

        nest_settings%write_live           = .true.  ! write out the physical live points?
        nest_settings%write_dead           = .true.  ! write out the dead points?
        nest_settings%write_stats          = .true.  ! write a stats file
        nest_settings%write_paramnames     = .false. ! don't write a paramnames file

        ! Whether to do clustering
        nest_settings%do_clustering =Ini%Read_Logical('do_clustering', .false.)

        nest_settings%equals        = .true.      ! turn on all posterior files
        nest_settings%posteriors    = .true.
        nest_settings%cluster_posteriors = Ini%Read_Logical('do_clustering', .false.) 

        nest_settings%generate_from_seed = Ini%Read_Logical('generate_from_seed', .false.) 
        nest_settings%write_prior = Ini%Read_Logical('write_prior', .true.) 
        nest_settings%nprior = Ini%Read_Int('nprior', -1) 
        nest_settings%nprior_repeat = Ini%Read_Int('nprior_repeat', -1) 

        ! Whether to do a maximisation once finished
        nest_settings%maximise =Ini%Read_Logical('nest_maximise', .false.)

        nest_settings%logzero = -logzero

        ! Read in the information from the ini file about priors

        ! Sorted uniform parameters:
        sorted_uniform_parameter_names = Ini%Read_String_Default('sorted_uniform_parameters','')
        sorted_uniform_max = Ini%Read_Double('sorted_uniform_max',1.d0) 
        sorted_uniform_min = Ini%Read_Double('sorted_uniform_min',0.d0) 


    end subroutine Initialise_PolyChord_Settings


    !-------------------------------------------------------------------------

    subroutine setup_polychord
        use params_module,only:param_type,add_parameter
        use ini_module, only: initialise_program
        use array_module, only: reallocate_1_d
        implicit none


        real(dp) lnew
        integer i,ix
        integer :: speed
        integer :: prior_type

        integer :: slow_index,semi_slow_index,semi_fast_index,fast_index
        integer :: num_speed

        type(param_type),dimension(:),allocatable :: nest_params         ! Parameter array
        type(param_type),dimension(:),allocatable :: nest_derived_params ! Derived parameter array

        real(mcp), allocatable :: derived(:) ! Pointer used to calculate derived parameters

        integer num_sorted_uniform
        integer sorted_uniform_indices(max_num_params) ! vector of parameter indices



        !retain CurParams throughout, for the code to be run in the only fast params mode
        CurParams%P(:num_params) = BaseParams%center(:num_params) 
        lnew = -Setup%LikeCalculator%GetLogLike(CurParams)


        ! ------- Initialise random number generator -------
        ! Initialise the random number generator with the system time
        ! (Provide an argument to this if you want to set a specific seed
        ! leave argumentless if you want to use the system time)
#ifdef MPI
        call initialise_random(MPI_COMM_WORLD)
#else
        call initialise_random(0)
#endif

        
        num_slow      = BaseParams%num_slow - BaseParams%num_semi_slow
        num_semi_slow = BaseParams%num_semi_slow
        num_semi_fast = BaseParams%num_semi_fast
        num_fast      = BaseParams%num_fast - BaseParams%num_semi_fast

        slow_index = 1
        semi_slow_index = slow_index+num_slow
        semi_fast_index = semi_slow_index+num_semi_slow
        fast_index = semi_fast_index+num_semi_fast


        ! Derived parameters
        call Setup%Config%Parameterization%CalcDerivedParams(CurParams,CurParams%Theory, derived) 
        call DataLikelihoods%addLikelihoodDerivedParams(CurParams%P, CurParams%Theory, derived, CurParams%Likelihoods, -lnew)

        if (allocated(derived)) then
            size_derived = size(derived)
        else
            size_derived = 0
        endif



        if (trim(sorted_uniform_parameter_names)/='') then
            ! Get the indices of the named variables contained on the line
            ! 'sorted_uniform_parameters' and store them in sorted_uniform_indices
            ! also store the number of indices in num_sorted_uniform (found by setting
            ! temp_number=-1)
            num_sorted_uniform = -1
            call BaseParams%NameMapping%ReadIndices(trim(sorted_uniform_parameter_names), sorted_uniform_indices, num_sorted_uniform)
        else
            num_sorted_uniform=0
        end if


        if (IsMainMPI()) open(1990,file=trim(prior_info_file(nest_settings)), action='write') 


        ! Allocate a temporary array of parameters
        allocate(nest_params(0),nest_derived_params(size_derived))

        do i=1,num_params_used
            ix = params_used(i)

            ! Set the speed appropriately
            if(i>=slow_index)      speed=1
            if(i>=semi_slow_index) speed=2
            if(i>=semi_fast_index) speed=3
            if(i>=fast_index)      speed=4
            

            if(BaseParams%GaussPriors%std(ix)/=0d0) then
                ! If this is a gaussian prior, then:
                call add_parameter(nest_params,&
                    trim(BaseParams%NameMapping%name(ix)),&
                    BaseParams%NameMapping%label(ix),&
                    speed,&
                    gaussian_type,&
                    gaussian_type,&
                    [ BaseParams%GaussPriors%mean(ix) , BaseParams%GaussPriors%std(ix) ] &
                    )
                ! Now set this to zero
                if (IsMainMPI()) write(1990,'("gaussian ", A, " ", 4E24.15E3)') trim(BaseParams%NameMapping%name(ix)), BaseParams%PMin(ix) , BaseParams%PMax(ix) ,BaseParams%GaussPriors%mean(ix) , BaseParams%GaussPriors%std(ix)
                BaseParams%GaussPriors%std(ix) = 0d0
            else if(any(ix==sorted_uniform_indices(:num_sorted_uniform))) then
                ! If this is a uniform prior, then:
                call add_parameter(nest_params,&
                    trim(BaseParams%NameMapping%name(ix)),&
                    BaseParams%NameMapping%label(ix),&
                    speed,&
                    sorted_uniform_type,&
                    sorted_uniform_type,&
                    [ sorted_uniform_min , sorted_uniform_max ] &
                    )
                if (IsMainMPI()) write(1990,'("sorted ", A, " ", 2E24.15E3)') trim(BaseParams%NameMapping%name(ix)), BaseParams%PMin(ix) , BaseParams%PMax(ix)
            else
                ! If this is a uniform prior, then:
                call add_parameter(nest_params,&
                    trim(BaseParams%NameMapping%name(ix)),&
                    BaseParams%NameMapping%label(ix),&
                    speed,&
                    uniform_type,&
                    uniform_type,&
                    [ BaseParams%PMin(ix) , BaseParams%PMax(ix) ] &
                    )
                if (IsMainMPI()) write(1990,'("uniform ", A, " ", 2E24.15E3)') trim(BaseParams%NameMapping%name(ix)), BaseParams%PMin(ix) , BaseParams%PMax(ix)
            end if

        end do

        if (IsMainMPI()) close(1990) 

        ! Initialise the program
        call initialise_program(nest_settings,priors,nest_params,nest_derived_params)

        if(allocated(nest_settings%grade_frac)) deallocate(nest_settings%grade_frac)
        allocate(nest_settings%grade_frac(0))
        num_speed=0
        if(num_slow/=0) then
            num_speed=num_speed+1
            call reallocate_1_d(nest_settings%grade_frac,num_speed)
            nest_settings%grade_frac(num_speed) = frac_slow
        endif
        if(num_semi_slow/=0) then
            num_speed=num_speed+1
            call reallocate_1_d(nest_settings%grade_frac,num_speed)
            nest_settings%grade_frac(num_speed) = frac_semi_slow
        endif
        if(num_semi_fast/=0) then
            num_speed=num_speed+1
            call reallocate_1_d(nest_settings%grade_frac,num_speed)
            nest_settings%grade_frac(num_speed) = frac_semi_fast
        endif
        if(num_fast/=0) then
            num_speed=num_speed+1
            call reallocate_1_d(nest_settings%grade_frac,num_speed)
            nest_settings%grade_frac(num_speed) = frac_fast
        endif


        allocate(nest_settings%seed_point(nest_settings%nTotal))
        nest_settings%seed_point(nest_settings%p0:nest_settings%p1) = CurParams%P(params_used)
        nest_settings%seed_point(nest_settings%d0:nest_settings%d1) = derived
        nest_settings%seed_point(nest_settings%l0) = lnew

        nest_settings%seed_point(nest_settings%h0:nest_settings%h1) = physical_to_hypercube(nest_settings%seed_point(nest_settings%p0:nest_settings%p1),priors)

    end subroutine setup_polychord

    !-------------------------------------------------------------------------

    function cosmomc_loglikelihood(theta,phi) result(loglikelihood)
        implicit none
        !> Input parameters
        real(dp), intent(in), dimension(:)   :: theta
        !> Output derived parameters
        real(dp), intent(out),  dimension(:) :: phi

        ! The return value
        real(dp) :: loglikelihood

        Type(ParamSet) TrialParams  !Temp storage for Params
        integer :: i

        real(mcp), allocatable :: derived(:) ! Pointer used to calculate derived parameters

        ! Set the local variable TrialParams equal to the last saved value
        TrialParams = CurParams

        ! Pass on the input parameters
        TrialParams%P(params_used) = theta

        ! Calculate the loglikelihood
        loglikelihood = -Setup%LikeCalculator%GetLogLike(TrialParams)

        ! Replace the variable CurParams with TrialPrams
        call CurParams%AcceptReject(TrialParams, loglikelihood>-logzero)

        if (loglikelihood>-logzero) then

            CurParams = TrialParams

            ! Derived parameters
            call Setup%Config%Parameterization%CalcDerivedParams(TrialParams,TrialParams%Theory, derived) 
            call DataLikelihoods%addLikelihoodDerivedParams(TrialParams%P, TrialParams%Theory, derived, TrialParams%Likelihoods, -loglikelihood)

            if (size_derived>0) phi(:size_derived) = derived

        end if

    end function cosmomc_loglikelihood



    subroutine nest_Sample(mpi_communicator)
        implicit none

        integer,intent(in) :: mpi_communicator

        ! Call the nested sampling algorithm on our chosen likelihood and priors
        output_info = NestedSampling(cosmomc_loglikelihood,prior_wrapper,dumper_wrapper,nest_settings,mpi_communicator)

        contains

        function prior_wrapper(cube) result(theta)
            implicit none
            real(dp), intent(in), dimension(:) :: cube
            real(dp), dimension(size(cube))    :: theta
            theta = hypercube_to_physical(cube,priors)
        end function prior_wrapper

        subroutine dumper_wrapper(live, dead, logweights, logZ, logZerr)
            real(dp), intent(in) :: live(:,:), dead(:,:), logweights(:)
            real(dp), intent(in) :: logZ, logZerr
        end subroutine dumper_wrapper

    end subroutine nest_Sample

    !-------------------------------------------------------------------------






end module nestwrap
