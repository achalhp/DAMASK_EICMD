! Copyright 2011-2022 Max-Planck-Institut für Eisenforschung GmbH
! 
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU Affero General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Affero General Public License for more details.
! 
! You should have received a copy of the GNU Affero General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) phenopowerlaw

  type :: tParameters
    real(pReal) :: &                                                                                !< declaring real number variables 
      dot_gamma_0_sl = 1.0_pReal, &                                                                 !< reference shear strain rate for slip
      dot_gamma_0_tw = 1.0_pReal, &                                                                 !< reference shear strain rate for twin
      n_sl           = 1.0_pReal, &                                                                 !< stress exponent for slip
      n_tw           = 1.0_pReal, &                                                                 !< stress exponent for twin
      f_sat_sl_tw    = 1.0_pReal, &                                                                 !< push-up factor for slip saturation due to twinning
      c_1            = 1.0_pReal, &                                                                 !< hardening parameters
      c_2            = 1.0_pReal, &
      c_3            = 1.0_pReal, &
      c_4            = 1.0_pReal, &
      h_0_sl_sl      = 1.0_pReal, &                                                                 !< reference hardening slip - slip
      h_0_tw_sl      = 1.0_pReal, &                                                                 !< reference hardening twin - slip
      h_0_tw_tw      = 1.0_pReal, &                                                                 !< reference hardening twin - twin
      h_0_tw_tw_nuc  = 1.0_pReal, &                                                                 !< Achal twin nucleation
      h_0_tw_tw_grt  = 1.0_pReal, &                                                                 !< Achal twin growth
      a_sl           = 1.0_pReal, &                                                                 !< non-Schmid Coefficient
      chkstep_nucl   = 1.0_pReal, &                                                                 !< Achal Monte Carlo sampling frequency (for twin nucleation)
      chkstep_grow   = 1.0_pReal, &                                                                 !< Achal Monte Carlo sampling frequency (for twin growth)
      chkgrowth_twin = 1.0_pReal, &                                                                 !< Achal flag for twin growth
      prefdecay_slip = 1.0_pReal, &                                                                 !< Achal pre-factor for the flag_slip decay
      twin_inclusion = 1.0_pReal                                                                    !< Achal flag_slip to introduce specific twin vol fraction to particular phase
      
    real(pReal),               allocatable, dimension(:) :: &                                       !< 1D array
      xi_inf_sl, &                                                                                  !< maximum critical shear stress for slip
      h_int, &                                                                                      !< per family hardening activity (optional)
      gamma_char                                                                                    !< characteristic shear for twins
    real(pReal),               allocatable, dimension(:,:) :: &                                     !< 2D array
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      h_sl_tw, &                                                                                    !< slip resistance from twin activity
      h_tw_sl, &                                                                                    !< twin resistance from slip activity
      h_tw_tw, &                                                                                    !< twin resistance from twin activity
      h_tw_tw_grow, &
      h_tw_tw_nucl
    real(pReal),               allocatable, dimension(:,:,:) :: &                                   !< 3D array
      P_sl, &                                                                                       !< Schmid slip
      P_tw, &                                                                                       !< Schmid twin
      P_nS_pos, &                                                                                   !< Non Schmid +ve
      P_nS_neg, &                                                                                   !< Non Schmid -ve
      Corrs_Matr                                                                                    !< Achal Correspondence Matrix
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      sum_N_tw                                                                                      !< total number of active twin systems
    logical :: &                                                                                    !< flag if non-Schmid slip is Active
      nonSchmidActive = .false.
    character(len=pStringLen), allocatable, dimension(:) :: &                                       !< Storing output messages
      output
    character(len=:),          allocatable, dimension(:) :: &                                       !< Labels of active slip and twin systems
      systems_sl, &
      systems_tw
  end type tParameters

  type :: tIndexDotState                                                                            !< Index to access variables in 2D tPhenopowerlawState
    integer, dimension(2) :: &
      xi_sl, &
      xi_tw, &
      gamma_sl, &
      gamma_tw, &
      xi_tw_nucl, &
      xi_tw_grow
  end type tIndexDotState

  type :: tPhenopowerlawState                                                                       !< state variables at material points
    real(pReal), pointer, dimension(:,:) :: &
      xi_sl, &                                                                                      !< critical shear stress for slip
      xi_tw, &                                                                                      !< critical shear stress for twin
      gamma_sl, &                                                                                   !< shear strain due to slip
      gamma_tw, &                                                                                   !< shear strain due to twin
      xi_tw_nucl, &
      xi_tw_grow, &
      f_tw_nucl, &
      f_tw_grow, &
      fmc_tw_nucl, &
      fmc_tw_grwo
    !real(pReal),  pointer, dimension(:) :: &
    !  variant_twin, &                                                                               !< flag used to assign twin variant
    !  frozen                                                                                        !< 0 to 1
      
  end type tPhenopowerlawState

!--------------------------------------------------------------------------------------------------
! containers for parameters, dot state index,  and state
! containers are arrays used to store parameters, state indices and state variables for multiple material points
  type(tParameters),         allocatable, dimension(:) :: param
  type(tIndexDotState),      allocatable, dimension(:) :: indexDotState
  type(tPhenopowerlawState), allocatable, dimension(:) :: state

contains

!< "contains" means the above types are used in below functions.

!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_phenopowerlaw_init() result(myPlasticity)                                   !< function "plastic_phenopowerlaw_init" returns myPlasticity

  logical, dimension(:), allocatable :: myPlasticity                                                !< logical array: phenopowerlaw is active or not
  integer :: &                                                                                      !< integer variables to iterate 
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &                                                       !< active number of slip and twin systems
    N_sl, N_tw
  real(pReal), dimension(:), allocatable :: &
    xi_0_sl, &                                                                                      !< initial critical shear stress for slip
    xi_0_tw, &                                                                                      !< initial critical shear stress for twin
    a                                                                                               !< non-Schmid coefficients
  character(len=pStringLen) :: &                                                                    !< to store error messages
    extmsg = ''
  class(tNode), pointer :: &                                                                        !< pointers variables, tNode objects
    phases, &
    phase, &
    mech, &                                                                                         !< mechanical model info
    pl                                                                                              !< plastic model info


  myPlasticity = plastic_active('phenopowerlaw')                                                    !< function defined in phase_mechanical_plastic
  if(count(myPlasticity) == 0) return                                                               !< if phenopowerlaw is active or not

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:phenopowerlaw init  -+>>>'
  print'(/,a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)                              !< print for which phases phenopowerlaw is active


  phases => config_material%get('phase')                                                            !< get "phase" objects from config_material object
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  

  do ph = 1, phases%length                                                                           !< looping for each phase
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), stt => state(ph), &                                                  !< parameters, state var, dot state indices of current phase
              idx_dot => indexDotState(ph))

    phase => phases%get(ph)                                                                          !< getting objects
    mech => phase%get('mechanical')
    pl => mech%get('plastic')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)                                   !< Read and Store no. of slip systems
    prm%sum_N_sl = sum(abs(N_sl))                                                                    !< total no. of slip systems
    slipActive: if (prm%sum_N_sl > 0) then                                                           !< if +ve, active slip system exists
      prm%systems_sl = lattice_labels_slip(N_sl,phase_lattice(ph))                                   !< 
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then                                                            !< cI: BCC
        a = pl%get_as1dFloat('a_nonSchmid',defaultVal=emptyRealArray)                                !< a: non-Schmid coefficients for projections
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%P_nS_pos = lattice_nonSchmidMatrix(N_sl,a,+1)                                            !< function returns non-schmid projections
        prm%P_nS_neg = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%P_nS_pos = prm%P_sl
        prm%P_nS_neg = prm%P_sl
      end if
      prm%h_sl_sl = lattice_interaction_SlipBySlip(N_sl,pl%get_as1dFloat('h_sl-sl'),phase_lattice(ph))

      xi_0_sl             = pl%get_as1dFloat('xi_0_sl',   requiredSize=size(N_sl))
      prm%xi_inf_sl       = pl%get_as1dFloat('xi_inf_sl', requiredSize=size(N_sl))
      prm%h_int           = pl%get_as1dFloat('h_int',     requiredSize=size(N_sl), &
                                            defaultVal=[(0.0_pReal,i=1,size(N_sl))])

      prm%dot_gamma_0_sl  = pl%get_asFloat('dot_gamma_0_sl')
      prm%n_sl            = pl%get_asFloat('n_sl')
      prm%a_sl            = pl%get_asFloat('a_sl')
      prm%h_0_sl_sl       = pl%get_asFloat('h_0_sl-sl')

      ! expand: family => system
      xi_0_sl             = math_expand(xi_0_sl,      N_sl)
      prm%xi_inf_sl       = math_expand(prm%xi_inf_sl,N_sl)
      prm%h_int           = math_expand(prm%h_int,    N_sl)

      ! sanity checks
      if (    prm%dot_gamma_0_sl  <= 0.0_pReal)      extmsg = trim(extmsg)//' dot_gamma_0_sl'
      if (    prm%a_sl            <= 0.0_pReal)      extmsg = trim(extmsg)//' a_sl'
      if (    prm%n_sl            <= 0.0_pReal)      extmsg = trim(extmsg)//' n_sl'
      if (any(xi_0_sl             <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_0_sl'
      if (any(prm%xi_inf_sl       <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_inf_sl'

    else slipActive
      xi_0_sl = emptyRealArray
      allocate(prm%xi_inf_sl,prm%h_int,source=emptyRealArray)
      allocate(prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    N_tw         = pl%get_as1dInt('N_tw', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(abs(N_tw))
    twinActive: if (prm%sum_N_tw > 0) then
      prm%systems_tw = lattice_labels_twin(N_tw,phase_lattice(ph))
      prm%P_tw     = lattice_SchmidMatrix_twin(N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%h_tw_tw  = lattice_interaction_TwinByTwin(N_tw,pl%get_as1dFloat('h_tw-tw'),phase_lattice(ph))
      prm%gamma_char = lattice_characteristicShear_twin(N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%h_tw_tw_nucl = lattice_interaction_TwinByTwin(N_tw,pl%get_as1dFloat('h_tw-tw'),phase_lattice(ph))
      prm%h_tw_tw_grow = lattice_interaction_TwinByTwin(N_tw,pl%get_as1dFloat('h_tw-tw'),phase_lattice(ph))
      !prm%chkstep_nucl = 
      !prm%chkstep_grow =
      !prm%chkgrowth_twin = 
      !prm%twin_inclusion =
      !prm%prefdecay_slip = 


      xi_0_tw             = pl%get_as1dFloat('xi_0_tw',requiredSize=size(N_tw))

      prm%c_1             = pl%get_asFloat('c_1',defaultVal=0.0_pReal)
      prm%c_2             = pl%get_asFloat('c_2',defaultVal=1.0_pReal)
      prm%c_3             = pl%get_asFloat('c_3',defaultVal=0.0_pReal)
      prm%c_4             = pl%get_asFloat('c_4',defaultVal=0.0_pReal)
      prm%dot_gamma_0_tw  = pl%get_asFloat('dot_gamma_0_tw')
      prm%n_tw            = pl%get_asFloat('n_tw')
      prm%f_sat_sl_tw     = pl%get_asFloat('f_sat_sl-tw')
      prm%h_0_tw_tw       = pl%get_asFloat('h_0_tw-tw')

      ! expand: family => system
      xi_0_tw       = math_expand(xi_0_tw,N_tw)

      ! sanity checks
      if (prm%dot_gamma_0_tw <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_gamma_0_tw'
      if (prm%n_tw           <= 0.0_pReal)  extmsg = trim(extmsg)//' n_tw'

    else twinActive
      xi_0_tw = emptyRealArray
      allocate(prm%gamma_char,source=emptyRealArray)
      allocate(prm%h_tw_tw(0,0))
      allocate(prm%h_tw_tw_nucl(0,0))
      allocate(prm%h_tw_tw_grow(0,0))
    end if twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
    slipAndTwinActive: if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
      prm%h_0_tw_sl  = pl%get_asFloat('h_0_tw-sl')
      prm%h_sl_tw    = lattice_interaction_SlipByTwin(N_sl,N_tw,pl%get_as1dFloat('h_sl-tw'), &
                                                      phase_lattice(ph))
      prm%h_tw_sl    = lattice_interaction_TwinBySlip(N_tw,N_sl,pl%get_as1dFloat('h_tw-sl'), &
                                                      phase_lattice(ph))
    else slipAndTwinActive
      allocate(prm%h_sl_tw(prm%sum_N_sl,prm%sum_N_tw))                                              ! at least one dimension is 0
      allocate(prm%h_tw_sl(prm%sum_N_tw,prm%sum_N_sl))                                              ! at least one dimension is 0
      prm%h_0_tw_sl = 0.0_pReal
    end if slipAndTwinActive

!--------------------------------------------------------------------------------------------------
!  output pararameters

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(pl)
#else
    prm%output = pl%get_as1dString('output',defaultVal=emptyStringArray)
#endif

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_phaseID == ph)
    sizeDotState = size(['xi_sl   ','gamma_sl']) * prm%sum_N_sl &
                 + size(['xi_tw   ','gamma_tw']) * prm%sum_N_tw !&
                 !+ size(['xi_tw_nucl','xi_tw_grow']) * prm%sum_N_tw &             ! Why not size(['xi_tw_nucl','gamma_tw'])?
                 !+ size(['f_tw_nucl','f_tw_grow']) * prm%sum_N_tw &
                 !+ size(['variant_twin','frozen']) * prm%sum_N_tw &
    sizeState = size(['xi_sl   ','gamma_sl']) * prm%sum_N_sl &
                 + size(['xi_tw   ','gamma_tw']) * prm%sum_N_tw !&
                 !+ size(['xi_tw_nucl','xi_tw_grow']) * prm%sum_N_tw &             ! Why not size(['xi_tw_nucl','gamma_tw'])?
                 !+ size(['f_tw_nucl','f_tw_grow']) * prm%sum_N_tw &
                 !+ size(['fmc_tw_nucl','fmc_tw_grow']) * prm%sum_N_tw &
                 !+ size(['variant_twin','frozen']) * prm%sum_N_tw &

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%xi_sl = [startIndex,endIndex]
    stt%xi_sl => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_sl =  spread(xi_0_sl, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%xi_tw = [startIndex,endIndex]
    stt%xi_tw => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_tw =  spread(xi_0_tw, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%gamma_tw = [startIndex,endIndex]
    stt%gamma_tw => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%xi_tw_nucl = [startIndex,endIndex]
    stt%xi_tw_nucl => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_tw_nucl =  spread(xi_0_tw, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%xi_tw_grow = [startIndex,endIndex]
    stt%xi_tw_grow => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_tw_grow =  spread(xi_0_tw, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    !startIndex = endIndex + 1
    !endIndex   = endIndex + prm%sum_N_tw
    !idx_dot%f_twin_nucl = [startIndex,endIndex]
    !stt%f_twin_nucl => plasticState(ph)%state(startIndex:endIndex,:)
    !plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    !startIndex = endIndex + 1
    !endIndex   = endIndex + prm%sum_N_tw
    !idx_dot%f_twin_grow = [startIndex,endIndex]
    !stt%f_twin_grow => plasticState(ph)%state(startIndex:endIndex,:)
    !plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    !startIndex = endIndex + 1
    !endIndex   = endIndex + prm%sum_N_tw
    !idx_dot%fmc_twin_nucl = [startIndex,endIndex]
    !stt%fmc_twin_nucl => plasticState(ph)%state(startIndex:endIndex,:)
    !plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    !startIndex = endIndex + 1
    !endIndex   = endIndex + prm%sum_N_tw
    !idx_dot%fmc_twin_grow = [startIndex,endIndex]
    !stt%fmc_twin_grow => plasticState(ph)%state(startIndex:endIndex,:)
    !plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    
    !startIndex = endIndex + 1
    !endIndex   = endIndex + prm%sum_N_tw
    !idx_dot%frozen = [startIndex,endIndex]
    !stt%frozen => plasticState(ph)%state(startIndex:endIndex,:)
    !plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)    

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!> @details asummes that deformation by dislocation glide affects twinned and untwinned volume
!  equally (Taylor assumption). Twinning happens only in untwinned volume
!--------------------------------------------------------------------------------------------------
module subroutine phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,               intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl_pos,dot_gamma_sl_neg, &
    ddot_gamma_dtau_sl_pos,ddot_gamma_dtau_sl_neg
    
  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw,ddot_gamma_dtau_tw 
  real(pReal), dimension(1) :: fdot_twin_nucl
  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(ph))

  call kinetics_sl(Mp,ph,en,dot_gamma_sl_pos,dot_gamma_sl_neg,ddot_gamma_dtau_sl_pos,ddot_gamma_dtau_sl_neg)
  slipSystems: do i = 1, prm%sum_N_sl
    Lp = Lp + (dot_gamma_sl_pos(i)+dot_gamma_sl_neg(i))*prm%P_sl(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_sl_pos(i) * prm%P_sl(k,l,i) * prm%P_nS_pos(m,n,i) &
                       + ddot_gamma_dtau_sl_neg(i) * prm%P_sl(k,l,i) * prm%P_nS_neg(m,n,i)
  end do slipSystems

  call kinetics_tw(Mp,ph,en,dot_gamma_tw,ddot_gamma_dtau_tw,fdot_twin_nucl)
  twinSystems: do i = 1, prm%sum_N_tw
    Lp = Lp + dot_gamma_tw(i)*prm%P_tw(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_tw(i)*prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
  end do twinSystems

  end associate

end subroutine phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function phenopowerlaw_dotState(Mp,ph,en) result(dotState)                                   !< Rate of change of state variables

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pReal) :: &
    xi_sl_sat_offset,&
    sumF
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl_pos,dot_gamma_sl_neg, &
    right_SlipSlip
  !real(pReal), dimension(param(ph)%sum_N_tw) :: &
    !fdot_twin_nucl, fdot_twin_grow    dot_xi_tw
  real(pReal), dimension(1) :: fdot_twin_nucl


  associate(prm => param(ph), stt => state(ph), &
            dot_xi_sl => dotState(indexDotState(ph)%xi_sl(1):indexDotState(ph)%xi_sl(2)), &
            dot_xi_tw => dotState(indexDotState(ph)%xi_tw(1):indexDotState(ph)%xi_tw(2)), &
            dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)), &
            dot_gamma_tw => dotState(indexDotState(ph)%gamma_tw(1):indexDotState(ph)%gamma_tw(2)), &
            dot_gamma_tw_nucl => dotState(indexDotState(ph)%gamma_tw(1):indexDotState(ph)%gamma_tw(2)), &
            dot_gamma_tw_grow => dotState(indexDotState(ph)%gamma_tw(1):indexDotState(ph)%gamma_tw(2)), &
            dot_xi_tw_nucl => dotState(indexDotState(ph)%xi_tw_nucl(1):indexDotState(ph)%xi_tw_nucl(2)), &
            dot_xi_tw_grow => dotState(indexDotState(ph)%xi_tw_grow(1):indexDotState(ph)%xi_tw_grow(2)))
            
    !sumGamma 
    !sumF_nucl
    !sumF_grow

    call kinetics_sl(Mp,ph,en,dot_gamma_sl_pos,dot_gamma_sl_neg)
    dot_gamma_sl = abs(dot_gamma_sl_pos+dot_gamma_sl_neg)
    call kinetics_tw(Mp,ph,en,dot_gamma_tw,fdot_twin_nucl)

    sumF = sum(stt%gamma_tw(:,en)/prm%gamma_char)
    xi_sl_sat_offset = prm%f_sat_sl_tw*sqrt(sumF)
    right_SlipSlip = sign(abs(1.0_pReal-stt%xi_sl(:,en) / (prm%xi_inf_sl+xi_sl_sat_offset))**prm%a_sl, &
                          1.0_pReal-stt%xi_sl(:,en) / (prm%xi_inf_sl+xi_sl_sat_offset))

    dot_xi_sl = prm%h_0_sl_sl * (1.0_pReal + prm%c_1*sumF** prm%c_2) * (1.0_pReal + prm%h_int) &
                * matmul(prm%h_sl_sl,dot_gamma_sl*right_SlipSlip) &
              + matmul(prm%h_sl_tw,dot_gamma_tw)
!< Rate of change of critical shear stress
    !dot_xi_tw = prm%h_0_tw_sl * sum(stt%gamma_sl(:,en))**prm%c_3 &
    !            * matmul(prm%h_tw_sl,dot_gamma_sl) &
    !          + prm%h_0_tw_tw * sumF**prm%c_4 * matmul(prm%h_tw_tw,dot_gamma_tw)

    dot_xi_tw_nucl = prm%h_0_tw_sl * sum(stt%gamma_sl(:,en))**prm%c_3 &
              * matmul(prm%h_tw_sl,dot_gamma_sl) &
            + prm%h_0_tw_tw_nuc * sumF**prm%c_4 * matmul(prm%h_tw_tw_nucl,dot_gamma_tw)   !sumF_grow
            
    dot_xi_tw_grow = prm%h_0_tw_sl * sum(stt%gamma_sl(:,en))**prm%c_3 &
              * matmul(prm%h_tw_sl,dot_gamma_sl) &
            + prm%h_0_tw_tw_grt * sumF**prm%c_4 * matmul(prm%h_tw_tw_grow,dot_gamma_tw)   !sumF_nucl

    !if(en==4) write(6,*) 'twin volume fraction_1', twin_volume_fraction
    if(en==4) write(6,*) 'dot_xi_tw_nucl_new', dot_xi_tw_grow

  end associate

end function phenopowerlaw_dotState

!--------------------------------------------------------------------------------------------------
!> @brief calculates instantaneous incremental change of kinematics and associated jump state
!--------------------------------------------------------------------------------------------------
! module subroutine plastic_kinematic_deltaFp(twin)
!   use prec, only: &
!    dNeq, &
!    dEq0
! #ifdef DEBUG
!  use debug, only: &
!    debug_level, &
!    debug_constitutive,&
!    debug_levelExtensive, &
!    debug_levelSelective
! #endif

!  use mesh, only: &
!    mesh_element, &
!    mesh_ipNeighborhood, &
!    mesh_ipCoordinates, &
!    mesh_ipVolume, &
!    mesh_ipAreaNormal, &
!    mesh_ipArea, &
!    FE_NipNeighbors, &
!    mesh_maxNipNeighbors, &
!    FE_geomtype, &
!    FE_celltype

!  use lattice
!  use math, only: &
!    math_I3

!  use material, only: &
!    phaseAt, phasememberAt, &
!    phase_plasticityInstance

!  implicit none
!  integer(pInt) :: &
!    ph, of, instance, & 
!    neighbor_el, &                                                  !< element number of neighboring material point
!    neighbor_ip, &                                                  !< integration point of neighboring material point
!    np, &                                                           !< neighbor phase
!    no, n                                                           !< nieghbor offset and index for loop at neighbor

!  integer(pInt),                intent(in)  :: &       
!    el, &                                                           !< element index
!    ip, &                                                           !< integration point index
!    ipc                                                             !< grain index
!  real(pReal), dimension(3,3),  intent(out) :: &
!    deltaFp                                                                                             
!  logical ,                     intent(out) :: &
!    twinJump
! !  real(pReal), dimension(3,3,param(instance)%totalNslip) :: &
! !    CorrespondanceMatrix 
!  integer(pInt), dimension(52) :: &
!    twin_el_incl         
!  real(pReal), dimension(6)     :: &
!    neighbor_stt                  
!  real(pReal)   :: &
!    random, random1
!  integer(pInt) :: &
!    i,j,var_growth,var_nucl
!  var_growth = 0_pInt
!  var_nucl   = 0_pInt
!  ph       = phaseAt(ipc, ip, el)
!  of       = phasememberAt(ipc, ip, el)
!  instance = phase_plasticityInstance(ph)     

!  associate(prm => param(instance), stt => state(instance), dlt => deltaState(instance))

!  twinJump = .false.
!  deltaFp  = math_I3
 
! ! for eshelby circular inclusion
!  twin_el_incl = (/ 10913,10914,10915,10916,10917,10918,10919,10920,10921,10922,10923,10924,10925, &
!                    10993,10994,10995,10996,10997,10998,10999,11000,11001,11002,11074,11075,11076, &
!                    11077,11078,11079,10751,10752,10753,10754,10755,10756,10757,10758,10759,10760, &
!                    10670,10671,10672,10673,10674,10675,10676,10677,10678,10679,10680,10681,10682 /)
! !  TwinLooptest: do i=1_pInt, prm%totalNtwin
! !     write(6,*)'CorrespondenceMatrix for system',i, prm%CorrespondanceMatrix(:,:,i)
! !  enddo TwinLooptest


! !Saving the neighbor information in an array
!  NeighborLoop1: do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))                ! only 4 neighbors for quasi 2D (1 element in z direction)           
!     neighbor_el = mesh_ipNeighborhood(1,n,ip,el)
!     neighbor_ip = mesh_ipNeighborhood(2,n,ip,el)
!     np = phaseAt(1,neighbor_ip,neighbor_el)
!     no = phasememberAt(1,neighbor_ip,neighbor_el)
!     neighbor_stt(n) = state(phase_plasticityInstance(np))%variant_twin(no)
!  enddo NeighborLoop1
 
! !checking if any of my neighbor is twinned if yes recognize the variant and exit
! !  NeighborLoop2: do n = 1_pInt,FE_NipNeighbors(FE_celltype(FE_geomtype(mesh_element(2,el))))                ! only 4 neighbors for quasi 2D (1 element in z direction)           
! !     neighbor_el = mesh_ipNeighborhood(1,n,ip,el)
! !     neighbor_ip = mesh_ipNeighborhood(2,n,ip,el)
! !     np = phaseAt(1,neighbor_ip,neighbor_el)
! ! !     if(of == 1) write(6,*)'phaseAt neighbor_ip of neighbor_el', np      
! !     no = phasememberAt(1,neighbor_ip,neighbor_el)
! ! !     if(of == 1) write(6,*)'phasememberAt at neighbor_ip of neighbor_el', no
! !     if (state(phase_plasticityInstance(np))%variant_twin(no) > 0_pInt) then
! !         var_growth = state(phase_plasticityInstance(np))%variant_twin(no)
! ! 
! !         exit NeighborLoop2
! !     endif       
! !  enddo NeighborLoop2

!  call RANDOM_NUMBER(random)
!  call RANDOM_NUMBER(random1)
! !   Sampling: if (var_growth > 0_pInt) then
! ! !     write(6,*)'I am sampling for growth with variant',var_growth
! !     Ability_Growth: if (stt%f_twin_grow(var_growth,of) > stt%fmc_twin_grow(var_growth,of) &
! !                                                           + prm%checkstep_grow) then
! !       stt%fmc_twin_grow(var_growth,of) = stt%fmc_twin_grow(var_growth,of) &
! !                                                           + prm%checkstep_grow     
! !       Success_Growth: if (random <= stt%f_twin_grow(var_growth,of) .or. ALL(neighbor_stt > 0_pReal)) then
! !           write(6,*)'growth sampling is successful for elem',el        
! !           twinJump                = .true.
! !           deltaFp                 =  prm%CorrespondanceMatrix(:,:,var_growth)  
! !           dlt%f_twin_grow(:,of)   =  0.0_pReal - stt%f_twin_grow(:,of)
! !           dlt%f_twin_nucl(:,of)   =  0.0_pReal - stt%f_twin_nucl(:,of)
! !           dlt%fmc_twin_grow(:,of) =  0.0_pReal - stt%fmc_twin_grow(:,of)
! !           dlt%fmc_twin_nucl(:,of) =  0.0_pReal - stt%fmc_twin_nucl(:,of)          
! !           dlt%frozen(of)          =  1.0_pReal - stt%frozen(of)              
! !           dlt%variant_twin(of)    =  var_growth - stt%variant_twin(of)                        
! !       endif Success_Growth
! !     endif Ability_Growth

! !   elseif (var_growth == 0_pInt .and. prm%checkgrowth_twin > 0_pReal ) then
!   if (var_growth == 0_pInt .and. prm%checkgrowth_twin > 0_pReal ) then
!     var_nucl   = maxloc(stt%f_twin_nucl(:,of), dim=1)
! !     write(6,*)'I am sampling for nucleation with variant',var_nucl,stt%f_twin_nucl(var_nucl,of)
!     Ability_Nucleation: if (stt%f_twin_nucl(var_nucl,of) > stt%fmc_twin_nucl(var_nucl,of) &
!                                                           + prm%checkstep_nucl) then
!       stt%fmc_twin_nucl(var_nucl,of) = stt%fmc_twin_nucl(var_nucl,of) &
!                                                           + prm%checkstep_nucl     
!       Success_Nucleation: if (random <= stt%f_twin_nucl(var_nucl,of) &
!                                                          .and. random1 <= 0.20) then 
!           write(6,*)'nucleation sampling is successful for elem',el                                                               
!           twinJump                = .true.
!           deltaFp                 =  prm%CorrespondanceMatrix(:,:,var_nucl)  
!           dlt%f_twin_nucl(:,of)   =  0.0_pReal - stt%f_twin_nucl(:,of)
!           dlt%f_twin_grow(:,of)   =  0.0_pReal - stt%f_twin_grow(:,of)
!           dlt%fmc_twin_nucl(:,of) =  0.0_pReal - stt%fmc_twin_nucl(:,of)
!           dlt%fmc_twin_grow(:,of) =  0.0_pReal - stt%fmc_twin_grow(:,of)          
!           dlt%frozen(of)          =  1.0_pReal - stt%frozen(of)              
!           dlt%variant_twin(of)    =  var_nucl - stt%variant_twin(of)                        
!       endif Success_Nucleation
!     endif Ability_Nucleation
!   endif
! !   endif Sampling
!  end associate
    
! end subroutine plastic_kinematic_deltaFp

!--------------------------------------------------------------------------------------------------
!> @brief calculates (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
!subroutine plastic_phenopowerlaw_deltaState(instance,of)
!  use prec, only: &
!    dNeq, &
!    dEq0
! #ifdef DEBUG
!  use debug, only: &
!    debug_level, &
!    debug_constitutive,&
!    debug_levelExtensive, &
!    debug_levelSelective
! #endif
 
!  implicit none
!  integer(pInt),                intent(in) :: &
!    instance, &
!    of
! 
!  associate(prm => param(instance), stt => state(instance), dlt => deltaState(instance))
 
! #ifdef DEBUG
!  if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
!             .and. (of == prm%of_debug &
!                    .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
!    write(6,'(a)') '======= phenopowerlaw delta state ======='
! !    write(6,*) sense,state(instance)%sense(:,of)
!  endif
! #endif
 
! !--------------------------------------------------------------------------------------------------
!    dlt%f_twin_nucl(:,of)       = 0.0_pReal
!    dlt%f_twin_grow(:,of)       = 0.0_pReal
!    dlt%fmc_twin_nucl(:,of)     = 0.0_pReal
!    dlt%fmc_twin_grow(:,of)     = 0.0_pReal
!    dlt%frozen(of)              = 0.0_pReal              
!    dlt%variant_twin(of)        = 0.0_pInt              
!  
!  end associate
 
!end subroutine plastic_phenopowerlaw_deltaState 

!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_phenopowerlaw_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('xi_sl')
          call results_writeDataset(stt%xi_sl,group,trim(prm%output(ou)), &
                                    'resistance against plastic slip','Pa',prm%systems_sl)
        case('gamma_sl')
          call results_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                    'plastic shear','1',prm%systems_sl)

        case('xi_tw')
          call results_writeDataset(stt%xi_tw,group,trim(prm%output(ou)), &
                                    'resistance against twinning','Pa',prm%systems_tw)
        case('gamma_tw')
          call results_writeDataset(stt%gamma_tw,group,trim(prm%output(ou)), &
                                    'twinning shear','1',prm%systems_tw)

      end select

    end do

  end associate

end subroutine plastic_phenopowerlaw_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
subroutine kinetics_sl(Mp,ph,en, &
                            dot_gamma_sl_pos,dot_gamma_sl_neg,ddot_gamma_dtau_sl_pos,ddot_gamma_dtau_sl_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal),                  intent(out), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl_pos, &
    dot_gamma_sl_neg
  real(pReal),                  intent(out), optional, dimension(param(ph)%sum_N_sl) :: &
    ddot_gamma_dtau_sl_pos, &
    ddot_gamma_dtau_sl_neg

  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau_sl_pos, &
    tau_sl_neg
  integer :: i

  associate(prm => param(ph), stt => state(ph))

    do i = 1, prm%sum_N_sl
      tau_sl_pos(i) =       math_tensordot(Mp,prm%P_nS_pos(1:3,1:3,i))
      tau_sl_neg(i) = merge(math_tensordot(Mp,prm%P_nS_neg(1:3,1:3,i)), &
                            0.0_pReal, prm%nonSchmidActive)
    end do

    where(dNeq0(tau_sl_pos))
      dot_gamma_sl_pos = prm%dot_gamma_0_sl * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &     ! 1/2 if non-Schmid active
                       * sign(abs(tau_sl_pos/stt%xi_sl(:,en))**prm%n_sl,  tau_sl_pos)
    else where
      dot_gamma_sl_pos = 0.0_pReal
    end where

    where(dNeq0(tau_sl_neg))
      dot_gamma_sl_neg = prm%dot_gamma_0_sl * 0.5_pReal &                                           ! only used if non-Schmid active, always 1/2
                       * sign(abs(tau_sl_neg/stt%xi_sl(:,en))**prm%n_sl,  tau_sl_neg)
    else where
      dot_gamma_sl_neg = 0.0_pReal
    end where

    if (present(ddot_gamma_dtau_sl_pos)) then
      where(dNeq0(dot_gamma_sl_pos))
        ddot_gamma_dtau_sl_pos = dot_gamma_sl_pos*prm%n_sl/tau_sl_pos
      else where
        ddot_gamma_dtau_sl_pos = 0.0_pReal
      end where
    end if
    if (present(ddot_gamma_dtau_sl_neg)) then
      where(dNeq0(dot_gamma_sl_neg))
        ddot_gamma_dtau_sl_neg = dot_gamma_sl_neg*prm%n_sl/tau_sl_neg
      else where
        ddot_gamma_dtau_sl_neg = 0.0_pReal
      end where
    end if

  end associate

end subroutine kinetics_sl


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved
!         stress. Twinning is assumed to take place only in untwinned volume.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
! Mp: Mandel Stress, ph: , en: , dot_gamma_tw: , ddot_gamma_dtau_tw:
!--------------------------------------------------------------------------------------------------
subroutine kinetics_tw(Mp,ph,en,&
                            dot_gamma_tw, fdot_twin_nucl,ddot_gamma_dtau_tw) !fdot_twin_nucl, fdot_twin_grow)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal), dimension(param(ph)%sum_N_tw), intent(out) :: &
    dot_gamma_tw, fdot_twin_nucl!, fdot_twin_grow
  real(pReal), dimension(param(ph)%sum_N_tw), intent(out), optional :: &
    ddot_gamma_dtau_tw

  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    tau_tw
  integer :: i


  associate(prm => param(ph), stt => state(ph))

    tau_tw = [(math_tensordot(Mp,prm%P_tw(1:3,1:3,i)),i=1,prm%sum_N_tw)]

    where(tau_tw > 0.0_pReal) !and stt%frozen(en) < 0.9_pReal)
      dot_gamma_tw = (1.0_pReal-sum(stt%gamma_tw(:,en)/prm%gamma_char)) &                           ! only twin in untwinned volume fraction
                   * prm%dot_gamma_0_tw*(abs(tau_tw)/stt%xi_tw(:,en))**prm%n_tw
      fdot_twin_nucl = max(sum(stt%gamma_tw(:,en)/prm%gamma_char),1.0_pReal)                        ! allocating volume fraction
      !fdot_twin_nucl = 
    else where
      dot_gamma_tw = 0.0_pReal
    end where

    if (present(ddot_gamma_dtau_tw)) then
      where(dNeq0(dot_gamma_tw))
        ddot_gamma_dtau_tw = dot_gamma_tw*prm%n_tw/tau_tw
      else where
        ddot_gamma_dtau_tw = 0.0_pReal
      end where
    end if

    !twin_volume_fraction = 1.0_pReal ! max(sum(stt%gamma_tw(:,en)/prm%gamma_char),1.0_pReal) !f_tw
    if(en==4) write(6,*) "twin volume fraction", sum(stt%gamma_tw(:,en)/prm%gamma_char) !,1.0_pReal)

  end associate

end subroutine kinetics_tw

end submodule phenopowerlaw

