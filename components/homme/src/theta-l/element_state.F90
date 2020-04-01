#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_state

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d

  implicit none
  private
#ifdef ARKODE
  integer, public, parameter :: timelevels = 50
#else
  ! integer, public, parameter :: timelevels = 3 ! commented out by SXM
  integer,   public, parameter :: timelevels = 4 ! tmLev4 is to store values at nm1 (actually n0 for nstep=0) for PM to avoid undersired situations, e.g., nm1 is used to store values at nm1+dt/5 for tstep_type=5 (ASXM)
#endif
  integer, public, parameter :: diagtimes = 6

  ! maximum number of Newton iterations taken for an IMEX-RK stage per time-step
  integer, public               :: max_itercnt=0
  real (kind=real_kind), public :: max_deltaerr=0
  real (kind=real_kind), public :: max_reserr=0

  ! pressure based TOM sponge layer
  real (kind=real_kind),public :: nu_scale_top(nlev)
  integer, public              :: nlev_tom


! =========== PRIMITIVE-EQUATION DATA-STRUCTURES =====================

  type, public :: elem_state_t

    ! prognostic variables for preqx solver

    ! prognostics must match those in prim_restart_mod.F90
    ! vertically-lagrangian code advects dp3d instead of ps_v
    ! tracers Q, Qdp always use 2 level time scheme

    real (kind=real_kind) :: v   (np,np,2,nlev,timelevels)        ! horizontal velocity 
    real (kind=real_kind) :: w_i (np,np,nlevp,timelevels)         ! vertical velocity at interfaces
    real (kind=real_kind) :: vtheta_dp(np,np,nlev,timelevels)     ! virtual potential temperature (mass)
    real (kind=real_kind) :: phinh_i(np,np,nlevp,timelevels)      ! geopotential used by NH model at interfaces
    real (kind=real_kind) :: dp3d(np,np,nlev,timelevels)          ! delta p on levels                  
    real (kind=real_kind) :: ps_v(np,np,timelevels)               ! surface pressure                   
    real (kind=real_kind) :: phis(np,np)                          ! surface geopotential (prescribed)  
    real (kind=real_kind) :: Q   (np,np,nlev,qsize_d)             ! Tracer concentration               

    ! real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,2)         ! Tracer mass (commented out by SXM)
    real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,3)           ! Tracer mass (tmLev3 is to store values at n0_qpd in previous nstep for PM, ASXM)

    ! ASXM (BEG)    
    real (kind=real_kind) :: Fv(np,np,2,nlev)                     ! F-forcing to p
    real (kind=real_kind) :: v0(np,np,2,nlev)                     ! v-basic state 
    real (kind=real_kind) :: Fw_i(np,np,nlevp)                    ! F-forcing to w_i
    real (kind=real_kind) :: w_i0(np,np,nlevp)                    ! w_i-basic state 
    ! real (kind=real_kind) :: Ftheta_dp_cp(np,np,nlev)           ! F-forcing to theta_dp_cp ! commented out by SXM (MSXM, GitHub Mar2020)
    real (kind=real_kind) :: Fvtheta_dp(np,np,nlev)               ! F-forcing to vtheta_dp   ! added by SXM (ASXM, GitHub Mar2020)
    ! real (kind=real_kind) :: theta_dp_cp0(np,np,nlev)           ! theta_dp_cp-basic state  ! commented out by SXM (MSXM, GitHub Mar2020)
    real (kind=real_kind) :: vtheta_dp0(np,np,nlev)               ! vtheta_dp-basic state    ! added by SXM (ASXM, GitHub Mar2020)
    real (kind=real_kind) :: Fphinh_i(np,np,nlevp)                ! F-forcing to phinh_i
    real (kind=real_kind) :: phinh_i0(np,np,nlevp)                ! phinh_i-basic state
    real (kind=real_kind) :: Fdp3d(np,np,nlev)                    ! F-forcing to dp3d
    real (kind=real_kind) :: dp3d0(np,np,nlev)                    ! dp3d-basic state 
    real (kind=real_kind) :: Fps_v(np,np)                         ! F-forcing to ps_v
    real (kind=real_kind) :: ps_v0(np,np)                         ! ps_v-basic state
    real (kind=real_kind) :: FQdp(np,np,nlev,qsize_d)             ! F-forcing to Qdp
    real (kind=real_kind) :: Qdp0(np,np,nlev,qsize_d)             ! Qdp-basic state 
    real (kind=real_kind) :: prvQ   (np,np,nlev,qsize_d)          ! Tracer concentration in previous nstep
    real (kind=real_kind) :: FQPM(np,np,nlev,qsize_d)             ! F-forcing to Q (not named as FQ to avoid the same name used in later part of this F90)
    real (kind=real_kind) :: Q0PM(np,np,nlev,qsize_d)             ! Q-basic state  (named to be consistent with FQPM)
    ! ASXM (END)
  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t

    ! diagnostic variables for preqx solver

    ! storage for subcycling tracers/dynamics

    real (kind=real_kind) :: vn0  (np,np,2,nlev)                      ! velocity for SE tracer advection
    real (kind=real_kind) :: vstar(np,np,2,nlev)                      ! velocity on Lagrangian surfaces
    real (kind=real_kind) :: dpdiss_biharmonic(np,np,nlev)            ! mean dp dissipation tendency, if nu_p>0
    real (kind=real_kind) :: dpdiss_ave(np,np,nlev)                   ! mean dp used to compute psdiss_tens

    ! diagnostics 
    real (kind=real_kind) :: omega_p(np,np,nlev)                      ! vertical tendency (derived)
    real (kind=real_kind) :: eta_dot_dpdn(np,np,nlevp)                ! mean vertical flux from dynamics
    real (kind=real_kind) :: eta_dot_dpdn_prescribed(np,np,nlevp)     ! prescribed wind test cases

    ! tracer advection fields used for consistency and limiters
    real (kind=real_kind) :: dp(np,np,nlev)                           ! for dp_tracers at physics timestep
    real (kind=real_kind) :: divdp(np,np,nlev)                        ! divergence of dp
    real (kind=real_kind) :: divdp_proj(np,np,nlev)                   ! DSSed divdp

    real (kind=real_kind) :: FQ(np,np,nlev,qsize_d)                ! tracer forcing
    real (kind=real_kind) :: FM(np,np,3,nlev)                      ! momentum forcing
    real (kind=real_kind) :: FT(np,np,nlev)                        ! temperature forcing
    real (kind=real_kind) :: FVTheta(np,np,nlev)                   ! potential temperature forcing
    real (kind=real_kind) :: FPHI(np,np,nlevp)                     ! PHI (NH) forcing
    real (kind=real_kind) :: FQps(np,np)                   ! forcing of FQ on ps_v

    real (kind=real_kind) :: gradphis(np,np,2)   ! grad phi at the surface, computed once in model initialization
    real (kind=real_kind) :: dp_ref(np,np,nlev)    ! ref states based on PHIS
    real (kind=real_kind) :: theta_ref(np,np,nlev)
    real (kind=real_kind) :: phi_ref(np,np,nlevp)  
  end type derived_state_t
  

  !___________________________________________________________________
  type, public :: elem_accum_t

#ifdef ENERGY_DIAGNOSTICS
    ! Energy equation:
    real (kind=real_kind) :: KEu_horiz1(np,np)
    real (kind=real_kind) :: KEu_horiz2(np,np)
    real (kind=real_kind) :: KEu_vert1(np,np)
    real (kind=real_kind) :: KEu_vert2(np,np)
    real (kind=real_kind) :: KEw_horiz1(np,np)  ! nonhydro only
    real (kind=real_kind) :: KEw_horiz2(np,np)  ! nonhydro only
    real (kind=real_kind) :: KEw_horiz3(np,np)  ! nonhydro only
    real (kind=real_kind) :: KEw_vert1(np,np)   ! nonhydro only
    real (kind=real_kind) :: KEw_vert2(np,np)   ! nonhydro only

    real (kind=real_kind) :: IEvert1(np,np)
    real (kind=real_kind) :: IEvert2(np,np)     ! nonhydro only
    real (kind=real_kind) :: PEvert1(np,np)
    real (kind=real_kind) :: PEvert2(np,np)
    real (kind=real_kind) :: PEhoriz1(np,np)
    real (kind=real_kind) :: PEhoriz2(np,np)

    real (kind=real_kind) :: T01(np,np)
    real (kind=real_kind) :: T2(np,np)
    real (kind=real_kind) :: S1(np,np)
    real (kind=real_kind) :: S2(np,np)
    real (kind=real_kind) :: P1(np,np)
    real (kind=real_kind) :: P2(np,np)
    real (kind=real_kind) :: T2_nlevp_term(np,np)

    real (kind=real_kind) :: CONV(np,np,2,nlev)                       ! dpdn u dot CONV = T1 + T2
#endif

    ! the "4" timelevels represents data computed at:
    !  1  t-.5
    !  2  t+.5   after dynamics
    !  3  t+.5   after forcing
    !  4  t+.5   after Robert
    ! after calling TimeLevelUpdate, all times above decrease by 1.0

    real (kind=real_kind) :: KEner(np,np,diagtimes)
    real (kind=real_kind) :: PEner(np,np,diagtimes)
    real (kind=real_kind) :: IEner(np,np,diagtimes)
    real (kind=real_kind) :: Qvar(np,np,qsize_d,diagtimes)                    ! Q variance at half time levels
    real (kind=real_kind) :: Qmass(np,np,qsize_d,diagtimes)                   ! Q mass at half time levels
    real (kind=real_kind) :: Q1mass(np,np,qsize_d)                    ! Q mass at full time levels

  end type elem_accum_t



contains

end module 
