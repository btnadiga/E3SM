#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module element_state

  use kinds,                  only: real_kind, long_kind, int_kind
  use dimensions_mod,         only: np, npsq, nlev, nlevp, qsize_d

  implicit none
  private
  ! integer, public, parameter :: timelevels = 3 ! commented out by SXM
  integer,   public, parameter :: timelevels = 4 ! to go through compilation for all cases when doing BW-PM (ASXM)


  type, public :: elem_state_t

    ! prognostic variables for shallow-water solver
     real (kind=real_kind) :: p(np,np,nlev,timelevels)
     real (kind=real_kind) :: ps(np,np)                               ! surface geopotential
     real (kind=real_kind) :: gradps(np,np,2)                         ! gradient of surface geopotential
     real (kind=real_kind) :: v(np,np,2,nlev,timelevels)              ! contravarient comp
     real (kind=real_kind) :: Fp(np,np)                               ! JRUB F-forcing to p
     real (kind=real_kind) :: p0(np,np)                               ! JRUB p-basic state 
     ! real (kind=real_kind) :: Fv(np,np,2)                           ! JRUB F-forcing to p ! commented out by SXM (to go through compilation for all cases when doing BW-PM, CSXM)
     ! real (kind=real_kind) :: v0(np,np,2)                           ! JRUB p-basic state  ! commented out by SXM (to go through compilation for all cases when doing BW-PM, CSXM)
     real (kind=real_kind) :: Qdp (np,np,nlev,qsize_d,3)              ! Tracer mass (to go through compilation for all cases when doing BW-PM, ASXM)

     ! to go through compilation for all cases when doing BW-PM (ASXM, BEG)
     real (kind=real_kind) :: w_i (np,np,nlevp,timelevels)             ! vertical velocity at interfaces
     ! real (kind=real_kind) :: theta_dp_cp(np,np,nlev,timelevels)     ! potential temperature         ! commented out by SXM (MSXM, GitHub Mar2020)
     real (kind=real_kind) :: vtheta_dp(np,np,nlev,timelevels)         ! virtual potential temperature ! added by SXM (ASXM, GitHub Mar2020)

     real (kind=real_kind) :: phinh_i(np,np,nlevp,timelevels)          ! geopotential used by NH model at interfaces
     real (kind=real_kind) :: Fv(np,np,2,nlev)                         ! F-forcing to p
     real (kind=real_kind) :: v0(np,np,2,nlev)                         ! v-basic state 
     real (kind=real_kind) :: Fw_i(np,np,nlevp)                        ! F-forcing to w_i
     real (kind=real_kind) :: w_i0(np,np,nlevp)                        ! w_i-basic state 
     ! real (kind=real_kind) :: Ftheta_dp_cp(np,np,nlev)               ! F-forcing to theta_dp_cp ! commented out by SXM (MSXM, GitHub Mar2020)
     real (kind=real_kind) :: Fvtheta_dp(np,np,nlev)                   ! F-forcing to vtheta_dp   ! added by SXM (ASXM, GitHub Mar2020)
     ! real (kind=real_kind) :: theta_dp_cp0(np,np,nlev)               ! theta_dp_cp-basic state  ! commented out by SXM (MSXM, GitHub Mar2020)
     real (kind=real_kind) :: vtheta_dp0(np,np,nlev)                   ! vtheta_dp-basic state    ! added by SXM (ASXM, GitHub Mar2020)
     real (kind=real_kind) :: Fphinh_i(np,np,nlevp)                    ! F-forcing to phinh_i
     real (kind=real_kind) :: phinh_i0(np,np,nlevp)                    ! phinh_i-basic state
     real (kind=real_kind) :: Fdp3d(np,np,nlev)                        ! F-forcing to dp3d
     real (kind=real_kind) :: dp3d0(np,np,nlev)                        ! dp3d-basic state 
     real (kind=real_kind) :: Fps_v(np,np)                             ! F-forcing to ps_v 
     real (kind=real_kind) :: ps_v0(np,np)                             ! ps_v-basic state 
     real (kind=real_kind) :: FQdp(np,np,nlev,qsize_d)                 ! F-forcing to Qdp
     real (kind=real_kind) :: Qdp0(np,np,nlev,qsize_d)                 ! Qdp-basic state 
     real (kind=real_kind) :: prvQ   (np,np,nlev,qsize_d)              ! Tracer concentration in previous nstep
     real (kind=real_kind) :: FQPM(np,np,nlev,qsize_d)                 ! F-forcing to Q (not named as FQ to avoid the same name used in later part of this F90)
     real (kind=real_kind) :: Q0PM(np,np,nlev,qsize_d)                 ! Q-basic state  (named to be consistent with FQPM)
     ! to go through compilation for all cases when doing BW-PM (ASXM, END)

  end type elem_state_t

  !___________________________________________________________________
  type, public :: derived_state_t
  end type derived_state_t

  type, public :: elem_accum_t
  end type elem_accum_t


contains
end module 
