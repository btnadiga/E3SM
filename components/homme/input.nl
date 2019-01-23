!
! theta: namelist for dcmip2016 test1: moist baroclinic wave
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test1"         ! test identifier
  ne                = 60                        ! number of elements per cube face
  qsize             = 6                         ! num tracer fields
  nmax             = 18001
  statefreq         = 10                        ! number of steps between screen dumps
  restartfreq       = 1800                        ! don't write restart files if < 0
  runtype           = 0                            ! 0 => new run
  restartfile	    = "./restart/R000296220"
  tstep             = 0.48                        ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5
  rsplit            = 1
  qsplit            = 1
  nu                = 13602352.55150193967211946123                     ! default= 1e15*(ne30/ne60)**3.2 = 1.1e14
  nu_s              = 13602352.55150193967211946123 
  nu_p              = 13602352.55150193967211946123 
  nu_top            = 0                         ! default = 2.5e5
  limiter_option    = 8
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  rearth            = 31880.00000000000000000000
  omega             = .01458400
  moisture          = 'dry'
  theta_hydrostatic_mode = .true.
  dcmip16_prec_type = -1                        ! 0=kessler physics
  dcmip16_pbl_type  = -1                        ! 0=reed-jablonowski pbl, -1 = none
/
&vert_nl
  vform             = "ccm"
  vfile_mid         = "../vcoord/camm-30.ascii"
  vfile_int         = "../vcoord/cami-30.ascii"
/
&analysis_nl
  output_prefix     = ""
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 0,                        ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 38                         ! every 6 hours
  output_varnames1  ='T','ps','pn','pnh','geo','u','v','w','omega','Th','Q','Q2','Q3','Q4','Q5','rho','precl','zeta','zeta_x','zeta_y','gradTh_x','gradTh_y','gradTh_z' ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16			! Lower the better for netcdf. Increase if run out of mem.
  interp_nlon       = 720
  interp_nlat       = 360
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/
