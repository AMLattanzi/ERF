# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 50000

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 1024 1024 1024

# PROBLEM SIZE & GEOMETRY
#geometry.prob_extent =  12.0     6.0    1.0
#amr.n_cell           =   384     192     32

#geometry.prob_extent =  12.0     1.0    1.0
#amr.n_cell           =   384      32     32

geometry.prob_extent =  1.0     1.0    1.0
amr.n_cell           =  32      32     32

# CONVERGENCE
#======================================
#geometry.prob_extent =  0.5     0.5    1.0
#amr.n_cell           =  32      32     64
#zhi.velocity = 0.5664029656909406 0 0
    
#geometry.prob_extent =  0.5     0.5    1.0
#amr.n_cell           =  16      16     32
#zhi.velocity = 0.5646234542509774 0 0

#geometry.prob_extent =  0.5     0.5    1.0
#amr.n_cell           =  8      8     16
#zhi.velocity = 0.5576738677495996 0 0

#geometry.prob_extent =  0.5     0.5    1.0
#amr.n_cell           =  4      4     8
#zhi.velocity = 0.5322727144609587 0 0
    
geometry.is_periodic = 1 1 0

# DEBUG ONLY
#zlo.velocity = 1 1 0 (fcompare w/ & w/o this passed) 
zlo.type = "NoSlipWall"
    
zhi.type = "SlipWall"
#zhi.type = "NoSlipWall"
 

# TIME STEP CONTROL
erf.no_substepping = 1
erf.incompressible = 1
erf.fixed_dt       = 0.01  # fixed time step depending on grid resolution
erf.horiz_spatial_order = 5
erf.vert_spatial_order  = 5

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
erf.check_file      = chk        # root name of checkpoint file
erf.check_int       = 10000       # number of timesteps between checkpoints

# PLOTFILES
erf.plot_file_1     = plt        # prefix of plotfile name
erf.plot_int_1      = 5000       # number of timesteps between plotfiles
erf.plot_vars_1     = density rhoadv_0 x_velocity y_velocity z_velocity pressure temp theta

# SOLVER CHOICE
erf.alpha_T = 0.0
erf.alpha_C = 1.0

# MATCH AMR-WIND SET UP
#=========================================
# Pressure gradient driven flow
erf.abl_driver_type   = "PressureGradient"
#erf.abl_pressure_grad = -0.003969 0 0    # Negative to force flow left -> right
erf.abl_pressure_grad = -0.001722 0 0    # Negative to force flow left -> right
    
# No gravity to get const rho field
erf.use_gravity = false

# Constant molecular viscosity and Smag model
erf.molec_diff_type  = "Constant"
erf.dynamicViscosity = 8.0e-6 #3.5e-3 # 3.5e-4
erf.les_type         = "Smagorinsky"
erf.Cs               = 0.1

# Specified Re_tau
erf.re_tau = 5200.0 #18.0 #180.0

# NOTE: useful relations
#-----------------------
# u_tau = mu * re_tau / rho
# u_tau = \sqrt( tau_wall / rho )
# (mu * re_tau)^2 = tau_wall = P_x * L_z

# Output for postprocessing
erf.data_log    = scalar_hist.log havg_profiles_mean.log havg_profiles_covar.log havg_profiles_sgs.log
erf.profile_int = 5000 # output interval in number of timesteps
#=========================================
       
# PROBLEM PARAMETERS
prob.rho_0 = 1.0
prob.A_0   = 1.0

prob.U_0 = 0.0
prob.V_0 = 0.0
prob.W_0 = 0.0
prob.T_0 = 300.0

# Higher values of perturbations lead to instability
# Instability seems to be coming from BC
prob.pert_factor    = 0.1 
prob.pert_periods_V = 1.0
prob.pert_periods_W = 4.0
