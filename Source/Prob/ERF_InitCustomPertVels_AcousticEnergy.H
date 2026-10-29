
    amrex::Gpu::streamSynchronize();
    Real U_0 = zero;
    Real V_0 = zero;

    ParmParse pp_for_pert_vels("prob");
    pp_for_pert_vels.query("U_0", U_0);
    pp_for_pert_vels.query("V_0", V_0);

    // Optional Gaussian pulse in x-velocity (same shape as the pressure pulse
    // set in ERF_InitCustomPert_AcousticEnergy.H):
    //    u'(x,z) = pert_u_amp * exp(-((x-pert_xc)^2 + (z-pert_zc)^2) / pert_radius^2)
    Real pert_u_amp  = zero;               pp_for_pert_vels.query("pert_u_amp" , pert_u_amp);
    Real pert_xc     = zero;               pp_for_pert_vels.query("pert_xc"    , pert_xc);
    Real pert_zc     = amrex::Real(3000.); pp_for_pert_vels.query("pert_zc"    , pert_zc);
    Real pert_radius = amrex::Real(1000.); pp_for_pert_vels.query("pert_radius", pert_radius);

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real du = zero;
        if (pert_u_amp != zero) {
            const auto prob_lo = geomdata.ProbLo();
            const auto dx      = geomdata.CellSize();
            const Real x = prob_lo[0] + i * dx[0];
            const Real z = (z_nd) ? fourth * ( z_nd(i,j,k  ) + z_nd(i,j+1,k  )
                                             + z_nd(i,j,k+1) + z_nd(i,j+1,k+1) )
                                  : prob_lo[2] + (k + myhalf) * dx[2];
            const Real r2 = ( (x-pert_xc)*(x-pert_xc) + (z-pert_zc)*(z-pert_zc) ) / (pert_radius*pert_radius);
            du = pert_u_amp * std::exp(-r2);
        }
        x_vel_pert(i, j, k) = U_0 + du;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = V_0;
    });

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = one / dx[0]; dxInv[1] = one / dx[1]; dxInv[2] = one / dx[2];

    // Set the z-velocity from impenetrable condition
    if (sc.terrain_type == TerrainType::StaticFittedMesh ||
        sc.terrain_type == TerrainType::MovingFittedMesh   ) {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = WFromOmega(i, j, k, zero,
                                             x_vel_pert, y_vel_pert,
                                             mf_u, mf_v, z_nd, dxInv);
        });
    }
