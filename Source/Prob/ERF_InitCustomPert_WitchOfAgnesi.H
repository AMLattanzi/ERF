    // Optional Gaussian acoustic pulse in pressure:
    //    p'(x,z) = pert_p_amp * exp(-((x-pert_xc)^2 + (z-pert_zc)^2) / pert_radius^2)
    // The pulse is isentropic, so theta is unchanged and rho' = (rho theta)' / theta
    ParmParse pp("prob");
    Real pert_p_amp  = zero;              pp.query("pert_p_amp" , pert_p_amp);
    Real pert_xc     = zero;              pp.query("pert_xc"    , pert_xc);
    Real pert_zc     = Real(3000.);       pp.query("pert_zc"    , pert_zc);
    Real pert_radius = Real(1000.);       pp.query("pert_radius", pert_radius);

    const bool const_rho = (sc.fixed_density[lev] == 1);
    const bool is_moist  = (sc.moisture_type != MoistureType::None);

    if (pert_p_amp != zero) {
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const auto prob_lo = geomdata.ProbLo();
            const auto dx      = geomdata.CellSize();

            const Real x = ( prob_lo[0] + (i + myhalf) * dx[0] ) / mf_m(i,j,0);
            const Real z = (z_cc) ? z_cc(i,j,k) : prob_lo[2] + (k + myhalf) * dx[2];

            const Real r2 = ( (x-pert_xc)*(x-pert_xc) + (z-pert_zc)*(z-pert_zc) ) / (pert_radius*pert_radius);
            const Real dp = pert_p_amp * std::exp(-r2);

            const Real rho      = state(i,j,k,Rho_comp);
            const Real rhotheta = state(i,j,k,RhoTheta_comp);
            const Real qv       = (is_moist) ? state(i,j,k,RhoQ1_comp)/rho : zero;

            const Real p            = getPgivenRTh(rhotheta, qv);
            const Real rhotheta_new = getRhoThetagivenP(p + dp, qv);

            state_pert(i,j,k,RhoTheta_comp) = rhotheta_new - rhotheta;
            if (!const_rho) {
                state_pert(i,j,k,Rho_comp) = rho * (rhotheta_new / rhotheta - one);
            }
        });
    }
