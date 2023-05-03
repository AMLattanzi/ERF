/**
 *  Wrapper for advancing the solution with the slow RHS in the absence of acoustic substepping
 */
    auto no_substep_fun = [&](Vector<MultiFab>& S_sum,
                              Vector<MultiFab>& S_old,
                              Vector<MultiFab>& F_slow,
                              const Real time_for_fp, const Real slow_dt)
    {
        BL_PROFILE("no_substep_fun");
        int n_data = IntVar::NumVars;

        const auto& dxInv = fine_geom.InvCellSizeArray();

        const amrex::GpuArray<int, IntVar::NumVars> scomp_fast = {0,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_fast = {2,1,1,1};

        // Update S_sum = S_stage only for the fast variables
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for ( MFIter mfi(S_sum[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.tilebox();
                Box tbx = surroundingNodes(mfi.tilebox(),0);
                Box tby = surroundingNodes(mfi.tilebox(),1);
                Box tbz = surroundingNodes(mfi.tilebox(),2);

                Vector<Array4<Real> >  ssum_h(n_data);
                Vector<Array4<Real> >  sold_h(n_data);
                Vector<Array4<Real> > fslow_h(n_data);

                for (int i = 0; i < n_data; ++i) {
                     ssum_h[i]  = S_sum[i].array(mfi);
                     sold_h[i]  = S_old[i].array(mfi);
                    fslow_h[i] = F_slow[i].array(mfi);
                }

                Gpu::AsyncVector<Array4<Real> >  sold_d(n_data);
                Gpu::AsyncVector<Array4<Real> >  ssum_d(n_data);
                Gpu::AsyncVector<Array4<Real> > fslow_d(n_data);

                Gpu::copy(Gpu::hostToDevice,  sold_h.begin(),  sold_h.end(),  sold_d.begin());
                Gpu::copy(Gpu::hostToDevice,  ssum_h.begin(),  ssum_h.end(),  ssum_d.begin());
                Gpu::copy(Gpu::hostToDevice, fslow_h.begin(), fslow_h.end(), fslow_d.begin());

                Array4<Real>*  sold =  sold_d.dataPtr();
                Array4<Real>*  ssum =  ssum_d.dataPtr();
                Array4<Real>* fslow = fslow_d.dataPtr();

                // Moving terrain
                if ( solverChoice.use_terrain && solverChoice.terrain_type == 1 )
                {
                    const Array4<const Real>& dJ_old =     detJ_cc[level]->const_array(mfi);
                    const Array4<const Real>& dJ_new = detJ_cc_new[level]->const_array(mfi);

                    const Array4<const Real>& z_nd_new = z_phys_nd_new[level]->const_array(mfi);

                    const Array4<Real      >& z_t_arr  =  z_t_rk[level]->array(mfi);

                    // We have already scaled the slow source term to have the extra factor of dJ
                    ParallelFor(bx, ncomp_fast[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) = dJ_old(i,j,k) *  sold[IntVar::cons](i,j,k,n)
                                                          + slow_dt * fslow[IntVar::cons](i,j,k,n);
                        ssum[IntVar::cons](i,j,k,n) /= dJ_new(i,j,k);
                    });

                    // We have already scaled the slow source term to have the extra factor of dJ
                    ParallelFor(tbx, tby,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        Real dJ_old_xface = 0.5 * (dJ_old(i,j,k) + dJ_old(i-1,j,k));
                        Real dJ_new_xface = 0.5 * (dJ_new(i,j,k) + dJ_new(i-1,j,k));
                        ssum[IntVar::xmom](i,j,k) = dJ_old_xface *  sold[IntVar::xmom](i,j,k)
                                                       + slow_dt * fslow[IntVar::xmom](i,j,k);
                        ssum[IntVar::xmom](i,j,k) /= dJ_new_xface;
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        Real dJ_old_yface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j-1,k));
                        Real dJ_new_yface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j-1,k));
                        ssum[IntVar::ymom](i,j,k) = dJ_old_yface *  sold[IntVar::ymom](i,j,k)
                                                       + slow_dt * fslow[IntVar::ymom](i,j,k);
                        ssum[IntVar::ymom](i,j,k) /= dJ_new_yface;
                    });
                    ParallelFor(tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (k == 0) {
                            Real rho_on_face = 0.5 * ( sold[IntVar::cons](i,j,k,Rho_comp) +sold[IntVar::cons](i,j,k-1,Rho_comp) );
                            ssum[IntVar::zmom](i,j,k) = WFromOmega(i,j,k,rho_on_face*z_t_arr(i,j,k),z_nd_new,dxInv);
                        } else {
                            Real dJ_old_zface = 0.5 * (dJ_old(i,j,k) + dJ_old(i,j,k-1));
                            Real dJ_new_zface = 0.5 * (dJ_new(i,j,k) + dJ_new(i,j,k-1));
                            ssum[IntVar::zmom](i,j,k) = dJ_old_zface *  sold[IntVar::zmom](i,j,k)
                                                           + slow_dt * fslow[IntVar::zmom](i,j,k);
                            ssum[IntVar::zmom](i,j,k) /= dJ_new_zface;
                        }
                    });

                } else { // Fixed or no terrain
                    ParallelFor(bx, ncomp_fast[IntVar::cons],
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
                        const int n = scomp_fast[IntVar::cons] + nn;
                        ssum[IntVar::cons](i,j,k,n) = sold[IntVar::cons](i,j,k,n) + slow_dt *
                           ( fslow[IntVar::cons](i,j,k,n) );
                    });
                    ParallelFor(tbx, tby, tbz,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::xmom](i,j,k) = sold[IntVar::xmom](i,j,k)
                                                  + slow_dt * fslow[IntVar::xmom](i,j,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::ymom](i,j,k) = sold[IntVar::ymom](i,j,k)
                                                  + slow_dt * fslow[IntVar::ymom](i,j,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        ssum[IntVar::zmom](i,j,k) = sold[IntVar::zmom](i,j,k)
                                                  + slow_dt * fslow[IntVar::zmom](i,j,k);
                    });
                }

                // OUTPUT L2 NORM (1 proc, no tiling)
                const Real* prob_lo = fine_geom.ProbLo();
                const Real* dx      = fine_geom.CellSize();
                Real Delta          = std::pow(dx[0]*dx[1]*dx[2],1.0/3.0);

                ParmParse pp("erf");
                Real mu; pp.query("dynamicViscosity",mu);
                Real Cs; pp.query("Cs",Cs);

                amrex::Vector<amrex::Real> abl_pressure_grad_in = {0.0, 0.0, 0.0};
                pp.queryarr("abl_pressure_grad",abl_pressure_grad_in);
                Real Px = abl_pressure_grad_in[0];

                Real Factor = mu / (Cs * Cs * Delta * Delta);

                Vector<Real> C0_arr = {605551.4603806001,  38115.76743991355, 2449.4399999999982, 169.90632344067848};
                Vector<Real> C1_arr = {3250.1208744082833, 811.5733178848386, 201.95839999999998, 49.63299783004268};
                Vector<Real> dz_arr = {1./64., 1./32., 1./16., 1./8.};

                Real C0, C1;
                for (int i(0); i<C0_arr.size(); ++i) {
                    if (std::abs(dx[2] - dz_arr[i]) < 1.0e-6){
                        C0 = C0_arr[i];
                        C1 = C1_arr[i];
                    }
                }

                Real sum(0.);
                ParallelFor(tbx, [=, &sum] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    Real z = prob_lo[2] + (k + 0.5) * dx[2];

                    Real Tmp1 = -Factor * z;
                    Real Tmp2 = (Cs * Cs * Delta * Delta) / (3.0 * Px);
                    Real Tmp3 = 2.0 / (Cs * Cs * Delta * Delta) * Px * z;
                    Real Tmp4 = std::pow(Factor,2.0);
                    Real u_ref = Tmp1 + Tmp2 * std::pow(Tmp3 + Tmp4 + C1,1.5) + C0;
                    Real u_sol = ssum[IntVar::xmom](i,j,k);
                    Real mag   = dx[0] * dx[1] * dx[2] * (u_sol - u_ref) * (u_sol - u_ref);

                    sum +=  mag;
                });
                const auto total_vol = fine_geom.ProbDomain().volume();
                amrex::Print() << "\n";
                amrex::Print() << "L2: " << std::sqrt(sum/total_vol) <<"\n";
                amrex::Print(amrex::ErrorStream()) << std::sqrt(sum/total_vol) << ',' <<"\n";
            }
        }

        // Even if we update all the conserved variables we don't need
        // to fillpatch the slow ones every acoustic substep
        int ng_cons = S_sum[IntVar::cons].nGrow();
        int ng_vel  = S_sum[IntVar::xmom].nGrow();
        bool fast_only          = true;
        bool vel_and_mom_synced = false;
        apply_bcs(S_sum, time_for_fp, ng_cons, ng_vel, fast_only, vel_and_mom_synced);
    };
