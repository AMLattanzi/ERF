#include <cstdio>

#include <AMReX_LevelBld.H>
#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>

#ifdef ERF_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "ERF.H"
#include "Derive.H"
#include "IndexDefines.H"
#include "prob.H"

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall, UserBC
//
static int scalar_bc[] = {INT_DIR,      EXT_DIR,      FOEXTRAP, REFLECT_EVEN,
                          REFLECT_EVEN, REFLECT_EVEN, EXT_DIR};

static int norm_vel_bc[] = {INT_DIR,     EXT_DIR,     FOEXTRAP, REFLECT_ODD,
                            REFLECT_ODD, REFLECT_ODD, EXT_DIR};

static int tang_vel_bc[] = {INT_DIR,      EXT_DIR,     FOEXTRAP, REFLECT_EVEN,
                            REFLECT_EVEN, REFLECT_ODD, EXT_DIR};

static void
set_scalar_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bc.setLo(dir, scalar_bc[lo_bc[dir]]);
    bc.setHi(dir, scalar_bc[hi_bc[dir]]);
  }
}

static void
set_x_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, norm_vel_bc[lo_bc[0]]); bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]););
}

static void
set_y_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, norm_vel_bc[lo_bc[1]]); bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    , bc.setLo(2, tang_vel_bc[lo_bc[2]]); bc.setHi(2, tang_vel_bc[hi_bc[2]]););
}

static void
set_z_vel_bc(amrex::BCRec& bc, const amrex::BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  AMREX_D_TERM(
    bc.setLo(0, tang_vel_bc[lo_bc[0]]); bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    , bc.setLo(1, tang_vel_bc[lo_bc[1]]); bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    , bc.setLo(2, norm_vel_bc[lo_bc[2]]); bc.setHi(2, norm_vel_bc[hi_bc[2]]););
}

void
ERF::variableSetUp()
{
  // ERF::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job

  // initialize the start time for our CPU-time tracker
  startCPUTime = amrex::ParallelDescriptor::second();

  // Output the git commit hashes used to build the executable.

  if (amrex::ParallelDescriptor::IOProcessor()) {
    const char* erf_hash = amrex::buildInfoGetGitHash(1);
    const char* amrex_hash = amrex::buildInfoGetGitHash(2);
    const char* buildgithash = amrex::buildInfoGetBuildGitHash();
    const char* buildgitname = amrex::buildInfoGetBuildGitName();

    if (strlen(erf_hash) > 0) {
      amrex::Print() << "\n"
                     << "ERF git hash: " << erf_hash << "\n";
    }
    if (strlen(amrex_hash) > 0) {
      amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
    }
    if (strlen(buildgithash) > 0) {
      amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
    }

    amrex::Print() << "\n";
  }

  AMREX_ASSERT(desc_lst.size() == 0);

  // Get options, set phys_bc
  read_params();

  init_transport();

  indxmap::init();

#ifdef ERF_USE_MASA
  if (do_mms) {
    init_mms();
  }
#endif

  //
  // Set number of state variables and pointers to components
  //

  int cnt = 0;
  Density = cnt++;
  Xmom = cnt++;
  Ymom = cnt++;
  Zmom = cnt++;
  Eden = cnt++;
  Eint = cnt++;
  Temp = cnt++;

#ifdef NUM_ADV
  NumAdv = NUM_ADV;
#else
  NumAdv = 0;
#endif

  if (NumAdv > 0) {
    FirstAdv = cnt;
    cnt += NumAdv;
  }

  int dm = AMREX_SPACEDIM;

  amrex::Vector<amrex::Real> center(AMREX_SPACEDIM, 0.0);
  amrex::ParmParse ppc("erf");
  ppc.queryarr("center", center, 0, AMREX_SPACEDIM);

  amrex::Interpolater* interp;

  if (state_interp_order == 0) {
    interp = &amrex::pc_interp;
  } else {
    if (lin_limit_state_interp == 1) {
      interp = &amrex::lincc_interp;
    } else {
      interp = &amrex::cell_cons_interp;
    }
  }

  // Note that the default is state_data_extrap = false,
  // store_in_checkpoint = true.  We only need to put these in
  // explicitly if we want to do something different,
  // like not store the state data in a checkpoint directory
  bool state_data_extrap = false;
  bool store_in_checkpoint;

  int ngrow_state = 2;

  store_in_checkpoint = true;
  desc_lst.addDescriptor(
    State_Type, amrex::IndexType::TheCellType(), amrex::StateDescriptor::Point,
    ngrow_state, NVAR, interp, state_data_extrap, store_in_checkpoint);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
  amrex::Vector<std::string> name(NVAR);

  amrex::BCRec bc;
  cnt = 0;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "density";
  cnt++;
  set_x_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "xmom";
  cnt++;
  set_y_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "ymom";
  cnt++;
  set_z_vel_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "zmom";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "rho_E";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "rho_e";
  cnt++;
  set_scalar_bc(bc, phys_bc);
  bcs[cnt] = bc;
  name[cnt] = "Temp";

  for (int i = 0; i < NumAdv; ++i) {
    char buf[64];
    sprintf(buf, "adv_%d", i);
    cnt++;
    set_scalar_bc(bc, phys_bc);
    bcs[cnt] = bc;
    name[cnt] = std::string(buf);
  }

  amrex::StateDescriptor::BndryFunc bndryfunc1(pc_bcfill_hyp);
  bndryfunc1.setRunOnGPU(true);

  desc_lst.setComponent(State_Type, Density, name, bcs, bndryfunc1);

  if (do_mol_load_balance) {
    desc_lst.addDescriptor(
      Work_Estimate_Type, amrex::IndexType::TheCellType(),
      amrex::StateDescriptor::Point, 0, 1, &amrex::pc_interp);
    // Because we use piecewise constant interpolation, we do not use bc and
    // BndryFunc.
    desc_lst.setComponent(
      Work_Estimate_Type, 0, "WorkEstimate", bc,
      amrex::StateDescriptor::BndryFunc(pc_nullfill));
  }

  // 
  // Create face-based StateData for velocities on each face,
  //     and for all conserved quantities as well
  //     (unnecessary right now but convenient)
  // 
  store_in_checkpoint = true;
  amrex::IndexType xface(amrex::IntVect(1,0,0));
  desc_lst.addDescriptor(X_Vel_Type, xface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  amrex::IndexType yface(amrex::IntVect(0,1,0));
  desc_lst.addDescriptor(Y_Vel_Type, yface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  amrex::IndexType zface(amrex::IntVect(0,0,1));
  desc_lst.addDescriptor(Z_Vel_Type, zface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  num_state_type = desc_lst.size();

  //
  // DEFINE DERIVED QUANTITIES
  //
  // Pressure
  //
  derive_lst.add(
    "pressure", amrex::IndexType::TheCellType(), 1, pc_derpres, the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, Density, NVAR);

  //
  // Sound speed (c)
  //
  derive_lst.add(
    "soundspeed", amrex::IndexType::TheCellType(), 1, pc_dersoundspeed,
    the_same_box);
  derive_lst.addComponent("soundspeed", desc_lst, State_Type, Density, NVAR);

  //
  // Mach number(M)
  //
  derive_lst.add(
    "MachNumber", amrex::IndexType::TheCellType(), 1, pc_dermachnumber,
    the_same_box);
  derive_lst.addComponent("MachNumber", desc_lst, State_Type, Density, NVAR);

  //
  // Vorticity
  //
  derive_lst.add(
    "magvort", amrex::IndexType::TheCellType(), 1, pc_dermagvort,
    grow_box_by_one);
  derive_lst.addComponent("magvort", desc_lst, State_Type, Density, NVAR);

  //
  // Div(u)
  //
  derive_lst.add(
    "divu", amrex::IndexType::TheCellType(), 1, pc_derdivu, grow_box_by_one);
  derive_lst.addComponent("divu", desc_lst, State_Type, Density, NVAR);

  //
  // Velocities
  //

  // This calculcates cell-centered x-velocity from x-face-centered values
  derive_lst.add(
    "x_velocity", amrex::IndexType::TheCellType(), 1, pc_dervelx, the_same_box);
  derive_lst.addComponent("x_velocity", desc_lst, State_Type, Density, NVAR);

  // This calculcates cell-centered y-velocity from y-face-centered values
  derive_lst.add(
    "y_velocity", amrex::IndexType::TheCellType(), 1, pc_dervely, the_same_box);
  derive_lst.addComponent("y_velocity", desc_lst, State_Type, Density, NVAR);

  // This calculcates cell-centered z-velocity from z-face-centered values
  derive_lst.add(
    "z_velocity", amrex::IndexType::TheCellType(), 1, pc_dervelz, the_same_box);
  derive_lst.addComponent("z_velocity", desc_lst, State_Type, Density, NVAR);

  derive_lst.add(
    "magvel", amrex::IndexType::TheCellType(), 1, pc_dermagvel, the_same_box);
  derive_lst.addComponent("magvel", desc_lst, State_Type, Density, NVAR);

  //
  // LES coefficients
  //
  if (do_les) {
    derive_lst.add(
      "C_s2", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);
    derive_lst.addComponent("C_s2", desc_lst, State_Type, Density, 1);

    derive_lst.add(
      "C_I", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);
    derive_lst.addComponent("C_I", desc_lst, State_Type, Density, 1);

    derive_lst.add(
      "Pr_T", amrex::IndexType::TheCellType(), 1, pc_dernull, the_same_box);
    derive_lst.addComponent("Pr_T", desc_lst, State_Type, Density, 1);
  }

  // MMS derives
#ifdef ERF_USE_MASA
  if (do_mms) {
    derive_lst.add(
      "rhommserror", amrex::IndexType::TheCellType(), 1, pc_derrhommserror,
      the_same_box);
    derive_lst.addComponent("rhommserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "ummserror", amrex::IndexType::TheCellType(), 1, pc_derummserror,
      the_same_box);
    derive_lst.addComponent("ummserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "vmmserror", amrex::IndexType::TheCellType(), 1, pc_dervmmserror,
      the_same_box);
    derive_lst.addComponent("vmmserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "wmmserror", amrex::IndexType::TheCellType(), 1, pc_derwmmserror,
      the_same_box);
    derive_lst.addComponent("wmmserror", desc_lst, State_Type, Density, NVAR);

    derive_lst.add(
      "pmmserror", amrex::IndexType::TheCellType(), 1, pc_derpmmserror,
      the_same_box);
    derive_lst.addComponent("pmmserror", desc_lst, State_Type, Density, NVAR);
  }
#endif

  // Problem-specific derives
  add_problem_derives<ProblemDerives>(derive_lst, desc_lst);

  // Set list of active sources
  set_active_sources();
}

void
ERF::set_active_sources()
{
  if (do_diffuse && !do_mol) {
    src_list.push_back(diff_src);
  }

  // optional external source
  if (add_ext_src == 1) {
    src_list.push_back(ext_src);
  }

  // optional forcing source
  if (add_forcing_src == 1) {
    src_list.push_back(forcing_src);
  }

#ifdef ERF_USE_MASA
  // optional MMS source
  if (do_mms) {
    src_list.push_back(mms_src);
  }
#endif
}
