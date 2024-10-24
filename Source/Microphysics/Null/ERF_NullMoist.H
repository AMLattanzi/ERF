#ifndef ERF_NULLMOIST_H
#define ERF_NULLMOIST_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>
#include <ERF_DataStruct.H>

class NullMoist {

public:
    NullMoist () {}

    virtual ~NullMoist () = default;

    virtual
    void
    Define (SolverChoice& /*sc*/) { }

    virtual
    void Init (const amrex::MultiFab& /*cons_in*/,
               const amrex::BoxArray& /*grids*/,
               const amrex::Geometry& /*geom*/,
               const amrex::Real& /*dt_advance*/,
               std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
               std::unique_ptr<amrex::MultiFab>& /*detJ_cc*/) { }

    virtual
    void
    Advance (const amrex::Real& /*dt_advance*/,
             const SolverChoice& /*solverChoce*/) { }

    virtual
    void
    Update_Micro_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Update_State_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_State_to_Micro (const amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_Micro_to_State (amrex::MultiFab& /*cons_in*/) { }

    virtual
    amrex::MultiFab*
    Qmoist_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    int
    Qmoist_Size () { return NullMoist::m_qmoist_size; }

    virtual
    int
    Qstate_Size () { return NullMoist::m_qstate_size; }

    virtual
    void
    Qmoist_Restart_Vars ( const SolverChoice&, std::vector<int>& a_idx, std::vector<std::string>& a_names) const
    {
        a_idx.clear();
        a_names.clear();
    }

private:
    int m_qmoist_size = 1;
    int m_qstate_size = 0;
};

#endif
