// Author: Daniel F Potter
// Date: 07-Apr-2010
// Place: Dutton Park, QLD, Oz

#ifndef COUPLING_COMPONENT_HH
#define COUPLING_COMPONENT_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/species-energy-modes.hh"

#include "reaction.hh"

class Coupling_component {
public:
    /// \brief Lua file constructor
    Coupling_component( lua_State * L, Reaction * r, std::string type, std::string mode, int idc );
    
    /// \brief Clone constructor
    Coupling_component( const Coupling_component &c );
    
    /// \brief Default destructor
    virtual ~Coupling_component();
    
    /// \brief Clone function
    virtual Coupling_component * clone() const;
    
    /// \brief Set the energy and number density from the beginning of the timestep
    void set_e_and_N_old( double e_old, double N_old )
    { e_old_ = e_old; N_old_ = N_old; }
    
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double compute_contribution( Gas_data &Q, std::vector<double> &delta_c )
    { return specific_compute_contribution( Q, delta_c ); }
    
    /// \brief Compute the the source term from this component in J
    double compute_source_term( Gas_data &Q, std::vector<double> &dcdt )
    { return specific_compute_source_term( Q, dcdt ); }
    
    int get_isp() { return isp_; }
    
    std::string get_mode() { return mode_; }
    
    std::string get_type() { return type_; }
    
protected:
    std::string type_;
    std::string mode_;
    int idc_;
    int isp_;
    double m_;
    int nu_;
    std::vector<Species_energy_mode*> sems_;
    int imode_;
    
    double e_old_;
    double N_old_;
    
protected:
    virtual double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c ) = 0;
    virtual double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt ) = 0;
};

void create_Coupling_components_for_reaction( lua_State * L, Reaction * r, int ir, std::vector<Coupling_component*> & ccs );

/*********************** Dissociation components **********************/

class Simple_dissociation_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Simple_dissociation_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Simple_dissociation_component( const Simple_dissociation_component &c );
    
     /// \brief Default destructor
    ~Simple_dissociation_component();
    
    /// \brief Clone function
    Simple_dissociation_component * clone() const;

private:
    double D_hat_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

class TreanorMarrone_dissociation_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    TreanorMarrone_dissociation_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    TreanorMarrone_dissociation_component( const TreanorMarrone_dissociation_component &c );
    
     /// \brief Default destructor
    ~TreanorMarrone_dissociation_component();
    
    /// \brief Clone function
    TreanorMarrone_dissociation_component * clone() const;

private:
    double U_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

class Park_dissociation_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Park_dissociation_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Park_dissociation_component( const Park_dissociation_component &c );
    
     /// \brief Default destructor
    ~Park_dissociation_component();
    
    /// \brief Clone function
    Park_dissociation_component * clone() const;

private:
    double n_;
    double D_;
    double s_v_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

class Macheret_dissociation_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Macheret_dissociation_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Macheret_dissociation_component( const Macheret_dissociation_component &c );
    
     /// \brief Default destructor
    ~Macheret_dissociation_component();
    
    /// \brief Clone function
    Macheret_dissociation_component * clone() const;

private:
    double A_;
    double n_;
    double theta_d_;
    
    double theta_v_;
    double alpha_;
    
    double theta_dstar_,delta_d_;
    double b_;
    double fac_;
    
    bool monatomic_collider_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/*********************** Recombination components **********************/

class Simple_recombination_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Simple_recombination_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Simple_recombination_component( const Simple_recombination_component &c );
    
     /// \brief Default destructor
    ~Simple_recombination_component();
    
    /// \brief Clone function
    Simple_recombination_component * clone() const;

private:
    double D_hat_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

class TreanorMarrone_recombination_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    TreanorMarrone_recombination_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    TreanorMarrone_recombination_component( const TreanorMarrone_recombination_component &c );
    
     /// \brief Default destructor
    ~TreanorMarrone_recombination_component();
    
    /// \brief Clone function
    TreanorMarrone_recombination_component * clone() const;
    
private:
    double U_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

class Park_recombination_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Park_recombination_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Park_recombination_component( const Park_recombination_component &c );
    
     /// \brief Default destructor
    ~Park_recombination_component();
    
    /// \brief Clone function
    Park_recombination_component * clone() const;
    
private:
    double n_;
    double D_;
    double s_v_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

class Macheret_recombination_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Macheret_recombination_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Macheret_recombination_component( const Macheret_recombination_component &c );
    
     /// \brief Default destructor
    ~Macheret_recombination_component();
    
    /// \brief Clone function
    Macheret_recombination_component * clone() const;
    
private:
    double n_;
    double theta_d_;
    
    double alpha_;
    
    double theta_dstar_,delta_d_;
    double b_;
    double fac_;
    
    bool monatomic_collider_;

private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/************************** Electron impact ionization ************************/

class Electron_impact_ionization_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Electron_impact_ionization_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Electron_impact_ionization_component( const Electron_impact_ionization_component &c );
    
     /// \brief Default destructor
    ~Electron_impact_ionization_component();
    
    /// \brief Clone function
    Electron_impact_ionization_component * clone() const;
    
private:
    double I_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/************************** Associative ionization ************************/

class Associative_ionization_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Associative_ionization_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    Associative_ionization_component( const Associative_ionization_component &c );
    
     /// \brief Default destructor
    ~Associative_ionization_component();
    
    /// \brief Clone function
    Associative_ionization_component * clone() const;
    
private:
    double alpha_;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/************************** Electron impact ionization-recombination ************************/

class EII_recombination_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    EII_recombination_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    EII_recombination_component( const Electron_impact_ionization_component &c );
    
     /// \brief Default destructor
    ~EII_recombination_component();
    
    /// \brief Clone function
    EII_recombination_component * clone() const;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/************************** Associative ionization-recombination ************************/

class AI_recombination_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    AI_recombination_component( lua_State * L, Reaction * r, int idc );
    
    /// \brief Clone constructor
    AI_recombination_component( const Associative_ionization_component &c );
    
     /// \brief Default destructor
    ~AI_recombination_component();
    
    /// \brief Clone function
    AI_recombination_component * clone() const;
    
private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );
    
    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/****************** Knab vanishing component ********************/

class Knab_vanishing_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Knab_vanishing_component( lua_State * L, Reaction * r, int idc );

    /// \brief Clone constructor
    Knab_vanishing_component( const Knab_vanishing_component &c );

     /// \brief Default destructor
    ~Knab_vanishing_component();

    /// \brief Clone function
    Knab_vanishing_component * clone() const;

private:
    double U0_;
    double U1_;
    double alpha_;
    double A_var_;

private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );

    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

/****************** Knab appearing component ********************/

class Knab_appearing_component : public Coupling_component {
public:
    /// \brief Lua file constructor
    Knab_appearing_component( lua_State * L, Reaction * r, int idc );

    /// \brief Clone constructor
    Knab_appearing_component( const Knab_appearing_component &c );

     /// \brief Default destructor
    ~Knab_appearing_component();

    /// \brief Clone function
    Knab_appearing_component * clone() const;

private:
    double U0_;
    double U1_;
    double alpha_;
    double A_var_;

private:
    /// \brief Compute the change in total average energy ( e_star - e_old ) * ( N_new - N_old )
    double specific_compute_contribution( Gas_data &Q, std::vector<double> &delta_c );

    /// \brief Compute the the source term from this component in J
    double specific_compute_source_term( Gas_data &Q, std::vector<double> &dcdt );
};

#endif
