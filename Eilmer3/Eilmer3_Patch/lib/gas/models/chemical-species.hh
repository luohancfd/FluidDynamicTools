// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef CHEMICAL_SPECIES_HH
#define CHEMICAL_SPECIES_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "gas_data.hh"
#include "species-energy-modes.hh"
#include "../../nm/source/segmented-functor.hh"

class Chemical_species;

Chemical_species * new_chemical_species_from_file( std::string name, std::string inFile );

class Chemical_species {
public:
    Chemical_species( std::string name, std::string type, int isp, double min_massf, lua_State * L  );
    virtual ~Chemical_species();
    
    Species_energy_mode* get_mode_pointer_from_type( std::string type );

    Electronic* get_electronic_mode_pointer();
    
    Multi_level_electronic* get_multi_level_electronic_mode_pointer();

    int get_element_count( std::string X );
    
    void partial_equilibrium_participants( std::vector<int> &betas, 
    					   std::vector<std::string> &E );
    
    // NOTE: following functions give values per kg of this species (ie unweighted)
    
    double eval_energy(const Gas_data &Q)
    { return s_eval_energy(Q); }
    
    double eval_modal_energy(int itm, const Gas_data &Q)
    { return s_eval_modal_energy(itm, Q); }

    double eval_enthalpy(const Gas_data &Q)
    { return s_eval_enthalpy(Q); }
    
    double eval_CEA_enthalpy(const Gas_data &Q)
    { return s_eval_CEA_enthalpy(Q); }
    
    double eval_modal_enthalpy(int itm, const Gas_data &Q)
    { return s_eval_modal_enthalpy(itm,Q); }
    
    double eval_entropy(const Gas_data &Q)
    { return s_eval_entropy(Q); }
    
    double eval_Cv(const Gas_data &Q)
    { return s_eval_Cv(Q); }
    
    virtual double eval_Cv_elec( const Gas_data &Q )
    { return modes_[1]->eval_Cv(Q); }	// electronic is always the second mode

    double eval_Cp(const Gas_data &Q)
    { return s_eval_Cp(Q); }
    
    double eval_gibbs_free_energy( double T )
    { return s_eval_gibbs_free_energy(T); }
    
    double eval_CEA_Gibbs_free_energy( double T )
    { return s_eval_CEA_Gibbs_free_energy(T); }
    
    double eval_partition_function( double T )
    { return s_eval_partition_function(T); }

    std::string get_name()
    { return name_; }
    
    std::string get_type()
    { return type_; }
    
    int get_n_modes()
    { return (int) modes_.size(); }
    
    Species_energy_mode* get_mode_pointer( int iem )
    { return modes_[iem]; }
    
    double get_h_f()
    { return h_f_; }
    
    int get_Z()
    { return Z_; }
    
    double get_R()
    { return R_; }
    
    double get_M()
    { return M_; }
    
    int get_isp()
    { return isp_; }
    
    double get_eps0()
    { return eps0_; }
    
    double get_sigma()
    { return sigma_; }
    
    int get_iT_trans()
    { return modes_[0]->get_iT(); }	// translation is always the first mode
    
    int get_iT_elec()
    { return modes_[1]->get_iT(); }	// electronic mode always second

protected:
    std::string name_;
    std::string type_;
    int isp_;
    double R_;
    double M_;
    double s_0_;
    double h_f_;
    double I_;
    int Z_;
    double eps0_;
    double sigma_;
    double min_massf_;
    std::vector<Species_energy_mode*> modes_;
    Segmented_functor * h_;
    Segmented_functor * s_;
    Segmented_functor * Cp_;
    
    double s_eval_energy(const Gas_data &Q);
    double s_eval_modal_energy(int itm, const Gas_data &Q);
    double s_eval_enthalpy(const Gas_data &Q);
    double s_eval_CEA_enthalpy(const Gas_data &Q);
    double s_eval_modal_enthalpy(int itm, const Gas_data &Q);
    virtual double s_eval_entropy(const Gas_data &Q);
    double s_eval_Cv(const Gas_data &Q);
    double s_eval_Cp(const Gas_data &Q);
    double s_eval_CEA_Gibbs_free_energy( double T );
    double s_eval_gibbs_free_energy(double T);
    virtual double s_eval_partition_function(double T);
};

class Atomic_species : public Chemical_species {
public:
    Atomic_species( std::string name, std::string type, int isp, double min_massf, lua_State * L );
    ~Atomic_species();
};

class Diatomic_species : public Chemical_species {
public:
    Diatomic_species( std::string name, std::string type, int isp, double min_massf, lua_State * L );
    ~Diatomic_species();
    
    int get_iT_elec()
    { return modes_[1]->get_iT(); }	// electronic mode always second
    
    int get_iT_rot()
    { return modes_[2]->get_iT(); }	// rotational mode always third
    
    int get_iT_vib()
    { return modes_[3]->get_iT(); }	// vibrational mode always fourth
    
    bool get_polar_flag()
    { return polar_flag_; }
    
    double get_r0()
    { return r0_; }
    
    double get_r_eq()
    { return r_eq_; }
    
    double get_f_m()
    { return f_m_; }
    
    double get_mu()
    { return mu_; }
    
    double get_alpha()
    { return alpha_; }
    
    double get_mu_B()
    { return mu_B_; }
    
    double get_theta_v()
    { return theta_v_; }
    
    double eval_Cv_rot( const Gas_data &Q )
    { return modes_[2]->eval_Cv(Q); }	// rotational mode always third
    
    double eval_Cv_vib( const Gas_data &Q )
    { return modes_[3]->eval_Cv(Q); }	// vibrational mode always fourth
    
private:
    std::string oscillator_type_;
    bool polar_flag_;
    double r0_;
    double r_eq_;
    double f_m_;
    double mu_;
    double alpha_;
    double mu_B_;
    double theta_v_;
};

class Fully_coupled_diatomic_species : public Chemical_species {
public:
    Fully_coupled_diatomic_species( std::string name, std::string type, int isp, double min_massf, lua_State * L );
    ~Fully_coupled_diatomic_species();
    
    void set_modal_temperature_indices();

    int get_iT_elec()
    { return modes_[1]->get_iT(); }	// electronic mode always second
    
    int get_iT_rot()
    { return modes_[2]->get_iT(); }	// rotational mode always third
    
    int get_iT_vib()
    { return modes_[3]->get_iT(); }	// vibrational mode always fourth
    
    bool get_polar_flag()
    { return polar_flag_; }
    
    double get_r0()
    { return r0_; }

    double get_r_eq()
    { return r_eq_; }

    double get_f_m()
    { return f_m_; }

    double get_mu()
    { return mu_; }
    
    double get_alpha()
    { return alpha_; }
    
    double get_mu_B()
    { return mu_B_; }
    
    Fully_coupled_diatom_internal * get_fcd_int_pointer()
    { return fcd_int_; }

#   if TABULATED_COUPLED_DIATOMIC_MODES==0
    int get_nlevs()
    { return elevs_.size(); }

    Diatom_electronic_level * get_elev_pointer(int ilev)
    { return elevs_[ilev]; }
#   else
    int get_nlevs()
    { return 0; }

    Diatom_electronic_level * get_elev_pointer(int ilev)
    { return NULL; }
#   endif

private:
    std::string oscillator_type_;
    bool polar_flag_;
    double r0_;
    double r_eq_;
    double f_m_;
    double mu_;
    double alpha_;
    double mu_B_;

#   if TABULATED_COUPLED_DIATOMIC_MODES==0
    std::vector<Diatom_electronic_level*> elevs_;
#   endif
    
    Fully_coupled_diatom_internal * fcd_int_;
    
    double s_eval_entropy(const Gas_data &Q);

    double s_eval_partition_function(double T);
};

class Polyatomic_species : public Chemical_species {
public:
    Polyatomic_species( std::string name, std::string type, int isp, double min_massf, lua_State * L );
    ~Polyatomic_species();
    
    int get_iT_elec()
    { return modes_[1]->get_iT(); }	// electronic mode always second
    
    int get_iT_rot()
    { return modes_[2]->get_iT(); }	// rotational mode always third
    
    int get_iT_vib()
    { return modes_[3]->get_iT(); }	// primary vibrational mode always fourth
    
    double get_theta_v()
    { return theta_v_; }
    
    bool get_linear_flag()
    { return linear_flag_; }

    double eval_Cv_rot( const Gas_data &Q )
    { return modes_[2]->eval_Cv(Q); }
    
    double eval_Cv_vib( const Gas_data &Q )
    { return s_eval_Cv_vib(Q); }
    
private:
    std::string oscillator_type_;
    bool linear_flag_;
    double theta_v_;
    
    double s_eval_Cv_vib( const Gas_data &Q );
};

class Fully_coupled_polyatomic_species : public Chemical_species {
public:
    Fully_coupled_polyatomic_species( std::string name, std::string type, int isp, double min_massf, lua_State * L );
    ~Fully_coupled_polyatomic_species();

    void set_modal_temperature_indices();

    int get_iT_elec()
    { return modes_[1]->get_iT(); }     // electronic mode always second

    int get_iT_rot()
    { return modes_[2]->get_iT(); }     // rotational mode always third

    int get_iT_vib()
    { return modes_[3]->get_iT(); }     // vibrational mode always fourth

    bool get_polar_flag()
    { return polar_flag_; }

    bool get_linear_flag()
    { return linear_flag_; }

    double eval_Cv_rot( const Gas_data &Q )
    { return 0.0; }

    double eval_Cv_vib( const Gas_data &Q )
    { return 0.0; }

    Fully_coupled_polyatom_internal * get_fcp_int_pointer()
    { return fcp_int_; }

private:
    std::string oscillator_type_;
    bool polar_flag_;
    bool linear_flag_;

    std::vector<Polyatom_electronic_level*> elevs_;

    Fully_coupled_polyatom_internal * fcp_int_;

    double s_eval_entropy(const Gas_data &Q)
    { return 0.0; }
};

class Free_electron_species : public Chemical_species {
public:
    Free_electron_species( std::string name, std::string type, int isp, double min_massf, lua_State * L );
    ~Free_electron_species();
    
    double eval_Cv_elec( const Gas_data &Q )
    { return 0.0; }
};


#endif
