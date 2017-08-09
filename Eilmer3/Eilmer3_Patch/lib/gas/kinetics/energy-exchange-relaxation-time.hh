// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_RELAXATION_TIME_HH
#define ENERGY_EXCHANGE_RELAXATION_TIME_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/gas-model.hh"

class Relaxation_time {
public:
    Relaxation_time();
    
    virtual ~Relaxation_time();

    double compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return specific_relaxation_time(Q,molef); }

    double compute_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return specific_transition_probability(Q,molef); }
protected:
    virtual double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef) = 0;
    virtual double specific_transition_probability(Gas_data &Q, std::vector<double> &molef) = 0;
};

class VT_MillikanWhite : public Relaxation_time {
public:
    VT_MillikanWhite(lua_State *L, int ip, int iq, int itrans);

    ~VT_MillikanWhite();
private:
    int ip_, iq_;
    double M_p_, M_q_;
    double theta_;
    int iT_;
    double mu_;
    double a_;
    double b_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class VT_LandauTeller_cf : public Relaxation_time {
public:
    VT_LandauTeller_cf(lua_State *L, int ip, int iq, int itrans);

    ~VT_LandauTeller_cf();
private:
    int ip_, iq_, iT_;
    double A_, B_, C_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class VT_Thivet_cf : public Relaxation_time {
public:
    VT_Thivet_cf(lua_State *L, int ip, int iq, int itrans);

    ~VT_Thivet_cf();
private:
    int ip_, iq_, iT_;
    double B_, C_;

    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class VT_high_temperature_cross_section {
public:
    VT_high_temperature_cross_section() {}

    virtual ~VT_high_temperature_cross_section() {}
    
    double eval_sigma( double T )
    { return s_eval_sigma(T); }
protected:
    virtual double s_eval_sigma( double T ) = 0;
};

class Park_VT_HTCS : public VT_high_temperature_cross_section {
public:
    Park_VT_HTCS(lua_State *L);
    
    ~Park_VT_HTCS();
private:
    double sigma_dash_;
    
    double s_eval_sigma( double T );
};

class Fujita_VT_HTCS : public VT_high_temperature_cross_section {
public:
    Fujita_VT_HTCS();
    
    ~Fujita_VT_HTCS();
private:
    double s_eval_sigma( double T );
};

VT_high_temperature_cross_section*
create_VT_high_temperature_cross_section_model( lua_State *L );


class VT_MillikanWhite_HTC : public Relaxation_time {
public:
    VT_MillikanWhite_HTC(lua_State *L, int ip, int iq, int itrans);

    ~VT_MillikanWhite_HTC();
private:
    int ip_, iq_, iT_;
    double m_av_;
    VT_MillikanWhite *VT_MW_;
    VT_high_temperature_cross_section *HTC_model_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ET_AppletonBray_Ion : public Relaxation_time {
public:
    ET_AppletonBray_Ion(lua_State *L, int ie, int ic);
    
    ~ET_AppletonBray_Ion();
private:
    int ie_;
    int ic_;
    int iTe_;
    double M_c_;
    double M_e_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ET_AppletonBray_Neutral : public Relaxation_time {
public:
    ET_AppletonBray_Neutral(lua_State *L, int ie, int ic);
    
    ~ET_AppletonBray_Neutral();
private:
    int ic_;
    double M_c_;
    int iTe_;
    std::vector<double> C_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ET_AppletonBray_TwoRangeNeutral : public Relaxation_time {
public:
    ET_AppletonBray_TwoRangeNeutral(lua_State *L, int ie, int ic);

    ~ET_AppletonBray_TwoRangeNeutral();
private:
    int ic_;
    double M_c_;
    int iTe_;
    std::vector<double> LT_C_, HT_C_;
    double T_switch_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ER_Abe_Ion : public Relaxation_time {
public:
    ER_Abe_Ion(lua_State *L, int ie, int ic);

    ~ER_Abe_Ion();
private:
    int ie_;
    int ic_;
    int iTe_;
    double M_c_;
    double M_e_;
    double g_rot_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ER_Abe_Neutral : public Relaxation_time {
public:
    ER_Abe_Neutral(lua_State *L, int ie, int ic);

    ~ER_Abe_Neutral();
private:
    int ic_;
    double M_c_;
    double g_rot_;
    int iTe_;
    std::vector<double> C_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};


// Helper functions for SSH calculations
double collider_distance_A(Diatomic_species &p, Diatomic_species &q);
double collider_distance_B(Diatomic_species &p, Diatomic_species &n);
double xi_correction(Diatomic_species &p, Diatomic_species &n);
double potential_well_A(Diatomic_species &p, Diatomic_species &q);
double potential_well_B(Diatomic_species &p, Diatomic_species &n);
double collision_frequency(double sigma, double mu, double T, double nd);
double SSH_beta(double eps0, double mu, double r0,
		double del_E, double T);
double SSH_D_star(double beta);
double SSH_r_c_star(double beta, double r0);
double SSH_r_c_star2(double D_star, double r0);
double SSH_del_star(double beta, double r0);
double SSH_del_star2(double D_star, double r0);
double SSH_alpha_pq(double mu, double del_E, double del_star);
double SSH_chi_pq(double alpha_pq, double T);
double SSH_Z_0(double del_star, double r_eq);
double SSH_Z_V(double f_m, double mu_pp, double mu_pq, double alpha_pq,
	       double theta_v, double del_E, int i);
double SSH_Z_T(double del_E, double alpha_pq, double T);
double SSH_Z_plus(double eps0, double chi, double T);
double SSH_A_factor(double r_c, double sigma);

class VT_SSH : public Relaxation_time {
public:
    VT_SSH(lua_State *L, int ip, int iq, int itrans);
    ~VT_SSH();
private:
    int ip_, iq_, iT_;
    double M_p_, r_eq_p_, f_m_p_, mu_p_, theta_v_p_;
    double M_q_, r_eq_q_, f_m_q_, mu_q_, theta_v_q_;
    double r_, eps_, mu_, delta_E_, sigma_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef);
};

class VV_SSH : public Relaxation_time {
public:
    VV_SSH(lua_State *L, int ip, int iq, int itrans);
    ~VV_SSH();
private:
    int ip_, iq_, iT_;
    double M_p_, r_eq_p_, f_m_p_, mu_p_, theta_v_p_;
    double M_q_, r_eq_q_, f_m_q_, mu_q_, theta_v_q_;
    double r_, eps_, mu_, delta_E_, sigma_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef);
};

class VV_MTLandauTeller : public Relaxation_time {
public:
    VV_MTLandauTeller(lua_State *L, int ip, int iq, int itrans);
    ~VV_MTLandauTeller();
private:
    int ip_, iq_, iT_;
    double A_, n_, B1_, B2_, B3_, beta_, theta_v_p_, theta_v_q_ ;
    double mu_, R0_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef);
    double rate_coefficient(double T);
};

class VV_from_eq : public Relaxation_time {
public:
    VV_from_eq(lua_State *L, int ip, int iq, int itrans);
    ~VV_from_eq();
private:
    int ip_, iq_, iT_;
    Relaxation_time* rt_;
    double theta_v_p_, theta_v_q_, R0_, mu_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef);
};

// class VV_Candler : public Relaxation_time {
// public:
//     VV_Candler(lua_State *L);
    
//     ~VV_Candler();
// private:
//     // Global data
//     int iT_;			// translational temperature
//     double sigma_;		// collision cross-section (dxd)
//     double mu_;			// reduced molecular weight
//     double P_;			// transition probability
    
//     // Species p data
//     int ip_;			// species index
    
//     // Species q data
//     int iq_;			// species index
//     double M_q_;		// molecular weight
      
//     double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
//     double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
//     { return 0.0; }
// };

class VE_Lee : public Relaxation_time {
public:
    VE_Lee(lua_State *L, int ie, int iv);
    
    ~VE_Lee();
private:
    // Electron data
    int ie_;
    int iTe_;
    
    // Vibrational species data
    int iv_;
    
    std::vector<double> T_switches_;
    std::vector< std::vector<double> > ptau_coeffs_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

// class RT_Parker : public Relaxation_time {
// public:
//     RT_Parker(lua_State *L);
    
//     ~RT_Parker();
// private:
//     // Global data
//     int iT_;					// translational temperature
//     std::vector<double> sigmas_;		// rotational collision cross-sections
//     std::vector<double> mus_;			// reduced molecular weights

//     // Rotator data
//     int ip_;					// species index
    
//     // Collider data
//     std::vector<int> iqs_;			// species indices
//     std::vector<double> M_qs_;			// molecular weights

//     double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
//     double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
//     { return 0.0; }
// };

Relaxation_time* create_new_relaxation_time(lua_State *L, int ip, int iq, int itrans);
void parse_input_for_rts(std::string cfile, Gas_model &g, lua_State *L);
Relaxation_time* get_rt_from_file(int irt, std::string cfile, Gas_model &g);
int get_no_rts_from_file(std::string cfile, Gas_model &g);


#endif
