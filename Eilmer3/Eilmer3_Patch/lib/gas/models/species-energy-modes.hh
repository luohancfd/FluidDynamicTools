// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#ifndef SPECIES_ERNEGY_MODES_HH
#define SPECIES_ERNEGY_MODES_HH

#define EVAL_ENTROPY_AT_1ATM 1
#define TABULATED_COUPLED_DIATOMIC_MODES 0

#include <string>

#include "gas_data.hh"
#include "physical_constants.hh"
#include "diatom-electronic-level.hh"
#include "coupled-diatom-LUT.hh"
#include "polyatom-electronic-level.hh"
#include "../../nm/source/segmented-functor.hh"

class Species_energy_mode {
public:
    Species_energy_mode( int isp=-1, double R=0.0, double min_massf_=1.0e-10,
    	                 std::string type="none", double theta=0.0, int iT=-1 );
    virtual ~Species_energy_mode() {}
    
    void set_iT(int iT)
    { iT_ = iT; }
    
    int get_iT()
    { return iT_; }
    
    int get_isp()
    { return isp_; }
    
    double eval_weighted_energy( const Gas_data &Q )
    { return Q.massf[isp_]*s_eval_energy( Q ); }
    
    double eval_energy( const Gas_data &Q )
    { return s_eval_energy( Q ); }
    
    double eval_energy_from_T( double T, double A=-1.0 )
    { return s_eval_energy_from_T( T, A ); }

    double eval_weighted_enthalpy( const Gas_data &Q )
    { return Q.massf[isp_]*s_eval_enthalpy( Q ); }
    
    double eval_enthalpy( const Gas_data &Q )
    { return s_eval_enthalpy( Q ); }
    
    double eval_enthalpy_from_T( double T, double A=-1.0 )
    { return s_eval_enthalpy_from_T( T, A ); }

    double eval_weighted_entropy( const Gas_data &Q )
    { return Q.massf[isp_]*s_eval_entropy(Q); }
    
    double eval_entropy( const Gas_data &Q )
    { return s_eval_entropy(Q); }
    
    double eval_entropy_from_T( double T )
    { return s_eval_entropy_from_T( T ); }

    double eval_weighted_Cv( Gas_data &Q )
    { return Q.massf[isp_]*s_eval_Cv( Q ); }
    
    double eval_Cv( const Gas_data &Q )
    { return s_eval_Cv( Q ); }
    
    double eval_Cv_from_T( double T )
    { return s_eval_Cv_from_T( T ); }

    double eval_weighted_Cp( Gas_data &Q )
    { return Q.massf[isp_]*s_eval_Cp( Q ); }
    
    double eval_Cp( const Gas_data &Q )
    { return s_eval_Cp( Q ); }
    
    double eval_Cp_from_T( double T )
    { return s_eval_Cp_from_T( T ); }

    double eval_Q_from_T( double T=0.0, double A=-1.0 )
    { return s_eval_Q_from_T(T,A); }
    
    std::string get_type()
    { return type_; }

    double get_theta()
    { return theta_; }
    
protected:
    int isp_;
    double R_;
    double min_massf_;
    std::string type_;
    double theta_;
    int iT_;
    
    virtual double s_eval_energy( const Gas_data &Q  ) = 0;
    virtual double s_eval_energy_from_T( double T, double A ) = 0;
    virtual double s_eval_enthalpy( const Gas_data &Q  ) = 0;
    virtual double s_eval_enthalpy_from_T( double T, double A ) = 0;
    virtual double s_eval_entropy( const Gas_data &Q ) = 0;
    virtual double s_eval_entropy_from_T( double T ) = 0;
    virtual double s_eval_Cv( const Gas_data &Q  ) = 0;
    virtual double s_eval_Cv_from_T( double T ) = 0;
    virtual double s_eval_Cp( const Gas_data &Q  ) = 0;
    virtual double s_eval_Cp_from_T( double T ) = 0;

    virtual double s_eval_Q_from_T( double T, double A ) = 0;
    
};

class Electronic : public Species_energy_mode {
public:
    Electronic( int isp, double R, double min_massf, double theta );
    ~Electronic() {}
    
protected:  
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_energy(Q); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_T(T,A); }
    double s_eval_Cp( const Gas_data &Q  ) { return s_eval_Cv(Q); }
    double s_eval_Cp_from_T( double T  ) { return s_eval_Cv_from_T(T); }
};

class One_level_electronic : public Electronic {
public:
    One_level_electronic( int isp, double R, double min_massf, int g, double theta );
    ~One_level_electronic() {}
    
private:
    int g_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

class Two_level_electronic : public Electronic {
public:
    Two_level_electronic( int isp, double R, double min_massf, int g0, 
    	                  double theta0, int g1, double theta1 );
    ~Two_level_electronic() {}
    
private:
    int g0_;
    int g1_;
    double delta_theta_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

class Multi_level_electronic : public Electronic {
public:
    Multi_level_electronic( int isp, double R, double min_massf, 
    			    std::vector<int> &g, std::vector<double> &theta);
    ~Multi_level_electronic() {}
    
    int get_n_levels()
    { return (int) g_vec_.size(); }

    int get_g( int ilev )
    { return g_vec_[ilev]; }

    double get_theta( int ilev )
    { return theta_vec_[ilev]; }

private:
    std::vector<int> g_vec_;
    std::vector<double> theta_vec_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

#if TABULATED_COUPLED_DIATOMIC_MODES==0
class Coupled_diatomic_electronic : public Electronic {
public:
    Coupled_diatomic_electronic( int isp, double R, double min_massf, int sigma_r, 
    	double fT_, std::vector<Diatom_electronic_level*> &elevs );
    ~Coupled_diatomic_electronic();
    
    void set_iTs( int iTe, int iTv, int iTr )
    { iTe_ = iTe; iTv_ = iTv; iTr_ = iTr; }

    Diatom_electronic_level * get_elev_pointer( int ilev )
    { return elevs_[ilev]; }

private:
    int iTe_;
    int iTv_;
    int iTr_;
    
    int sigma_r_;
    double m_;
    double fT_;
    std::vector<Diatom_electronic_level*> elevs_;

    double eval_total_internal_energy( double T_el, double T_vib, double T_rot );
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_energy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_entropy_from_T( double T ) { return s_eval_entropy_from_Ts(T,T,T); }
    double s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cv_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#else
class Coupled_diatomic_electronic : public Electronic {
public:
    Coupled_diatomic_electronic( int isp, double R, double min_massf, double theta, std::string lut_fname );
    ~Coupled_diatomic_electronic();
    
    void set_iTs( int iTe, int iTv, int iTr )
    { iTe_ = iTe; iTv_ = iTv; iTr_ = iTr; }

private:
    int iTe_;
    int iTv_;
    int iTr_;
    
    double m_;
    NoneqCoupledDiatomicLUT * e_LUT_;
    NoneqCoupledDiatomicLUT * Cv_LUT_;

    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_energy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_entropy_from_T( double T ) { return s_eval_entropy_from_Ts(T,T,T); }
    double s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cv_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#endif

class Fully_excited_translation : public Species_energy_mode {
public:
    Fully_excited_translation( int isp, double R, double min_massf );
    ~Fully_excited_translation() {}
    
private:
    double Cv_;
    double Cp_;
    double entropy_constant_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_enthalpy_from_T(Q.T[iT_],-1.0); }
    double s_eval_enthalpy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q );
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Cp( const Gas_data &Q ) { return s_eval_Cp_from_T(Q.T[iT_]); }
    double s_eval_Cp_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

class Rotation : public Species_energy_mode {
public:
    Rotation( int isp, double R, double min_massf, double theta );
    ~Rotation() {}
    
protected:   
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_energy(Q); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_T(T,A); }
    double s_eval_Cp( const Gas_data &Q  ) { return s_eval_Cv(Q); }
    double s_eval_Cp_from_T( double T  ) { return s_eval_Cv_from_T(T); }
};

class Fully_excited_rotation : public Rotation {
public:
    // NOTE: sigma = 1 or 2 for hetero/homonuclear molecules
    Fully_excited_rotation( int isp, double R, double min_massf, double theta, int sigma );
    ~Fully_excited_rotation() {}
    
private:
    double Cv_;
    int sigma_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

class Fully_excited_nonlinear_rotation : public Rotation {
public:
    Fully_excited_nonlinear_rotation( int isp, double R, double min_massf, 
    				      double theta_A0, double theta_B0, double theta_C0,
    				      int sigma );
    ~Fully_excited_nonlinear_rotation() {}
    
private:
    double Cv_;
    double theta_A0_;
    double theta_B0_;
    double theta_C0_;
    int sigma_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

#if TABULATED_COUPLED_DIATOMIC_MODES==0
class Coupled_diatomic_rotation : public Rotation {
public:
    Coupled_diatomic_rotation( int isp, double R, double min_massf, int sigma_r,
    	double fT_, std::vector<Diatom_electronic_level*> &elevs );
    ~Coupled_diatomic_rotation();
    
    void set_iTs( int iTe, int iTv, int iTr )
    { iTe_ = iTe; iTv_ = iTv; iTr_ = iTr; }

private:
    int iTe_;
    int iTv_;
    int iTr_;
    
    int sigma_r_;
    double m_;
    double fT_;
    std::vector<Diatom_electronic_level*> elevs_;

    double eval_total_internal_energy( double T_el, double T_vib, double T_rot );
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_energy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_entropy_from_T( double T ) { return s_eval_entropy_from_Ts(T,T,T); }
    double s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cv_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#else
class Coupled_diatomic_rotation : public Rotation {
public:
    Coupled_diatomic_rotation( int isp, double R, double min_massf, double theta_r, std::string lut_fname );
    ~Coupled_diatomic_rotation();
    
    void set_iTs( int iTe, int iTv, int iTr )
    { iTe_ = iTe; iTv_ = iTv; iTr_ = iTr; }

private:
    int iTe_;
    int iTv_;
    int iTr_;
    
    double m_;
    NoneqCoupledDiatomicLUT * e_LUT_;
    NoneqCoupledDiatomicLUT * Cv_LUT_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_energy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_enthalpy_from_Ts( double T_el, double T_vib, double T_rot ) { return s_eval_energy_from_Ts(T_el,T_vib,T_rot); }
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_entropy_from_T( double T ) { return s_eval_entropy_from_Ts(T,T,T); }
    double s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cv_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cp( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cp_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cp_from_Ts( double T_el, double T_vib, double T_rot ) { return s_eval_Cv_from_Ts(T_el,T_vib,T_rot); }
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#endif

class Vibration : public Species_energy_mode {
public:
    Vibration( int isp, double R, double min_massf, double theta );
    ~Vibration() {}
    
protected:   
    double s_eval_enthalpy( const Gas_data &Q  ) { return s_eval_energy(Q); }
    double s_eval_enthalpy_from_T( double T, double A  ) { return s_eval_energy_from_T(T,A); }
    double s_eval_Cp( const Gas_data &Q  ) { return s_eval_Cv(Q); }
    double s_eval_Cp_from_T( double T  ) { return s_eval_Cv_from_T(T); }
};

class Harmonic_vibration : public Vibration {
public:
    Harmonic_vibration( int isp, double R, double min_massf, double theta );
    ~Harmonic_vibration() {}
    
private:
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};

class Anharmonic_vibration : public Vibration {
public:
    Anharmonic_vibration( int isp, double R, double min_massf, double theta, 
        double thetaR, int sigma,
        Segmented_functor* h, Segmented_functor* s, Segmented_functor* Cp_);
    ~Anharmonic_vibration() {}
    
private:
    Segmented_functor * h_;
    Segmented_functor * s_;
    Segmented_functor * Cp_;
    double thetaR_;
    int sigma_;
    double Tt_;
    double p_;
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q );
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Q_from_T( double T, double A );
};


class Truncated_harmonic_vibration : public Vibration {
public:
    Truncated_harmonic_vibration( int isp, double R, double min_massf, double theta, double theta_D_ );
    ~Truncated_harmonic_vibration() {}
    
    double eval_HO_energy_from_T( double T )
    { return s_eval_HO_energy( T ); }
    
private:
    double theta_D_;
    
private:
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_HO_energy( double T );
    double s_eval_Q_from_T( double T, double A );
};

#if TABULATED_COUPLED_DIATOMIC_MODES==0
class Coupled_diatomic_vibration : public Vibration {
public:
    Coupled_diatomic_vibration( int isp, double R, double min_massf, int sigma_r, 
    	double fT, std::vector<Diatom_electronic_level*> &elevs );
    ~Coupled_diatomic_vibration();
    
    void set_iTs( int iTe, int iTv, int iTr )
    { iTe_ = iTe; iTv_ = iTv; iTr_ = iTr; }

private:
    int iTe_;
    int iTv_;
    int iTr_;
    
    int sigma_r_;
    double m_;
    double fT_;
    std::vector<Diatom_electronic_level*> elevs_;
    
    double eval_total_internal_energy( double T_el, double T_vib, double T_rot );
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_energy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_enthalpy_from_Ts( double T_el, double T_vib, double T_rot ) { return s_eval_energy_from_Ts(T_el,T_vib,T_rot); }
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_entropy_from_T( double T ) { return s_eval_entropy_from_Ts(T,T,T); }
    double s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cv_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cp( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cp_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cp_from_Ts( double T_el, double T_vib, double T_rot ) { return s_eval_Cv_from_Ts(T_el,T_vib,T_rot); }
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#else
class Coupled_diatomic_vibration : public Vibration {
public:
    Coupled_diatomic_vibration( int isp, double R, double min_massf, double theta_v, std::string lut_fname );
    ~Coupled_diatomic_vibration();
    
    void set_iTs( int iTe, int iTv, int iTr )
    { iTe_ = iTe; iTv_ = iTv; iTr_ = iTr; }

private:
    int iTe_;
    int iTv_;
    int iTr_;
    
    double m_;
    NoneqCoupledDiatomicLUT * e_LUT_;
    NoneqCoupledDiatomicLUT * Cv_LUT_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_energy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_energy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_Ts(T,T,T); }
    double s_eval_enthalpy_from_Ts( double T_el, double T_vib, double T_rot ) { return s_eval_energy_from_Ts(T_el,T_vib,T_rot); }
    double s_eval_entropy( const Gas_data &Q ) { return s_eval_entropy_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_entropy_from_T( double T ) { return s_eval_entropy_from_Ts(T,T,T); }
    double s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cv_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot );
    double s_eval_Cp( const Gas_data &Q ) { return s_eval_Cv_from_Ts(Q.T[iTe_],Q.T[iTv_],Q.T[iTr_]); }
    double s_eval_Cp_from_T( double T ) { return s_eval_Cv_from_Ts(T,T,T); }
    double s_eval_Cp_from_Ts( double T_el, double T_vib, double T_rot ) { return s_eval_Cv_from_Ts(T_el,T_vib,T_rot); }
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#endif

#if TABULATED_COUPLED_DIATOMIC_MODES==0
class Fully_coupled_diatom_internal : public Species_energy_mode {
public:
    Fully_coupled_diatom_internal( int isp, double R, double min_massf, int sigma_r, 
    	double fT, std::vector<Diatom_electronic_level*> &elevs );
    ~Fully_coupled_diatom_internal();
    
    Diatom_electronic_level * get_elev_pointer( int ilev );

    int get_nlevs() { return (int) elevs_.size(); }

private:
    int sigma_r_;
    double m_;
    double fT_;
    std::vector<Diatom_electronic_level*> elevs_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_enthalpy( const Gas_data &Q  ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_T(T,A); }
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Cp( const Gas_data &Q  ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cp_from_T( double T  ) { return s_eval_Cv_from_T(T); }
    double s_eval_Q_from_T( double T, double A );
};
#else
class Fully_coupled_diatom_internal : public Species_energy_mode {
public:
    Fully_coupled_diatom_internal( int isp, double R, double min_massf, std::string lut_fname );
    ~Fully_coupled_diatom_internal();

private:
    double m_;
    EqCoupledDiatomicLUT * e_LUT_;
    EqCoupledDiatomicLUT * Cv_LUT_;
    EqCoupledDiatomicLUT * s_LUT_;
    
    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_enthalpy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_T(T,A); }
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Cp( const Gas_data &Q  ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cp_from_T( double T  ) { return s_eval_Cv_from_T(T); }
    double s_eval_Q_from_T( double T, double A ) { return 0.0; }
};
#endif

class Fully_coupled_polyatom_internal : public Species_energy_mode {
public:
    Fully_coupled_polyatom_internal( int isp, double R, double min_massf,
        double fT, std::vector<Polyatom_electronic_level*> &elevs );
    ~Fully_coupled_polyatom_internal();

    Polyatom_electronic_level * get_elev_pointer( int ilev );

    int get_nlevs() { return (int) elevs_.size(); }

private:
    double m_;
    double fT_;
    std::vector<Polyatom_electronic_level*> elevs_;

    double s_eval_energy( const Gas_data &Q ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_energy_from_T( double T, double A );
    double s_eval_enthalpy( const Gas_data &Q  ) { return s_eval_energy_from_T(Q.T[iT_],-1.0); }
    double s_eval_enthalpy_from_T( double T, double A ) { return s_eval_energy_from_T(T,A); }
    double s_eval_entropy( const Gas_data &Q  ) { return s_eval_entropy_from_T(Q.T[iT_]); }
    double s_eval_entropy_from_T( double T );
    double s_eval_Cv( const Gas_data &Q ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cv_from_T( double T );
    double s_eval_Cp( const Gas_data &Q  ) { return s_eval_Cv_from_T(Q.T[iT_]); }
    double s_eval_Cp_from_T( double T  ) { return s_eval_Cv_from_T(T); }
    double s_eval_Q_from_T( double T, double A );
};

class Energy_of_formation : public Species_energy_mode {
public:
    Energy_of_formation(int isp, double R, double min_massf, double h_f)
	: Species_energy_mode(isp, R, min_massf, "enthalpy-of-formation", 0), h_f_(h_f) {}
    ~Energy_of_formation() {}

private:
    double h_f_;
    double s_eval_energy(const Gas_data &Q)
    { return h_f_; }
    double s_eval_energy_from_T(double T,double A)
    { return h_f_; }
    double s_eval_enthalpy(const Gas_data &Q)
    { return 0.0; }
    double s_eval_enthalpy_from_T(double T,double A)
    { return 0.0; }
    double s_eval_entropy(const Gas_data &Q)
    { return 0.0; }
    double s_eval_entropy_from_T(double T)
    { return 0.0; }
    double s_eval_Cv(const Gas_data &Q)
    { return 0.0; }
    double s_eval_Cv_from_T(double T)
    { return 0.0; }
    double s_eval_Cp(const Gas_data &Q)
    { return 0.0; }
    double s_eval_Cp_from_T(double T)
    { return 0.0; }
    double s_eval_Q_from_T(double T, double A)
    { return 0.0; }

};

#endif
