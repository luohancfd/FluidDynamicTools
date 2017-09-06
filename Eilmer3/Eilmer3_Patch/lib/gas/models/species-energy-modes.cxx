// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "../../util/source/useful.h"

#include "species-energy-modes.hh"
#include "physical_constants.hh"

// Forward declaration
std::string get_library_species_name( int isp );

using namespace std;

Species_energy_mode::
Species_energy_mode( int isp, double R, double min_massf, string type, 
                     double theta, int iT )
 : isp_( isp ), R_( R ), min_massf_( min_massf ), type_( type ), 
theta_( theta ), iT_( iT ) {}
 
/* ------- Generic electronic --------- */

Electronic::
Electronic( int isp, double R, double min_massf, double theta )
: Species_energy_mode( isp, R, min_massf, "electronic", theta )
{}
 
/* ------- One level electronic species energy mode ------- */
// NOTES: - This class uses the Multi_level_electronic expressions, simplified
//          for one level
//        - The primary use is for electrons and electronic state pseudo-species
//          where theta=0

One_level_electronic::
One_level_electronic( int isp, double R, double min_massf, int g, double theta )
 : Electronic( isp, R, min_massf, theta ), g_( g )
{}
 
double
One_level_electronic::
s_eval_energy_from_T( double T, double A )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.31
    return R_ * theta_;
}

double
One_level_electronic::
s_eval_entropy_from_T( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.36
    return R_*log( g_*exp( - theta_ / T ) ) + (R_/T)*theta_;
}

double
One_level_electronic::
s_eval_Cv_from_T( double T )
{
    return 0.0;
}

double
One_level_electronic::
s_eval_Q_from_T( double T, double A )
{
    UNUSED_VARIABLE(A);
    return g_*exp( - theta_ / T );
}

/* ------- Two level electronic species energy mode ------- */

Two_level_electronic::
Two_level_electronic( int isp, double R, double min_massf, int g0, 
                      double theta0, int g1, double theta1 )
 : Electronic( isp, R, min_massf, theta1 ), g0_( g0 ), 
   g1_( g1 ), delta_theta_( theta1 - theta0 )
{}
 
double
Two_level_electronic::
s_eval_energy_from_T( double T, double A )
{
    // Ref: Vincenti and Kruger (1975) Eq. 11.3 p 131
    double tmp = double(g1_)/double(g0_) * exp( - delta_theta_ / T );
    return R_ * delta_theta_ * tmp / ( 1.0 + tmp );
}

double
Two_level_electronic::
s_eval_entropy_from_T( double T )
{
    // Ref: My own derivation from the definition of entropy as a function of
    //      the partition function and temperature (see Eq 3.2 in PhD thesis)
    double tmp = double(g1_) * exp( -delta_theta_ / T );
    double Q = double(g0_) + tmp;

    return R_ * ( log(Q) + delta_theta_ / T / ( double(g0_) / tmp + 1.0 ) );
}

double
Two_level_electronic::
s_eval_Cv_from_T( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 11.4 p 131
    double tmp = double(g1_)/double(g0_) * exp( - delta_theta_ / T );
    return R_ * tmp / pow( T / delta_theta_ * ( 1.0 + tmp ), 2 );
    // return R_ * pow( delta_theta_ / T, 2 ) * tmp / pow( 1.0 + tmp, 2 );
}

double
Two_level_electronic::
s_eval_Q_from_T( double T, double A )
{
    double Q = double(g0_) + double(g1_) * exp( - delta_theta_ / T );

    return Q;
}

/* ------- Multi level electronic species energy mode ------- */

Multi_level_electronic::
Multi_level_electronic( int isp, double R, double min_massf, 
    			vector<int> &g, vector<double> &theta )
 : Electronic( isp, R, min_massf, theta[1] ), g_vec_( g ), theta_vec_( theta )
{}
 
double
Multi_level_electronic::
s_eval_energy_from_T( double T, double A )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.31
    double numerator = 0.0;		// sum of weighted energies
    double denominator = 0.0;		// total partition function
    double tmp;				// level partition function
    for ( size_t i=0; i<theta_vec_.size(); ++i ) {
    	tmp = double(g_vec_[i]) * exp( - theta_vec_[i] / T );
    	numerator += theta_vec_[i] * tmp;
    	denominator += tmp;
    }

    return R_ * numerator / denominator;
}

double
Multi_level_electronic::
s_eval_entropy_from_T( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.36
    double tmpA = 0.0;	// sum of weighted energies
    double tmpB = 0.0;	// total partition function
    double tmpC;	// level partition function
    for ( size_t i=0; i<theta_vec_.size(); ++i ) {
    	tmpC = double(g_vec_[i]) * exp( - theta_vec_[i] / T );
    	tmpA += theta_vec_[i] * tmpC;
    	tmpB += tmpC;
    }
    return R_*log( tmpB) + R_/T*tmpA/tmpB;
}

double
Multi_level_electronic::
s_eval_Cv_from_T( double T )
{
    // Ref: Bottin, B. "Aerothermodynamic model of an Inductively-Coupled Plasma
    //                  Wind Tunnel" VKI PhD Thesis 1999, Eq. A.45
    double tmp;           // level partition function
    double u = 0.0;	  // sum of weighted energies
    double v = 0.0;	  // total partition function
    double u_dash_star = 0.0;  // derivative of u without 1/T**2
    double v_dash_star = 0.0;  // derivative of v without 1/T**2
    for ( size_t i=0; i<theta_vec_.size(); ++i ) {
    	tmp = double(g_vec_[i]) * exp( - theta_vec_[i] / T );
    	u += theta_vec_[i] * tmp;
    	v += tmp;
    	u_dash_star += tmp * theta_vec_[i] * theta_vec_[i];
    	v_dash_star += theta_vec_[i] * tmp;
    }

    return R_/(T*T)*( u_dash_star * v - u * v_dash_star ) / ( v * v );
}

double
Multi_level_electronic::
s_eval_Q_from_T( double T, double A )
{
    double Q = 0.0;
    for ( size_t i=0; i<theta_vec_.size(); ++i ) {
        Q += double(g_vec_[i]) * exp( - theta_vec_[i] / T );
    }
    return Q;
}

/* ------- Coupled diatomic electronic mode ------- */

#if TABULATED_COUPLED_DIATOMIC_MODES==0
Coupled_diatomic_electronic::
Coupled_diatomic_electronic( int isp, double R, double min_massf, int sigma_r,
    double fT, std::vector<Diatom_electronic_level*> &elevs )
 : Electronic( isp, R, min_massf, elevs[0]->E / PC_k_SI ), sigma_r_( sigma_r ),
 m_( PC_R_u / R_ / PC_Avogadro ), fT_( fT )
{
    // Set the elev pointers
    for ( size_t ilev=0; ilev < elevs.size(); ++ilev ) {
    	elevs_.push_back( elevs[ilev] );
    }
    UNUSED_VARIABLE(sigma_r_);
}

Coupled_diatomic_electronic::
~Coupled_diatomic_electronic()
{
    // elevs are just pointers now
}

double
Coupled_diatomic_electronic::
eval_total_internal_energy( double T_el, double T_vib, double T_rot  )
{
    // Calculate TOTAL (all modes) internal energy
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
                double E_rot_dash = E_el + E_vib + E_rot;
                double Q_rot_dash = Q_el*Q_vib*Q_rot;
    	    	if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
                E_weighted += E_rot_dash * Q_rot_dash;
                Q_total += Q_rot_dash;
            }
        }
    }

    return E_weighted / Q_total / m_;   // convert J/particle -> J/kg
}

double
Coupled_diatomic_electronic::
s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    // Loop over all rovibronic levels, sum up contributions from definition
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        double Q_vib_sum = 0.0;
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            double Q_rot_sum = 0.0;
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
       	    	if ( std::isnan(Q_rot) || std::isinf( Q_rot ) ) continue;
                Q_rot_sum += Q_rot;
            }
            Q_vib_sum += Q_vib*Q_rot_sum;
        }
        Q_total += Q_el * Q_vib_sum;
        E_weighted += E_el * Q_el * Q_vib_sum;
    }

    return E_weighted / Q_total / m_;   // convert J/particle -> J/kg
}

double
Coupled_diatomic_electronic::
s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot  )
{
#   if 0
    // Loop over all rovibronic levels, sum up contributions from definition
    // NOTE: not sure that this formulation is correct
    double E_weighted = 0.0;
    double Q_total = 0.0;
    double Q_el_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        double Q_vib_sum = 0.0;
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            double Q_rot_sum = 0.0;
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
       	    	if ( isnan(Q_rot) || isinf( Q_rot ) ) continue;
                Q_rot_sum += Q_rot;
            }
            Q_vib_sum += Q_vib*Q_rot_sum;
        }
        Q_total += Q_el * Q_vib_sum;
        Q_el_total += Q_el;
        E_weighted += E_el * Q_el * Q_vib_sum;
    }
    
    double e = E_weighted / Q_total / m_;   // convert J/particle -> J/kg
    
    return R_*log(Q_el_total) + e/T_el;		// units of J/kg/K
#   else
    cout << "Coupled_diatomic_electronic::s_eval_entropy_from_Ts()" << endl
         << "No expression available." << endl;
    exit( FAILURE );
#   endif
}

double
Coupled_diatomic_electronic::
s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot  )
{
    // Step size
    double h = T_el * fT_;

    // five-point stencil method
    double fxp2h = eval_total_internal_energy(T_el + 2.0 * h,T_vib,T_rot);
    double fxph  = eval_total_internal_energy(T_el + h,       T_vib,T_rot);
    double fxmh  = eval_total_internal_energy(T_el - h,       T_vib,T_rot);
    double fxm2h = eval_total_internal_energy(T_el - 2.0 * h,T_vib,T_rot);

    return ( - fxp2h + 8.0 * fxph - 8.0 * fxmh + fxm2h ) / 12.0 / h;
}
#else
Coupled_diatomic_electronic::
Coupled_diatomic_electronic( int isp, double R, double min_massf, double theta, string lut_fname )
 : Electronic( isp, R, min_massf, theta ), m_( PC_R_u / R_ / PC_Avogadro )
{
    // Make the look-up-tables
    e_LUT_ = new NoneqCoupledDiatomicLUT( lut_fname, 5 );
    Cv_LUT_ = new NoneqCoupledDiatomicLUT( lut_fname, 9 );
}

Coupled_diatomic_electronic::
~Coupled_diatomic_electronic()
{
    delete e_LUT_;
    delete Cv_LUT_;
}

double
Coupled_diatomic_electronic::
s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    return e_LUT_->eval(T_el,T_vib,T_rot);
}

double
Coupled_diatomic_electronic::
s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    cout << "Coupled_diatomic_vibration::s_eval_entropy_from_Ts()" << endl
         << "No expression available." << endl;
    exit( FAILURE );
}

double
Coupled_diatomic_electronic::
s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot  )
{
    return Cv_LUT_->eval(T_el,T_vib,T_rot);
}
#endif


/* ------- Fully excited translation ------- */

Fully_excited_translation::
Fully_excited_translation( int isp, double R, double min_massf )
 : Species_energy_mode( isp, R, min_massf, "translation" ), Cv_( 1.5 * R ), 
   Cp_( 2.5 * R )
{
    // calculate constant expression for entropy calculation
    double m = PC_R_u/R_/PC_Avogadro;
    double tmp = pow( 2.0 * M_PI * m * PC_k_SI / PC_h_SI / PC_h_SI, 1.5 ) * PC_k_SI;
    entropy_constant_ = R_ * ( log( tmp ) + 2.5 );
}

double
Fully_excited_translation::
s_eval_energy_from_T( double T, double A )
{
    return Cv_*T;
}

double
Fully_excited_translation::
s_eval_enthalpy_from_T( double T, double A )
{
    return Cp_*T;
}

double
Fully_excited_translation::
s_eval_entropy( const Gas_data &Q )
{
    // Ref: Vincenti and Kruger (1975) Eq. 9.5b p 123
    // NOTE: - check use of SI units in 'tmp'
    //       - if 'p' is zero then returned value will be '-inf'
    //       - constant part of expression has been pre-calculated in constructor
    // double m = PC_R_u/R_/PC_Avogadro;
    // double tmp = pow( 2.0 * M_PI * m * PC_k_SI / PC_h_SI / PC_h_SI, 1.5 ) * PC_k_SI;
    double T = Q.T[iT_];
    double p = Q.p;
    return Cp_*log(T) - R_*log(p) + entropy_constant_;
}

double
Fully_excited_translation::
s_eval_entropy_from_T( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 9.5b p 123
    // NOTE: - check use of SI units in 'tmp'
    //       - if 'p' is zero then returned value will be '-inf'
    //       - constant part of expression has been pre-calculated in constructor
    // double m = PC_R_u/R_/PC_Avogadro;
    // double tmp = pow( 2.0 * M_PI * m * PC_k_SI / PC_h_SI / PC_h_SI, 1.5 ) * PC_k_SI;

    return Cp_*log(T) - R_*log(PC_P_atm) + entropy_constant_;
}

double
Fully_excited_translation::
s_eval_Cv_from_T( double T  )
{
    return Cv_;
}

double
Fully_excited_translation::
s_eval_Cp_from_T( double T  )
{
    return Cp_;
}

double
Fully_excited_translation::
s_eval_Q_from_T( double T, double A )
{
    // NOTE: must divide by avogadro's number to get units of 1/m3 (I think)
    return pow( 2.0 * M_PI * PC_R_u / R_ / PC_Avogadro * PC_k_SI * T / PC_h_SI / PC_h_SI, 1.5 ) / PC_Avogadro;
}

/* ------- Generic rotation --------- */

Rotation::
Rotation( int isp, double R, double min_massf, double theta )
: Species_energy_mode( isp, R, min_massf, "rotation", theta )
{}

/* ------- Fully excited rotation ------- */

Fully_excited_rotation::
Fully_excited_rotation( int isp, double R, double min_massf, double theta, int sigma )
 : Rotation( isp, R, min_massf, theta ), sigma_( sigma )
{
    Cv_ = R_;
}

double
Fully_excited_rotation::
s_eval_energy_from_T( double T, double A  )
{
    return Cv_*T;
}

double
Fully_excited_rotation::
s_eval_entropy_from_T( double T  )
{
    // Ref: Vincenti and Kruger (1975) Ex 12.1 p 138
    return R_ * ( log( T / double(sigma_) / theta_ ) + 1.0 );
}

double
Fully_excited_rotation::
s_eval_Cv_from_T( double T  )
{
    return Cv_;
}

double
Fully_excited_rotation::
s_eval_Q_from_T( double T, double A )
{
    return T / sigma_ / theta_ + 1.0 / ( 3.0 * sigma_);
}

/* ------- Fully excited nonlinear rotation ------- */

Fully_excited_nonlinear_rotation::
Fully_excited_nonlinear_rotation( int isp, double R, double min_massf, 
    				  double theta_A0, double theta_B0, double theta_C0,
    				  int sigma )
 : Rotation( isp, R, min_massf, theta_A0 ), theta_A0_( theta_A0 ), 
   theta_B0_( theta_B0 ), theta_C0_( theta_C0 ), sigma_( sigma )
{
    Cv_ = 1.5 * R_;
    UNUSED_VARIABLE(sigma_);
}

double
Fully_excited_nonlinear_rotation::
s_eval_energy_from_T( double T, double A  )
{
    return Cv_*T;
}

double
Fully_excited_nonlinear_rotation::
s_eval_entropy_from_T( double T  )
{
    // Ref: My derivation from Capitelli ESA STR 246 p21 Q_rot expression
    double PF = sqrt( T*T*T * M_PI / ( theta_A0_ * theta_B0_ * theta_C0_ ) );
    
    return R_ * ( log( PF ) + 1.5 );
}

double
Fully_excited_nonlinear_rotation::
s_eval_Cv_from_T( double T  )
{
    return Cv_;
}

double
Fully_excited_nonlinear_rotation::
s_eval_Q_from_T( double T, double A )
{
    return  sqrt( T*T*T * M_PI / ( theta_A0_ * theta_B0_ * theta_C0_ ) );
}

/* ------- Coupled diatomic rotation ------- */

#if TABULATED_COUPLED_DIATOMIC_MODES==0
Coupled_diatomic_rotation::
Coupled_diatomic_rotation( int isp, double R, double min_massf, int sigma_r,
    double fT, std::vector<Diatom_electronic_level*> &elevs )
 : Rotation( isp, R, min_massf, elevs[0]->B_e / PC_k_SI ), sigma_r_( sigma_r ),
 m_( PC_R_u / R_ / PC_Avogadro ), fT_( fT )
{
    // Set the elev pointers
    for ( size_t ilev=0; ilev < elevs.size(); ++ilev ) {
    	elevs_.push_back( elevs[ilev] );
    }
    UNUSED_VARIABLE(sigma_r_);
}

Coupled_diatomic_rotation::
~Coupled_diatomic_rotation()
{
    // elevs are just pointers now
}

double
Coupled_diatomic_rotation::
eval_total_internal_energy( double T_el, double T_vib, double T_rot  )
{
    // Calculate TOTAL (all modes) internal energy
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
                double E_rot_dash = E_el + E_vib + E_rot;
                double Q_rot_dash = Q_el*Q_vib*Q_rot;
    	    	if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
                E_weighted += E_rot_dash * Q_rot_dash;
                Q_total += Q_rot_dash;
            }
        }
    }

    return E_weighted / Q_total / m_;   // convert J/particle -> J/kg
}

double
Coupled_diatomic_rotation::
s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    // Loop over all rovibronic levels, sum up contributions from definition
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        double Q_vib_sum = 0.0;
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            double Q_rot_sum = 0.0;
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
       	    	if ( std::isnan(Q_rot) || std::isinf( Q_rot ) ) continue;
                Q_rot_sum += Q_rot;
                E_weighted += E_rot * Q_el * Q_vib * Q_rot;
            }
            Q_vib_sum += Q_vib*Q_rot_sum;
        }
        Q_total += Q_el * Q_vib_sum;
    }

    return E_weighted / Q_total / m_;   // convert J/particle -> J/kg
}

double
Coupled_diatomic_rotation::
s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot  )
{
#   if 0
    // Loop over all rovibronic levels, sum up contributions from definition
    // NOTE: not sure that this formulation is correct
    double E_weighted = 0.0;
    double Q_total = 0.0;
    double Q_rot_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        double Q_vib_sum = 0.0;
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            double Q_rot_sum = 0.0;
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
       	    	if ( isnan(Q_rot) || isinf( Q_rot ) ) continue;
                Q_rot_sum += Q_rot;
                E_weighted += E_rot * Q_el * Q_vib * Q_rot;
            }
            Q_vib_sum += Q_vib*Q_rot_sum;
            Q_rot_total += Q_rot_sum;
        }
        Q_total += Q_el * Q_vib_sum;
    }

    double e = E_weighted / Q_total / m_;   // convert J/particle -> J/kg
    
    return R_*log(Q_rot_total) + e/T_rot;		// units of J/kg/K
#   else
    cout << "Coupled_diatomic_rotation::s_eval_entropy_from_Ts()" << endl
         << "No expression available." << endl;
    exit( FAILURE );
#   endif
}

double
Coupled_diatomic_rotation::
s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot  )
{
    // Step size
    double h = T_rot * fT_;

    // five-point stencil method
    double fxp2h = eval_total_internal_energy(T_el,T_vib,T_rot + 2.0 * h);
    double fxph  = eval_total_internal_energy(T_el,T_vib,T_rot + h);
    double fxmh  = eval_total_internal_energy(T_el,T_vib,T_rot - h);
    double fxm2h = eval_total_internal_energy(T_el,T_vib,T_rot - 2.0 * h);

    return ( - fxp2h + 8.0 * fxph - 8.0 * fxmh + fxm2h ) / 12.0 / h;
}
#else
Coupled_diatomic_rotation::
Coupled_diatomic_rotation( int isp, double R, double min_massf, double theta_r, string lut_fname )
 : Rotation( isp, R, min_massf, theta_r ), m_( PC_R_u / R_ / PC_Avogadro )
{
    // Make the look-up-tables
    e_LUT_ = new NoneqCoupledDiatomicLUT( lut_fname, 7 );
    Cv_LUT_ = new NoneqCoupledDiatomicLUT( lut_fname, 11 );
    
}

Coupled_diatomic_rotation::
~Coupled_diatomic_rotation()
{
    delete e_LUT_;
    delete Cv_LUT_;
}

double
Coupled_diatomic_rotation::
s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    return e_LUT_->eval(T_el,T_vib,T_rot);
}

double
Coupled_diatomic_rotation::
s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    cout << "Coupled_diatomic_vibration::s_eval_entropy_from_Ts()" << endl
         << "No expression available." << endl;
    exit( FAILURE );
}

double
Coupled_diatomic_rotation::
s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot  )
{
    return Cv_LUT_->eval(T_el,T_vib,T_rot);
}
#endif

/* ------- Generic vibration ------- */

Vibration::
Vibration(  int isp, double R, double min_massf, double theta )
: Species_energy_mode( isp, R, min_massf, "vibration", theta ) {}

/* ------- Un-truncated harmonic vibration ------- */

Harmonic_vibration::
Harmonic_vibration(  int isp, double R, double min_massf, double theta )
: Vibration( isp, R, min_massf, theta ) {}
double
Harmonic_vibration::
s_eval_energy_from_T( double T, double A  )
{
    // Ref: Vincenti and Kruger (1975) Eq. 12.12 p 135
    return R_ * theta_ / ( exp( theta_ / T ) - 1.0 );
}

double
Harmonic_vibration::
s_eval_entropy_from_T( double T  )
{
    // Ref: Vincenti and Kruger (1975) Ex 12.2 p 138
    return R_ * ( - log ( 1 - exp( - theta_ / T ) ) + ( theta_ / T ) / ( exp( theta_ / T ) - 1.0 ) );
}

double
Harmonic_vibration::
s_eval_Cv_from_T( double T  )
{
    // Ref: Vincenti and Kruger (1975) Eq 12.13 p 135
    double theta_2T = theta_ / ( 2.0 * T );
    return R_ * pow( theta_2T / sinh( theta_2T ), 2 );
}

double
Harmonic_vibration::
s_eval_Q_from_T( double T, double A )
{
    // Ref: Vincenti and Kruger (1975) p 135
    return 1.0 / ( 1.0 - exp( - theta_ / T ) );
}

/* ------- Anharmonic vibration ------- */

Anharmonic_vibration::
Anharmonic_vibration(  int isp, double R, double min_massf, double theta, 
    double thetaR, int sigma,
    Segmented_functor* h, Segmented_functor* s, Segmented_functor* Cp )
: Vibration( isp, R, min_massf, theta ),h_(h),s_(s),Cp_(Cp),thetaR_(thetaR),sigma_(sigma) {}

double
Anharmonic_vibration::
s_eval_energy_from_T( double Tv, double A)
{
    return (*h_)(Tv) - 3.5*R_*Tv;
}

double 
Anharmonic_vibration::
s_eval_entropy( const Gas_data &Q)
{
    Tt_ = Q.T[0];
    p_ = Q.p;
    return s_eval_entropy_from_T(Q.T[iT_]);
}

double
Anharmonic_vibration::
s_eval_entropy_from_T( double Tv )
{
    // Evaluate entropy at Tt for translational mode
    double m = PC_R_u/R_/PC_Avogadro;
    double tmp = pow( 2.0 * M_PI * m * PC_k_SI / PC_h_SI / PC_h_SI, 1.5 ) * PC_k_SI;
    double entropy_constant_ = R_ * ( log( tmp ) + 2.5 );
    double s_t =  2.5*R_*log(Tt_) - R_*log(p_) + entropy_constant_;
    // Evaluate entropy at Tv for rotational mode
    double s_rot =  R_ * ( log( Tv / double(sigma_) / thetaR_ ) + 1.0 );
    
    return (*s_)(Tv) - s_t - s_rot;
}

double
Anharmonic_vibration::
s_eval_Cv_from_T( double T  )
{
    // Ref: Vincenti and Kruger (1975) Eq 12.13 p 135
    return (*Cp_)(T)-3.5*R_;
}

double
Anharmonic_vibration::
s_eval_Q_from_T( double T, double A )
{
    // Ref: Vincenti and Kruger (1975) p 135
    return 1.0 / ( 1.0 - exp( - theta_ / T ) );
}


/* ------- Truncated harmonic vibration ------- */

Truncated_harmonic_vibration::
Truncated_harmonic_vibration(  int isp, double R, double min_massf, double theta, double theta_D )
: Vibration( isp, R, min_massf, theta ), theta_D_( theta_D ) {}

double
Truncated_harmonic_vibration::
s_eval_energy_from_T( double T, double A  )
{
    // Ref: Gollan, R.G. PhD Thesis June 2008
    // NOTE: 'A', if positive, is the truncation temperature
    double theta_lim = A;
    if ( theta_lim < 0.0 ) theta_lim = theta_D_;
    double exp_theta_D_T = exp( theta_lim / T );
    if ( std::isnan(exp_theta_D_T) ) {
    	// use harmonic oscillator solution
    	return R_ * ( theta_ / ( exp( theta_ / T ) - 1.0 ) );
    }
    return R_ * ( theta_ / ( exp( theta_ / T ) - 1.0 ) - theta_lim / ( exp_theta_D_T - 1.0 ) );
}

double
Truncated_harmonic_vibration::
s_eval_entropy_from_T( double T  )
{
    // Ref: My own derivation from the definition of entropy as a function of
    //      the partition function and temperature (see Eq 3.2 in PhD thesis)
    double exp_theta_v_T = exp( - theta_ / T );
    double exp_theta_D_T = exp( - theta_D_ / T );
    
    return R_ * ( - log ( ( 1.0 - exp_theta_v_T ) / ( 1.0 - exp_theta_D_T ) ) + \
    		 ( 1.0 / T ) * ( ( theta_ * exp_theta_v_T ) / ( 1.0 - exp_theta_v_T ) - \
    		     		   theta_D_ * exp_theta_D_T / ( 1.0 - exp_theta_D_T ) ) );
}

double
Truncated_harmonic_vibration::
s_eval_Cv_from_T( double T  )
{
    // Ref: My own derivation from the definition of Cv as the derivative of
    //      energy w.r.t temperature
    double tmpA = theta_ / T;
    double tmpB = theta_D_ / T;
    double exp_tmpA = exp( tmpA );
    double exp_tmpB = exp( tmpB );
    double tmpC = tmpA * tmpA * exp_tmpA /  ( (exp_tmpA - 1.0) * (exp_tmpA - 1.0) );
    double tmpD = tmpB * tmpB * exp_tmpB /  ( (exp_tmpB - 1.0) * (exp_tmpB - 1.0) );
    
    // if tmpD is nan then it is converging to zero
    if ( std::isnan(tmpD) ) tmpD = 0.0;
    
    return R_ * ( tmpC - tmpD );
}

double
Truncated_harmonic_vibration::
s_eval_Q_from_T( double T, double A )
{
    // NOTE: 'A', if positive, is the truncation temperature
    double theta_lim = A;
    if ( theta_lim < 0.0 ) theta_lim = theta_D_;
    // Ref: DFP thesis
    return ( 1.0 - exp( - theta_lim / T ) ) / ( 1.0 - exp( - theta_ / T ) );
}

double
Truncated_harmonic_vibration::
s_eval_HO_energy( double T )
{
    // Ref: Vincenti and Kruger (1975) Eq. 12.12 p 135
    return R_ * theta_ / ( exp( theta_ / T ) - 1.0 );
}

/* ------- Coupled diatomic vibration mode ------- */

#if TABULATED_COUPLED_DIATOMIC_MODES==0
Coupled_diatomic_vibration::
Coupled_diatomic_vibration( int isp, double R, double min_massf, int sigma_r,
    double fT, std::vector<Diatom_electronic_level*> &elevs )
 : Vibration( isp, R, min_massf, elevs[0]->omega_e / PC_k_SI ), sigma_r_( sigma_r ),
 m_( PC_R_u / R_ / PC_Avogadro ), fT_( fT )
{
    // Set the elev pointers
    for ( size_t ilev=0; ilev < elevs.size(); ++ilev ) {
    	elevs_.push_back( elevs[ilev] );
    }
    UNUSED_VARIABLE(sigma_r_);
}

Coupled_diatomic_vibration::
~Coupled_diatomic_vibration()
{
    // elevs are just pointers now
}

double
Coupled_diatomic_vibration::
eval_total_internal_energy( double T_el, double T_vib, double T_rot  )
{
    // Calculate TOTAL (all modes) internal energy
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
                double E_rot_dash = E_el + E_vib + E_rot;
                double Q_rot_dash = Q_el*Q_vib*Q_rot;
    	    	if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
                E_weighted += E_rot_dash * Q_rot_dash;
                Q_total += Q_rot_dash;
            }
        }
    }

    return E_weighted / Q_total / m_;   // convert J/particle -> J/kg
}

double
Coupled_diatomic_vibration::
s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    // Loop over all rovibronic levels, sum up contributions from definition
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        double Q_vib_sum = 0.0;
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            double Q_rot_sum = 0.0;
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
       	    	if ( std::isnan(Q_rot) || std::isinf( Q_rot ) ) continue;
                Q_rot_sum += Q_rot;
            }
            Q_vib_sum += Q_vib*Q_rot_sum;
            E_weighted += E_vib * Q_el * Q_vib * Q_rot_sum;
        }
        Q_total += Q_el * Q_vib_sum;
    }
    
    return E_weighted / Q_total / m_;   // convert J/particle -> J/kg
}

double
Coupled_diatomic_vibration::
s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot  )
{
#   if 0
    // Loop over all rovibronic levels, sum up contributions from definition
    // NOTE: not sure that this formulation is correct
    double E_weighted = 0.0;
    double Q_total = 0.0;
    double Q_vib_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Diatom_electronic_level * elev = elevs_[ilev];
        double E_el = elev->E;
        double Q_el = elev->g * exp( - E_el / PC_k_SI / T_el );
        double Q_vib_sum = 0.0;
        for ( int iV=0; iV<elev->V_max; ++iV ) {
            double E_rot_0 = elev->eval_E_rot(iV,0);
            double E_vib = elev->eval_E_vib( iV ) + E_rot_0;
            // cout << "ilev = " << ilev << ", iV = " << iV << ", E_el_0 = " << E_el_0 << ", E_vib_0 = " << E_vib_0 << ", E_rot_0 = " << E_rot_0 << endl;
            double Q_vib = exp( - E_vib / PC_k_SI / T_vib );
            double Q_rot_sum = 0.0;
            for ( int iJ=0; iJ<elev->J_max[iV]; ++iJ ){
                double E_rot = elev->eval_E_rot(iV,iJ) - E_rot_0;
                double Q_rot = double(2*iJ+1) * exp( - E_rot / PC_k_SI / T_rot );
       	    	if ( isnan(Q_rot) || isinf( Q_rot ) ) continue;
                Q_rot_sum += Q_rot;
            }
            Q_vib_sum += Q_vib*Q_rot_sum;
            E_weighted += E_vib * Q_el * Q_vib * Q_rot_sum;
            Q_vib_total += Q_vib;
        }
        Q_total += Q_el * Q_vib_sum;
    }

    double e = E_weighted / Q_total / m_;   // convert J/particle -> J/kg
    
    return R_*log(Q_vib_total) + e/T_vib;		// units of J/kg/K
#   else
    cout << "Coupled_diatomic_vibration::s_eval_entropy_from_Ts()" << endl
         << "No expression available." << endl;
    exit( FAILURE );
#   endif
}

double
Coupled_diatomic_vibration::
s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot  )
{
    // Step size
    double h = T_vib * fT_;

    // five-point stencil method
    double fxp2h = eval_total_internal_energy(T_el,T_vib + 2.0 * h,T_rot);
    double fxph  = eval_total_internal_energy(T_el,T_vib + h,T_rot);
    double fxmh  = eval_total_internal_energy(T_el,T_vib - h,T_rot);
    double fxm2h = eval_total_internal_energy(T_el,T_vib - 2.0 * h,T_rot);

    return ( - fxp2h + 8.0 * fxph - 8.0 * fxmh + fxm2h ) / 12.0 / h;
}
#else
Coupled_diatomic_vibration::
Coupled_diatomic_vibration( int isp, double R, double min_massf, double theta_v, string lut_fname )
 : Vibration( isp, R, min_massf, theta_v ), m_( PC_R_u / R_ / PC_Avogadro )
{
    // Make the look-up-tables
    e_LUT_ = new NoneqCoupledDiatomicLUT( lut_fname, 6 );
    Cv_LUT_ = new NoneqCoupledDiatomicLUT( lut_fname, 10 );
    
}

Coupled_diatomic_vibration::
~Coupled_diatomic_vibration()
{
    delete e_LUT_;
    delete Cv_LUT_;
}

double
Coupled_diatomic_vibration::
s_eval_energy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    return e_LUT_->eval(T_el,T_vib,T_rot);
}

double
Coupled_diatomic_vibration::
s_eval_entropy_from_Ts( double T_el, double T_vib, double T_rot  )
{
    cout << "Coupled_diatomic_vibration::s_eval_entropy_from_Ts()" << endl
         << "No expression available." << endl;
    exit( FAILURE );
}

double
Coupled_diatomic_vibration::
s_eval_Cv_from_Ts( double T_el, double T_vib, double T_rot  )
{
    return Cv_LUT_->eval(T_el,T_vib,T_rot);
}
#endif

/* ------- Fully coupled diatomic internal mode ------- */

#if TABULATED_COUPLED_DIATOMIC_MODES==0
Fully_coupled_diatom_internal::
Fully_coupled_diatom_internal( int isp, double R, double min_massf, int sigma_r,
    double fT, std::vector<Diatom_electronic_level*> &elevs )
 : Species_energy_mode( isp, R, min_massf, "internal" ), sigma_r_( sigma_r ),
 m_( PC_R_u / R_ / PC_Avogadro ), fT_( fT )
{
    // Set the elev pointers
    for ( size_t ilev=0; ilev < elevs.size(); ++ilev ) {
    	elevs_.push_back( elevs[ilev] );
    }
    UNUSED_VARIABLE(sigma_r_);
    UNUSED_VARIABLE(fT_);
}

Fully_coupled_diatom_internal::
~Fully_coupled_diatom_internal()
{
    // elevs are just pointers
}

Diatom_electronic_level * Fully_coupled_diatom_internal::get_elev_pointer( int ilev )
{
    if ( ilev >= int(elevs_.size()) ) {
        cout << "Fully_coupled_diatom_internal::get_elev_pointer()" << endl
             << "Requested electronic level does not exist." << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }

    return elevs_[ilev];
}

double
Fully_coupled_diatom_internal::
s_eval_energy_from_T( double T, double A  )
{
    // Loop over all rovibronic levels, sum up contributions from definition
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
    	double E_el = elevs_[ilev]->E;
    	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
            double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
    	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) + E_rot_0;
    	    for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ){
    	    	double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
    	    	double E_rot_dash = E_el + E_vib + E_rot;
    	    	double Q_rot_dash = elevs_[ilev]->g * double(2*iJ+1) * exp( - E_rot_dash / PC_k_SI / T );
    	    	if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
    	    	E_weighted += E_rot_dash * Q_rot_dash;
    	    	Q_total += Q_rot_dash;
    	    }
    	}
    }

    return E_weighted / Q_total / m_;	// convert J/particle -> J/kg
}

double
Fully_coupled_diatom_internal::
s_eval_entropy_from_T( double T  )
{
    // Loop over all rovibronic levels, sum up contributions from definition:
    // S = N k ln(Q) + E/T; s = R ln(Q) + e/T
    double E_weighted = 0.0;
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
    	double E_el = elevs_[ilev]->E;
    	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
            double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
    	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) + E_rot_0;
    	    for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ){
    	    	double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
    	    	double E_rot_dash = E_el + E_vib + E_rot;
    	    	double Q_rot_dash = elevs_[ilev]->g * double(2*iJ+1) * exp( - E_rot_dash / PC_k_SI / T );
    	    	if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
    	    	E_weighted += E_rot_dash * Q_rot_dash;
    	    	Q_total += Q_rot_dash;
    	    }
    	}
    }
    
    double e = E_weighted / Q_total / m_;	// convert J/particle -> J/kg
    
    return R_*log(Q_total) + e/T;		// units of J/kg/K
}

double
Fully_coupled_diatom_internal::
s_eval_Cv_from_T( double T  )
{
    // differentiation of energy wrt T using the quotient rule
    double u = 0.0, u_dash = 0.0;
    double v = 0.0, v_dash = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
    	double E_el = elevs_[ilev]->E;
    	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
            double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
    	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) + E_rot_0;
            for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ) {
    	    	double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
    	    	double E_rot_dash = E_el + E_vib + E_rot;
    	    	double g = elevs_[ilev]->g * double(2*iJ+1);
    	    	double Q_rot_dash = g * exp( - E_rot_dash / PC_k_SI / T );
    	    	if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
    	    	u += E_rot_dash * Q_rot_dash;
    	    	u_dash += E_rot_dash / PC_k_SI / T / T * E_rot_dash * Q_rot_dash;
    	    	v += Q_rot_dash;
    	    	v_dash += E_rot_dash / PC_k_SI / T / T * Q_rot_dash;
    	    }
    	}
    }
    
    return ( u_dash * v - u * v_dash ) / ( v * v ) / m_;	// convert J/particle/K -> J/kg/K
}

double
Fully_coupled_diatom_internal::
s_eval_Q_from_T( double T, double A )
{
    // Loop over all rovibronic levels, sum up contributions from definition
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
	double E_el = elevs_[ilev]->E;
	for ( int iV=0; iV<elevs_[ilev]->V_max; ++iV ) {
	    double E_rot_0 = elevs_[ilev]->eval_E_rot(iV,0);
	    double E_vib = elevs_[ilev]->eval_E_vib( iV ) + E_rot_0;
	    for ( int iJ=0; iJ<elevs_[ilev]->J_max[iV]; ++iJ ){
		double E_rot = elevs_[ilev]->eval_E_rot(iV,iJ) - E_rot_0;
		double E_rot_dash = E_el + E_vib + E_rot;
		double Q_rot_dash = elevs_[ilev]->g * double(2*iJ+1) * exp( - E_rot_dash / PC_k_SI / T );
		if ( std::isnan(Q_rot_dash) || std::isinf( Q_rot_dash ) ) continue;
		Q_total += Q_rot_dash;
	    }
	}
    }

    return Q_total;
}

#else
Fully_coupled_diatom_internal::
Fully_coupled_diatom_internal( int isp, double R, double min_massf, string lut_fname )
 : Species_energy_mode( isp, R, min_massf, "internal" ), m_( PC_R_u / R_ / PC_Avogadro )
{
    // Make the look-up-tables
    e_LUT_ = new EqCoupledDiatomicLUT( lut_fname, 4 );
    Cv_LUT_ = new EqCoupledDiatomicLUT( lut_fname, 8 );
    s_LUT_ = new EqCoupledDiatomicLUT( lut_fname, 12 );
}

Fully_coupled_diatom_internal::
~Fully_coupled_diatom_internal()
{
    delete e_LUT_;
    delete Cv_LUT_;
    delete s_LUT_;
}

double
Fully_coupled_diatom_internal::
s_eval_energy_from_T( double T, double A  )
{
    return e_LUT_->eval(T);
}

double
Fully_coupled_diatom_internal::
s_eval_entropy_from_T( double T  )
{
    return s_LUT_->eval(T);
}

double
Fully_coupled_diatom_internal::
s_eval_Cv_from_T( double T  )
{
    return Cv_LUT_->eval(T);
}
#endif

/* ------- Fully coupled polyatomic internal mode ------- */

Fully_coupled_polyatom_internal::
Fully_coupled_polyatom_internal( int isp, double R, double min_massf,
    double fT, std::vector<Polyatom_electronic_level*> &elevs )
 : Species_energy_mode( isp, R, min_massf, "internal" ),
 m_( PC_R_u / R_ / PC_Avogadro ), fT_( fT )
{
    // Set the elev pointers
    for ( size_t ilev=0; ilev < elevs.size(); ++ilev ) {
        elevs_.push_back( elevs[ilev] );
    }
    UNUSED_VARIABLE(fT_);
    UNUSED_VARIABLE(m_);
}

Fully_coupled_polyatom_internal::
~Fully_coupled_polyatom_internal()
{
    // elevs are just pointers
}

Polyatom_electronic_level * Fully_coupled_polyatom_internal::get_elev_pointer( int ilev )
{
    if ( ilev >= int(elevs_.size()) ) {
        cout << "Fully_coupled_polyatom_internal::get_elev_pointer()" << endl
             << "Requested electronic level does not exist." << endl
             << "Exiting program." << endl;
        exit( FAILURE );
    }

    return elevs_[ilev];
}

double
Fully_coupled_polyatom_internal::
s_eval_energy_from_T( double T, double A  )
{
    return 0.0;
}

double
Fully_coupled_polyatom_internal::
s_eval_entropy_from_T( double T  )
{
    return 0.0;
}

double
Fully_coupled_polyatom_internal::
s_eval_Cv_from_T( double T  )
{
    return 0.0;
}

double
Fully_coupled_polyatom_internal::
s_eval_Q_from_T( double T, double A )
{
    double Q_total = 0.0;
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev ) {
        Q_total += elevs_[ilev]->eval_Q_from_T(T);
    }

    return Q_total;
}
