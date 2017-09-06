// Author: Daniel F. Potter
// Date: 18-Nov-2009


#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "../models/gas_data.hh"
#include "../models/physical_constants.hh"
#include "../models/chemical-species.hh"
#include "../models/chemical-species-library.hh"
#include "../models/species-energy-modes.hh"

#include "energy-exchange-relaxation-time.hh"

using namespace std;

Relaxation_time::
Relaxation_time() {}

Relaxation_time::
~Relaxation_time() {}

VT_MillikanWhite::
VT_MillikanWhite(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{

    Chemical_species *p = get_library_species_pointer(ip);
    Chemical_species *q = get_library_species_pointer(iq);

    M_p_ = p->get_M();
    M_q_ = q->get_M();
    if ( p->get_type().find("diatomic")!=string::npos ) {
        Diatomic_species *P = dynamic_cast<Diatomic_species*>(p);
        theta_ = P->get_theta_v();
    }
    else if ( p->get_type().find("polyatomic")!=string::npos ) {
        Polyatomic_species *P = dynamic_cast<Polyatomic_species*>(p);
        theta_ = P->get_theta_v();
    }
    else {
        cout << "VT_MillikanWhite::VT_MillikanWhite()" << endl
             << "Species: " << p->get_name() << " is not declared as a molecule!" << endl;
        exit( BAD_INPUT_ERROR );
    }

    mu_ = ((M_p_ * M_q_) / (M_p_ + M_q_))*1000.0;

    lua_getfield(L, -1, "a");
    if ( lua_isnil(L, -1) ) {
	// No 'a' value supplied. Compute according to Millikan and White.
	a_ = 1.16e-3*sqrt(mu_)*pow(theta_, 4.0/3.0);
    }
    else {
	// Take directly from table.
	a_ = luaL_checknumber(L, -1);
    }
    lua_pop(L, 1);

    lua_getfield(L, -1, "b");
    if ( lua_isnil(L, -1) ) {
	// No 'b' value supplied. Compute according to Millikan and White.
	b_ = 0.015*pow(mu_, 1.0/4.0);
    }
    else {
	// Take directly from table.
	b_ = luaL_checknumber(L, -1);
    }
    lua_pop(L, 1);

    UNUSED_VARIABLE(ip_);
}

VT_MillikanWhite::
~VT_MillikanWhite() {}

double
VT_MillikanWhite::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    // When the 'bath' pressure is very small,
    // practically no relaxation occurs.
    if ( molef[iq_] <= DEFAULT_MIN_MASS_FRACTION ) {
	return -1.0;
    }
    double p_bath;
    double tau;
    // Set the 'bath' pressure as that of the 'q' colliders
    // and compute in atm for use in Millikan-White expression
    p_bath = molef[iq_]*Q.p/PC_P_atm;
    tau = (1.0/p_bath)*exp(a_*(pow(Q.T[iT_], -1.0/3.0) - b_) - 18.42);
    //    cout << "tau= " << tau << endl;
    return tau;
}

VT_PolyFit::
VT_PolyFit(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{

    Chemical_species *p = get_library_species_pointer(ip);
    Chemical_species *q = get_library_species_pointer(iq);

    M_p_ = p->get_M();
    M_q_ = q->get_M();
    if ( p->get_type().find("diatomic")!=string::npos ) {
        Diatomic_species *P = dynamic_cast<Diatomic_species*>(p);
        theta_ = P->get_theta_v();
    }
    else if ( p->get_type().find("polyatomic")!=string::npos ) {
        Polyatomic_species *P = dynamic_cast<Polyatomic_species*>(p);
        theta_ = P->get_theta_v();
    }
    else {
        cout << "VT_Polyfit::VT_PolyFit()" << endl
             << "Species: " << p->get_name() << " is not declared as a molecule!" << endl;
        exit( BAD_INPUT_ERROR );
    }

    lua_getfield(L, -1, "thetav");
    if ( ! lua_isnil(L, -1) ) {
       // Take directly from table.
       theta_ = luaL_checknumber(L, -1);
    }
    lua_pop(L, 1);

    mu_ = ((M_p_ * M_q_) / (M_p_ + M_q_))*1000.0;

    type_ = get_number(L,-1,"type");
    if (type_ > 2 || type_ < 0){
       ostringstream ost;
       ost << "VT_PolyFit::VT_PolyFit()\n";
       ost << "Incorrect value for 'type', it should be between 0 and 2\n";
       input_error(ost);
    }

    /* type = 0:
     *     p*tau = a * sqrt(T) / (1-exp(thetav/T))*exp(b/T^(-1/3))
     * type = 1:
     *     polybase < 0:
     *           p*tau = scale*exp(polyval[coeffs,T/Tnorm^(T_pow)])
     *     polybase > 0:
     *           p*tau = scale*pow(polybase,polyval[coeffs,T/Tnorm^(T_pow)])
     * type = 2:
     *     T < Thigh
     *        use type 1
     *     T > Thigh
     *        use type 0*/

    if (type_ == 2) {
       T_high_ = get_positive_number(L,-1,"T_high");
    }

    if ( type_ != 1 ) {
       a_ = get_number(L,-1,"a");
       b_ = get_number(L,-1,"b");
    }

    if ( type_ != 0 ){
       polybase_ = get_number(L,-1,"polybase");

       lua_getfield(L, -1, "coeff");
       if ( !lua_istable(L, -1) ) {
          ostringstream ost;
          ost << "VT_PolyFit::VT_PolyFit()\n";
          ost << "A 'coeffs' table is expected with the 'parameters' table, but it's not found.\n";
          input_error(ost);
       }

       n_poly_ = lua_objlen(L, -1);
       if ( n_poly_ < 1 ) {
          cout << "VT_PolyFit::VT_PolyFit()\n";
          cout << " At least one coeffs should be provided. \n";
          exit(BAD_INPUT_ERROR);
       }

       c_.resize(n_poly_);
       for ( size_t i = 0; i < c_.size(); ++i ) {
          lua_rawgeti(L, -1, i+1);
          c_[i] = luaL_checknumber(L, -1);
          lua_pop(L, 1);
       }
       lua_pop(L,1); //pop coeffs

       lua_getfield(L, -1, "T_norm");
       if ( lua_isnil(L, -1) ) {
          // No 'T_norm' value supplied. Compute according to Millikan and White.
          T_norm_ = 1000.0;
       }
       else {
          // Take directly from table.
          T_norm_ = luaL_checknumber(L, -1);
       }
       lua_pop(L, 1);

       lua_getfield(L, -1, "T_pow");
       if ( lua_isnil(L, -1) ) {
          // No 'T_pow' value supplied. Compute according to Millikan and White.
          T_pow_ = 1.0;
       }
       else {
          // Take directly from table.
          T_pow_ = luaL_checknumber(L, -1);
       }
       lua_pop(L, 1);

       lua_getfield(L, -1, "scale");
       if ( lua_isnil(L, -1) ) {
          // No 'T_pow' value supplied. Compute according to Millikan and White.
          prescale_ = 1.0;
       }
       else {
          // Take directly from table.
          prescale_ = luaL_checknumber(L, -1);
       }
       lua_pop(L, 1);
    }
}

VT_PolyFit::
~VT_PolyFit() {}

double
VT_PolyFit::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    // When the 'bath' pressure is very small,
    // practically no relaxation occurs.
    if ( molef[iq_] <= DEFAULT_MIN_MASS_FRACTION ) {
	return -1.0;
    }
    double ptau;
    double T = Q.T[iT_];
    // Set the 'bath' pressure as that of the 'q' colliders
    // and compute in atm for use in Millikan-White expression
    double p_bath = molef[iq_]*Q.p/PC_P_atm;

    if ( type_ == 0 || (type_ == 2 && T <= T_high_))
       ptau = a_ * sqrt(T) / (1.0 - exp(-theta_/ T))*exp(b_*pow(T,-1.0/3.0));
    else if (type_ == 1 || (type_ == 2 && T > T_high_)){
       ptau = 0.0;
       double T2 = pow(T/T_norm_,T_pow_);
       for (size_t i = 0; i < n_poly_; i++){
          ptau += c_[i]*pow(T2,n_poly_-i-1);
       }
       if (polybase_ < 0.0) ptau =  prescale_ * exp(ptau);
       else  ptau = prescale_*pow(polybase_,ptau);
    }
    //    cout << "tau= " << tau << endl;
    return ptau/p_bath;
}


VT_LandauTeller_cf::
VT_LandauTeller_cf(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    A_ = get_number(L, -1, "A");
    B_ = get_number(L, -1, "B");
    C_ = get_number(L, -1, "C");
    UNUSED_VARIABLE(ip_);
}

VT_LandauTeller_cf::
~VT_LandauTeller_cf() {}

double
VT_LandauTeller_cf::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double p_bath;
    double tau;
    double T = Q.T[iT_];
    // Set the 'bath' pressure as that of the 'q' colliders
    // and compute in atm for use in Millikan-White expression
    p_bath = molef[iq_]*Q.p/PC_P_atm;
    //p_bath = Q.p/PC_P_atm;
    tau = (A_/p_bath) * exp((B_/pow(T, 1.0/3.0)) + C_);
    return tau;
}

VT_Thivet_cf::
VT_Thivet_cf(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    B_ = get_number(L, -1, "B");
    C_ = get_number(L, -1, "C");
    UNUSED_VARIABLE(ip_);
}

VT_Thivet_cf::
~VT_Thivet_cf() {}

double
VT_Thivet_cf::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double p_bath;
    double tau;
    double T = Q.T[iT_];
    if ( T > 5.0e4 ) T = 5.0e4;
    // Set the 'bath' pressure as that of the 'q' colliders
    // and compute in atm for use in Millikan-White expression
    p_bath = molef[iq_]*Q.p/PC_P_atm;
    tau = (1.0/p_bath) * pow(T, 2.0/3.0) * exp((B_/pow(T, 1.0/3.0)) - C_);
    return tau;
}

VT_high_temperature_cross_section *
create_VT_high_temperature_cross_section_model(lua_State *L)
{
    string type = get_string(L, -1, "type");
    if ( type=="Park" ) return new Park_VT_HTCS(L);
    else if ( type=="Fujita" ) return new Fujita_VT_HTCS();
    else {
	cout << "create_VT_high_temperature_correction_model():\n";
	cout << "HTC model type: " << type << " was not recognised.\n";
	cout << "Bailing out!\n";
	exit( BAD_INPUT_ERROR );
    }
}

Park_VT_HTCS::Park_VT_HTCS(lua_State *L)
{
    sigma_dash_ = get_positive_number(L, -1, "sigma_dash" );
}

Park_VT_HTCS::~Park_VT_HTCS() {}

double
Park_VT_HTCS::s_eval_sigma(double T)
{
    // if ( T > 5.0e4 ) T = 5.0e4;
    return sigma_dash_ * ( 50000.0 / T ) * ( 50000.0 / T );		// cm**2
}

Fujita_VT_HTCS::Fujita_VT_HTCS() {}

Fujita_VT_HTCS::~Fujita_VT_HTCS() {}

double
Fujita_VT_HTCS::s_eval_sigma( double T )
{
    return ( 1.8e-14 ) / sqrt( T ) + ( 1.0e-4 ) / pow( T, 3.0 );	// cm**2
}


VT_MillikanWhite_HTC::
VT_MillikanWhite_HTC(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    // 1. Set Millikan-White relaxation time model
    VT_MW_ = new VT_MillikanWhite(L, ip, iq, itrans);
    Chemical_species *p = get_library_species_pointer(ip);
    double M_p = p->get_M();
    if ( ip == iq ) {
	m_av_ = M_p / PC_Avogadro * 1.0e3;
    }
    else {
	Chemical_species *q = get_library_species_pointer(iq);
	double M_q = q->get_M();
	m_av_ = 0.5*(M_p + M_q)/PC_Avogadro*1.0e3;
    }
    // 2. Pick up the high-temperature correction cross-section model
    lua_getfield(L,-1,"HTCS");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VT_MillikanWhite_HTC::VT_MillikanWhite_HTCS():\n";
	ost << "Error in the declaration of HTCS model: a table is expected.\n";
	input_error(ost);
    }
    HTC_model_ = create_VT_high_temperature_cross_section_model(L);
    lua_pop(L, 1);	// pop HTCS_model

    UNUSED_VARIABLE(ip_);
}

VT_MillikanWhite_HTC::
~VT_MillikanWhite_HTC()
{
    delete HTC_model_;
}

double
VT_MillikanWhite_HTC::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double p_bath = molef[iq_]*Q.p/PC_P_atm;
    double n_col = p_bath * PC_P_atm / PC_k_SI / Q.T[iT_] * 1.0e-6;
    double tau = VT_MW_->compute_relaxation_time(Q, molef);
    // high-temperature correction
    double sigma = HTC_model_->eval_sigma(Q.T[iT_]);
    tau += 1.0 / (n_col * sqrt(8.0*PC_k_CGS*Q.T[iT_] / (M_PI*m_av_)) * sigma);
    return tau;
}


ET_AppletonBray_Ion::
ET_AppletonBray_Ion( lua_State * L, int ie, int ic )
    : Relaxation_time(), ie_ ( ie ), ic_ ( ic )
{
    // 1. Colliding species data
    Chemical_species * X = get_library_species_pointer( ic );
    if ( X->get_Z()<1 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_Ion::ET_AppletonBray_Ion():\n";
	ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an ion.\n";
	input_error(ost);
    }
    M_c_ = X->get_M();

    // 2. Electron data
    X = get_library_species_pointer( ie );
    if ( X->get_Z()!=-1 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_Ion::ET_AppletonBray_Ion():\n";
	ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an electron.\n";
	input_error(ost);
    }
    iTe_ = X->get_iT_trans();
    M_e_ = X->get_M();
}

ET_AppletonBray_Ion::
~ET_AppletonBray_Ion() {}

double
ET_AppletonBray_Ion::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)
    double n_e = Q.massf[ie_] * Q.rho / M_e_ * PC_Avogadro * 1.0e-6;	// Number density of electrons (in cm**-3)

    // If there are no participating species set a negative relaxation time
    if ( n_e==0.0 || n_c==0.0 ) return -1.0;

    double tmpA = 8.0 / 3.0 * sqrt( M_PI / PC_m_CGS );
    double tmpB = n_c * pow ( PC_e_CGS, 4.0 ) / pow ( 2.0 * PC_k_CGS * Q.T[iTe_], 1.5 );
    double tmpC = log( pow( PC_k_CGS * Q.T[iTe_], 3.0 ) / ( M_PI * n_e * pow ( PC_e_CGS, 6.0 ) ) );
    double nu_ec = tmpA * tmpB * tmpC;

    double tau_ec = M_c_ / nu_ec;

    return tau_ec;

}

ET_AppletonBray_Neutral::
ET_AppletonBray_Neutral( lua_State * L, int ie, int ic )
    : Relaxation_time(), ic_( ic )
{
    int nCs;

    // 1. Colliding species data
    Chemical_species * X = get_library_species_pointer( ic );
    if ( X->get_Z()!=0 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_Neutral::ET_AppletonBray_Neutral():\n";
	ost << "Error in the declaration of colliding species: " << X->get_name() << " is not a neutral.\n";
	input_error(ost);
    }
    M_c_ = X->get_M();

    // 2. Electron data
    X = get_library_species_pointer( ie );
    if ( X->get_Z()!=-1 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_Neutral::ET_AppletonBray_Neutral():\n";
	ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an electron.\n";
	input_error(ost);
    }
    iTe_ = X->get_iT_trans();

    // 3.  Sigma quadratic curve fit coefficients
    lua_getfield(L, -1, "sigma_coefficients" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ET_AppletonBray_Neutral::ET_AppletonBray_Neutral():\n";
	ost << "Error in the declaration of sigma coefficients: a table is expected.\n";
	input_error(ost);
    }
    nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_Neutral::ET_AppletonBray_Neutral():\n";
	ost << "Error in the declaration of sigma coefficients: 3 coefficients are expected.\n";
	input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	C_.push_back( luaL_checknumber(L,-1) );
    	lua_pop(L,1);
    }
    lua_pop(L,1);	// pop 'sigma_coefficients'
}

ET_AppletonBray_Neutral::
~ET_AppletonBray_Neutral()
{
    C_.resize(0);
}

double
ET_AppletonBray_Neutral::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)

    // If there are no colliding species set a negative relaxation time
    if ( n_c==0.0 ) return -1.0;

    // calculate sigma using the appropriate curve fit for the electron temperature
    double T = Q.T[iTe_];
    double sigma_ec = 0.0;
    sigma_ec = C_[0] + C_[1] * T + C_[2] * T * T;
    /* Convert sigma_ec from m**2 -> cm**2 */
    sigma_ec *= 1.0e4;
    double nu_ec = n_c * sigma_ec * sqrt( 8.0 * PC_k_CGS * T / ( M_PI * PC_m_CGS ) );
    double tau_ec = M_c_ / nu_ec;

    return tau_ec;
}

ET_AppletonBray_TwoRangeNeutral::
ET_AppletonBray_TwoRangeNeutral( lua_State * L, int ie, int ic )
    : Relaxation_time(), ic_( ic )
{
    int nCs;

    // 1. Colliding species data
    Chemical_species * X = get_library_species_pointer( ic );
    if ( X->get_Z()!=0 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_TwoRangeNeutral::ET_AppletonBray_TwoRangeNeutral():\n";
	ost << "Error in the declaration of colliding species: " << X->get_name() << " is not a neutral.\n";
	input_error(ost);
    }
    M_c_ = X->get_M();

    // 2. Electron data
    X = get_library_species_pointer( ie );
    if ( X->get_Z()!=-1 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_TwoRangeNeutral::ET_AppletonBray_TwoRangeNeutral():\n";
	ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an electron.\n";
	input_error(ost);
    }
    iTe_ = X->get_iT_trans();

    // 3.  Sigma quadratic curve fit coefficients
    // 3a. Temperature switch
    T_switch_ = get_positive_number(L, -1, "T_switch");
    // 3b. Low temperature range
    lua_getfield(L, -1, "sigma_low_T" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "ET_AppletonBray_TwoRangeNeutral::ET_AppletonBray_TwoRangeNeutral():\n";
	ost << "Error in the declaration of low temperature sigma coefficients: a table is expected.\n";
	input_error(ost);
    }
    nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
	ostringstream ost;
	ost << "ET_AppletonBray_TwoRangeNeutral::ET_AppletonBray_TwoRangeNeutral():\n";
	ost << "Error in the declaration of low temperature sigma coefficients: 3 coefficients are expected.\n";
	input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
    	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
    	LT_C_.push_back( luaL_checknumber(L,-1) );
    	lua_pop(L,1);
    }
    lua_pop(L,1);	// pop 'LT_sigma_coefficients'
    // 3c. High temperature range
    lua_getfield(L, -1, "sigma_high_T" );
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "ET_AppletonBray_TwoRangeNeutral::ET_AppletonBray_TwoRangeNeutral():\n";
        ost << "Error in the declaration of high temperature sigma coefficients: a table is expected.\n";
        input_error(ost);
    }
    nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
        ostringstream ost;
        ost << "ET_AppletonBray_TwoRangeNeutral::ET_AppletonBray_TwoRangeNeutral():\n";
        ost << "Error in the declaration of high temperature sigma coefficients: 3 coefficients are expected.\n";
        input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
        lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
        HT_C_.push_back( luaL_checknumber(L,-1) );
        lua_pop(L,1);
    }
    lua_pop(L,1);       // pop 'HT_sigma_coefficients'
}

ET_AppletonBray_TwoRangeNeutral::
~ET_AppletonBray_TwoRangeNeutral()
{
    LT_C_.resize(0);
    HT_C_.resize(0);
}

double
ET_AppletonBray_TwoRangeNeutral::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;	// Number density of colliders (in cm**-3)

    // If there are no colliding species set a negative relaxation time
    if ( n_c==0.0 ) return -1.0;

    // calculate sigma using the appropriate curve fit for the electron temperature
    double T = Q.T[iTe_];
    double sigma_ec = 0.0;
    if ( T < T_switch_ )
        sigma_ec = LT_C_[0] + LT_C_[1] * T + LT_C_[2] * T * T;
    else
        sigma_ec = HT_C_[0] + HT_C_[1] * T + HT_C_[2] * T * T;
    /* Convert sigma_ec from m**2 -> cm**2 */
    sigma_ec *= 1.0e4;
    double nu_ec = n_c * sigma_ec * sqrt( 8.0 * PC_k_CGS * T / ( M_PI * PC_m_CGS ) );
    double tau_ec = M_c_ / nu_ec;

    return tau_ec;
}

ER_Abe_Ion::
ER_Abe_Ion( lua_State * L, int ie, int ic )
    : Relaxation_time(), ie_ ( ie ), ic_ ( ic )
{
    // 1. Colliding species data
    Chemical_species * X = get_library_species_pointer( ic );
    if ( X->get_Z()<1 ) {
        ostringstream ost;
        ost << "ER_Abe_Ion::ER_Abe_Ion():\n";
        ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an ion.\n";
        input_error(ost);
    }
    // Check that it is a diatomic molecule
    if ( X->get_type().find("diatomic")==string::npos ) {
        ostringstream ost;
        ost << "ER_Abe_Ion::ER_Abe_Ion():\n";
        ost << "Error in the declaration of colliding species: " << X->get_name() << " is not a diatomic molecule.\n";
        input_error(ost);
    }

    M_c_ = X->get_M();

    // 2. Electron data
    X = get_library_species_pointer( ie );
    if ( X->get_Z()!=-1 ) {
        ostringstream ost;
        ost << "ER_Abe_Ion::ER_Abe_Ion():\n";
        ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an electron.\n";
        input_error(ost);
    }
    iTe_ = X->get_iT_trans();
    M_e_ = X->get_M();

    // 3. Rotational scaling factor ( from Abe et al 2002: sigma_ER = sigma_ET * g_rot )
    g_rot_ = get_positive_number( L, -1, "g_rot" );
}

ER_Abe_Ion::
~ER_Abe_Ion() {}

double
ER_Abe_Ion::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;    // Number density of colliders (in cm**-3)
    double n_e = Q.massf[ie_] * Q.rho / M_e_ * PC_Avogadro * 1.0e-6;    // Number density of electrons (in cm**-3)

    // If there are no participating species set a negative relaxation time
    if ( n_e==0.0 || n_c==0.0 ) return -1.0;

    double tmpA = 8.0 / 3.0 * sqrt( M_PI / PC_m_CGS );
    double tmpB = n_c * pow ( PC_e_CGS, 4.0 ) / pow ( 2.0 * PC_k_CGS * Q.T[iTe_], 1.5 );
    double tmpC = log( pow( PC_k_CGS * Q.T[iTe_], 3.0 ) / ( M_PI * n_e * pow ( PC_e_CGS, 6.0 ) ) );
    double nu_ec = tmpA * tmpB * tmpC;

    double tau_ec = M_c_ / nu_ec / g_rot_;

    return tau_ec;

}

ER_Abe_Neutral::
ER_Abe_Neutral( lua_State * L, int ie, int ic )
    : Relaxation_time(), ic_( ic )
{
    int nCs;

    // 1. Colliding species data
    Chemical_species * X = get_library_species_pointer( ic );
    if ( X->get_Z()!=0 ) {
        ostringstream ost;
        ost << "ER_Abe_Neutral::ER_Abe_Neutral():\n";
        ost << "Error in the declaration of colliding species: " << X->get_name() << " is not a neutral.\n";
        input_error(ost);
    }
    // Check that it is a diatomic molecule
    if ( X->get_type().find("diatomic")==string::npos ) {
        ostringstream ost;
        ost << "ER_Abe_Neutral::ER_Abe_Neutral():\n";
        ost << "Error in the declaration of colliding species: " << X->get_name() << " is not a diatomic molecule.\n";
        input_error(ost);
    }
    M_c_ = X->get_M();

    // 2. Electron data
    X = get_library_species_pointer( ie );
    if ( X->get_Z()!=-1 ) {
        ostringstream ost;
        ost << "ER_Abe_Neutral::ER_Abe_Neutral():\n";
        ost << "Error in the declaration of colliding species: " << X->get_name() << " is not an electron.\n";
        input_error(ost);
    }
    iTe_ = X->get_iT_trans();

    // 3.  Sigma quadratic curve fit coefficients
    lua_getfield(L, -1, "sigma_coefficients" );
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "ER_Abe_Neutral::ER_Abe_Neutral():\n";
        ost << "Error in the declaration of sigma coefficients: a table is expected.\n";
        input_error(ost);
    }
    nCs = lua_objlen(L, -1);
    if ( nCs!=3 ) {
        ostringstream ost;
        ost << "ER_Abe_Neutral::ER_Abe_Neutral():\n";
        ost << "Error in the declaration of sigma coefficients: 3 coefficients are expected.\n";
        input_error(ost);
    }
    for ( int i=0; i<nCs; ++i ) {
        lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
        C_.push_back( luaL_checknumber(L,-1) );
        lua_pop(L,1);
    }
    lua_pop(L,1);       // pop 'sigma_coefficients'

    // 3. Rotational scaling factor ( from Abe et al 2002: sigma_ER = sigma_ET * g_rot )
    g_rot_ = get_positive_number( L, -1, "g_rot" );
}

ER_Abe_Neutral::
~ER_Abe_Neutral()
{
    C_.resize(0);
}

double
ER_Abe_Neutral::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double n_c = Q.massf[ic_] * Q.rho / M_c_ * PC_Avogadro * 1.0e-6;    // Number density of colliders (in cm**-3)

    // If there are no colliding species set a negative relaxation time
    if ( n_c==0.0 ) return -1.0;

    // calculate sigma using the appropriate curve fit for the electron temperature
    double T = Q.T[iTe_];
    double sigma_ec = 0.0;
    sigma_ec = C_[0] + C_[1] * T + C_[2] * T * T;
    /* Convert sigma_ec from m**2 -> cm**2 */
    sigma_ec *= 1.0e4;
    double nu_ec = n_c * sigma_ec * sqrt( 8.0 * PC_k_CGS * T / ( M_PI * PC_m_CGS ) );
    double tau_ec = M_c_ / nu_ec / g_rot_;

    return tau_ec;
}

/**
 * \brief Calculate the collider distance for like colliders.
 *
 * This function (A) is applicable for like colliders, ie.
 * polar-polar or nonpolar-nonpolar.
 * See Hirschfelder. p. 222
 *
 * \[ r^0_{pq} = \frac{1}{2} \left( r^0_{pp} + r^{0}_{qq} \right) \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double collider_distance_A(Diatomic_species &p, Diatomic_species &q)
{
    return (0.5 * (p.get_r0() + q.get_r0()));
}

/**
 * \brief Calculate the collider distance for unlike colliders.
 *
 * This function (B) is applicable when the colliders are unalike,
 * ie. polar-nonpolar
 *
 * \[ r^{0}_{pn} = \frac{1}{2} \left( r^{0}_{pp} + r^{0}_{nn} \right) \xi^{-1/6} \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double collider_distance_B(Diatomic_species &p, Diatomic_species &n)
{
    double r_pn, xi;

    xi = xi_correction(p, n);
    r_pn = 0.5 * (p.get_r0() + n.get_r0()) * pow(xi, -1.0/6.0);

    return (r_pn);
}

/**
 * \brief Calculate xi (a correction factor)
 *
 * This function calculates a correction factor needed for the
 * the empirical combining laws when there is a collision
 * between polar and nonpolar molecules.
 *
 * \[ \xi = \left\[ 1 + \frac{\alpha_n \mu^{*2}_p}{4 \r^3_n}
 *                      \sqrt{\frac{\epsilon_p}{\epsilon_n}} \right\] \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double xi_correction(Diatomic_species &p, Diatomic_species &n)
{
    double mu_star, tmp_a, tmp_b, xi;
    // Require eps in K (convert from J).
    double eps_p = p.get_eps0()/PC_k_SI;
    double eps_n = n.get_eps0()/PC_k_SI;
    // Require r0 in Angstrom
    double r_p = p.get_r0() * 1.0e10;
    double r_n = n.get_r0() * 1.0e10;

    mu_star = p.get_mu_B() / sqrt( eps_p * pow(r_p, 3.0));
    tmp_a = (n.get_alpha() * mu_star*mu_star) / pow(r_n, 3.0);
    tmp_b = sqrt(eps_p / eps_n);

    xi = 1.0 + 0.25*tmp_a*tmp_b;

    return xi;
}

/**
 * \brief Calculate eps0_pq: potential well for like colliders.
 *
 * This function (A) calculates the potential well in a Lennard-Jones
 * type collision for two alike molecules, ie. polar-polar or
 * nonpolar-nonpolar.
 *
 * \[ \epsilon^0_{pq} = \sqrt{ \epsilon^0_{pp} \epsilon^0_{qq} } \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double potential_well_A(Diatomic_species &p, Diatomic_species &q)
{
    return ( sqrt(p.get_eps0() * q.get_eps0()) );
}

/**
 * \brief Calculate eps0_pq: potential well for unlike colliders.
 *
 * This function (B) calculates the potential well for collisions
 * between unalike molecules, ie. polar-nonpolar.
 *
 * \[ \epsilon^0_{pn} = \sqrt{\epsilon^0_{pp} \epsilon^0_{qq}} \xi^2 \]
 *
 * \author Rowan J Gollan
 * \date 10-Mar-05
 **/

double potential_well_B(Diatomic_species &p, Diatomic_species &n)
{
    double eps_pq, xi;
    xi = xi_correction(p, n);
    eps_pq = sqrt(p.get_eps0()*n.get_eps0())*xi*xi;

    return eps_pq;
}

/**
 * \brief Calculate the collision frequency based on the hard-sphere model.
 *
 * \[ Z^{pq}_{\mbox{coll}} = 2\sqrt{2} nB \sigma^2_{pq} \sqrt{\pi k_b T / \mu_{pq}} \]
 *
 * Note that this gives collision frequency: the average number of collisions
 * for any one molecule per unit time per unit volume. This is not to be confused
 * with total number of collisions (sometimes called collision rate).
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double collision_frequency(double sigma, double mu, double T, double nB)
{
    double tmp;
    double Z;

    tmp = sqrt(2.0 * M_PI * PC_k_SI * T / mu);
    Z = 2.0 * nB * sigma * sigma * tmp;

    return Z;
}

/**
 * \brief Calculate expansion parameter beta.
 *
 * \[ \beta = \{ \frac{1}{2} (2 \epsilon^0_{pq} / \mu_{pq})^9
 *              (3h\mu_{pq}/\pi^2 r^0_{pq} k T \Delta E)^6 \}^{1/19} \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_beta(double eps0, double mu, double r0,
		double del_E, double T)
{

    double tmp_a, tmp_b;
    double beta;

    tmp_a = 2.0 * eps0 / mu;
    tmp_b = 3.0 * PC_h_SI * mu /
	(M_PI * M_PI * r0 * PC_k_SI * T * del_E);

    beta = pow( 0.5 * pow(tmp_a, 9) * pow(tmp_b, 6), 1.0/19.0);

    return (beta);
}

/**
 * \brief Calculate D_star using a first-order expression.
 *
 * Implements:
 * \[ D^* = (1/\beta^2)( 1 - \frac{14}{19} \beta ) \]
 *
 * \author Rowan J Gollan
 * \date 01-March-2005
 **/

double SSH_D_star(double beta)
{
    return ((1.0 / (beta * beta)) * (1.0 - (14.0/19.0) * beta));
}

/**
 * \brief Calculate r_c_star.
 *
 * r_c_star is the classical distance of closest
 * approach based on kinetic theory.
 * Implements:
 * \[ r^{c*}_{pq} = r^0_{pq}/ ( (1/(2\beta)^{1/6})
 *                              (1 + 2/19 \beta) ) \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_r_c_star(double beta, double r0)
{
    double tmp_a, tmp_b;

    tmp_a = 1.0 / pow(2.0 * beta, 1.0/6.0);
    tmp_b = 1.0 + 2.0/19.0 * beta;

    return ( r0 / (tmp_a * tmp_b) );
}

double SSH_r_c_star2(double D_star, double r0)
{

    return (r0 * pow( (0.5 + 0.5*sqrt(D_star)), 6.0));

}

/**
 * \brief Calculate del_star.
 *
 * del_star is a parameter of the exponential potential.
 * Implements:
 * \[ \delta^*_{pq} = \frac{1}{r^0_{pq}} ( 12/(2\beta)^{1/6}) (1 + 2/19 \beta) \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_del_star(double beta, double r0)
{
    double tmp_a, tmp_b;

    tmp_a = 12.0 / pow(2.0 * beta, 1.0/6.0);
    tmp_b = 1.0 + 21.0 / 19.0 * beta;

    return ( (tmp_a * tmp_b) / r0 );
}

double SSH_del_star2(double D_star, double r0)
{
    double tmp_a, tmp_b;

    tmp_a = sqrt(D_star);
    tmp_b = pow((0.5 + 0.5*tmp_a), 1.0/6.0);
    return ( (12.0/r0) * tmp_b * (1.0 + 1.0/tmp_a) );
}

/**
 * \brief Calculate alpha_pq.
 *
 * alpha_pq is the main parameter which describes
 * the relaxation tims'e dependence on temperature.
 * Implements:
 * \[ \alpha_{pq} = \frac{16 \pi^2 \mu_{pq} \Delta E^2}
 *                       {\delta^{*2}_{pq} h^2 k}  \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_alpha_pq(double mu, double del_E, double del_star)
{
    double tmp_a, tmp_b;

    tmp_a = 16.0 * pow(M_PI, 4) * mu * del_E * del_E;
    tmp_b = del_star * del_star * PC_h_SI * PC_h_SI * PC_k_SI;

    return ( tmp_a/tmp_b );
}

/**
 * \brief Calculate chi_pq.
 *
 * chi_pq is an interaction parameter of molecules
 * p and q.
 * Implements:
 * \[ \chi^*_{pq} = \frac{1}{2} \left( \frac{alpha_{pq}}{T} \right)^{1/3} \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_chi_pq(double alpha_pq, double T)
{
    return ( 0.5 * pow(alpha_pq/T, 1.0/3.0) );
}

/**
 * \brief Calculate Z_0.
 *
 * Z_0 is a steric factor reqpresenting the orientation
 * of collisions.
 * Implements:
 * \[ Z^p_0 = \delta^*_{pq} * r^p_{eq} +
 *            \frac{5}{2} \left[1/(\delta^*_{pq} r^p_{eq})^2 \right] \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_0(double del_star, double r_eq)
{
    double tmp_a;

    tmp_a = del_star * r_eq;

    return ( tmp_a + 5.0/2.0 * ( 1.0 / (tmp_a * tmp_a)) );
}

/**
 * \brief Calculate Z_V.
 *
 * Z_V is the vibrational factor.
 * Implements:
 * \[ Z^p_V = \frac{f^p_m}{\pi^2} \frac{\mu_{pp}}{\mu_{pq}}
 *            \frac{\alpha_{pq}}{\theta^p_v}
 *            \left( \frac{k \theta^p_v}{\Delta E} \right)^2 \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_V(double f_m, double mu_pp, double mu_pq, double alpha_pq,
	       double theta_v, double del_E, int i)
{
    double tmp_a, tmp_b, tmp_c, tmp_d;
    tmp_a = f_m / ((i+1.0)*M_PI*M_PI);
    tmp_b = mu_pp / mu_pq;
    tmp_c = alpha_pq / theta_v;
    tmp_d = PC_k_SI * theta_v / del_E;
    return ( tmp_a * tmp_b * tmp_c * tmp_d*tmp_d );
}

/**
 * \brief Calculate Z_T
 *
 * Z_T is the translational factor.
 * Implements:
 * \[ Z^{pq}_T = \pi^2 \left( \frac{3}{2\pi} \right)^{1/2}
 *               \left( \frac{\Delta E}{k \alpha_{pq}} \right)^2
 *               \left( \frac{T}{\alpha_{pq}} \right)^{1/6}
 *               \times \exp \left[ \frac{3}{2}
 *                                  \left( \frac{\alpha_{pq}}{T} \right)^{1/3}
 *                                   - \frac{\Delta E}{2kT} \right] \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_T(double del_E, double alpha_pq, double T)
{
    double tmp_a, tmp_b, tmp_c;

    tmp_a = del_E / (PC_k_SI * alpha_pq);
    tmp_b = 3.0/2.0 * pow(alpha_pq/T, 1.0/3.0) - (del_E / (2.0*PC_k_SI*T));
    tmp_c = sqrt(3.0/(2.0*M_PI)) * tmp_a * tmp_a * pow(T/alpha_pq, 1.0/6.0);

    return (M_PI*M_PI * tmp_c * exp(tmp_b));
}

/**
 * \brief Calculate Z_plus.
 *
 * Z_plus represents the contribution from attractive forces.
 * This function implements:
 * \[ Z^{pq}_{+} = \exp \left[ - \frac{4}{\pi} \left(\frac{\epsilon^0_{pq} \chi^{*}_{pq}}{kT} \right)^{1/2}
 *                             - \frac{16}{3\pi^2} \frac{ \epsilon^0_{pq}}{kT} \right] \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-2005
 **/

double SSH_Z_plus(double eps0, double chi, double T)
{
    double tmp_a, tmp_b;

    tmp_a = (-4.0 / M_PI) * sqrt(eps0*chi / (PC_k_SI*T));
    tmp_b = (-16.0*eps0) / (3.0*M_PI*M_PI*PC_k_SI*T);

    return ( exp(tmp_a + tmp_b) );
}

/**
 * \brief Calculate A.
 *
 * The scripted A represents a collision cross-reference
 * factor.
 * This function implements:
 * \[ \cal{A} = \left( \frac{r_{c*}_{pq}}{\sigma_{pq}} \right)^2 \]
 *
 * \author Rowan J Gollan
 * \date 01-Mar-05
 **/

double SSH_A_factor(double r_c, double sigma)
{
    return ( (r_c*r_c) / (sigma*sigma) );
}

VT_SSH::
VT_SSH(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    // 1. Setup values related to 'p' species
    Diatomic_species *p = get_library_diatom_pointer(ip);
    M_p_ = p->get_M();
    r_eq_p_ = p->get_r_eq();
    f_m_p_ = p->get_f_m();
    mu_p_ = p->get_mu();
    Truncated_harmonic_vibration *p_vib = dynamic_cast<Truncated_harmonic_vibration*>(p->get_mode_pointer_from_type("vibration"));
    if ( p_vib == 0 ) {
	cout << "The vibrating species " << p->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_p_ = p_vib->get_theta();

    // 2. Setup values related to 'q' species
    Diatomic_species *q = get_library_diatom_pointer(iq);
    M_q_ = q->get_M();
    r_eq_q_ = q->get_r_eq();
    f_m_q_ = q->get_f_m();
    mu_q_ = q->get_mu();
    Truncated_harmonic_vibration *q_vib = dynamic_cast<Truncated_harmonic_vibration*>(q->get_mode_pointer_from_type("vibration"));
    if ( q_vib == 0 ) {
	cout << "The vibrating species " << q->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_q_ = q_vib->get_theta();

    // 3. Derived parameters
    if ( p->get_polar_flag() ) {
	if ( q->get_polar_flag() ) {
	    r_ = collider_distance_A(*p, *q);
	    eps_ = potential_well_A(*p, *q);
	}
	else {
	    r_ = collider_distance_B(*p, *q);
	    eps_ = potential_well_B(*p, *q);
	}
    }
    else {
	if ( q->get_polar_flag() ) {
	    r_ = collider_distance_B(*q, *p);
	    eps_ = potential_well_B(*q, *p);
	}
	else {
	    r_ = collider_distance_A(*p, *q);
	    eps_ = potential_well_A(*p, *q);
	}
    }

    mu_ = ((M_p_ * M_q_) / (M_p_ + M_q_)) / PC_Avogadro;
    delta_E_ = PC_k_SI * theta_v_p_;
    sigma_ = r_;
}

VT_SSH::
~VT_SSH() {}

double
VT_SSH::
specific_transition_probability(Gas_data &Q, vector<double> &molef)
{
    // Values that change with temperature
    double T = Q.T[iT_];
    double beta = SSH_beta(eps_, mu_, r_, delta_E_, T);
    double r_c_star = SSH_r_c_star(beta, r_);
    double del_star = SSH_del_star(beta, r_);
    double alpha_pq = SSH_alpha_pq(mu_, delta_E_, del_star);
    double chi_pq = SSH_chi_pq(alpha_pq, T);

    // The steric factors
    double Z_0 = SSH_Z_0(del_star, r_eq_p_);
    double Z_V = SSH_Z_V(f_m_p_, mu_p_, mu_, alpha_pq, theta_v_p_, delta_E_, 0);
    double Z_T = SSH_Z_T(delta_E_, alpha_pq, T);
    double Z_plus = SSH_Z_plus(eps_, chi_pq, T);

    // Collision cross-reference factor
    double A = SSH_A_factor(r_c_star, sigma_);
    double P = A/(Z_0*Z_V*Z_T*Z_plus);

    return P;
}

double
VT_SSH::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double T = Q.T[iT_];

    // If either collider has very little concentration,
    // then a relaxation time has no meaning
    if ( molef[ip_] <= DEFAULT_MIN_MOLE_FRACTION ||
	 molef[iq_] <= DEFAULT_MIN_MOLE_FRACTION ) {
	return -1.0;
    }

    // The computed collision frequency is independent
    // of whether the collisions are between:
    //     like colliders : ip == iq; or
    //     unlike colliders : ip != iq
    // See pp 87,88 of Chapman and Cowling
    //
    // Chapman and Cowling (1970)
    // The Mathematical Theory of Non-Uniform Gases, Third edition
    // Cambridge University Press, London
    //

    double n_q = molef[iq_]*Q.p/(PC_k_SI*T);
    double Z = collision_frequency(sigma_, mu_, T, n_q);
    double P = specific_transition_probability(Q, molef);
    double tau = 1.0 / (Z*P*(1 - exp(-1.0*theta_v_p_/T)));
    return tau;
}

VV_SSH::
VV_SSH(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    // 1. Setup values related to 'p' species
    Diatomic_species *p = get_library_diatom_pointer(ip);
    M_p_ = p->get_M();
    r_eq_p_ = p->get_r_eq();
    f_m_p_ = p->get_f_m();
    mu_p_ = p->get_mu();
    Truncated_harmonic_vibration *p_vib = dynamic_cast<Truncated_harmonic_vibration*>(p->get_mode_pointer_from_type("vibration"));
    if ( p_vib == 0 ) {
	cout << "The vibrating species " << p->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_p_ = p_vib->get_theta();

    // 2. Setup values related to 'q' species
    Diatomic_species *q = get_library_diatom_pointer(iq);
    M_q_ = q->get_M();
    r_eq_q_ = q->get_r_eq();
    f_m_q_ = q->get_f_m();
    mu_q_ = q->get_mu();
    Truncated_harmonic_vibration *q_vib = dynamic_cast<Truncated_harmonic_vibration*>(q->get_mode_pointer_from_type("vibration"));
    if ( q_vib == 0 ) {
	cout << "The vibrating species " << q->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_q_ = q_vib->get_theta();

    // 3. Derived parameters
    if ( p->get_polar_flag() ) {
	if ( q->get_polar_flag() ) {
	    r_ = collider_distance_A(*p, *q);
	    eps_ = potential_well_A(*p, *q);
	}
	else {
	    r_ = collider_distance_B(*p, *q);
	    eps_ = potential_well_B(*p, *q);
	}
    }
    else {
	if ( q->get_polar_flag() ) {
	    r_ = collider_distance_B(*q, *p);
	    eps_ = potential_well_B(*q, *p);
	}
	else {
	    r_ = collider_distance_A(*p, *q);
	    eps_ = potential_well_A(*p, *q);
	}
    }

    if ( p->get_name() == q->get_name() ) {
	cout << "VV_SSH::VV_SSH(): " << endl;
	cout << "Problem detected setting up VV exchange rate." << endl;
	cout << "The two different collider species are specified as the SAME species." << endl;
	cout << "VV exchanges model the transfer of vibrational energy between" << endl;
	cout << "two different species of molecules. " << endl;
	cout << "The ONE species specified is:" << endl;
	cout << "species name= " << p->get_name() << endl;
	cout << "This has no physical meaning. Any transfer of vibrational energy" << endl;
	cout << "between molecules of the SAME type has no overall effect on" << endl;
	cout << "the average vibrational energy of that species as a collection." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }

    mu_ = ((M_p_ * M_q_) / (M_p_ + M_q_)) / PC_Avogadro;
    delta_E_ = PC_k_SI * (theta_v_p_ - theta_v_q_);
    sigma_ = r_;
}

VV_SSH::
~VV_SSH() {}

double
VV_SSH::
specific_transition_probability(Gas_data &Q, vector<double> &molef)
{
    // Values that change with temperature
    double T = Q.T[iT_];
    double beta = SSH_beta(eps_, mu_, r_, delta_E_, T);
    double r_c_star = SSH_r_c_star(beta, r_);
    double del_star = SSH_del_star(beta, r_);
    double alpha_pq = SSH_alpha_pq(mu_, delta_E_, del_star);
    double chi_pq = SSH_chi_pq(alpha_pq, T);

    // The steric factors
    double Z_0_p = SSH_Z_0(del_star, r_eq_p_);
    double Z_0_q = SSH_Z_0(del_star, r_eq_q_);
    double Z_V_p = SSH_Z_V(f_m_p_, mu_p_, mu_, alpha_pq, theta_v_p_, delta_E_, 2);
    double Z_V_q = SSH_Z_V(f_m_q_, mu_q_, mu_, alpha_pq, theta_v_q_, delta_E_, 1);
    double Z_T = SSH_Z_T(delta_E_, alpha_pq, T);
    double Z_plus = SSH_Z_plus(eps_, chi_pq, T);

    // Collision cross-reference factor
    double A = SSH_A_factor(r_c_star, sigma_);
    double P = A/(Z_0_p*Z_0_q*Z_V_p*Z_V_q*Z_T*Z_plus);

    return P;
}

double
VV_SSH::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double T = Q.T[iT_];

    // If either collider has very small concentration
    // then a relaxation time has no meaning
    if ( molef[ip_] <= DEFAULT_MIN_MOLE_FRACTION ||
	 molef[iq_] <= DEFAULT_MIN_MOLE_FRACTION ) {
 	return -1.0;
    }

    // The computed collision frequency is independent
    // of whether the collisions are between:
    //     like colliders : ip == iq; or
    //     unlike colliders : ip != iq
    // See pp 87,88 of Chapman and Cowling
    //
    // Chapman and Cowling (1970)
    // The Mathematical Theory of Non-Uniform Gases, Third edition
    // Cambridge University Press, London
    //
    double n_q = molef[iq_]*Q.p/(PC_k_SI*T);
    double Z = collision_frequency(sigma_, mu_, T, n_q);
    double P = specific_transition_probability(Q, molef);
    double tau = (1.0/(Z*P))*exp(-1.0*(theta_v_p_ - theta_v_q_)/T);

    return tau;
}

// VV_Candler::
// VV_Candler( lua_State * L )
//     : Relaxation_time()
// {
//     // 1. Vibrating species 'p'
//     string p_name = get_string(L, -1, "p_name");
//     Diatomic_species * P = get_library_diatom_pointer_from_name( p_name );
//     ip_ = P->get_isp();
//     iT_ = P->get_mode_pointer_from_type("translation")->get_iT();
//     double M_p = P->get_M();
//     double r0_p = P->get_r0();

//     // 2. Vibrating species 'q'
//     string q_name = get_string(L, -1, "q_name");
//     Diatomic_species * Q = get_library_diatom_pointer_from_name( q_name );
//     iq_ = Q->get_isp();
//     M_q_ = Q->get_M();
//     double r0_q = Q->get_r0();

//     // 3. derived parameters
//     mu_ = ((M_p * M_q_) / (M_p + M_q_)) / PC_Avogadro;
//     // CHECKME: Candler uses d_p*dq here. is there a difference between r0 and d?
//     sigma_ = r0_p * r0_q;	// collision cross-section

//     // 4. Transition probability (use fixed value as recommended in Candler thesis)
//     P_ = 1.0e-2;
// }

// VV_Candler::
// ~VV_Candler()
// {}

// double
// VV_Candler::
// specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
// {
//     // If either collider has zero concentration
//     // then a relaxation time has no meaning
//     if( molef[ip_] <= 0.0 || molef[iq_] <= 0.0 ) {
// 	return -1.0;
//     }

//     // Calculate collision frequency x average vibrational energy
//     double Z_eps = PC_Avogadro * sigma_ * sqrt( 8.0 * PC_R_u * Q.T[iT_] / ( M_PI * mu_ ) )
//     			* Q.rho * Q.massf[iq_] / M_q_;

//     // Calculate relaxation time
//     double tau = 1.0 / ( P_ * Z_eps );

//     return tau;
// }


VV_MTLandauTeller::
VV_MTLandauTeller(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    // 1. Setup values related to 'p' species
    Diatomic_species *p = get_library_diatom_pointer(ip);
    double M_p = p->get_M();
    Truncated_harmonic_vibration *p_vib = dynamic_cast<Truncated_harmonic_vibration*>(p->get_mode_pointer_from_type("vibration"));
    if ( p_vib == 0 ) {
	cout << "The vibrating species " << p->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_p_ = p_vib->get_theta();

    // 2. Setup values related to 'q' species
    Diatomic_species *q = get_library_diatom_pointer(iq);
    double M_q = q->get_M();
    Truncated_harmonic_vibration *q_vib = dynamic_cast<Truncated_harmonic_vibration*>(q->get_mode_pointer_from_type("vibration"));
    if ( q_vib == 0 ) {
	cout << "The vibrating species " << q->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_q_ = q_vib->get_theta();

    // 3. Pull model parameters from Lua table
    A_ = get_number(L, -1, "A");
    n_ = get_number(L, -1, "n");
    B1_ = get_number(L, -1, "B1");
    B2_ = get_number(L, -1, "B2");
    B3_ = get_number(L, -1, "B3");
    beta_ = get_number(L, -1, "beta");

    // 4. Derived parameters
    if ( p->get_polar_flag() ) {
	if ( q->get_polar_flag() ) {
	    R0_ = collider_distance_A(*p, *q);
	}
	else {
	    R0_ = collider_distance_B(*q, *p);
	}
    }
    else {
	if ( q->get_polar_flag() ) {
	    R0_ = collider_distance_B(*p, *q);
	}
	else {
	    R0_ = collider_distance_A(*p, *q);
	}
    }

    mu_ = ((M_p * M_q) / (M_p + M_q)) / PC_Avogadro;

    UNUSED_VARIABLE(ip_);
}

VV_MTLandauTeller::
~VV_MTLandauTeller() {}

double
VV_MTLandauTeller::
specific_relaxation_time(Gas_data &Q, vector<double> &molef)
{
    double T = Q.T[iT_];
    double k_10 = rate_coefficient(T);
    double p_in_atm = Q.p/PC_P_atm;
    double tau = (1.0/p_in_atm)*(PC_k_SI*T)/(k_10*(1.0 - exp(-theta_v_p_/T)));
    return tau;
}

double
VV_MTLandauTeller::
specific_transition_probability(Gas_data &Q, vector<double> &molef)
{
    double T = Q.T[iT_];
    double tau = compute_relaxation_time(Q, molef);
    double n_q = molef[iq_]*Q.p/(PC_k_SI*T);
    double Z = collision_frequency(R0_, mu_, T, n_q);
    double P_10 = (1.0/(Z*tau))*exp(-1.0*(theta_v_p_ - theta_v_q_)/T);
    return P_10;
}

double
VV_MTLandauTeller::
rate_coefficient(double T)
{
    double k_10 = A_*pow(T, n_)*exp(B1_*pow(T, -1.0/3.0) + B2_*pow(T, -2.0/3.0) + B3_*pow(T, beta_));
    return k_10;
}

VV_from_eq::
VV_from_eq(lua_State *L, int ip, int iq, int itrans)
    : Relaxation_time(), ip_(ip), iq_(iq), iT_(itrans)
{
    // 1. Setup values related to 'p' species
    Diatomic_species *p = get_library_diatom_pointer(ip);
    double M_p = p->get_M();
    Truncated_harmonic_vibration *p_vib = dynamic_cast<Truncated_harmonic_vibration*>(p->get_mode_pointer_from_type("vibration"));
    if ( p_vib == 0 ) {
	cout << "The vibrating species " << p->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_p_ = p_vib->get_theta();

    // 2. Setup values related to 'q' species
    Diatomic_species *q = get_library_diatom_pointer(iq);
    double M_q = q->get_M();
    Truncated_harmonic_vibration *q_vib = dynamic_cast<Truncated_harmonic_vibration*>(q->get_mode_pointer_from_type("vibration"));
    if ( q_vib == 0 ) {
	cout << "The vibrating species " << q->get_name() << " could not be cast" << endl;
	cout << "as a truncated harmonic oscillator." << endl;
	cout << "Bailing out!" << endl;
	exit(BAD_INPUT_ERROR);
    }
    theta_v_q_ = q_vib->get_theta();

    // 3. Some derived parameters
    if ( p->get_polar_flag() ) {
	if ( q->get_polar_flag() ) {
	    R0_ = collider_distance_A(*p, *q);
	}
	else {
	    R0_ = collider_distance_B(*q, *p);
	}
    }
    else {
	if ( q->get_polar_flag() ) {
	    R0_ = collider_distance_B(*p, *q);
	}
	else {
	    R0_ = collider_distance_A(*p, *q);
	}
    }

    mu_ = ((M_p * M_q) / (M_p + M_q)) / PC_Avogadro;

    // 4. Get forward relaxation time model to TOS
    lua_getfield(L, -1, "forward_rt");
    lua_rawgeti(L, -1, 1);
    string model = luaL_checkstring(L, -1);
    lua_pop(L, 1);
    if ( model == "Landau-Teller-VV" ) {
	// Swap order of iq and ip because using relaxation time in other direction
	rt_ = new VV_MTLandauTeller(L, iq, ip, itrans);
    }
    else if ( model == "SSH-VV" ) {
	// Swap order of iq and ip because using relaction time in other direction
	rt_ = new VV_SSH(L, iq, ip, itrans);
    }
    else {
	cout << "Error initialising VV_from_eq().\n";
	cout << "This model can only be used with the following models for forward_rt:\n";
	cout << "   Landau-Teller-VV\n";
	cout << "   SSH-VV\n";
	cout << "Bailing out!\n";
	exit(1);
    }
    lua_pop(L, 1);

    UNUSED_VARIABLE(ip_);
}

VV_from_eq::
~VV_from_eq()
{
    delete rt_;
}

double
VV_from_eq::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    double T = Q.T[iT_];
    double P_10 = compute_transition_probability(Q, molef);
    double n_q = molef[iq_]*Q.p/(PC_k_SI*T);
    double Z = collision_frequency(R0_, mu_, T, n_q);
    double tau = (1.0/(Z*P_10))*exp(-1.0*(theta_v_p_ - theta_v_q_)/T);
    return tau;
}

double
VV_from_eq::
specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
{
    double T = Q.T[iT_];
    double P_01 = rt_->compute_transition_probability(Q, molef);
    double P_10 = P_01*exp(theta_v_p_/T)/exp(theta_v_q_/T);
    return P_10;
}

VE_Lee::
VE_Lee( lua_State * L, int ie, int iv )
    : Relaxation_time(), ie_( ie ), iv_( iv )
{
    // 1. Free electron and vibrational species index (and iTe)
    Chemical_species * E = get_library_species_pointer( ie );
    // Check that it is an electron
    if ( E->get_type().find("electron")==string::npos ) {
        ostringstream ost;
        ost << "VE_Lee::VE_Lee():\n";
        ost << "Error in the declaration of colliding species: " << E->get_name() << " is not an electron.\n";
        input_error(ost);
    }
    iTe_ = E->get_mode_pointer_from_type("translation")->get_iT();

    // 2. Temperature switches
    lua_getfield(L, -1, "T_switches" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "VE_Lee::VE_Lee():\n";
        ost << "Error in the declaration of T switches: a table is expected.\n";
 	input_error(ost);
     }
     for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
     	lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
     	T_switches_.push_back( luaL_checknumber(L,-1) );
     	lua_pop(L,1);
     }
     lua_pop(L,1);	// pop 'T_switches'

     // 3. p.tau coefficients
     ptau_coeffs_.resize( T_switches_.size() + 1 );
     for ( size_t ic=0; ic<T_switches_.size()+1; ++ic ) {
     	ostringstream ptau_label;
     	ptau_label << "ptau_coefficients_" << ic;
 	lua_getfield(L, -1, ptau_label.str().c_str() );
 	if ( !lua_istable(L, -1) ) {
 	    ostringstream ost;
 	    ost << "VE_Lee::VE_Lee():\n";
 	    ost << "Error in the declaration of " << ptau_label.str() << ": a table is expected.\n";
 	    input_error(ost);
 	}
 	int nCs = lua_objlen(L, -1);
 	if ( nCs!=3 ) {
 	    ostringstream ost;
 	    ost << "VE_Lee::VE_Lee():\n";
 	    ost << "Error in the declaration of p.tau coefficients: 3 coefficients are expected.\n";
 	    input_error(ost);
 	}
 	for ( int i=0; i<nCs; ++i ) {
 	    lua_rawgeti(L, -1, i+1); // A Lua list is offset one from the C++ vector index
 	    ptau_coeffs_[ic].push_back( luaL_checknumber(L,-1) );
 	    lua_pop(L,1);
 	}
 	lua_pop(L,1);	// pop 'ptau_coefficients'
    }
    UNUSED_VARIABLE(ie_);
    UNUSED_VARIABLE(iv_);
}

VE_Lee::
~VE_Lee()
{
    T_switches_.resize(0);
    ptau_coeffs_.resize(0);
}

double
VE_Lee::
specific_relaxation_time(Gas_data &Q, std::vector<double> &molef)
{
    // Select the ptau coefficients to use from Te
    double Te = Q.T[iTe_];
    size_t iptc;
    for ( iptc=0; iptc<T_switches_.size(); ++iptc ) {
     	if ( Te < T_switches_[iptc] ) break;
    }

    // Calculate tau
    // NOTE: pe_atm should be > 0 as molef[ie_] > 0
    double pe_atm = Q.p_e / PC_P_atm;
    double log_ptau = ptau_coeffs_[iptc][0] * log10(Te) * log10(Te) + ptau_coeffs_[iptc][1] * log10(Te) + ptau_coeffs_[iptc][2];
    double tau = pow( 10.0, log_ptau ) / pe_atm;

    return tau;
}

Relaxation_time* create_new_relaxation_time(lua_State *L, int ip, int iq, int itrans)
{

    // 1. get name of relaxation time type to be created
    string relaxation_time = get_string(L, -1, "model");

    // 2. create the relaxation time instance
    if ( relaxation_time == "Millikan-White" ) {
	return new VT_MillikanWhite(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "PolyFit" ){
        return new VT_PolyFit(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "Landau-Teller-cf" ) {
	return new VT_LandauTeller_cf(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "Thivet-cf" ) {
        return new VT_Thivet_cf(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "SSH-VT" ) {
	return new VT_SSH(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "Millikan-White:HTC" ) {
     	return new VT_MillikanWhite_HTC(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "SSH-VV" ) {
    	return new VV_SSH(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "Landau-Teller-VV" ) {
	return new VV_MTLandauTeller(L, ip, iq, itrans);
    }
    else if ( relaxation_time == "from-eq" ) {
	return new VV_from_eq(L, ip, iq, itrans);
    }
    else if( relaxation_time == "Appleton-Bray:Ion" ) {
	return new ET_AppletonBray_Ion(L, ip, iq);
    }
    else if( relaxation_time == "Appleton-Bray:Neutral" ) {
	return new ET_AppletonBray_Neutral(L, ip, iq);
    }
    else if( relaxation_time == "Appleton-Bray:TwoRangeNeutral" ) {
	return new ET_AppletonBray_TwoRangeNeutral(L, ip, iq);
    }
    else if( relaxation_time == "Abe-ER:Ion" ) {
        return new ER_Abe_Ion(L, ip, iq);
    }
    else if( relaxation_time == "Abe-ER:Neutral" ) {
        return new ER_Abe_Neutral(L, ip, iq);
    }
    else if( relaxation_time == "Lee-VE" ) {
        return new VE_Lee(L, ip, iq);
    }
    // else if( relaxation_time == "VV_Candler" ) {
    // 	return new VV_Candler(L);
    // }
    // else if( relaxation_time == "RT_Parker" ) {
    // 	return new RT_Parker(L);
    // }
    else {
    	cout << "create_new_relaxation_time()\n";
    	cout << "The relaxation_time: " << relaxation_time << " is not known.\n";
	exit(BAD_INPUT_ERROR);
    }
}

void parse_input_for_rts(string cfile, Gas_model &g, lua_State *L)
{
    // Set up a species table
    lua_newtable(L);
    for ( int isp = 0; isp < g.get_number_of_species(); ++isp ) {
	// This species table maps to C++ indices, because
	// it is used to setup the integer maps for
	// the reaction coefficients.
	lua_pushinteger(L, isp);
	lua_setfield(L, -2, g.species_name(isp).c_str());
	// Also add the reverse lookup
	lua_pushinteger(L, isp);
	lua_pushstring(L, g.species_name(isp).c_str());
	lua_settable(L, -3);
    }
    // Plus add a field 'size': no of species
    lua_pushinteger(L, g.get_number_of_species());
    lua_setfield(L, -2, "size");
    lua_setglobal(L, "species");

    // Setup a table of thermal modes
    lua_newtable(L);
    for ( int imode = 0; imode < g.get_number_of_modes(); ++imode ) {
	lua_newtable(L);
	for ( int ic = 0; ic < g.mode_no_components(imode); ++ic ) {
	    lua_pushinteger(L, ic);
	    lua_setfield(L, -2, g.mode_component_name(imode, ic).c_str());
	}
	lua_setfield(L, -2, g.mode_name(imode).c_str());
    }
    lua_setglobal(L, "modes");

    // Setup a table to find index of a given mode name
    lua_newtable(L);
    for ( int imode = 0; imode < g.get_number_of_modes(); ++imode) {
	lua_pushinteger(L, imode);
	lua_setfield(L, -2, g.mode_name(imode).c_str());
    }
    lua_setglobal(L, "mode_idx");

    // Path to reaction parsing script
    char *e3bin = getenv("E3BIN");
    string home;
    if ( e3bin == NULL ) {
	// Assume default location of $HOME/e3bin
	home.append(getenv("HOME")); home.append("/e3bin");
    }
    else {
	home.append(e3bin);
    }
    string script_file(home);
    script_file.append("/energy_exchange_parser.lua");

    if ( luaL_dofile(L, script_file.c_str()) != 0 ) {
	ostringstream ost;
	ost << "create_energy_exchange_update():\n";
	ost << "Error in loading script file: " << script_file << endl;
	ost << lua_tostring(L, -1) << endl;
	input_error(ost);
    }

    // Parse the input file...
    lua_getglobal(L, "main");
    lua_pushstring(L, cfile.c_str());
    if ( lua_pcall(L, 1, 0, 0) != 0 ) {
	ostringstream ost;
	ost << "parse_input_for_rts():\n";
	ost << "Error trying to load/parse energy exchange input file: " << cfile << endl;
	ost << lua_tostring(L, -1) << endl;
	input_error(ost);
    }
}

Relaxation_time* get_rt_from_file(int irt, string cfile, Gas_model &g)
{
    // Setup lua_State for parsing as per energy_exchange_update.cxx
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    parse_input_for_rts(cfile, g, L);

    lua_getglobal(L, "mechs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "get_rt_from_file():\n";
	ost << "Error finding 'mechs' table.\n";
	input_error(ost);
    }

    int nmechs = lua_objlen(L, -1);
    if ( irt >= nmechs ) {
	ostringstream ost;
	ost << "get_rt_from_file():\n";
	ost << "Error: requested relaxation time index " << irt << endl;
	ost << "       is greater than the number of available mechanisms.\n";
	input_error(ost);
    }

    lua_rawgeti(L, -1, irt+1); // user supplies counting from '0'
    int ip = get_int(L, -1, "ip");
    int iq = get_int(L, -1, "iq");
    int itrans = get_int(L, -1, "itrans");

    lua_getfield(L,-1,"relaxation_time");
    Relaxation_time *rt = create_new_relaxation_time(L, ip, iq, itrans);
    lua_close(L);

    return rt;
}

int get_no_rts_from_file(string cfile, Gas_model &g)
{
    // Setup lua_State for parsing as per energy_exchange_update.cxx
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    parse_input_for_rts(cfile, g, L);

    lua_getglobal(L, "mechs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "get_rt_from_file():\n";
	ost << "Error finding 'mechs' table.\n";
	input_error(ost);
    }

    int nmechs = lua_objlen(L, -1);
    lua_close(L);
    return nmechs;
}


