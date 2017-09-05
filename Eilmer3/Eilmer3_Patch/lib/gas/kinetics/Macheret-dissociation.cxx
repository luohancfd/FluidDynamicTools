// Author: Daniel F Potter
// Date: 07-Dec-2009

#include <cmath>
#include <sstream>
#include <iostream>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "Macheret-dissociation.hh"
#include "../models/physical_constants.hh"
#include "../models/chemical-species-library.hh"

using namespace std;

Macheret_dissociation::
Macheret_dissociation(lua_State *L, Gas_model &g, double T_upper, double T_lower)
    : Generalised_Arrhenius(L, g, T_upper, T_lower)
{
    // 0. GA data
    A_ = Generalised_Arrhenius::get_A();
    n_ = Generalised_Arrhenius::get_n();
    E_a_ = Generalised_Arrhenius::get_E_a();
  
    lua_getfield(L, -1, "khigh");
    if ( lua_isnil(L, -1) ) khigh_ =0; 
    else khigh_ = luaL_checknumber(L, -1);
    lua_pop(L, 1);
    cout<<"khigh = "<<khigh_<<endl;
 
    // 1. Vibrating species data
    string v_name = get_string(L,-1,"v_name");
    Diatomic_species * D = get_library_diatom_pointer_from_name( v_name );
    iTv_ = D->get_mode_pointer_from_type("vibration")->get_iT();
    theta_v_ = dynamic_cast<Vibration*>(D->get_mode_pointer_from_type("vibration"))->get_theta();
    theta_d_ = E_a_ / PC_k_SI;	// use reaction energy as dissociation temperature
    double M_v = D->get_M();

    // 2. Colliding species data
    string c_name = get_string(L,-1,"c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    iT_ = X->get_mode_pointer_from_type("translation")->get_iT();
    if ( X->get_type()=="monatomic" ) monatomic_collider_ = true;
    else monatomic_collider_ = false;
    double M_c = X->get_M();

    // 3. Pre-calculate alpha
    if (monatomic_collider_)
       alpha_ = pow( ( M_v*0.5  / ( M_v*0.5 + M_c ) ), 2.0 );
    else
       alpha_ = pow( ( M_v / (M_v + M_c)),2.0);
    cout<<"alpha for MF model"<<alpha_<<endl;

    b_ = 2.0;
    delta_d_ = 3.0*b_*alpha_*alpha_*theta_d_;
    theta_dstar_ = theta_d_ - delta_d_;

    // 3. Set the reaction rate coefficient type
    type_ = "dissociation";

    UNUSED_VARIABLE(ic_);
}

Macheret_dissociation::
Macheret_dissociation(double A, double n, double E_a, double T_upper, double T_lower, string v_name, string c_name, int khigh)
    : Generalised_Arrhenius(A, n, E_a, T_upper, T_lower)
{
    // 0. GA data
    A_ = A;
    n_ = n;
    E_a_ = E_a;

    // 1. Vibrating species data
    Diatomic_species * D = get_library_diatom_pointer_from_name( v_name );
    iTv_ = D->get_mode_pointer_from_type("vibration")->get_iT();
    theta_v_ = dynamic_cast<Vibration*>(D->get_mode_pointer_from_type("vibration"))->get_theta();
    theta_d_ = E_a_ / PC_k_SI;	// use reaction energy as dissociation temperature
    double M_v = D->get_M();

    // 2. Colliding species data
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    iT_ = X->get_mode_pointer_from_type("translation")->get_iT();
    if ( X->get_type()=="monatomic" ) monatomic_collider_ = true;
    else monatomic_collider_ = false;
    double M_c = X->get_M();

    khigh_ = khigh;
    // 3. Pre-calculate alpha
    if (monatomic_collider_)
       alpha_ = pow( ( M_v*0.5  / ( M_v*0.5 + M_c ) ), 2.0 );
    else
       alpha_ = pow( ( M_v / (M_v + M_c)),2.0);
    cout<<"alpha for MF model"<<alpha_<<endl;

    b_ = 2.0;
    delta_d_ = 3.0*b_*alpha_*alpha_*theta_d_;
    theta_dstar_ = theta_d_ - delta_d_;

    // 3. Set the reaction rate coefficient type
    type_ = "dissociation";
}

Macheret_dissociation::
~Macheret_dissociation() {}

int
Macheret_dissociation::
s_eval(const Gas_data &Q)
{
    // 1. Eval k_ as Arrhenius rate
    double T = Q.T[iT_];
    if ( T > T_upper_ )
	T = T_upper_;
    if ( T < T_lower_ )
	T = T_lower_;
    Generalised_Arrhenius::eval_from_T(T); // This call sets 'k_'
                                           // 'k_' is then modified by 'Z' below

    // 2. Get the appropriate temperatures
    double Tv = Q.T[iTv_];
    T = Q.T[iT_]; // Reset T to T_trans even if it is above limits.
                  // because now we are computing the nonequilibrium factor.
    double Ta = alpha_ * Tv + ( 1.0 - alpha_ ) * T;
    fac_ = (1.0-exp(-theta_v_/Tv)) / (1.0-exp(-theta_v_/T));

    // 3. Calculate nonequilibrium factor
    double L;
    if ( monatomic_collider_ ) {
      //  L = 9.0 * sqrt( M_PI * ( 1.0 - alpha_ ) ) / 64.0 * pow( T / theta_d_, 1.0 - n_ ) * \
      //	    ( 1.0 + 5.0 * ( 1.0 - alpha_) * T / ( 2.0 * theta_d_ ));
        L = sqrt(1.0 - alpha_)/pow(M_PI,1.5)*(1.0+5.0*(1.0-alpha_)*T/2.0/theta_dstar_) * pow( T / theta_d_, 1.0 - n_ ) * \
	    sqrt(theta_d_/theta_dstar_)*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
    }
    else {
	//L = 2.0 * ( 1.0 - alpha_ ) / ( M_PI * M_PI * pow( alpha_, 0.75 ) ) * \
	    //pow( T / theta_d_, 1.5 - n_ ) * ( 1.0 + 7.0 * ( 1.0 - alpha_) * ( 1.0 + sqrt(alpha_) ) * T / ( 2.0 * theta_d_ ));
        L = 2.0 * ( 1.0 - alpha_ ) / ( M_PI * M_PI * pow( alpha_, 0.75 ) ) * theta_d_/theta_dstar_ * \
            pow( T / theta_d_, 1.5 - n_ ) * ( 1.0 + 7.0 * ( 1.0 - alpha_) * ( 1.0 + sqrt(alpha_) ) * T / ( 2.0 * theta_dstar_ ));
        L = L*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
    }

    //double Z = ( ( 1.0 - exp( - theta_v_ / Tv ) ) / ( 1.0 - exp( - theta_v_ / T ) ) ) * ( 1.0 - L ) * \
                    //exp( - theta_d_ * ( 1.0 / Tv - 1.0 / T ) ) + L * exp( - theta_d_ * ( 1.0 / Ta - 1.0 / T ) );
    double Zl,Zh;
    Zl = L * exp(-theta_d_*(1.0/Ta-1.0/T) + delta_d_*(1.0/Ta - 1.0/T));
    Zh = (1.0 - L)*fac_*exp(-(theta_d_-khigh_*delta_d_)*(1.0/Tv-1/T));
    
    // Augment Arrhenius rate by noneq factor
    k_ *= Zl+Zh;

    return SUCCESS;
}

Reaction_rate_coefficient* create_Macheret_dissociation_coefficient(lua_State *L, Gas_model &g, double T_upper, double T_lower)
{
    return new Macheret_dissociation(L, g, T_upper, T_lower);
}
