// Author: Daniel F Potter
// Date: 07-Apr-2010
// Place: Dutton Park, QLD, Oz
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"

#include "../models/chemical-species-library.hh"

#include "coupling-component.hh"

using namespace std;

static const double min_mass_frac = 1.0e-50;

/***************************** Coupling_component ****************************/

Coupling_component::Coupling_component( lua_State *L, Reaction *r, string type, string mode, int idc )
: type_( type ), mode_( mode ), idc_( idc )
{
    string species_name = get_string(L,-1,"species");
    Chemical_species * X = get_library_species_pointer_from_name( species_name );
    isp_ = X->get_isp();
    m_ = X->get_M() / PC_Avogadro;
    nu_ = r->get_nu(isp_);

    // Search for the corresponding energy modes
    bool found = false;
    for ( int itm=0; itm<X->get_n_modes(); ++itm ) {
    	Species_energy_mode * sem = X->get_mode_pointer(itm);
    	if ( sem->get_type()==mode_ ) {
    	    found = true;
    	    sems_.push_back( sem );
    	}
    }
    if ( !found ) {
    	cout << "Coupling_component::Coupling_component()" << endl
    	     << "mode: " << mode_ << " not found for species: "
    	     << X->get_name() << endl
    	     << "Exiting." << endl;
    	exit( BAD_INPUT_ERROR );
    }

    imode_ = sems_[0]->get_iT();
}

Coupling_component::Coupling_component( const Coupling_component &c )
: type_( c.type_ ), mode_( c.mode_ ), idc_( c.idc_ ), isp_( c.isp_ ), m_( c.m_ ),
  nu_( c.nu_ ), sems_( c.sems_ ), imode_( c.imode_ ) {}

Coupling_component::~Coupling_component() {}

Coupling_component*
Coupling_component::clone() const
{
    return 0;
}

/***************************** Creation function ****************************/

void create_Coupling_components_for_reaction( lua_State * L, Reaction * r, int ir, vector<Coupling_component*> & ccs )
{
    lua_getfield( L, -1, "chemistry_energy_coupling" );
    if ( !lua_isnil(L,-1) ) {
    	for ( size_t i=1; i<=lua_objlen(L,-1); ++i ) {
    	    lua_rawgeti(L, -1, i);
    	    string mode = get_string(L,-1,"mode");
    	    string model = get_string(L,-1,"model");
    	    if ( mode=="vibration" ) {
    	    	if ( model=="Simple" ) {
    	    	    ccs.push_back( new Simple_dissociation_component( L, r, 2*ir ) );
    	    	    ccs.push_back( new Simple_recombination_component( L, r, 2*ir+1 ) );
    	    	}
    	    	else if ( model=="TreanorMarrone" ) {
    	    	    ccs.push_back( new TreanorMarrone_dissociation_component( L, r, 2*ir ) );
    	    	    ccs.push_back( new TreanorMarrone_recombination_component( L, r, 2*ir+1 ) );
    	    	}
    	    	else if ( model=="Park" ) {
    	    	    ccs.push_back( new Park_dissociation_component( L, r, 2*ir ) );
    	    	    ccs.push_back( new Park_recombination_component( L, r, 2*ir+1 ) );
    	    	}
    	    	else if ( model=="Macheret" ) {
    	    	    ccs.push_back( new Macheret_dissociation_component( L, r, 2*ir ) );
    	    	    ccs.push_back( new Macheret_recombination_component( L, r, 2*ir+1 ) );
    	    	}
		else if ( model=="Knab_et_al" ) {
		    ccs.push_back( new Knab_vanishing_component( L, r, 2*ir ) );
		    ccs.push_back( new Knab_appearing_component( L, r, 2*ir+1 ) );
		}
    	    	else {
		    ostringstream oss;
		    oss << "create_Coupling_components_for_reaction()" << endl
		        << "The requested model: " << model
		        << " is not available for chemistry-vibration coupling." << endl;
		    input_error( oss );
		}
    	    }
    	    else if ( mode=="translation" ) {
    	    	if ( model=="electron impact ionization" ) {
    	    	    ccs.push_back( new Electron_impact_ionization_component( L, r, 2*ir ) );
    	    	    // ccs.push_back( new EII_recombination_component( L, r, 2*ir +1 ) );
    	    	}
    	    	else if ( model=="associative ionization" ) {
    	    	    ccs.push_back( new Associative_ionization_component( L, r, 2*ir ) );
    	    	    // ccs.push_back( new AI_recombination_component( L, r, 2*ir + 1 ) );
    	    	}
    	    }
    	    else {
    	    	ostringstream oss;
    	    	oss << "create_Coupling_components_for_reaction()" << endl
    	    	    << "The reqested mode: " << mode
    	    	    << " does not have any coupling models available" << endl;
    	    	input_error( oss );
    	    }
    	    lua_pop(L,1);	// pop the coupling component object
    	}
    }
    lua_pop(L,1);	// pop the 'chemistry_energy_coupling' field

    return;
}

/****************** Simple_dissociation_component ********************/

Simple_dissociation_component::
Simple_dissociation_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Simple_dissociation_component","vibration",idc)
{
    D_hat_ = get_positive_number(L,-1,"T_D") * PC_k_SI;
}

Simple_dissociation_component::
Simple_dissociation_component( const Simple_dissociation_component &c )
: Coupling_component( c ), D_hat_( c.D_hat_ ) {}

Simple_dissociation_component::
~Simple_dissociation_component()
{}

Simple_dissociation_component*
Simple_dissociation_component::clone() const
{
    return new Simple_dissociation_component(*this);
}

double
Simple_dissociation_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Set the vanishing vibration energy
    double e_va = D_hat_;

    // 2. Calculate delta_E
    double delta_N = nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_va - e_old_ ) * delta_N;

    // cout << "Simple_dissociation_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Simple_dissociation_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Set the vanishing vibration energy
    double e_va = D_hat_;

    // 2. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = ( e_va - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Simple_dissociation_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/****************** TreanorMarrone_dissociation_component ********************/

TreanorMarrone_dissociation_component::
TreanorMarrone_dissociation_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"TreanorMarrone_dissociation_component","vibration",idc)
{
    U_ = get_positive_number(L,-1,"U");
}

TreanorMarrone_dissociation_component::
TreanorMarrone_dissociation_component( const TreanorMarrone_dissociation_component &c )
: Coupling_component( c ), U_( c.U_ ) {}

TreanorMarrone_dissociation_component::
~TreanorMarrone_dissociation_component()
{}

TreanorMarrone_dissociation_component*
TreanorMarrone_dissociation_component::clone() const
{
    return new TreanorMarrone_dissociation_component(*this);
}

double
TreanorMarrone_dissociation_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[imode_];

    // 1. Calculate pseudo-temperature gamma
    double gamma_inv = 1.0/Tv - 1.0/T - 1.0/U_;
    double gamma = 1.0 / gamma_inv;

    // 2. Evaluate the vibrational energy at T=gamma
    double e_va = 0.0;
    for ( size_t i = 0; i < sems_.size(); ++i )
    	e_va += sems_[i]->eval_energy_from_T( gamma ) * m_;

    // 3. Calculate delta_E
    double delta_N = nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_va - e_old_ ) * delta_N;

    //    cout << "TreanorMarrone_dissociation_component::compute_contribution()" << endl
    //	 << "gamma = " << gamma << ", Tv = " << Tv << ", T = " << T << endl
    //	 << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //	 << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
TreanorMarrone_dissociation_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[imode_];

    // 1. Calculate pseudo-temperature gamma
    double gamma_inv = 1.0/Tv - 1.0/T - 1.0/U_;
    double gamma = 1.0 / gamma_inv;

    // 2. Evaluate the vibrational energy at T=gamma
    double e_va = 0.0;
    for ( size_t i = 0; i < sems_.size(); ++i )
    	e_va += sems_[i]->eval_energy_from_T( gamma ) * m_;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = ( e_va - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "TreanorMarrone_dissociation_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/****************** Park_dissociation_component ********************/

Park_dissociation_component::
Park_dissociation_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Park_dissociation_component","vibration",idc)
{
    n_ = get_number(L,-1,"n");
    D_ = get_positive_number(L,-1,"T_d") * PC_k_SI;
    s_v_ = get_positive_number(L,-1,"s_v");
}

Park_dissociation_component::
Park_dissociation_component( const Park_dissociation_component &c )
: Coupling_component( c ), n_( c.n_ ), D_( c.D_ ), s_v_( c.s_v_ ) {}

Park_dissociation_component::
~Park_dissociation_component()
{}

Park_dissociation_component*
Park_dissociation_component::clone() const
{
    return new Park_dissociation_component(*this);
}

double
Park_dissociation_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[imode_];

    // 2. Evaluate the vanishing vibrational energy (J/particle)
    //    NOTE: the '+ e_old_' is required, see AMOD TN 3.2
    double e_va = (1.0 - s_v_) * ( n_ * PC_k_SI * Tv + D_ * pow( Tv/T, s_v_ ) ); // + e_old_;

    // 3. Calculate delta_E
    double delta_N = nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_va - e_old_ ) * delta_N;

    // cout << "Park_dissociation_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Park_dissociation_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[imode_];

    // 2. Evaluate the vanishing vibrational energy (J/particle)
    // 2a. Energy addition
    double e_va = (1.0 - s_v_) * ( n_ * PC_k_SI * Tv + D_ * pow( Tv/T, s_v_ ) );
    // 2b. Average energy
    // for ( size_t i=0; i<sems_.size(); ++i )
    //	e_va += sems_[i]->eval_energy_from_T(Tv) * m_;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = ( e_va - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Park_dissociation_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/****************** Macheret_dissociation_component ********************/

Macheret_dissociation_component::
Macheret_dissociation_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Macheret_dissociation_component","vibration",idc)
{
    // 0. Check for just one mode
    if ( sems_.size() != 1 ) {
    	ostringstream oss;
    	oss << "Macheret_dissociation_component::Macheret_dissociation_component" << endl
    	    << sems_.size() << " energy modes found, just one expcted." << endl;
    	input_error( oss );
    }

    // 1. GA data
    A_ = get_number(L,-1,"A") * 1.0e-6;		// Convert cm**3/mole-s -> m**3/mole-s
    n_ = get_number(L,-1,"n");
    theta_d_ = get_number(L,-1,"T_d");

    // 2. Dissociating species data
    theta_v_ = dynamic_cast<Vibration*>(sems_[0])->get_theta();
    double M_v = m_ * PC_Avogadro;

    // 2. Colliding species data
    string c_name = get_string(L,-1,"c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    if ( X->get_type()=="monatomic" ) monatomic_collider_ = true;
    else monatomic_collider_ = false;
    double M_c = X->get_M();

    lua_getfield(L, -1, "khigh");
    if ( lua_isnil(L, -1) ) khigh_ = 0.0;
    else khigh_ = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    // 3. Pre-calculate alpha
    if (monatomic_collider_)
       alpha_ = pow( ( M_v*0.5  / ( M_v*0.5 + M_c ) ), 2.0 );
    else
       alpha_ = pow( ( M_v / (M_v + M_c)),2.0);

    b_ = 2.0;
    delta_d_ = 3.0*b_*alpha_*alpha_*theta_d_;
    theta_dstar_ = theta_d_ - delta_d_;
}

Macheret_dissociation_component::
Macheret_dissociation_component( const Macheret_dissociation_component &c )
: Coupling_component( c ), A_( c.A_ ), n_( c.n_ ), theta_d_( c.theta_d_ ),
theta_v_( c.theta_v_ ), alpha_( c.alpha_ ), monatomic_collider_( c.monatomic_collider_ ),khigh_(c.khigh_),
b_(c.b_), delta_d_(c.delta_d_), theta_dstar_(c.theta_dstar_){}

Macheret_dissociation_component::
~Macheret_dissociation_component()
{}

Macheret_dissociation_component*
Macheret_dissociation_component::clone() const
{
    return new Macheret_dissociation_component(*this);
}

double
Macheret_dissociation_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Get the appropriate temperatures
    double Tv = Q.T[imode_];
    double T = Q.T[0];
    double Ta = alpha_ * Tv + ( 1.0 - alpha_ ) * T;
    fac_ = (1.0-exp(-theta_v_/Tv)) / (1.0-exp(-theta_v_/T));

    // 2. Evaluate the vanishing vibrational energy (J/particle)
    // 2a. Eval Arrhenius rate
    double k_f = A_ * pow(T, n_) * exp(-theta_d_ / T );

    // 2b. Calculate nonequilibrium factor
    double L;
    if ( monatomic_collider_ ) {
        //L = 9.0 * sqrt( M_PI * ( 1.0 - alpha_ ) ) / 64.0 * pow( T / theta_d_, 1.0 - n_ ) * \
	    //( 1.0 + 5.0 * ( 1.0 - alpha_) * T / ( 2.0 * theta_d_ ));
        L = sqrt(1.0 - alpha_)/ pow(M_PI,1.5)*(1.0+5.0*(1.0-alpha_)*T/2.0/theta_dstar_) * pow( T / theta_d_, 1.0 - n_ ) * \
	    sqrt(theta_d_/theta_dstar_);
        L = L*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
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
    double Zl, Zh;
    Zl = L * exp(-theta_d_*(1.0/Ta-1.0/T) + delta_d_*(1.0/Ta - 1.0/T));
    Zh = (1.0 - L)*fac_*exp(-(theta_d_-khigh_*delta_d_)*(1.0/Tv-1.0/T) );

    // 2d. Eval low and high reaction rates
    //double k_f_l = ( 1.0 - L ) * A_ * pow(T,n_) * exp( - theta_d_ / Ta );
    //double k_f_h =  ( 1.0 - exp( - theta_v_ / Tv ) ) / ( 1.0 - exp( - theta_v_ / T ) ) * ( L * A_ * pow(T,n_) * exp( - theta_d_ / Tv ) );
    double k_f_l = Zl*k_f;
    double k_f_h =  Zh*k_f;

    // 2c. Augment Arrhenius rate by noneq factor
    k_f = k_f_l + k_f_h;

    // 2e. Eval the vanishing energy per particle (thus we use k instead of R_i)
    //double e_va = ( alpha_ * PC_k_SI * theta_d_ * pow( T / Ta, 2 ) * k_f_l + PC_k_SI * theta_d_ * k_f_h ) / k_f;
    double e_va = ( alpha_ * PC_k_SI * theta_dstar_ * pow( Tv / Ta, 2 ) * k_f_l + PC_k_SI * (theta_d_ - khigh_*delta_d_) * k_f_h ) / k_f;

    // 3. Calculate delta_E
    double delta_N = nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_va - e_old_ ) * delta_N;

    // cout << "Macheret_dissociation_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Macheret_dissociation_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Get the appropriate temperatures
    double Tv = Q.T[imode_];
    double T = Q.T[0];
    double Ta = alpha_ * Tv + ( 1.0 - alpha_ ) * T;
    fac_ = (1.0-exp(-theta_v_/Tv)) / (1.0-exp(-theta_v_/T));

    // 2. Evaluate the vanishing vibrational energy (J/particle)
    // 2a. Eval Arrhenius rate
    double k_f = A_ * pow(T, n_) * exp(-theta_d_ / T );

    // 2b. Calculate nonequilibrium factor
    double L;
    if ( monatomic_collider_ ) {
        //L = 9.0 * sqrt( M_PI * ( 1.0 - alpha_ ) ) / 64.0 * pow( T / theta_d_, 1.0 - n_ ) * \
	    //( 1.0 + 5.0 * ( 1.0 - alpha_) * T / ( 2.0 * theta_d_ ));
        L = sqrt(1.0 - alpha_ ) / pow(M_PI,1.5) *(1.0+5.0*(1.0-alpha_)*T/2.0/theta_dstar_) * pow( T / theta_d_, 1.0 - n_ );
        L = L*sqrt(theta_d_/theta_dstar_)*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
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
    Zh = (1.0 - L)*fac_*exp(-(theta_d_-khigh_*delta_d_)*(1.0/Tv-1.0/T) );

    // 2d. Eval low and high reaction rates
    //double k_f_l = ( 1.0 - L ) * A_ * pow(T,n_) * exp( - theta_d_ / Ta );
    //double k_f_h =  ( 1.0 - exp( - theta_v_ / Tv ) ) / ( 1.0 - exp( - theta_v_ / T ) ) * ( L * A_ * pow(T,n_) * exp( - theta_d_ / Tv ) );
    double k_f_l = Zl*k_f;
    double k_f_h =  Zh*k_f;

    // 2c. Augment Arrhenius rate by noneq factor
    k_f = k_f_l + k_f_h;

    // 2e. Eval the vanishing energy per particle (thus we use k instead of R_i)
    //double e_va = ( alpha_ * PC_k_SI * theta_d_ * pow( T / Ta, 2 ) * k_f_l + PC_k_SI * theta_d_ * k_f_h ) / k_f;
    double e_va = ( alpha_ * PC_k_SI * theta_dstar_ * pow( Tv / Ta, 2 ) * k_f_l + PC_k_SI * (theta_d_-khigh_*delta_d_) * k_f_h ) / k_f;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = ( e_va - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Macheret_dissociation_component::specific_compute_source_term()" << endl
        // << "dcdt[idc_] = " << dcdt[idc_] << endl
        // << "dEdt = " << dEdt << endl
        // << "k_f = " << k_f << endl
        // << "L = " << L << endl
        // << "Z = " << Z << endl
        // << "k_f_l = " << k_f_l << endl
        // << "k_f_h = " << k_f_h << endl
        // << "Tv = " << Tv << endl
        // << "Ta = " << Ta << endl;
    //
    // cout << "Macheret_dissociation_component::specific_compute_source_term()" << endl
         // << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl
         // << "e_va = " << e_va << ", e_average = " << sems_[0]->eval_energy_from_T(Tv)*m_ << endl;

    return dEdt;
}

/****************** Simple_recombination_component ********************/

Simple_recombination_component::
Simple_recombination_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Simple_recombination_component","vibration",idc)
{
    D_hat_ = get_positive_number(L,-1,"T_D") * PC_k_SI;

}

Simple_recombination_component::
Simple_recombination_component( const Simple_recombination_component &c )
: Coupling_component( c ), D_hat_( c.D_hat_ ) {}

Simple_recombination_component::
~Simple_recombination_component()
{}

Simple_recombination_component*
Simple_recombination_component::clone() const
{
    return new Simple_recombination_component(*this);
}

double
Simple_recombination_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Set the appearing vibration energy
    double e_app = D_hat_;

    // 2. Calculate delta_E
    double delta_N = - nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_app - e_old_ ) * delta_N;

    // cout << "Simple_recombination_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Simple_recombination_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Set the appearing vibration energy
    double e_app = D_hat_;

    // 2. Calculate delta_E
    double dEdt = - ( e_app - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Simple_recombination_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}


/****************** TreanorMarrone_recombination_component ********************/

TreanorMarrone_recombination_component::
TreanorMarrone_recombination_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"TreanorMarrone_recombination_component","vibration",idc)
{
    U_ = get_positive_number(L,-1,"U");
}

TreanorMarrone_recombination_component::
TreanorMarrone_recombination_component( const TreanorMarrone_recombination_component &c )
: Coupling_component( c ), U_( c.U_ ) {}

TreanorMarrone_recombination_component::
~TreanorMarrone_recombination_component()
{}

TreanorMarrone_recombination_component*
TreanorMarrone_recombination_component::clone() const
{
    return new TreanorMarrone_recombination_component(*this);
}

double
TreanorMarrone_recombination_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Calculate pseudo-temperature gamma
    //    NOTE: Since Tv=T for the recombination reaction, gamma = -U
    //          The negative temperature indicates formation in high lying vibrational levels
    double gamma = - U_;

    // 2. Evaluate the vibrational energy at T=gamma
    double e_app = 0.0;
    for ( size_t i = 0; i < sems_.size(); ++i )
    	e_app += sems_[i]->eval_energy_from_T( gamma ) * m_;

    // 3. Calculate delta_E
    double delta_N = - nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_app - e_old_ ) * delta_N;

    // cout << "TreanorMarrone_recombination_component::compute_contribution()" << endl
    //      << "gamma = " << gamma << ", Tv = " << Q.T[imode_] << ", T = " << Q.T[0] << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
TreanorMarrone_recombination_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Calculate pseudo-temperature gamma
    //    NOTE: Since Tv=T for the recombination reaction, gamma = -U
    //          The negative temperature indicates formation in high lying vibrational levels
    double gamma = - U_;

    // 2. Evaluate the vibrational energy at T=gamma
    double e_app = 0.0;
    for ( size_t i = 0; i < sems_.size(); ++i )
    	e_app += sems_[i]->eval_energy_from_T( gamma ) * m_;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = - ( e_app - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "TreanorMarrone_recombination_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/****************** Park_recombination_component ********************/

Park_recombination_component::
Park_recombination_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Park_recombination_component","vibration",idc)
{
    n_ = get_number(L,-1,"n");
    D_ = get_positive_number(L,-1,"T_d") * PC_k_SI;
    s_v_ = get_positive_number(L,-1,"s_v");
}

Park_recombination_component::
Park_recombination_component( const Park_recombination_component &c )
: Coupling_component( c ), n_( c.n_ ), D_( c.D_ ), s_v_( c.s_v_ ) {}

Park_recombination_component::
~Park_recombination_component()
{}

Park_recombination_component*
Park_recombination_component::clone() const
{
    return new Park_recombination_component(*this);
}

double
Park_recombination_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    //    NOTE: Tv=T for the recombination reaction
    double T = Q.T[0];
    double Tv = T;

    // 2. Evaluate the appearing vibrational energy (J/particle)
    //    NOTE: The '+ e_old_' is required, see AMOD TN 3.2
    //          Omitting 'pow( Tv/T, s_v_ )' which will be unity
    double e_app = (1.0 - s_v_) * ( n_ * PC_k_SI * Tv + D_ ); //  + e_old_;

    // 3. Calculate delta_E
    double delta_N = - nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_app - e_old_ ) * delta_N;

    // cout << "Park_recombination_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Park_recombination_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    //    NOTE: Tv=T for the recombination reaction
    double T = Q.T[0];
    double Tv = T;

    // 2. Evaluate the appearing vibrational energy (J/particle)
    double e_app = (1.0 - s_v_) * ( n_ * PC_k_SI * Tv + D_ );
    // 2b. Average energy
    // for ( size_t i=0; i<sems_.size(); ++i )
    // 	e_app += sems_[i]->eval_energy_from_T(Tv) * m_;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = - ( e_app - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Park_recombination_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/****************** Macheret_recombination_component ********************/

Macheret_recombination_component::
Macheret_recombination_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Macheret_recombination_component","vibration",idc)
{
    // 0. Check for just one mode
    if ( sems_.size()!=1 ) {
    	ostringstream oss;
    	oss << "Macheret_dissociation_component::Macheret_dissociation_component" << endl
    	    << sems_.size() << " energy modes found, just one expcted." << endl;
    	input_error( oss );
    }

    // 1. GA data
    n_ = get_number(L,-1,"n");
    theta_d_ = get_number(L,-1,"T_d");

    // 2. Dissociating species data
    double M_v = m_ * PC_Avogadro;

    // 2. Colliding species data
    string c_name = get_string(L,-1,"c_name");
    Chemical_species * X = get_library_species_pointer_from_name( c_name );
    if ( X->get_type()=="monatomic" ) monatomic_collider_ = true;
    else monatomic_collider_ = false;
    double M_c = X->get_M();

    lua_getfield(L, -1, "khigh");
    if ( lua_isnil(L, -1) ) khigh_ =0;
    else khigh_ = luaL_checknumber(L, -1);
    lua_pop(L, 1);

    // 3. Pre-calculate alpha
    if (monatomic_collider_)
       alpha_ = pow( ( M_v*0.5  / ( M_v*0.5 + M_c ) ), 2.0 );
    else
       alpha_ = pow( ( M_v / (M_v + M_c)),2.0);

    b_ = 2.0;
    delta_d_ = 3.0*b_*alpha_*alpha_*theta_d_;
    theta_dstar_ = theta_d_ - delta_d_;
}

Macheret_recombination_component::
Macheret_recombination_component( const Macheret_recombination_component &c )
: Coupling_component( c ), n_( c.n_ ), theta_d_( c.theta_d_ ), alpha_( c.alpha_ ),
monatomic_collider_( c.monatomic_collider_ ),khigh_(c.khigh_),
b_(c.b_), delta_d_(c.delta_d_), theta_dstar_(c.theta_dstar_){}

Macheret_recombination_component::
~Macheret_recombination_component()
{}

Macheret_recombination_component*
Macheret_recombination_component::clone() const
{
    return new Macheret_recombination_component(*this);
}

double
Macheret_recombination_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: Tv=T for the recombination reaction
    double T = Q.T[0];

    // 2. Evaluate the appearing vibrational energy (J/particle)
    double L;
    if ( monatomic_collider_ ) {
        L = sqrt(1.0 - alpha_ )/pow(M_PI,1.5) *(1.0+5.0*(1.0-alpha_)*T/2.0/theta_dstar_) * pow( T / theta_d_, 1.0 - n_ ) * \
	    sqrt(theta_d_/theta_dstar_)*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
    }
    else {
        L = 2.0 * ( 1.0 - alpha_ ) / ( M_PI * M_PI * pow( alpha_, 0.75 ) ) * theta_d_/theta_dstar_ * \
            pow( T / theta_d_, 1.5 - n_ ) * ( 1.0 + 7.0 * ( 1.0 - alpha_) * ( 1.0 + sqrt(alpha_) ) * T / ( 2.0 * theta_dstar_ ));
        L = L*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
    }

    // Ta = T
    double e_app = PC_k_SI * (theta_dstar_*(alpha_)*(1.0-L) + (theta_d_ - delta_d_*khigh_)*L);

    // 3. Calculate delta_E
    double delta_N = - nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_app - e_old_ ) * delta_N;

    // cout << "Macheret_recombination_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Macheret_recombination_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: Tv=T for the recombination reaction
    double T = Q.T[0];

    // 2. Evaluate the appearing vibrational energy (J/particle)
    double L;
    if ( monatomic_collider_ ) {
        L = 9.0 * sqrt( M_PI * ( 1.0 - alpha_ ) ) / 64.0 *(1.0+5.0*(1.0-alpha_)*T/2.0/theta_dstar_) * pow( T / theta_d_, 1.0 - n_ ) * \
	    sqrt(theta_d_/theta_dstar_)*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
    }
    else {
        L = 2.0 * ( 1.0 - alpha_ ) / ( M_PI * M_PI * pow( alpha_, 0.75 ) ) * theta_d_/theta_dstar_ * \
            pow( T / theta_d_, 1.5 - n_ ) * ( 1.0 + 7.0 * ( 1.0 - alpha_) * ( 1.0 + sqrt(alpha_) ) * T / ( 2.0 * theta_dstar_ ));
        L = L*sqrt(12.0*M_PI*b_*alpha_*(1.0-alpha_)*theta_d_/T);
    }
    double e_app = PC_k_SI * (theta_dstar_*(alpha_)*(1.0-L) + (theta_d_ - delta_d_*khigh_)*L);

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = - ( e_app - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Macheret_recombination_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl
    //      << "e_app = " << e_app << ", e_average = " << sems_[0]->eval_energy_from_T(Q.T[1])*m_ << endl;

    return dEdt;
}

/********************* Electron_impact_ionization_component *******************/

Electron_impact_ionization_component::
Electron_impact_ionization_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Electron_impact_ionization_component","translation",idc)
{
    I_ = get_positive_number(L,-1,"T_I") * PC_k_SI;
}

Electron_impact_ionization_component::
Electron_impact_ionization_component( const Electron_impact_ionization_component &c )
: Coupling_component( c ), I_( c.I_ ) {}

Electron_impact_ionization_component::
~Electron_impact_ionization_component()
{}

Electron_impact_ionization_component*
Electron_impact_ionization_component::clone() const
{
    return new Electron_impact_ionization_component(*this);
}

double
Electron_impact_ionization_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Calculate delta_E
    double delta_N = nu_ * ( delta_c[idc_] - delta_c[idc_+1] ) * PC_Avogadro;
    double delta_E = - ( I_ - e_old_ ) * delta_N;

    // cout << "Electron_impact_ionization_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Electron_impact_ionization_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Calculate dEdt:
    // NOTE: this is a special case as it represents energy lost in the reaction
    //       so we don't have to include e_old_
    double dEdt = - I_ * nu_ * ( dcdt[idc_] - dcdt[idc_+1] ) * PC_Avogadro;

    // cout << "Electron_impact_ionization_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/********************* Associative ionization *******************/

Associative_ionization_component::
Associative_ionization_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Associative_ionization_component","translation",idc)
{
    alpha_ = get_positive_number(L,-1,"alpha");
}

Associative_ionization_component::
Associative_ionization_component( const Associative_ionization_component &c )
: Coupling_component( c ), alpha_( c.alpha_ ) {}

Associative_ionization_component::
~Associative_ionization_component()
{}

Associative_ionization_component*
Associative_ionization_component::clone() const
{
    return new Associative_ionization_component(*this);
}

double
Associative_ionization_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Pull out the translational and free-electron temperatures
    double T = Q.T[0];
    double Te = Q.T[imode_];

    // 2. Model the free electrons as being formed at some temperature between
    //    the free-electron and translation temperatures
    double T_star =  Te + (T - Te)*alpha_;
    double e_app = 1.5 * T_star * PC_k_SI;

    // 1. Calculate delta_E
    double delta_N = nu_ * ( delta_c[idc_] - delta_c[idc_+1] ) * PC_Avogadro;
    double delta_E = ( e_app - e_old_ ) * delta_N;

    // cout << "Associative_ionization_component::compute_contribution()" << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl
    //      << "e_app = " << e_app << ", e_old_ = " << e_old_ << endl;

    return delta_E;
}

double
Associative_ionization_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 1. Pull out the translational and free-electron temperatures
    double T = Q.T[0];
    double Te = Q.T[imode_];

    // 2. Model the free electrons as being formed at some temperature between
    //    the free-electron and translation temperatures
    double T_star =  Te + (T - Te)*alpha_;
    double e_app = 1.5 * T_star * PC_k_SI;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = ( e_app - e_old_ ) * nu_ * ( dcdt[idc_] - dcdt[idc_+1] ) * PC_Avogadro;

    // cout << "Associative_ionization_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;


    return dEdt;
}

/********************* Electron_impact_ionization_component *******************/

EII_recombination_component::
EII_recombination_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"EII_recombination_component","translation",idc)
{}

// EII_recombination_component::
// EII_recombination_component( const EII_recombination_component &c )
// : Coupling_component( c ) {}

EII_recombination_component::
~EII_recombination_component()
{}

EII_recombination_component*
EII_recombination_component::clone() const
{
    return new EII_recombination_component(*this);
}

// NOTE: forward and backwards components are the same, so the recombination
//       energy is done in the ionisation function

double
EII_recombination_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    return 0.0;
}

double
EII_recombination_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    // Assume recombination at average electron energy

    return 0.0;
}

/********************* Associative ionization *******************/

AI_recombination_component::
AI_recombination_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"AI_recombination_component","translation",idc)
{}

// AI_recombination_component::
// AI_recombination_component( const AI_recombination_component &c )
// : Coupling_component( c )
// {}

AI_recombination_component::
~AI_recombination_component()
{}

AI_recombination_component*
AI_recombination_component::clone() const
{
    return new AI_recombination_component(*this);
}

// NOTE: forward and backwards components are the same, so the recombination
//       energy is done in the ionisation function

double
AI_recombination_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    return 0.0;
}

double
AI_recombination_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    return 0.0;
}

/****************** Knab components ********************/

double calculate_Knab_energy( std::vector<Species_energy_mode*> &sems,
			      double U0, double U1, double alpha, double A_var, double T, double Tv)
{
    // 1. Calculate pseudo-temperatures
    double U = U0 + U1*T;
    double gamma, T0, T_ast;
    if ( U > 0.0 ) {
	gamma = 1.0 / (1.0/Tv - 1.0/T - 1.0/U);
	T0 = 1.0 / (1.0/Tv - 1.0/U);
	T_ast = 1.0 / (1.0/T - 1.0/U);
    }
    else {
	// This is the special case of U=inf
	gamma = 1.0 / (1.0/Tv - 1.0/T);
	T0 = Tv;
	T_ast = T;
    }

    // 2. Calculate energies and partition functions
    double L_d_T0 = 0.0;
    double L_a_gamma = 0.0, L_a_T0 = 0.0;
    double Q_d_T = 1.0, Q_d_Tv = 1.0, Q_d_T0 = 1.0, Q_d_T_ast = 1.0;
    double Q_a_gamma = 1.0, Q_a_U = 1.0, Q_a_T0 = 1.0, Q_a_T_ast = 1.0;

    for ( size_t i=0; i<sems.size(); ++i ) {
	L_d_T0 += sems[i]->eval_energy_from_T(T0);
	L_a_gamma += sems[i]->eval_energy_from_T(gamma,alpha*A_var);
	L_a_T0 += sems[i]->eval_energy_from_T(T0,alpha*A_var);
	Q_d_T *= sems[i]->eval_Q_from_T(T);
	Q_d_Tv *= sems[i]->eval_Q_from_T(Tv);
	Q_d_T0 *= sems[i]->eval_Q_from_T(T0);
	Q_d_T_ast *= sems[i]->eval_Q_from_T(T_ast);
	Q_a_gamma *= sems[i]->eval_Q_from_T(gamma,alpha*A_var);
	Q_a_T0 *= sems[i]->eval_Q_from_T(T0,alpha*A_var);
	Q_a_T_ast *= sems[i]->eval_Q_from_T(T_ast,alpha*A_var);
	if ( U > 0.0 ) {
	    Q_a_U *= sems[i]->eval_Q_from_T(-U,alpha*A_var);
	}
	else {
	    // Special case of U = inf
	    // I showed this numerically by plotting
	    // partition function of truncated harmonic
	    // osicallator for large values of T
	    Q_a_U *= alpha*A_var/sems[i]->get_theta();
	}
    }

    // 3. Calculate the energy of the vanishing or appearing molecules
    double tmpA = exp( - alpha * A_var / T );
    double tmpB = tmpA * Q_a_gamma * L_a_gamma + Q_d_T0 * L_d_T0 - Q_a_T0 * L_a_T0;
    double tmpC = tmpA * Q_a_gamma + Q_d_T0 - Q_a_T0;

    double val = tmpB/tmpC;

    // At small values of temperature, the values of the partition functions
    // are essentially the same (to machine precision). At the same time,
    // the exponential term is essentially zero. This leads to a divide
    // by zero problem.
    //
    // Instead, we'll set the average vibrational energy
    // of either appearing or vanishing molecules based on
    // the local translational.

    if ( isnan(val) || isinf(val) ) {
	// Just return the energy from the first mode
	// assuming that it's dominant.
	return sems[0]->eval_energy_from_T(T);
    }

    return val;	// energy in J/kg
}

/****************** Knab vanishing component ********************/

Knab_vanishing_component::
Knab_vanishing_component(lua_State *L, Reaction *r, int idc )
: Coupling_component(L,r,"Knab_vanishing_component","vibration",idc)
{
    U0_ = get_number(L, -1, "U0");
    U1_ = get_number(L, -1, "U1");
    alpha_ = get_number(L, -1, "alpha");
    A_var_ = get_number(L, -1, "A");
}

Knab_vanishing_component::
Knab_vanishing_component(const Knab_vanishing_component &c)
    : Coupling_component(c), U0_(c.U0_), U1_(c.U1_), alpha_(c.alpha_), A_var_(c.A_var_) {}

Knab_vanishing_component::
~Knab_vanishing_component()
{}

Knab_vanishing_component*
Knab_vanishing_component::clone() const
{
    return new Knab_vanishing_component(*this);
}

double
Knab_vanishing_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[imode_];

    // 1. Calculate energy of vanishing molecules
    double e_va = calculate_Knab_energy( sems_, U0_, U1_, alpha_, A_var_, T, Tv );
    // Convert to J/particle
    e_va *= m_;

    // 2. Calculate delta_E
    double delta_N = nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_va - e_old_ ) * delta_N;

    // cout << "Knab_vanishing_component::compute_contribution()" << endl
    //      << "gamma = " << gamma << ", Tv = " << Tv << ", T = " << T << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Knab_vanishing_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: assuming mode '0' is translation, as it always should be
    double T = Q.T[0];
    double Tv = Q.T[imode_];

    // 1. Calculate energy of vanishing molecules
    double e_va = calculate_Knab_energy( sems_, U0_, U1_, alpha_, A_var_, T, Tv );
    // Convert to J/particle
    e_va *= m_;

    // 2. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = ( e_va - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Knab_vanishing_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}

/****************** Knab appearing component ********************/

Knab_appearing_component::
Knab_appearing_component(lua_State *L, Reaction *r, int idc )
  : Coupling_component(L,r,"Knab_appearing_component","vibration",idc)
{
    U0_ = get_number(L,-1,"U0");
    U1_ = get_number(L,-1,"U1");
    alpha_ = get_number(L,-1,"alpha");
    A_var_ = get_number(L,-1,"A");
}

Knab_appearing_component::
Knab_appearing_component( const Knab_appearing_component &c )
    : Coupling_component(c), U0_(c.U0_), U1_(c.U1_), alpha_(c.alpha_), A_var_(c.A_var_) {}

Knab_appearing_component::
~Knab_appearing_component()
{}

Knab_appearing_component*
Knab_appearing_component::clone() const
{
    return new Knab_appearing_component(*this);
}

double
Knab_appearing_component::
specific_compute_contribution( Gas_data &Q, vector<double> &delta_c )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: - assuming mode '0' is translation, as it always should be
    //          - assuming Tv=T for recombination processes
    double T = Q.T[0];
    double Tv = T;

    // 1. Calculate energy of appearing molecules
    double e_app = calculate_Knab_energy( sems_, U0_, U1_, alpha_, A_var_, T, Tv );
    // Convert to J/particle
    e_app *= m_;

    // 3. Calculate delta_E
    double delta_N = - nu_ * delta_c[idc_] * PC_Avogadro;
    double delta_E = ( e_app - e_old_ ) * delta_N;

    // cout << "Knab_appearing_component::compute_contribution()" << endl
    //      << "gamma = " << gamma << ", Tv = " << Q.T[imode_] << ", T = " << Q.T[0] << endl
    //      << "delta_c[idc_] = " << delta_c[idc_] << ", delta_N = " << delta_N << endl
    //      << "delta_E = " << delta_E << endl;

    return delta_E;
}

double
Knab_appearing_component::
specific_compute_source_term( Gas_data &Q, vector<double> &dcdt )
{
    if( Q.massf[isp_] < min_mass_frac ) return 0.0;

    // 0. Pull out translational and vibrational temperatures
    //    NOTE: - assuming mode '0' is translation, as it always should be
    //          - assuming Tv=T for recombination processes
    double T = Q.T[0];
    double Tv = T;

    // 1. Calculate energy of appearing molecules
    double e_app = calculate_Knab_energy( sems_, U0_, U1_, alpha_, A_var_, T, Tv );
    // Convert to J/particle
    e_app *= m_;

    // 3. Calculate dEdt
    // NOTE: minus e_old_ is required as average energy flux due to reactions is taken into account separately
    double dEdt = - ( e_app - e_old_ ) * nu_ * dcdt[idc_] * PC_Avogadro;

    // cout << "Knab_appearing_component::specific_compute_source_term()" << endl
    //      << "dEdt = " << dEdt << ", nu_ = " << nu_ << ", dcdt = " << dcdt[idc_] << endl;

    return dEdt;
}
