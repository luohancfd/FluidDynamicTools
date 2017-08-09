// Author: Daniel F Potter
// Version: 21-Sep-2009
//          Initial coding.
//
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"

#include "chemical-species.hh"
#include "physical_constants.hh"
#include "gas-model.hh"
#include "CEA-Cp-functor.hh"
#include "CEA-h-functor.hh"
#include "CEA-s-functor.hh"

using namespace std;

/* ------- Chemical species ------- */

Chemical_species * new_chemical_species_from_file( string name, string inFile )
{
    lua_State *L = initialise_lua_State();
    
    if( do_gzfile(L, inFile) != 0 ) {
	ostringstream ost;
	ost << "new_chemical_species_from_file():\n";
	ost << "Error in input file: " << inFile << endl;
	input_error(ost);
    }
    
    lua_getglobal(L, name.c_str());
    string species_type = get_string( L, -1, "species_type" );
    cout << "- Creating " << name << " as a new " << species_type << " species" << endl;
    Chemical_species * X = 0;
    if ( species_type.find("monatomic")!=string::npos )
	X = new Atomic_species( string(name), species_type, 0, 0.0, L );
    else if ( species_type.find("fully coupled diatomic")!=string::npos )
        X = new Fully_coupled_diatomic_species( string(name), species_type, 0, 0.0, L );
    else if ( species_type.find("diatomic")!=string::npos )
	X = new Diatomic_species( string(name), species_type, 0, 0.0, L );
    else if ( species_type.find("fully coupled polyatomic")!=string::npos )
        X = new Fully_coupled_polyatomic_species( string(name), species_type, 0, 0.0, L );
    else if ( species_type.find("polyatomic")!=string::npos )
        X = new Polyatomic_species( string(name), species_type, 0, 0.0, L );
    else if ( species_type.find("free electron")!=string::npos )
	X = new Free_electron_species( string(name), species_type, 0, 0.0, L );
    else {
	ostringstream ost;
	ost << "new_chemical_species_from_file():\n";
	ost << "Could not decode type label for species " << name << ": " << species_type << endl;
	input_error(ost);
    }
    
    return X;
}

Chemical_species::Chemical_species( string name, string type, int isp, double min_massf, lua_State * L )
 : name_( name ), type_( type), isp_( isp ), min_massf_( min_massf )
{
    M_ = get_positive_value( L, -1, "M" );	// now kg/mol
    R_   = PC_R_u/M_;
    s_0_ = get_value( L, -1, "s_0" );
    h_f_ = get_value( L, -1, "h_f" );
    I_   = get_value( L, -1, "I" );
    Z_   = (int) get_value( L, -1, "Z" );
    // LJ parameters
    // NOTE: using eps0 to test presence of this data
    lua_getfield(L, -1, "eps0");
    if ( !lua_istable(L, -1) ) {
	cout << "Chemical_species::Chemical_species()\n";
	cout << "Could not find 'eps0' table" << endl;
	cout << "Initialising all LJ parameters to zero" << endl;
	eps0_ = 0.0; sigma_ = 0.0;
	lua_pop(L,1);
    }
    else {
    	lua_pop(L,1);
	eps0_ = get_positive_value( L, -1, "eps0" );
	sigma_ = get_positive_value( L, -1, "sigma" );
    }
    
    // Create the CEA h and s segmented functors
    vector<Univariate_functor*> Cp;
    vector<Univariate_functor*> h;
    vector<Univariate_functor*> s;
    vector<double> breaks;
    double T_low, T_high;

    lua_getfield(L, -1, "CEA_coeffs");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Chemical_species::Chemical_species():\n";
	ost << "Error locating 'CEA_coeffs' table for species: " << name_ << endl;
	input_error(ost);
    }
    for ( size_t i = 1; i <= lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i);
	T_low = get_positive_number(L, -1, "T_low");
	T_high = get_positive_number(L, -1, "T_high");
	Cp.push_back(new CEA_Cp_functor(L, R_));
	h.push_back(new CEA_h_functor(L, (*Cp.back())(T_low), (*Cp.back())(T_high), R_));
	s.push_back(new CEA_s_functor(L, (*Cp.back())(T_low), (*Cp.back())(T_high), R_));
	breaks.push_back(T_low);
	lua_pop(L, 1);
    }
    breaks.push_back(T_high);
    lua_pop(L, 1); // pop coeffs
    h_ = new Segmented_functor(h, breaks);
    s_ = new Segmented_functor(s, breaks);
    Cp_ = new Segmented_functor(Cp, breaks);
    
    for ( size_t i = 0; i < Cp.size(); ++i ) {
	delete Cp[i];
	delete h[i];
	delete s[i];
    } 
}

Chemical_species::~Chemical_species()
{
    for ( size_t iem=0; iem<modes_.size(); ++iem )
    	delete modes_[iem];
    
    delete h_;
    delete s_;
    delete Cp_;
}

Species_energy_mode*
Chemical_species::get_mode_pointer_from_type( string type )
{
    for ( size_t iem=0; iem<modes_.size(); ++iem ) {
    	if ( modes_[iem]->get_type()==type ) return modes_[iem];
    }
    
    cout << "Chemical_species::get_mode_pointer_from_type()" << endl
    	 << "mode type: " << type << " not found." << endl
    	 << "Bailing out!" << endl;
    exit( BAD_INPUT_ERROR );
}

Electronic*
Chemical_species::get_electronic_mode_pointer()
{
    // First test that there is mode in the 1 position
    if ( modes_.size()<2 ) {
        cout << "Chemical_species::get_electronic_mode_pointer()" << endl
             << "No energy mode found in the expected position!" << endl;
        exit( FAILURE );
    }

    // Now ensure that this is the electronic mode
    // FIXME: somehow need to ensure this is a Multi_level_electronic instance
    if ( modes_[1]->get_type()!="electronic") {
        cout << "Chemical_species::get_multi_level_electronic_mode_pointer()" << endl
             << "Electronic mode not found in the expected position!" << endl;
        exit( FAILURE );
    }

    // Now can safely do the dynamic cast and return the pointer
    return dynamic_cast<Electronic*>(modes_[1]);
}

Multi_level_electronic*
Chemical_species::get_multi_level_electronic_mode_pointer()
{
    // First test that there is mode in the 1 position
    if ( modes_.size()<2 ) {
        cout << "Chemical_species::get_multi_level_electronic_mode_pointer()" << endl
             << "No energy mode found in the expected position!" << endl;
        exit( FAILURE );
    }

    // Now ensure that this is the electronic mode
    // FIXME: somehow need to ensure this is a Multi_level_electronic instance
    if ( modes_[1]->get_type()!="electronic") {
        cout << "Chemical_species::get_multi_level_electronic_mode_pointer()" << endl
             << "Electronic mode not found in the expected position!" << endl;
        exit( FAILURE );
    }

    // Now can safely do the dynamic cast and return the pointer
    return dynamic_cast<Multi_level_electronic*>(modes_[1]);
}

int
Chemical_species::
get_element_count( string X )
{
    // 0. Create neutral name for base-element searching
    string neutral_name = name_;
    if      ( Z_== 1 ) neutral_name.assign( name_, 0, name_.size()-5 );
    else if ( Z_==-1 ) neutral_name.assign( name_, 0, name_.size()-6 );
    
    // 1. Search for the requested base element
    int beta = 0;
    string::size_type loc = neutral_name.find( X );
    if ( loc!=string::npos ) {
    	// it is present, but how many?
    	if ( loc!=(neutral_name.size()-X.size()) ) {
    	    if ( neutral_name.compare(loc+X.size(),1,"2")==0 )
    	    	beta = 2;
    	    else if ( neutral_name.compare(loc+X.size(),1,"3")==0 )
    	    	beta = 3;
    	    else if ( neutral_name.compare(loc+X.size(),1,"4")==0 )
    	    	beta = 4;
    	    // else -> beta already set to zero
    	}
    	else {
    	    beta = 1;
    	}
    }
    // else -> beta already set to zero
    
    return beta;
}

const string base_elements[] = { "Ar", "He", "C", "H", "N", "O" };
const int nbe = 6;

void
Chemical_species::
partial_equilibrium_participants( vector<int> &betas,
    				  vector<string> &participants )
{
    // If this is an electron we don't need to do anything
    if ( name_=="e_minus" ) return;
    
    // 0. Create neutral name for base-element searching
    string neutral_name = name_;
    if      ( Z_== 1 ) neutral_name.assign( name_, 0, name_.size()-5 );
    else if ( Z_==-1 ) neutral_name.assign( name_, 0, name_.size()-6 );
    
    // 1. Search for base-elements and create as reactants
    while ( neutral_name.size()>0 ) {
        for ( int ib=0; ib<nbe; ++ib ) {
            string::size_type loc = neutral_name.find( base_elements[ib] );
            if ( loc!=string::npos ) {
                participants.push_back( base_elements[ib] );
                betas.push_back( -1 );
                if ( loc!=(neutral_name.size()-base_elements[ib].size()) ) {
                    if ( neutral_name.compare(loc+base_elements[ib].size(),1,"2")==0 ) {
                        betas.back() = -2;
                        neutral_name.erase( loc+base_elements[ib].size(), 1 );
                    }
                    else if ( neutral_name.compare(loc+base_elements[ib].size(),1,"3")==0 ) {
                        betas.back() = -3;
                        neutral_name.erase( loc+base_elements[ib].size(), 1 );
                    }
                    else if ( neutral_name.compare(loc+base_elements[ib].size(),1,"4")==0 ) {
                        betas.back() = -4;
                        neutral_name.erase( loc+base_elements[ib].size(), 1 );
                    }
                    // else -> beta already set to -1
                }
                neutral_name.erase( loc, base_elements[ib].size() );
            }
        }
    }
    
    // 2. Create products from this species and electrons if not neutral
    betas.push_back(1);
    participants.push_back( name_ );
    if ( Z_!=0 ) {
    	betas.push_back( Z_ );
    	participants.push_back( "e_minus" );
    }
    
#   if 0
    // FIXME: Remove these print statements after testing
    cout << "species: " << name_ << " has betas: ";
    for ( size_t ip = 0; ip < betas.size(); ++ip ) 
	cout << betas[ip] << ", ";
    cout << "\n";
    cout << "species: " << name_ << " has participants: ";
    for ( size_t ip = 0; ip < participants.size(); ++ip ) 
	cout << participants[ip] << ", ";
    cout << "\n";
#   endif

    return;
}

double Chemical_species::s_eval_energy(const Gas_data &Q)
{
    // 1. Heat of formation
    double e = h_f_;
    
    // 2. Thermal energy modes
    for ( size_t iem=0; iem<modes_.size(); ++iem )
    	e += modes_[iem]->eval_energy(Q);

    return e;
}

double Chemical_species::s_eval_modal_energy(int itm, const Gas_data &Q)
{
    // 1. Initialise h
    //    CHECK-ME: Omitt Heat of formation?
    double h = 0.0;
    
    // 2. Look for itm in thermal energy modes and add contribution
    //    FIXME: need a more efficient method eventually, but only used for
    //           diffusion at present
    for ( size_t iem=0; iem<modes_.size(); ++iem ) {
    	if ( modes_[iem]->get_iT()==itm )
    	    h += modes_[iem]->eval_energy(Q);
    }
    
    return h;
}

double Chemical_species::s_eval_enthalpy(const Gas_data &Q)
{
    // 1. Heat of formation
    double h = h_f_;
    
    // 2. Thermal energy modes
    for ( size_t iem=0; iem<modes_.size(); ++iem )
    	h += modes_[iem]->eval_enthalpy(Q);
    
    return h;
}

double
Chemical_species::
s_eval_CEA_enthalpy(const Gas_data &Q)
{
    // NOTE: returning enthalpy in units of J/kg-of-this-species
    double h = (*h_)(Q.T[0]);	// eval the CEA h functor
    
    return h;
}

double Chemical_species::s_eval_modal_enthalpy(int itm, const Gas_data &Q)
{
    // 1. Initialise h
    //    CHECK-ME: Omitt Heat of formation?
    double h = 0.0;
    
    // 2. Look for itm in thermal energy modes and add contribution
    //    FIXME: need a more efficient method eventually, but only used for
    //           diffusion at present
    for ( size_t iem=0; iem<modes_.size(); ++iem ) {
    	if ( modes_[iem]->get_iT()==itm )
    	    h += modes_[iem]->eval_enthalpy(Q);
    }
    
    return h;
}

double Chemical_species::s_eval_entropy(const Gas_data &Q)
{
    // NOTE: Standard state entropy (s_0_) is included implicitly in calculation
    double s = 0.0;
    
    for ( size_t iem=0; iem<modes_.size(); ++iem ) {
    	s += modes_[iem]->eval_entropy(Q);
    }
    
    return s;
}

double Chemical_species::s_eval_Cv(const Gas_data &Q)
{
    double Cv = 0.0;
    for ( size_t iem=0; iem<modes_.size(); ++iem )
    	Cv += modes_[iem]->eval_Cv(Q);
    
    return Cv;
}

double Chemical_species::s_eval_Cp(const Gas_data &Q)
{
    double Cp = 0.0;
    for ( size_t iem=0; iem<modes_.size(); ++iem )
    	Cp += modes_[iem]->eval_Cp(Q);
    
    return Cp;
}

double
Chemical_species::
s_eval_gibbs_free_energy( double T )
{
    // NOTE: returning Gibbs free energy in units of J/kg-of-this-species
    double g = h_f_;
    for ( size_t iem=0; iem<modes_.size(); ++iem )
    	g += modes_[iem]->eval_enthalpy_from_T(T) - T*modes_[iem]->eval_entropy_from_T(T);
    
    return g;
}

double
Chemical_species::
s_eval_CEA_Gibbs_free_energy( double T )
{
    // NOTE: returning Gibbs free energy in units of J/kg-of-this-species
    double h = (*h_)(T);	// eval the CEA h functor
    double s = (*s_)(T);	// eval the CEA s functor
    double g = h - T*s;
    
    return g;
}

double
Chemical_species::
s_eval_partition_function( double T )
{
    double Q = exp(-h_f_/R_/T);
    for ( size_t iem=0; iem<modes_.size(); ++iem )
        Q *= modes_[iem]->eval_Q_from_T(T);

    return Q;
}


/* ------- Atomic species ------- */

Atomic_species::Atomic_species( string name, string type, int isp, double min_massf, lua_State * L )
 : Chemical_species( name, type, isp, min_massf, L )
{
    lua_getfield(L, -1, "electronic_levels");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Atomic_species::Atomic_species()\n";
	ost << "Error locating 'level_data' table" << endl;
	input_error(ost);
    }
    
    int nlevs = get_int(L, -1, "n_levels");
    if ( nlevs < 1 ) {
    	ostringstream ost;
    	ost << "Atomic_species::Atomic_species()\n";
    	ost << "Require at least 1 electronic level - only " << nlevs << " present.\n";
    	input_error(ost);
    }
    
    // Get electronic level data
    vector<double> lev_data, theta_vec;
    vector<int> g_vec;
    
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	ostringstream lev_oss;
    	lev_oss << "ilev_" << ilev;
    	lua_getfield(L, -1, lev_oss.str().c_str());
    	lev_data.clear();
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    lev_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ilev
	g_vec.push_back( int(lev_data[2]) );
	theta_vec.push_back( lev_data[1]*PC_c*PC_h_SI/PC_k_SI );	// cm-1 -> K
    }
    
    lua_pop(L,1);	// pop 'electronic_levels'
    
    // Create energy modes
    // [0] Translation
    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
    // [1] Electronic
    if ( nlevs==1 ) 
    	modes_.push_back( new One_level_electronic( isp_, R_, min_massf_, g_vec[0], theta_vec[0] ) );
    else if ( nlevs==2 )
    	modes_.push_back( new Two_level_electronic( isp_, R_, min_massf_, g_vec[0], theta_vec[0], g_vec[1], theta_vec[1] ) );
    else 
    	modes_.push_back( new Multi_level_electronic( isp_, R_, min_massf_, g_vec, theta_vec ) );
}

Atomic_species::~Atomic_species()
{
    // NOTE: energy modes are cleared by ~Chemical_species()
}

/* ------- Diatomic species ------- */

Diatomic_species::Diatomic_species( string name, string type, int isp, double min_massf, lua_State * L )
 : Chemical_species( name, type, isp, min_massf, L )
{
    // polarity flag
    if ( type.find("nonpolar")!=string::npos )
    	polar_flag_ = false;
    else
    	polar_flag_ = true;
    
    // SSH theory parameters
    // NOTE: using r0 to test presence of this data
    lua_getfield(L, -1, "r0");
    if ( !lua_istable(L, -1) ) {
	cout << "Diatomic_species::Diatomic_species()\n";
	cout << "Could not find 'r0' table for species: " << name << endl;
	cout << "Initialising all SSH related parameters to zero" << endl;
	r0_ = 0.0; r_eq_ = 0.0; f_m_ = 0.0;
	mu_ = 0.0; alpha_ = 0.0; mu_B_ = 0.0;
	lua_pop(L,1);
    }
    else {
    	lua_pop(L,1);
	r0_ = get_positive_value( L, -1, "r0" );
	r_eq_ = get_positive_value( L, -1, "r_eq" );
	f_m_ = get_positive_value( L, -1, "f_m" );
	mu_ = get_positive_value( L, -1, "mu" );
	if ( polar_flag_ ) {
	    // We might need the dipole moment
	    mu_B_ = get_positive_value( L, -1, "mu_B" );
	}
	else {
	    // We might need the electric polarizability
	    alpha_ = get_positive_value( L, -1, "alpha" );
	}
    }
    
    

    lua_getfield(L, -1, "electronic_levels");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Diatomic_species::Diatomic_species()\n";
	ost << "Error locating 'level_data' table" << endl;
	input_error(ost);
    }
    
    int nlevs = get_int(L, -1, "n_levels");
    if ( nlevs < 1 ) {
    	ostringstream ost;
    	ost << "Diatomic_species::Diatomic_species()\n";
    	ost << "Require at least 1 electronic level - only " << nlevs << " present.\n";
    	input_error(ost);
    }
    
    // Get electronic level data
    vector<double> lev_data, theta_el_vec;
    vector<int> g_el_vec;
    
    // Ground state - outside of loop to get vibration and rotation constants
    string lev_str = "ilev_0";
    lua_getfield(L, -1, lev_str.c_str());
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	lev_data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L,1);	// pop ilev0
    g_el_vec.push_back( int(lev_data[2]) );
    theta_el_vec.push_back( lev_data[0]*PC_c*PC_h_SI/PC_k_SI );	// cm-1 -> K
    double theta_v = lev_data[4]*PC_c*PC_h_SI/PC_k_SI;
    double theta_r = lev_data[8]*PC_c*PC_h_SI/PC_k_SI;
    double theta_D = lev_data[3]*PC_c*PC_h_SI/PC_k_SI;
    // Check that a usable dissociation energy is present
    if ( theta_D < 0.0 && lev_data[5] > 0.0 ) {
    	// use D ~ omega_e**2 / ( 4*xomega_e )
    	theta_D = ( lev_data[4]*lev_data[4] / ( 4*lev_data[5] ) ) * PC_c*PC_h_SI/PC_k_SI;
    }
    else if ( theta_D < 0.0 ) {
        // set to 164,000K (just above H2 dissociation limit)
        theta_D = 1.64e5;
    }
    // Set the characteristic vibrational temperature in kelvin
    theta_v_ = theta_v;
    
    // Excited states - electron degeneracy and energy only
    for ( int ilev=1; ilev<nlevs; ++ilev ) {
    	ostringstream lev_oss;
    	lev_oss << "ilev_" << ilev;
    	lua_getfield(L, -1, lev_oss.str().c_str());
    	lev_data.clear();
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    lev_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ilev
	g_el_vec.push_back( int(lev_data[2]) );
	theta_el_vec.push_back( lev_data[0]*PC_c*PC_h_SI/PC_k_SI );	// cm-1 -> K
    }
    
    // Rotational homonuclear factor
    int sigma = 1;
    // sigma is 2 for homonuclear diatoms
    if ( name.find("2")!=string::npos ) sigma = 2;
    
    lua_pop(L,1);	// pop 'electronic_levels'
    
    // Vibrational oscillator type
    oscillator_type_ = get_string( L, -1, "oscillator_type" );

    // We already set characteristic vibrational temperature based 
    // on the energy levels, BUT we'll pick up an explicit value for
    // theta_v if given.

    lua_getfield(L, -1, "theta_v");
    if ( !lua_istable(L, -1) ) {
	lua_pop(L, 1);
    }
    else {
	lua_pop(L, 1);
	theta_v = get_positive_value(L, -1, "theta_v");
	theta_v_ = theta_v;
    }
    // Create energy modes
    // [0] Translation
    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
    // [1] Electronic
    if ( nlevs==1 ) 
    	modes_.push_back( new One_level_electronic( isp_, R_, min_massf_, g_el_vec[0], theta_el_vec[0] ) );
    else if ( nlevs==2 )
    	modes_.push_back( new Two_level_electronic( isp_, R_, min_massf_, g_el_vec[0], theta_el_vec[0], g_el_vec[1], theta_el_vec[1] ) );
    else 
    	modes_.push_back( new Multi_level_electronic( isp_, R_, min_massf_, g_el_vec, theta_el_vec ) );
    // [2] Rotation
    modes_.push_back( new Fully_excited_rotation( isp_, R_, min_massf_, theta_r, sigma ) );
    // [3] Vibration
    if ( oscillator_type_=="harmonic" ) 
    	modes_.push_back( new Harmonic_vibration( isp_, R_, min_massf_, theta_v ) );
    else if ( oscillator_type_=="truncated harmonic" ) 
    	modes_.push_back( new Truncated_harmonic_vibration( isp_, R_, min_massf_, theta_v, theta_D ) );
    else if ( oscillator_type_ == "anharmonic")
        modes_.push_back( new Anharmonic_vibration(isp_,R_,min_massf_,theta_v,theta_r,sigma,h_,s_,Cp_));
    else {
    	ostringstream oss;
    	oss << "Diatomic_species::Diatomic_species()" << endl
    	    << name_ << " oscillator_type: " << oscillator_type_ << " not available" << endl;
    	input_error( oss );
    }
}

Diatomic_species::~Diatomic_species()
{
    // NOTE: energy modes are cleared by ~Chemical_species()
}

/* ------- Diatomic species with fully coupled internal modes ------- */

Fully_coupled_diatomic_species::
Fully_coupled_diatomic_species( string name, string type, int isp, double min_massf, lua_State * L )
 : Chemical_species( name, type, isp, min_massf, L )
{
    // polarity flag
    if ( type.find("nonpolar")!=string::npos )
    	polar_flag_ = false;
    else
    	polar_flag_ = true;
    
    // SSH theory parameters
    // NOTE: using r0 to test presence of this data
    lua_getfield(L, -1, "r0");
    if ( !lua_istable(L, -1) ) {
	cout << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()\n";
	cout << "Could not find 'r0' table" << endl;
	cout << "Initialising all SSH related parameters to zero" << endl;
	r0_ = 0.0; r_eq_ = 0.0; f_m_ = 0.0;
	mu_ = 0.0; alpha_ = 0.0; mu_B_ = 0.0;
	lua_pop(L,1);
    }
    else {
    	lua_pop(L,1);
	r0_ = get_positive_value( L, -1, "r0" );
	r_eq_ = get_positive_value( L, -1, "r_eq" );
	f_m_ = get_positive_value( L, -1, "f_m" );
	mu_ = get_positive_value( L, -1, "mu" );
	if ( polar_flag_ ) {
	    // We might need the dipole moment
	    mu_B_ = get_positive_value( L, -1, "mu_B" );
	}
	else {
	    // We might need the electric polarizability
	    alpha_ = get_positive_value( L, -1, "alpha" );
	}
    }
    
#   if TABULATED_COUPLED_DIATOMIC_MODES==0
    // Set the temperature perturbation factor (used for finding derivatives numerically)
    double fT = 1.0e-6;
    
    lua_getfield(L, -1, "electronic_levels");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()\n";
	ost << "Error locating 'level_data' table" << endl;
	input_error(ost);
    }
    
    int nlevs = get_int(L, -1, "n_levels");
    if ( nlevs < 1 ) {
    	ostringstream ost;
    	ost << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()\n";
    	ost << "Require at least 1 electronic level - only " << nlevs << " present.\n";
    	input_error(ost);
    }
    
    // Create the electronic levels
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
    	vector<double> elev_data;
    	ostringstream lev_oss;
    	lev_oss << "ilev_" << ilev;
    	lua_getfield(L, -1, lev_oss.str().c_str());
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    elev_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
    	cout << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()" << endl
    	<< "Attempting to create ilev = " << ilev << " for species: " << name << endl;
    	elevs_.push_back( new Diatom_electronic_level( elev_data ) );
	lua_pop(L,1);	// pop ilev
    }
    
    // Rotational modes, temperature and sigma
    // int r_modes = 2;		// always 2 rot modes for diatoms
    int sigma_r = 1;
    // sigma is 2 for homonuclear diatoms
    if ( name.find("2")!=string::npos ) sigma_r = 2;
    
    lua_pop(L,1);	// pop 'electronic_levels'
    
    // Check the vibrational oscillator type
    oscillator_type_ = get_string( L, -1, "oscillator_type" );
    if ( oscillator_type_ != "truncated anharmonic" ) {
    	ostringstream ost;
    	ost << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()\n";
    	ost << "The oscillator type is required to be 'truncated anharmonic' for this species.\n";
    	input_error(ost);
    }
    
    // Create energy modes
    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
    modes_.push_back( new Coupled_diatomic_electronic( isp_, R_, min_massf_, sigma_r, fT, elevs_ ) );
    modes_.push_back( new Coupled_diatomic_rotation( isp_, R_, min_massf_, sigma_r, fT, elevs_ ) );
    modes_.push_back( new Coupled_diatomic_vibration( isp_, R_, min_massf_, sigma_r, fT, elevs_ ) );
    
    // Create a single temperature fully coupled internal mode
    fcd_int_ = new Fully_coupled_diatom_internal( isp_, R_, min_massf_, sigma_r, fT, elevs_ );
#   else
    lua_getfield(L, -1, "electronic_levels");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()\n";
	ost << "Error locating 'level_data' table" << endl;
	input_error(ost);
    }
    
    int nlevs = get_int(L, -1, "n_levels");
    if ( nlevs < 1 ) {
    	ostringstream ost;
    	ost << "Fully_coupled_diatomic_species::Fully_coupled_diatomic_species()\n";
    	ost << "Require at least 1 electronic level - only " << nlevs << " present.\n";
    	input_error(ost);
    }
    
    // Get the ground state spectroscopic parameters
    vector<double> lev_data;
    
    // Ground state - outside of loop to get vibration and rotation constants
    string lev_str = "ilev_0";
    lua_getfield(L, -1, lev_str.c_str());
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	lev_data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L,1);	// pop ilev0
    double theta_v = lev_data[4]*PC_c*PC_h_SI/PC_k_SI;
    double theta_r = lev_data[8]*PC_c*PC_h_SI/PC_k_SI;
    
    // Get the first excited state electronic energy
    lev_str = "ilev_1";
    lua_getfield(L, -1, lev_str.c_str());
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	lev_data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L,1);	// pop ilev1
    double theta_e = lev_data[0]*PC_c*PC_h_SI/PC_k_SI;
    
    lua_pop(L,1);   // pop 'electronic_levels'

    /* FIXME: assuming data is located in the current working directory */
    ostringstream LUT_fname;
    LUT_fname << "FCD_tables/" << name_ << ".data";

    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
    modes_.push_back( new Coupled_diatomic_electronic( isp_, R_, min_massf_, theta_e, LUT_fname.str()  ) );
    modes_.push_back( new Coupled_diatomic_rotation( isp_, R_, min_massf_, theta_r, LUT_fname.str() ) );
    modes_.push_back( new Coupled_diatomic_vibration( isp_, R_, min_massf_, theta_v, LUT_fname.str() ) );
    
    fcd_int_ = new Fully_coupled_diatom_internal( isp_, R_, min_massf_, LUT_fname.str() );
#   endif

}

Fully_coupled_diatomic_species::~Fully_coupled_diatomic_species()
{
#   if TABULATED_COUPLED_DIATOMIC_MODES==0
    // Clear the electronic levels
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev )
    	delete elevs_[ilev];
#   endif
    
    // Delete the fully coupled diatomic mode
    delete fcd_int_;
}

void Fully_coupled_diatomic_species::set_modal_temperature_indices()
{
    int iTt = modes_[0]->get_iT();
    int iTe = modes_[1]->get_iT();
    int iTr = modes_[2]->get_iT();
    int iTv = modes_[3]->get_iT();
    
    // cout << "Fully_coupled_diatomic_species::set_modal_temperature_indices()" << endl
    //      << "iTe = " << iTe << ", iTr = " << iTr << ", iTv = " << iTv << endl;
    
    dynamic_cast<Coupled_diatomic_electronic*>(modes_[1])->set_iTs(iTe,iTv,iTr);
    dynamic_cast<Coupled_diatomic_rotation*>(modes_[2])->set_iTs(iTe,iTv,iTr);
    dynamic_cast<Coupled_diatomic_vibration*>(modes_[3])->set_iTs(iTe,iTv,iTr);
    fcd_int_->set_iT(iTt);
}

double Fully_coupled_diatomic_species::s_eval_entropy(const Gas_data &Q)
{
    // NOTE: This species needs its own expression for entropy as the coupled 
    //       modes cannot presently calculate entropy.
    //       The calculation here assumes thermal equilibrium at T_trans.
    double s_trans = modes_[0]->eval_entropy(Q);
    double s_int = fcd_int_->eval_entropy(Q);
    
    return s_trans + s_int;
}

double Fully_coupled_diatomic_species::s_eval_partition_function( double T )
{
    // formation energy contribution
    double Q = exp(-h_f_/R_/T);
    // translational contribution
    Q *= modes_[0]->eval_Q_from_T(T);
    // internal contribution
    Q *= fcd_int_->eval_Q_from_T(T);

    return Q;
}

/* ------- Polyatomic species ------- */

Polyatomic_species::Polyatomic_species( string name, string type, int isp, double min_massf, lua_State * L )
 : Chemical_species( name, type, isp, min_massf, L )
{
    // Now explicitly setting the characteristic vibrational temperature
    theta_v_ = get_positive_value( L, -1, "theta_v" );

    lua_getfield(L, -1, "electronic_levels");
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Polyatomic_species::Polyatomic_species()\n";
	ost << "Error locating 'level_data' table" << endl;
	input_error(ost);
    }
    
    int nlevs = get_int(L, -1, "n_levels");
    if ( nlevs < 1 ) {
    	ostringstream ost;
    	ost << "Polyatomic_species::Polyatomic_species()\n";
    	ost << "Require at least 1 electronic level - only " << nlevs << " present.\n";
    	input_error(ost);
    }
    
    // Get electronic level data
    vector<double> lev_data, theta_el_vec, theta_v_vec;
    vector<int> g_el_vec;
    
    // Ground state - outside of loop to get vibration and rotation constants
    string lev_str = "ilev_0";
    lua_getfield(L, -1, lev_str.c_str());
    for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	lua_rawgeti(L, -1, i+1);
	lev_data.push_back( luaL_checknumber(L, -1) );
	lua_pop(L, 1 );
    }
    lua_pop(L,1);	// pop ilev0
    theta_el_vec.push_back( lev_data[0]*PC_c*PC_h_SI/PC_k_SI );	// cm-1 -> K
    g_el_vec.push_back( int(lev_data[2]) );
    double theta_D = lev_data[3]*PC_c*PC_h_SI/PC_k_SI;
    double theta_A0 = lev_data[4]*PC_c*PC_h_SI/PC_k_SI;
    double theta_B0 = lev_data[5]*PC_c*PC_h_SI/PC_k_SI;
    double theta_C0 = lev_data[6]*PC_c*PC_h_SI/PC_k_SI;
    int sigma = int(lev_data[7]);
    int sigma_r = int(lev_data[8]);
    for ( size_t i=9; i<lev_data.size(); ++i )
    	theta_v_vec.push_back( lev_data[i]*PC_c*PC_h_SI/PC_k_SI );
    // Currently sigma_r is not being used
    UNUSED_VARIABLE(sigma_r);
    
    // Excited states - electron degeneracy and energy only
    for ( int ilev=1; ilev<nlevs; ++ilev ) {
    	ostringstream lev_oss;
    	lev_oss << "ilev_" << ilev;
    	lua_getfield(L, -1, lev_oss.str().c_str());
    	lev_data.clear();
	for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
	    lua_rawgeti(L, -1, i+1);
	    lev_data.push_back( luaL_checknumber(L, -1) );
	    lua_pop(L, 1 );
	}
	lua_pop(L,1);	// pop ilev
	theta_el_vec.push_back( lev_data[0]*PC_c*PC_h_SI/PC_k_SI );	// cm-1 -> K
	g_el_vec.push_back( int(lev_data[2]) );
    }
    
    lua_pop(L,1);	// pop 'electronic_levels'
    
    // Create energy modes
    // [0] Translation
    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
    // [1] Electronic
    if ( nlevs==1 ) 
    	modes_.push_back( new One_level_electronic( isp_, R_, min_massf_, g_el_vec[0], theta_el_vec[0] ) );
    else if ( nlevs==2 )
    	modes_.push_back( new Two_level_electronic( isp_, R_, min_massf_, g_el_vec[0], theta_el_vec[0], g_el_vec[1], theta_el_vec[1] ) );
    else 
    	modes_.push_back( new Multi_level_electronic( isp_, R_, min_massf_, g_el_vec, theta_el_vec ) );
    // [2] Rotation
    if ( type.find("nonlinear")!=string::npos ) {
        linear_flag_ = false;
    	// Different s, e and c_v expressions for non-linear molecules
    	modes_.push_back( new Fully_excited_nonlinear_rotation( isp_, R_, min_massf_, theta_A0, theta_B0, theta_C0, sigma ) );
    }
    else {
        linear_flag_ = true;
    	modes_.push_back( new Fully_excited_rotation( isp_, R_, min_massf_, theta_B0, sigma ) );
    }
    // [3:] Vibration
    // Vibrational oscillator type
    oscillator_type_ = get_string( L, -1, "oscillator_type" );
    if ( oscillator_type_=="harmonic" ) {
    	// create multiple vibrational modes
    	for ( size_t iv=0; iv<theta_v_vec.size(); ++iv )
    	    modes_.push_back( new Harmonic_vibration( isp_, R_, min_massf_, theta_v_vec[iv] ) );
    }
    else if ( oscillator_type_=="truncated harmonic" ) {
    	// create multiple vibrational modes
    	// NOTE: setting equal parts of the dissociation energy to each vibrational mode
    	double T_d = theta_D;
    	if ( theta_v_vec.size() > 0 )
    	    T_d = theta_D / double(theta_v_vec.size());
    	for ( size_t iv=0; iv<theta_v_vec.size(); ++iv )
    	    modes_.push_back( new Truncated_harmonic_vibration( isp_, R_, min_massf_, theta_v_vec[iv], T_d) );
    }
    else {
    	ostringstream oss;
    	oss << "Polyatomic_species::Polyatomic_species()" << endl
    	    << name_ << " oscillator_type: " << oscillator_type_ << " not available" << endl;
    	input_error( oss );
    }

}

Polyatomic_species::~Polyatomic_species()
{
    // NOTE: energy modes are cleared by ~Chemical_species()
}

double
Polyatomic_species::
s_eval_Cv_vib( const Gas_data &Q )
{
    // Sum contributions from all vibrational modes
    double Cv_vib = 0.0;
    for ( size_t i=3; i<modes_.size(); ++i )
    	Cv_vib += modes_[i]->eval_Cv(Q);
    
    return Cv_vib;
}

/* ------- Fully coupled polyatomic species -------- */

Fully_coupled_polyatomic_species::
Fully_coupled_polyatomic_species( string name, string type, int isp, double min_massf, lua_State * L )
 : Chemical_species( name, type, isp, min_massf, L )
{
    // polarity flag
    if ( type.find("nonpolar")!=string::npos )
        polar_flag_ = false;
    else
        polar_flag_ = true;

    // linearity flag
    if ( type.find("nonlinear")!=string::npos )
        linear_flag_ = false;
    else
        linear_flag_ = true;

    // Set the temperature perturbation factor (used for finding derivatives numerically)
    double fT = 1.0e-6;

    lua_getfield(L, -1, "electronic_levels");
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "Fully_coupled_polyatomic_species::Fully_coupled_polyatomic_species()\n";
        ost << "Error locating 'level_data' table" << endl;
        input_error(ost);
    }

    int nlevs = get_int(L, -1, "n_levels");
    if ( nlevs < 1 ) {
        ostringstream ost;
        ost << "Fully_coupled_polyatomic_species::Fully_coupled_polyatomic_species()\n";
        ost << "Require at least 1 electronic level - only " << nlevs << " present.\n";
        input_error(ost);
    }

    // Create the electronic levels
    for ( int ilev=0; ilev<nlevs; ++ilev ) {
        vector<double> elev_data;
        ostringstream lev_oss;
        lev_oss << "ilev_" << ilev;
        lua_getfield(L, -1, lev_oss.str().c_str());
        for ( size_t i=0; i<lua_objlen(L, -1); ++i ) {
            lua_rawgeti(L, -1, i+1);
            elev_data.push_back( luaL_checknumber(L, -1) );
            lua_pop(L, 1 );
        }
        cout << "Fully_coupled_polyatomic_species::Fully_coupled_polyatomic_species()" << endl
             << "Attempting to create ilev = " << ilev << " for species: " << name << endl;
        int sigma_rot = elev_data[8];
        if ( sigma_rot==0 )
            elevs_.push_back( new Spherical_top_polyatom_electronic_level( elev_data ) );
        else {
            // Need to decide between symmetrical and asymmetric top...
            // FIXME: assuming asymmetrical top as most non-linear levels have a non-zero C0
            elevs_.push_back( new Asymmetric_top_polyatom_electronic_level( elev_data ) );
        }
        lua_pop(L,1);   // pop ilev
    }

    lua_pop(L,1);       // pop 'electronic_levels'

    // Create a single temperature fully coupled internal mode
    fcp_int_ = new Fully_coupled_polyatom_internal( isp_, R_, min_massf_, fT, elevs_ );

    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
}

Fully_coupled_polyatomic_species::~Fully_coupled_polyatomic_species()
{
    // Clear the electronic levels
    for ( size_t ilev=0; ilev<elevs_.size(); ++ilev )
        delete elevs_[ilev];

    // Delete the fully coupled polyatomic mode
    delete fcp_int_;
}

void Fully_coupled_polyatomic_species::set_modal_temperature_indices()
{
    return;
}

/* ------- Free electron species ------- */

Free_electron_species::Free_electron_species( string name, string type, int isp, double min_massf, lua_State * L )
 : Chemical_species( name, type, isp, min_massf, L )
{
    // Create translational energy mode
    modes_.push_back( new Fully_excited_translation( isp_, R_, min_massf_ ) );
    // Create electronic energy mode (necessary for correct entropy calculation)
    modes_.push_back( new One_level_electronic( isp_, R_, min_massf_, 2, 0.0 ) );
}

Free_electron_species::~Free_electron_species()
{
    // NOTE: energy modes are cleared by ~Chemical_species()
}

