/** \file species_pieces.cxx
 *  \brief Definition of SpeciesPieces class.
 *
 *  \author Rowan J Gollan
 *  \date 20-May-2006
 *
 **/

#include <vector>
#include "../../nm/source/no_fuss_linear_algebra.hh"
#include "species-pieces.hh"

using namespace std;

Species_pieces::Species_pieces( int isp, vector<int> reac_index, vector<int> nu )
    : reac_index_( reac_index ), nu_( nu ), init_conc_( 0.0 ) {}

Species_pieces::Species_pieces( const Species_pieces &s ) 
    : reac_index_( s.reac_index_ ), nu_( s.nu_ ), init_conc_( s.init_conc_ ) {}

Species_pieces::~Species_pieces() {}

Species_pieces*
Species_pieces::clone()
{
    return new Species_pieces(*this);
}

double
Species_pieces::eval_conc( const vector<double> &w ) 
{
    double conc = init_conc_;
    int ir;
    for( size_t i = 0; i < reac_index_.size(); ++i ) {
	ir = reac_index_[i];
	conc += nu_[i] * w[ir];
    }
    
    return conc;

}

