/** \file species_pieces.hh
 *  \brief Class to organise the data for the concentration of each species.
 *
 *  \author Rowan J Gollan
 *  \date 30-May-2006
 *
 **/

#ifndef S_PIECES_HH
#define S_PIECES_HH

#include <vector>
#include "../../nm/source/no_fuss_linear_algebra.hh"

using namespace std;

class Species_pieces {
public: 
    Species_pieces( int isp, std::vector<int> reac_index, std::vector<int> nu );
    Species_pieces( const Species_pieces &s );
    virtual ~Species_pieces();
    Species_pieces* clone();

    double eval_conc( const std::vector<double> &w );
    void set_init_conc(double init_conc) { init_conc_ = init_conc; }
    double get_init_conc() { return init_conc_; }
    std::vector<int> reac_index_;
    std::vector<int> nu_;

private:
    double init_conc_;
};


#endif
