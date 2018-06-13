// Author: Daniel F Potter
// Date: 07-Dec-2009

#ifndef MACHERET_DISSOCIATION_HH
#define MACHERET_DISSOCIATION_HH

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "reaction-rate-coeff.hh"
#include "../models/gas_data.hh"
#include "generalised-Arrhenius.hh"

class Macheret_dissociation : public Generalised_Arrhenius {
public:
    Macheret_dissociation(lua_State *L, Gas_model &g, double T_upper, double T_lower);
    Macheret_dissociation(double A, double n, double E_a, double T_upper, double T_lower,
			  std::string v_name, std::string c_name, double khigh );
    ~Macheret_dissociation();

private:
    int s_eval(const Gas_data &Q);
    
private:
    double A_;
    double n_;
    double E_a_;
    int iTv_;
    int iT_;
    int ic_;
    bool monatomic_collider_;
    double alpha_;
    double theta_v_;
    double theta_d_;
    double theta_dstar_,delta_d_;
    double b_;
    double fac_;

    double khigh_;
};

Reaction_rate_coefficient* create_Macheret_dissociation_coefficient(lua_State *L, Gas_model &g, double T_upper, double T_lower);

#endif
