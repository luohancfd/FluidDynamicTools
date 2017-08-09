scheme_t = {
    update = "energy exchange ODE",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

ode_t = {
    step_routine = 'rkf',
    max_step_attempts = 4,
    max_increase_factor = 1.15,
    max_decrease_factor = 0.01,
    decrease_factor = 0.333
}

-- all VT exchange mechanisms identified by Park (1993)
-- all ET exchange mechanisms from Gnoffo (1989)

mechanism{
   'O2 ~~ ( O ) : V-T',
   rt={'Poly-Fit', type=1, polybase=-10.0, coeffs={1.92, -5.119, -14.19},T_norm=1000.0,T_pow=-1.0/3.0 }
-- Ref: Borges Sebasti?o, Israel, Marat Kulakhmetov, and Alina Alexeenko. "DSMC study of oxygen shockwaves based on high-fidelity vibrational relaxation and dissociation models." Physics of Fluids 29.1 (2017): 017102.
}

mechanism{
   'O2 ~~ ( O2 ) : V-T',
   rt={'Poly-Fit', a=8.8e-14, b=172.7, thetav=2238.0, type=2, T_high=6000.0, polybase=10.0, coeffs={-31307.9, 6600.9, -407.1, 1.11},T_norm=1.0,T_pow=-1.0/3.0 }
-- Ref: Ibraguimova, L. B., et al. "Investigation of oxygen dissociation and vibrational relaxation at temperatures 4000¨C10 800 K." The Journal of chemical physics 139.3 (2013): 034317.
}
