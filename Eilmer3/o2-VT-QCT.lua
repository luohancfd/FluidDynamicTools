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
--   rt={'Millikan-White:HTC', a=47.7, b=0.059, HTCS = { type = "Park", sigma_dash = 3.0e-17 } }
   --   Marat O2+O
   --rt={'PolyFit', type=1, polybase=-10.0, coeff={1.92, -5.119, -14.19},T_norm=1000.0,T_pow=-1.0/3.0 }
   -- Ibraguimova
   rt={'PolyFit', type=0, a = 1.5e-12, b = 86.4, thetav = 2238.0} 
}

mechanism{
   'O2 ~~ ( O2 ) : V-T',
   -- Ibraguimova O2+O2
   rt={'PolyFit', type=2, a = 8.8e-14, b = 172.7, thetav = 2238.0, T_high=6000.0, polybase=10.0, coeff={-31307.9, 6600.9, -407.1, 1.11},T_norm=1.0,T_pow=-1.0/3.0 }
}
