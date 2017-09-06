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
   rt={'Millikan-White:HTC', a=47.7, b=0.059, HTCS = { type = "Park", sigma_dash = 3.0e-17 } }
}

mechanism{
   'O2 ~~ ( O2 ) : V-T',
   rt={'Millikan-White:HTC', a=138, b=0.03, HTCS = { type = "Park", sigma_dash = 3.0e-17 } }
}
