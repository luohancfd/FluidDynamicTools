-- Baulch reaction rates

scheme_t = {
    update = "chemical kinetic ODE MC",
    temperature_limits = {
        lower = 20.0,
        upper = 100000.0
    },
    error_tolerance = 0.000001
}

-- Dissociation reactions
-- Note by Han: Ta = p_name(p_mode)^s_p * q_name(q_mode)^s_q
-- Equilibrium contant calculate by T[iT]
-- iT = 0 translational
-- iT = 1 vibrational
-- iT = 2 electronic 
-- This is for three temperature gas, check python file for detail
reaction{
   'O2 + O <=> O + O + O',
   fr={'Park', A=2.5e18, n=-0.565, T_a=60491.819, p_name='O2', p_mode='vibration', s_p=0.5, q_name='O', q_mode='translation'},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + O2 <=> O + O + O2',
   fr={'Park', A=1.63019353e21, n=-1.319, T_a=63050, p_name='O2', p_mode='vibration', s_p=0.5, q_name='O2', q_mode='translation'},
   ec={model='from CEA curves',iT=0}
}


