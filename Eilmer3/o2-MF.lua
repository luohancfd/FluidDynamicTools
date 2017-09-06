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
   -- Marat QCT
   fr={'Macheret', A=2.5e18, n=-0.565, T_a=60491.819, v_name='O2', c_name='O'},
   chemistry_energy_coupling={{species='O2', mode='vibration', model='Macheret',A=2.5e18,n=-0.565,T_d=60491.819,c_name='O'}},
   ec={model='from CEA curves',iT=0}
}

reaction{
   'O2 + O2 <=> O + O + O2',
   -- Andrienko and Boyd
   --fr={'Macheret', A=1.63019353e21, n=-1.319, T_a=63050, v_name='O2', c_name='O2'},
   --chemistry_energy_coupling={{species='O2', mode='vibration', model='Macheret', A=1.63019353e21, n=-1.319, T_d=63050,c_name='O2'}},
   -- Park
   fr={'Macheret', A=2.0e21, n=-1.5, T_a=59500, v_name='O2', c_name='O2'},
   chemistry_energy_coupling={{species='O2', mode='vibration', model='Macheret', A=2.0e21, n=-1.5, T_d=59500,c_name='O2'}},
   ec={model='from CEA curves',iT=0}
}




