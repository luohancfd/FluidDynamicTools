-- Author: Rowan J. Gollan
-- Date: 13-Mar-2009 (Friday the 13th)
-- Place: NIA, Hampton, Virginia, USA

module(..., package.seeall)

local k_Boltz = 1.3806505e-23

function transform_species_str(sp)
  if string.match(sp, '+') then
    return string.gsub(sp, '+', '_plus')
  end
  if sp == 'e-' then
    return 'e_minus'
  end
  -- In all other cases return string unaltered
  return sp
end

function transform_rate_model(t, participants, third_body)
   local nc = 0
   for _,p in ipairs(participants) do
      nc = nc + p[2]
   end
   
   if third_body then nc = nc + 1 end

   local base = 1.0e-6
   local conv_factor = base^(nc-1)

   m = {}
   m.model = t[1]

   if m.model == "Arrhenius" then
      m.A = t.A * conv_factor
      m.n = t.n
      m.E_a = t.T_a * k_Boltz
   elseif m.model == "Park" then
      m.A = t.A * conv_factor
      m.n = t.n
      m.E_a = t.T_a * k_Boltz
      m.p_name = t.p_name
      m.p_mode = t.p_mode
      m.s_p = t.s_p
      m.q_name = t.q_name
      m.q_mode = t.q_mode
   elseif m.model == "Macheret" then
      m.A = t.A * conv_factor
      m.n = t.n
      m.E_a = t.T_a * k_Boltz
      m.v_name = transform_species_str(t.v_name)
      m.c_name = transform_species_str(t.c_name)
      m.khigh = t.khigh
   elseif m.model == "pressure dependent" then
      m.k_inf = {}
      m.k_inf.A = t.k_inf.A * conv_factor
      m.k_inf.n = t.k_inf.n
      m.k_inf.E_a = t.k_inf.T_a * k_Boltz
      
      conv_factor_low = base^nc
      m.k_0 = {}
      m.k_0.A = t.k_0.A * conv_factor_low
      m.k_0.n = t.k_0.n
      m.k_0.E_a = t.k_0.T_a * k_Boltz

      m.Troe = t.Troe
   elseif m.model == "MarroneTreanor" then
      m.A = t.A * conv_factor
      m.n = t.n
      m.E_a = t.T_a * k_Boltz
      m.v_name = t.v_name
      m.U = t.U
   elseif m.model == "Knab_et_al" then
      m.A = t.A * conv_factor
      m.n = t.n
      m.E_a = t.T_a * k_Boltz
      m.v_name = t.v_name
      m.U0 = t.U0
      m.U1 = t.U1
      m.alpha = t.alpha
      m.Z_limit = t.Z_limit
   else
      print("Reaction rate coefficient model: ", t[1])
      print("is unknown.")
      print("Bailing out!")
      os.exit(1)
   end

   return m
end
