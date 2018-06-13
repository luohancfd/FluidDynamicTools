-- Author: Rowan J. Gollan
-- Date: 13-Mar-2009 (Friday the 13th)
-- Place: NIA, Hampton, Virginia, USA
--
-- History:
--   19-Mar-2009 :- added checking of mass balance
--                  and charge balance

module(..., package.seeall)

require 'reaction_rate'
require 'lex_elems'

local transform_rate_model = reaction_rate.transform_rate_model

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

-- lexical elements for parsing the whole reaction string
-- get common elements from lex_elems.lua
for k,v in pairs(lex_elems) do
  _G[k] = v
end
-- module-specific elements
local PressureDependent = Open * Plus * "M" * Space * Close
local function pdstring() return "pressure dependent" end

-- Grammar
local Participant = lpeg.V"Participant"
local Reaction = lpeg.V"Reaction"
local Mechanism = lpeg.V"Mechanism"

G = lpeg.P{ Mechanism,
Mechanism = lpeg.Ct(Reaction * ( RArrow + FArrow ) * Reaction),
Reaction = lpeg.Ct(Participant * (Plus * Participant)^0 * (PressureDependent / pdstring)^0 ) * Space,
Participant = lpeg.Ct(lpeg.C(Number^0) * Space * Species * Space)
   }

G = Space * G * -1

SpeciesG = lpeg.P{ Species }

function parse_reaction_string(s)
  t = lpeg.match(G, s)
  return t
end

function validate_reaction(t)
  if type(t[1]) ~= 'string' then
    print("There was an error when parsing reaction number: ", #reactions+1)
    print("The first entry should be a string denoting the reaction mechanism.")
    print("Bailing out!")
    os.exit(1)
  end

  reac = parse_reaction_string(t[1])
  if reac == nil then
    print("There was an error parsing the reaction string for reaction number: ", #reactions+1)
    print("It seems the string is badly formed.  The given string is: ")
    print(t[1])
    print("Bailing out!")
    os.exit(1)
  end

  mass, charge = check_equation_balances(reac, #reactions+1)

  if not mass then
    print("The mass does not balance.")
    print("Bailing out!")
    os.exit(1)
  end

  if not charge then
    print("The charge does not balance.")
    print("Bailing out!")
    os.exit(1)
  end

  return true
end

local Species2 = lpeg.Ct(lpeg.Ct((lpeg.C(Element) * lpeg.C(Number^0) * lpeg.C(PM^0) * Solid^0))^1)
Species2G = lpeg.P{ Species2 }

function check_equation_balances(r, rnumber)
  -- Checks both mass and charge balance
  local elem = {}
  local charge = 0

  for _,p in ipairs(r[1]) do
    if type(p) == 'table' then
      local coeff = tonumber(p[1]) or 1
      if p[2] ~= "M" then
        -- break species into elements
        local ets = lpeg.match(Species2, p[2])
        for _,et in ipairs(ets) do
          local e = et[1]
          local ne = tonumber(et[2]) or 1
          -- collate the element counts
          if e ~= "e" then -- don't add electron to mass balance
            if elem[e] then
              elem[e] = elem[e] + coeff*ne
            else
              elem[e] = coeff*ne
            end
          end
          -- collate the charge counts
          if et[3] == "+" then
            charge = charge + coeff*1
          elseif et[3] == "-" then
            charge = charge - coeff*1
          end
        end
      end
    end
  end

  -- now check on the other side of the reaction
  for _,p in ipairs(r[3]) do
    if type(p) == 'table' then
      local coeff = tonumber(p[1]) or 1
      if p[2] ~= "M" then
        -- break species into elements
        local ets = lpeg.match(Species2, p[2])
        for _,et in ipairs(ets) do
          local e = et[1]
          local ne = tonumber(et[2]) or 1
          -- collate the element counts
          if e ~= "e" then -- don't add mass of electron to mass balance
            if elem[e] then
              elem[e] = elem[e] - coeff*ne
            else
              elem[e] = -coeff*ne
            end
          end		  
          -- collate the charge counts
          if et[3] == "+" then
            charge = charge - coeff*1
          elseif et[3] == "-" then
            charge = charge + coeff*1
          end
        end
      end
    end
  end

  -- mass check
  mass_balances = true
  for k,v in pairs(elem) do
    if v ~= 0 then
      mass_balances = false
      print("There is a problem with the mass balance for reaction: ", rnumber)
      print("In particular, the element: ", k)
      print("does not balance on both sides of the reaction equation.")
    end
  end

  -- charge check
  charge_balances = true
  if charge ~= 0 then
    charge_balances = false
    print("There is a problem with the charge balance for reaction: ", rnumber)
  end

  return mass_balances, charge_balances

end

function transform_reaction(t, species, suppress_warnings)
  r = {}
  r.equation = t[1]
  r.type = "normal reaction"

  rs = parse_reaction_string(t[1])
  third_body = false
  -- deal with forward elements

  f_coeffs = {}
  for _,p in ipairs(rs[1]) do
    if type(p) == 'table' then
      p[1] = tonumber(p[1]) or 1
      p[2] = transform_species_str(p[2])
      if p[2] == "M" then
        third_body = true
      else 
        sp_index = species[p[2]]
        if sp_index == nil then
          print("The following species has been declared in a reaction: ", p[2])
          print("but is not part of the declared gas model.")
          print("This occurred for reaction number: ", t.number)
          if t.label then
            print("label: ", t.label)
          end
          print("Bailing out!")
          os.exit(1)
        end
        -- check this is not already in f_coeffs
        found = false
        for _,e in ipairs(f_coeffs) do
          if e[1] == sp_index then
            found = true
            e[2] = e[2] + p[1]
          end
        end

        if not found then
          f_coeffs[#f_coeffs+1] = {sp_index, p[1]}
        end
      end
    end
  end

  b_coeffs = {}
  if ( rs[2] == "<=>" ) then
    -- We have a forward and reverse reaction
    -- do the same as above for b_coeffs
    for _,p in ipairs(rs[3]) do
      if type(p) == 'table' then
        p[1] = tonumber(p[1]) or 1
        p[2] = transform_species_str(p[2])
        if p[2] == "M" then
          third_body = true
        else 
          sp_index = species[p[2]]
          if sp_index == nil then
            print("The following species has been declared in a reaction: ", p[2])
            print("but is not part of the declared gas model.")
            print("This occurred for reaction number: ", t.number)
            if t.label then
              print("label: ", t.label)
            end
            print("Bailing out!")
            os.exit(1)
          end
          -- check this is not already in f_coeffs
          found = false
          for _,e in ipairs(b_coeffs) do
            if e[1] == sp_index then
              found = true
              e[2] = e[2] + p[1]
            end
          end

          if not found then
            b_coeffs[#b_coeffs+1] = {sp_index, p[1]}
          end
        end
      end
    end
  end

  r.f_coeffs = f_coeffs
  r.b_coeffs = b_coeffs
  r.third_body = third_body

  if t.fr then
    r.frc = transform_rate_model(t.fr, f_coeffs, third_body)
  end
  if t.br then
    r.brc = transform_rate_model(t.br, b_coeffs, third_body)
  end
  if t.ec then
    r.ec = t.ec
  else
    -- By default
    r.ec = {model="from thermo",iT=0}
  end

  -- Look for presence of pressure dependent reaction
  pressure_dependent = false
  for _,p in ipairs(rs[1]) do
    if type(p) == 'string' and p == 'pressure dependent' then
      pressure_dependent = true
    end
  end
  for _,p in ipairs(rs[3]) do
    if type(p) == 'string' and p == 'pressure dependent' then
      pressure_dependent = true
    end
  end

  -- Deal with any third body efficiencies (if needed)
  if third_body or pressure_dependent then
    if third_body then
      r.type = "third body reaction"
    end
    -- All efficiencies are set to 1.0
    r.efficiencies = {}
    for i=1,species.size do
      -- the { , } table needs to have C++ indices in it
      r.efficiencies[i] = {i-1, 1.0}
    end
    -- Next look at the special cases
    if t.efficiencies then
      for k,v in pairs(t.efficiencies) do
        if not species[k] then
          if not suppress_warnings then
            print("WARNING: One of the species given in the efficiencies list")
            print("is NOT one of the species in the gas model.")
            print("The efficiency for: ", k, " will be skipped.")
            print("This occurred for reaction number: ", t.number)
            if t.label then
              print("label: ", t.label)
            end
          end
        else
          -- species[k] gives back a C++ index,
          -- need to insert it at +1
          r.efficiencies[species[k]+1] = {species[k], v}
        end
      end
    end
    -- For efficient execution in inner loops,
    -- remove values whose efficiency is zero.
    local essentially_zero = 1.0e-9
    for i=#r.efficiencies,1,-1 do
      if r.efficiencies[i][2] <= essentially_zero then
        table.remove(r.efficiencies, i)
      end
    end

    if pressure_dependent then
      if r.frc then
        if r.frc.model == "pressure dependent" then
          r.frc.efficiencies = r.efficiencies
        end
      end
      if r.brc then
        if r.brc.model == "pressure dependent" then
          r.brc.efficiencies = r.efficiencies
        end
      end
    end

  end

  -- Look for chemistry_energy_coupling field
  if t.chemistry_energy_coupling then
    r.chemistry_energy_coupling = t.chemistry_energy_coupling
    if t.fr[1] == 'Macheret' then
      for i,p in ipairs(t.chemistry_energy_coupling) do
        if p['c_name'] then
          r.chemistry_energy_coupling[i]['c_name'] = transform_species_str(p['c_name'])
        end
      end
    end
  end

  return r

end
