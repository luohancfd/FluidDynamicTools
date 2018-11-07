-- a little script to evaluate gamma at 300.0 K
-- based on CEA curves.

-- N2_plus
coeffs = {  0.000000000e+00,  0.000000000e+00,  2.500000000e+00,
                0.000000000e+00,  0.000000000e+00,  0.000000000e+00,
                0.000000000e+00, -7.453750000e+02, -1.172081224e+01
    }
    
R = 8.31451/ 0.000548579903e-3
print("R= ", R)
T = 300.0

function eval_Cp(a, T)
   val = a[1]/(T^2) + a[2]/T + a[3] + a[4]*T
   val = val + a[5]*T^2 + a[6]*T^3 + a[7]*T^4
   return val*R
end

Cp = eval_Cp(coeffs, T)

gamma = Cp / (Cp - R)

print("gamma= ", gamma)
