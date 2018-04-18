function Z = SpecVParFun( T,E )
%data is spec data
Z = 0;
kb = 8.617328149741e-5;    % eV/K

for i = 1:length(E)
    Z_e = E(i).ne*exp(-E(i).Te/T);
    
    NV = E(i).Vmax+1;

    Z_v = zeros(1,NV);
    for l=1:NV
        Z_v(l) =exp(-E(i).Erv(l,1)/kb/T);
    end

    Z = Z + Z_e*sum(Z_v);
end