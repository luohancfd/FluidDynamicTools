function Z = SpecParFun( T,E )
%data is spec data
Z = 0;
kb = 8.617328149741e-5;    % eV/K

for i = 1:length(E)
    Z_e = E(i).ne*exp(-E(i).Te/T);
    
    NV = E(i).Vmax+1;
    NJM = E(i).Jmax(1)+1;
    
    Jmax = E(i).Jmax;
    Z_rv = zeros(NV,NJM);
    
    indexm = false(NV,NJM);
    for l=1:NV
        for j=1:Jmax(l)+1
            indexm(l,j) = true;
        end
    end
    
    ns = repmat(reshape(E(i).ns,1,[]),[NV,ceil(NJM/2)]);
    ns = ns(:,1:NJM);
    ns(~indexm) = 0;   %nuclear degeneracy
    rotg = repmat(0:Jmax(1),[NV,1]);
    rotg = 2*rotg+1;
    ns = ns.*rotg;
    
    Z_rv= ns.*exp(-E(i).Erv/kb/T);
    
    Z_rv(~indexm)=0;

    Z = Z + Z_e*sum(sum(Z_rv));
end