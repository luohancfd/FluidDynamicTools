function [ Kc,Kp ] = LewisEqConstant(T,spl,spr,data)
%% calculate equilibrium constant for at temperature T
% spl - > spr
R = 8.3144598/1e3; %J/mol/K
kb = 1.38064852e-23;
Kp = zeros(1,length(T));
Kc = zeros(1,length(T));
for i = 1:length(T)
    G0 = 0;
    eta = 0;
    [~,mspl] = size(spl);
    Gl = zeros(1,length(mspl));
    for ii = 1:mspl
        [ ~,~,Gl(ii),~ ] = Lewis(spl{1,ii},T(i),data);
        G0 = G0 - Gl(ii)*spl{2,ii};
        eta = eta - spl{2,ii};
    end
    
    [~,mspr] = size(spr);
    Gr = zeros(1,length(mspr));
    for ii = 1:mspr
        [ ~,~,Gr(ii),~ ] = Lewis(spr{1,ii},T(i),data);
        G0 = G0 + Gr(ii)*spr{2,ii};
        eta = eta + spr{2,ii};
    end    
    
    Kp(i) = exp(-G0/R/T(i)); %pressure in bar
    Kc(i) = Kp(i)/(kb*T(i)/1e5)^eta; % number density in m^(-3)
end    


end

