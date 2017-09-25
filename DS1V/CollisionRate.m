function [ v,v2 ] = CollisionRate(T,sp1,sp2)
%% calculate VHS collision rates
dat = struct('sp',{'N','O','O2','N2'},...
           'd' ,{3e-10, 3.458e-10, 3.985e-10, 4.17e-10},...
           'm', {14, 16, 32, 28},...
           'w' ,{0.71, 0.76, 0.71, 0.74});
splist = {'N','O','O2','N2'};
% m = find(strcmp(splist,sp1));
% n = find(strcmp(splist,sp2));
m=1;n=1;
for i = 1:length(splist)
    if strcmp(splist{i},sp1) == 1
        m = i;
        break
    end
end
for i = 1:length(splist)
    if strcmp(splist{i},sp2) == 1
        n = i;
        break
    end
end

d = 0.5*(dat(m).d + dat(n).d);
wave = (dat(m).w+dat(n).w)/2;

m1 = dat(m).m;
m2 = dat(n).m;
mr = m1*m2/(m1+m2)*1.660539040e-27;
kb = 1.38064852e-23;
Tref = 273.15;
v = 2*sqrt(pi)*d^2*(T/Tref).^(1-wave)*(2*kb*Tref/mr)^0.5;  %m^3/s
v2 = v*1e6*6.023e23;  %cm^3/mol/s
end

