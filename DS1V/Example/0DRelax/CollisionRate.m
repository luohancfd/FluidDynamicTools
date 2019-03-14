function [ v2,rv,v,A,eta,Ta ] = CollisionRate(T,sp1,sp2,wave,d)
%% calculate VHS collision rates
Na = 6.022140857e23;
dat = struct('sp',{'N','O','O2','N2'},...
           'd' ,{3e-10, 3.458e-10, 3.985e-10, 4.17e-10},...
           'm', {14, 16, 32, 28},...
           'w' ,{0.71, 0.76, 0.71, 0.74});
react = struct('sp',{'N2+O','N2+N','O2+O','N2+N2','O2+O2'},...
               'A',{8.934e-7*Na,6.6831e-6*Na,2.5e18,4.5e-6*Na,8.132E-10*Na*10},...
               'eta',{-0.4807,-0.6996,-0.565,-0.675,-0.131},...
               'Ta',{113950,1139200,60491.819,117000,59380});
%% make diatom be the first
if ~contains(sp1,'2')
    sp3 = sp1;
    sp1 = sp2;
    sp2 = sp3;
end
%% get VHS parameter
splist = {dat.sp};
m = find(strcmp(splist,sp1));
n = find(strcmp(splist,sp2));
if nargin == 3
    d = 0.5*(dat(m).d + dat(n).d);
    wave = (dat(m).w+dat(n).w)/2;
    if (strcmp(sp1,'O2') && strcmp(sp2,'O'))
        wave = 0.75; d = 3.7215e-10;
    elseif (strcmp(sp1,'O2') && strcmp(sp2,'O2'))
        wave = 0.7318;
        d = 4.1515e-10;
    end
end
% elseif (strcmp(sp1,'N2') && strcmp(sp2,'O'))
%     wave = 0.7468;
%     d = 3.63858e-10;
% end
%% Get reaction parameter
spcombine = [sp1,'+',sp2];
l = find(strcmp({react.sp},spcombine));

m1 = dat(m).m;
m2 = dat(n).m;
mr = m1*m2/(m1+m2)*1.660539040e-27;
kb = 1.38064852e-23;
Tref = 273.15;
v = 2*sqrt(pi)*d^2*(T/Tref).^(1-wave)*(2*kb*Tref/mr)^0.5;  %m^3/s
v2 = v*1e6*6.023e23;  %cm^3/mol/s

rv = react(l).A*T.^react(l).eta.*exp(-react(l).Ta./T);
A = react(l).A;
eta = react(l).eta;
Ta = react(l).Ta;
end

