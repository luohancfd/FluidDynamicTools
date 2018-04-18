function [ mu,k,D,omega11,omega12,omega13,omega22 ] = Capitelli(T,p,sp1,sp2,data)
%% calculate visocity and diffusion coefficient from Esposito's paper 
% JTHT April Jun 2000
cal2J = 4.1868;
% T in K, p in atm

data_len = length(data);
find_sp = false;
for i = 1:data_len
    if strcmp(sp1,data(i).sp1) && strcmp(sp2,data(i).sp2)
        find_sp = true;
        break
    elseif strcmp(sp2,data(i).sp1) && strcmp(sp2,data(i).sp2)
        find_sp = true;
        break
    end
end
if find_sp
    index = i;
else
    error('There is no data for %s-%s',sp1,sp2);
end
a = data(index).a;
m1 = data(index).m1;
m2 = data(index).m2;
mu = m1*m2/(m1+m2);
CollInt= @(T,b) (b(1)+b(2)*T.^b(3))./(b(4)+b(5)*T.^(b(6)));
%% Collision integral, with unit Angstrom^2, already divided by pi
omega11 = CollInt(T,a(1,:));
omega12 = CollInt(T,a(2,:));
omega13 = CollInt(T,a(3,:));
omega22 = CollInt(T,a(4,:));
% thermal conductivity: in unit cal/cm/sec/K
k = 1989.1e-7*sqrt(T/mu/2)./omega22;
% diffusion coefficient: in cm^2/s
D = 0.002628*sqrt(T.^3/mu/2)./p./omega11;
% viscosity: in unit g/cm/s
mu = 266.93e-7*sqrt(T*mu*2)./omega22;

end

