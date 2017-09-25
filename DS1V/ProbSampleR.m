function [ prob,Er ] = ProbSampleR(Ev,Et,N,T,D,alpha)
%% sample prob
% sample to get probability of reaction for given Ev and Et, Er is
% obtainted from temperature
% Ev in eV
% Et in kT
Ev = Ev*11604.319869319335;
Et = Et*T;
Er = -log(rand(N,1))*T;  % in K

theta = rand(N,1)*pi;
phi = rand(N,1)*pi;
gam1 = rand(N,1)*pi;
% gam1 = acos(1-2*rand(N,1));
gam2 = rand(N,1)*pi;
Dstar = D - Er + 2*Er.^(1.5)/(3*sqrt(6*D));


a = 1-sqrt(alpha);
b = 1+sqrt(alpha);
c = sqrt(Dstar - Ev*sin(phi).^2) + sqrt(Ev)*cos(phi);
F = a/b./cos(gam1).^2./cos(gam2).^2.*(c./cos(theta)/a - sqrt(Ev)*cos(phi).*cos(theta)+1).^2;
for i = 1:N
    if (~isreal(F(i)))
        F(i) = 0;
    end
end
prob = sum(Et > F)/N;


end

