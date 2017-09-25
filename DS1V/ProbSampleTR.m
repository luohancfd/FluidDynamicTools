function [ prob ] = ProbSampleTR(Ev,N,T,D,alpha,omega,Evlevel2,Tv)
%% sample prob
% sample to get probability of reaction for given Ev
% Er and Et is obtainted from temperature
% Ev in eV
% Et in kT
ev2k = 11604.319869319335;
Ev = Ev*ev2k;
Er = -log(rand(N,1))*T;  % in K
Et = random('Gamma',(2.5-omega)*ones(N,1),1)*T; %in K
if nargin == 6
    dia = 0;
else
    dia = 1;
    vibpar2 = exp(-Evlevel2*ev2k/Tv);
    vibpar2 = vibpar2/sum(vibpar2);
    vibpar2 = cumsum(vibpar2);
    X = rand(N,1);
    Y = discretize(X,[0 vibpar2],'IncludedEdge','right');
%     Y = discretize_bin(X,[0,vibpar2]);
    Ev2 = Evlevel2(Y)*ev2k;
    Ev2 = Ev2';
    Er2 = -log(rand(N,1))*T;  % in K
end

theta = rand(N,1)*pi;
phi = rand(N,1)*pi;
gam1 = rand(N,1)*pi;
% gam1 = acos(1-2*rand(N,1));
gam2 = rand(N,1)*pi;
beta1 = zeros(N,1); phi2 = zeros(N,1); delta = zeros(N,1);  beta2 = zeros(N,1);
if dia == 1
    beta1 = rand(N,1)*pi;
    phi2 = rand(N,1)*pi;
    delta = rand(N,1)*2*pi;
    theta = 2*theta;
    beta2 = rand(N,1)*2*pi;
end
F = zeros(N,1);

Dstar = D - Er + 2*Er.^(1.5)/(3*sqrt(6*D));

a = 1-sqrt(alpha);
b = 1+sqrt(alpha);
c = Dstar - Ev*sin(phi).^2;

i2 = find(c > 0);
c(i2) = sqrt(c(i2))+sqrt(Ev).*cos(phi(i2));
F(i2) = c(i2)./cos(theta(i2))/a - sqrt(Ev).*cos(phi(i2)).*cos(theta(i2));
if dia == 1
    d = sqrt(Er2(i2)).*(cos(delta(i2)).*cos(beta2(i2)).*sin(beta1(i2))+sin(beta2(i2)).*sin(delta(i2)));
    e = sqrt(Ev2(i2)).*cos(phi2(i2)).*cos(beta1(i2)).*cos(beta2(i2));
    F(i2) = F(i2) - sqrt(sqrt(alpha)/a)*(d-e);
end
F(i2) = a*F(i2).^2./cos(gam1(i2)).^2./cos(gam2(i2)).^2;
if dia == 0
    F(i2) = F(i2)/b;
end

prob = sum(Et > F)/N;


end

