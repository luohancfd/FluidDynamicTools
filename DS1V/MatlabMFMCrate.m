function [ Eqrate,EvRem,NEvRem,ProbM ] = MatlabMFMCrate(T,Tv,sp1,sp2,N,D,method,SHO)
%% sample prob
persistent N2 O2 N2thetaV O2thetaV
ev2K = 11604.319869319335;
% calculate rate from Matlab sampling
N2 = [0.06357786, 0.35378524, 0.63986163, 0.92182002, 1.19967443, ...
    1.47343384, 1.74310826, 2.00870469, 2.27022813, 2.52768458, ...
    2.78107604, 3.03040450, 3.27566898, 3.51686846, 3.75400095, ...
    3.98706045, 4.21604296, 4.44094048, 4.66174501, 4.87844654, ...
    5.09103309, 5.29949364, 5.50381220, 5.70397377, 5.89996135, ...
    6.09175694, 6.27933854, 6.46268514, 6.64177276, 6.81657639, ...
    6.98706902, 7.15322066, 7.31499932, 7.47237298, 7.62530565, ...
    7.77375833, 7.91769102, 8.05705973, 8.19181844, 8.32191716, ...
    8.44730289, 8.56791863, 8.68370238, 8.79458914, 8.90050692, ...
    9.00137970, 9.09712390, 9.18764890, 9.27285572, 9.35263555, ...
    9.42686849, 9.49542104, 9.55814401, 9.61486829, 9.66540048, ...
    9.70951598, 9.74694820, 9.77737265, 9.80037704, 9.81540013];
N2thetaV = 3352;
O2 = [0.09744979, 0.29064938, 0.48124897, 0.66944856, 0.85494817, ...
    1.03784777, 1.21794739, 1.39524701, 1.56964663, 1.74104627, ...
    1.90934591, 2.07464555, 2.23674520, 2.39554486, 2.55104453, ...
    2.70314420, 2.85174389, 2.99674358, 3.13814327, 3.27584298, ...
    3.40974269, 3.53984241, 3.66584214, 3.78784188, 3.90564163, ...
    4.01924138, 4.12854115, 4.23335092, 4.33361071, 4.42920050, ...
    4.52000031, 4.60589013, 4.68672995, 4.76237979, 4.83269964, ...
    4.89751950, 4.95666937, 5.00996926, 5.05723916, 5.09828907, ...
    5.13292100, 5.16098794, 5.18241589, 5.19734186, 5.20636404, ...
    5.21081273, 5.21244607];
O2thetaV = 2273.54;

if nargin == 7
    SHO =0;
end

if SHO == 1 && method == 2
    error("unsupported options");
end

dia = 1;
Evlevel2 = 0; thetaV2=0;
if strcmp(sp2,'N2')
    Evlevel2 = N2;
    thetaV2 = N2thetaV;
elseif strcmp(sp2,'O2')
    Evlevel2 = O2;
    thetaV2 = O2thetaV;
else
    dia = 0;
end

if strcmp(sp1,'N2')
    Evlevel = N2;    thetaV = N2thetaV;
else
    Evlevel = O2;    thetaV = O2thetaV;
end

colldat = struct('sp',{'N','O','O2','N2'},...
           'd' ,{3e-10, 3.458e-10, 3.985e-10, 4.17e-10},...
           'm', {14, 16, 32, 28},...
           'm2',{14,16,16,14},...
           'w' ,{0.71, 0.76, 0.71, 0.74});
splist = {'N','O','O2','N2'};
omega = (colldat(strcmp(splist,sp1)).w+colldat(strcmp(splist,sp2)).w)/2;
% omega = (colldat(m).w+colldat(n).w)/2;

if (strcmp(sp1,'O2') && strcmp(sp2,'O'))
    omega = 0.75;
elseif (strcmp(sp1,'O2') && strcmp(sp2,'O2'))
    omega = 0.75;
end
alpha = (colldat(strcmp(splist,sp1)).m2/(colldat(strcmp(splist,sp1)).m2+colldat(strcmp(splist,sp2)).m2))^2;
% alpha = (colldat(m).m2/(colldat(m).m2+colldat(n).m2))^2;
ProbM = zeros(1,length(Evlevel));
%% method 1, sample vibrational energy from Boltzmann distribution
if method == 1
    % sample vibrational energy, rotational and translational
    % level for molecule 1
    if SHO == 0
        vibpar = exp(-Evlevel*ev2K/Tv);
        vibpar = vibpar/sum(vibpar);
        vibpar = cumsum(vibpar);
        X = rand(N,1);
        Y = discretize(X,[0 vibpar],'IncludedEdge','right');
        %     Y = discretize_bin(X,[0,vibpar]);
        Ev = Evlevel(Y)*ev2K;
        Ev = Ev';
    elseif SHO == 1
        Ev = -log(rand(N,1))*T;  % in K
        Ev = floor(Ev/thetaV)*thetaV;
    else
        Ev = -log(rand(N,1))*T;  % in K
    end
    Er = -log(rand(N,1))*T;  % in K
    
    Ev2 = zeros(N,1); Er2 = Ev2;
    if dia == 1
        if SHO == 0
            vibpar2 = exp(-Evlevel2*ev2K/Tv);
            vibpar2 = vibpar2/sum(vibpar2);
            vibpar2 = cumsum(vibpar2);
            X = rand(N,1);
            Y = discretize(X,[0 vibpar2],'IncludedEdge','right');
            %         Y = discretize_bin(X,[0,vibpar2]);
            Ev2 = Evlevel2(Y)*ev2K;
            Ev2 = Ev2';
        elseif SHO ==1
            Ev2 = -log(rand(N,1))*T;  % in K
            Ev2 = floor(Ev2/thetaV2)*thetaV2;
        else
            Ev2 = -log(rand(N,1))*T;  % in K
        end
        Er2 = -log(rand(N,1))*T;  % in K
    end
    Et = random('Gamma',(2.5-omega)*ones(N,1),1)*T; %in K
    % sample angles
    theta = rand(N,1)*pi;
    phi = rand(N,1)*pi;
    gam1 = rand(N,1)*pi;    
    gam2 = rand(N,1)*pi;
    % gam1 = acos(1-2*rand(N,1));
    beta1 = zeros(N,1); phi2 = zeros(N,1); delta = zeros(N,1);  beta2 = zeros(N,1);
    if dia == 1
        beta1 = rand(N,1)*pi;
        phi2 = rand(N,1)*pi;
        delta = rand(N,1)*2*pi;
        theta = 2*theta;
        beta2 = rand(N,1)*2*pi;
    end
    Dstar = D - Er + 2*Er.^(1.5)/(3*sqrt(6*D));
    a = 1-sqrt(alpha);
    b = 1+sqrt(alpha);
    c = Dstar - Ev.*sin(phi).^2;
    F = zeros(N,1);
    i2 = find(c > 0);

    c(i2) = sqrt(c(i2))+sqrt(Ev(i2)).*cos(phi(i2));
    F(i2) = c(i2)./cos(theta(i2))/a - sqrt(Ev(i2)).*cos(phi(i2)).*cos(theta(i2));
    if dia == 1
        d = sqrt(Er2(i2)).*(cos(delta(i2)).*cos(beta2(i2)).*sin(beta1(i2))+sin(beta2(i2)).*sin(delta(i2)));
        e = sqrt(Ev2(i2)).*cos(phi2(i2)).*cos(beta1(i2)).*cos(beta2(i2));
        F(i2) = F(i2) - sqrt(sqrt(alpha)/a)*(d-e);
    end
    F(i2) = a*F(i2).^2./cos(gam1(i2)).^2./cos(gam2(i2)).^2;
    if dia == 0
        F(i2) = F(i2)/b;
    end
    index = Et>F;
    prob = sum(index)/N;
    EvRem = sum(Ev(index))/ev2K;
else
    parfor i = 1:length(Evlevel)
        if Evlevel(i)<alpha*D
            Nsample = 1e6;
        else
            Nsample = 1e5;
        end
        if dia == 0
            ProbM(i) = ProbSampleTR(Evlevel(i),Nsample,T,D,alpha,omega);
        else
            ProbM(i) = ProbSampleTR(Evlevel(i),Nsample,T,D,alpha,omega,Evlevel2,Tv);
        end
    end
    vibpar = exp(-Evlevel*ev2K/Tv);
    vibpar = vibpar/sum(vibpar);
    prob = ProbM*vibpar';
    EvRem = Evlevel*(ProbM'.*vibpar')/prob;
    % calculate rate for each vibrational level
end
%% collision rate
[~,crate] = CollisionRate(T,sp1,sp2);
Eqrate = crate*prob;
NEvRem = EvRem*ev2K/D;
end