clear
clc
load Parameter.mat Erv Ev prop
Trange = 300:50:6000;
Na = 6.022140857e23;
kb = 8.617328149741e-5;    % eV/K
ThetaV2 = prop.we(1)/kb;
ThetaR = 2.88;  % N2
R0 = Na*prop.k; %8.3144598; %J mol-1 K -1
R = R0*1000/28;  %J kg-1 K-1
TD = 78740*1.42879; %dissociation from eilmer3

conver = [8065.479 349.755 83.5935 0.695028792]; % 1 eV kcal/mol kJ/mol K Hato cm-1
Cprv = zeros(1,length(Trange));    %AHO+RR
Cprv2 = zeros(1,length(Trange));   %Rvcouple
Cprv3 = zeros(1,length(Trange));   %SHO+RR
Cprv4 = zeros(1,length(Trange));   %Spectro
Cprv5 = zeros(1,length(Trange));   %Spectro+RR


%% QCT Vibrational level with RR rotational
Prv = zeros(3,length(Trange));
dT  = 1;
tdata =[];
tdata.varnames = {'Temperature (K)','C<sub>p</sub>/R'};
for i = 1:length(Trange)
    Prv(2,i) = MolParFun(Trange(i),1,Erv,prop.Vmax,prop.Jmax,1);
    Prv(1,i) = MolParFun(Trange(i)-dT,1,Erv,prop.Vmax,prop.Jmax,1);
    Prv(3,i) = MolParFun(Trange(i)+dT,1,Erv,prop.Vmax,prop.Jmax,1);
    Cprv(i) = R0*((2*Trange(i)* (log(Prv(3,i)/Prv(1,i)))/(2*dT)+ ...
        Trange(i)^2*log(Prv(3,i)/Prv(2,i)^2*Prv(1,i))/(dT*dT))+1+3/2+1);
end
plot(Trange,Cprv/R0,'DisplayName','Rigid Rotator')
tdata.lines(1).x = Trange;
tdata.lines(1).y = Cprv/R0;
tdata.lines(1).zonename = 'AHO(Gamallo PES)+RR';
hold on

%% QCT Rovibrational level
Prv = zeros(3,length(Trange));
for i = 1:length(Trange)
    Prv(2,i) = MolParFun(Trange(i),1,Erv,prop.Vmax,prop.Jmax,0);
    Prv(1,i) = MolParFun(Trange(i)-dT,1,Erv,prop.Vmax,prop.Jmax,0);
    Prv(3,i) = MolParFun(Trange(i)+dT,1,Erv,prop.Vmax,prop.Jmax,0);
    Cprv2(i) = R0*((2*Trange(i)* (log(Prv(3,i)/Prv(1,i)))/(2*dT)+ ...
        Trange(i)^2*log(Prv(3,i)/Prv(2,i)^2*Prv(1,i))/(dT*dT))+1+3/2);
end
plot(Trange,Cprv2/R0,'DisplayName','Ro-vib coupling')
tdata.lines(2).x = Trange;
tdata.lines(2).y = Cprv2/R0;
tdata.lines(2).zonename = 'Gamallo PES';

%% QCT Simple harmonic oscillator + RR
Tv = ThetaV2./Trange;
parelec = 0;
Cprv3= R0*(Tv.^2.*exp(Tv)./(exp(Tv)-1).^2+3/2+1+1+parelec);
plot(Trange,Cprv3/R0,'DisplayName','SHO+Rigid Rotator')
tdata.lines(3).x = Trange;
tdata.lines(3).y = Cprv3/R0;
tdata.lines(3).zonename = 'SHO+RR';

%% QCT Eilmer3 SHO
TVeilmer = 2358.569*1.42879;
Tv = TVeilmer./Trange;
parelec = 0;
Cprv7= R0*(Tv.^2.*exp(Tv)./(exp(Tv)-1).^2+3/2+1+1+parelec);
plot(Trange,Cprv7/R0,'DisplayName','ESHO+Rigid Rotator')
% tdata.lines(3).x = Trange;
% tdata.lines(3).y = Cprv7/R0;
% tdata.lines(3).zonename = 'ESHO+RR';
hold on

%% QCT Eilmer3 truncated SHO
TVeilmer = 2358.569*1.42879;
Tv = TVeilmer./Trange;
parelec = 0;
Cprv6= R0*(Tv.^2.*exp(Tv)./(exp(Tv)-1).^2+3/2+1+1+parelec);
Cprv6= Cprv6-R0*(Tv.^2.*exp(TD./Trange)./(exp(TD./Trange)-1).^2);


plot(Trange,Cprv6/R0,'DisplayName','ETSHO+Rigid Rotator')
% tdata.lines(3).x = Trange;
% tdata.lines(3).y = Cprv3/R0;
% tdata.lines(3).zonename = 'SHO+RR';

%% Spectropic data
% TermSymbol Te | we wexe weye | Be alphae gammae  | De betae
Prv = zeros(3,length(Trange));
data = ReadNIST('NIST.dat');
Spec = find(ismember({data.Particle},'N2'));
Spec = data(Spec).data;
diss = 9.759*conver(1); % in cm-1
E    = SpecTableEnergy(Spec,diss);
E(1).ns = [6 3];
E(2).ns = [3 6];
E(3).ns = [6 3];
for i = 1:length(Trange)
    Prv(2,i) = SpecParFun(Trange(i),E);
    Prv(1,i) = SpecParFun(Trange(i)-dT,E);
    Prv(3,i) = SpecParFun(Trange(i)+dT,E);
    Cprv4(i) = R0*((2*Trange(i)* (log(Prv(3,i)/Prv(1,i)))/(2*dT)+ ...
        Trange(i)^2*log(Prv(3,i)/Prv(2,i)^2*Prv(1,i))/(dT*dT))+1+3/2);
end
plot(Trange,Cprv4/R0,'DisplayName','Spectroscopy')
tdata.lines(4).x = Trange;
tdata.lines(4).y = Cprv4/R0;
tdata.lines(4).zonename = 'Spectroscopy';
%% Spectropic data with RR
Prv = zeros(3,length(Trange));
for i = 1:length(Trange)
    Prv(2,i) = SpecVParFun(Trange(i),E);
    Prv(1,i) = SpecVParFun(Trange(i)-dT,E);
    Prv(3,i) = SpecVParFun(Trange(i)+dT,E);
    Cprv5(i) = R0*((2*Trange(i)* (log(Prv(3,i)/Prv(1,i)))/(2*dT)+ ...
        Trange(i)^2*log(Prv(3,i)/Prv(2,i)^2*Prv(1,i))/(dT*dT))+1+3/2+1);
end
plot(Trange,Cprv5/R0,'DisplayName','Spectroscopy+RR')
tdata.lines(5).x = Trange;
tdata.lines(5).y = Cprv5/R0;
tdata.lines(5).zonename = 'Spectroscop+RR';
 
%% NIST fit
Trange = 300:50:6000;
plot(Trange,NistN2(Trange)/R0,'k--','DisplayName','NIST')
tdata.lines(6).x = Trange;
tdata.lines(6).y = NistN2(Trange)/R0;
tdata.lines(6).zonename = 'NIST-JANAF';

xlabel('Temperature (K)')
ylabel('C_p/R')
legend show
mat2tecplot(tdata,'N2.plt')