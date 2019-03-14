clear
clc
folder = 'ZV_Corrected';
VHS_DS1V = [4.1515 0.7318];
Tv0 = 500;
T_run = [600 1000 1500 2000 4000 6000 8000 10000 15000 1750 20000];
nd = 1e21;
tauV = zeros(1,length(T_run));
tauColl = zeros(1,length(T_run));

result = repmat(struct('time',[]),1,length(T_run));
Zv = zeros(1,length(T_run));
for i = 1:length(T_run)
    result(i).file = sprintf('./%s/RELAX_%03d.DAT',folder,i-1);
    data = load(result(i).file);
    result(i).time = data(:,1);
    result(i).T = data(:,2);
    result(i).Tv = data(:,3);
    result(i).Ev = data(:,38);
    result(i).EvEq = data(:,39);
    result(i).tauColl = data(:,end);
    result(i).tauCollTheory = 1./(CollisionRate(result(i).T,'O2','O2',VHS_DS1V(2),VHS_DS1V(1)*1e-10)/6.023e23/1e6*nd);
%     plot(time,tauColl,time,tauCollTheory);
%     plot(time, (Tv-Tv0)./(T_run(i)-Tv0));
    phi = (result(i).Ev - result(i).EvEq)./(result(i).Ev(1) - result(i).EvEq);
%     semilogx(result(i).time, phi);
    result(i).phi = phi;
    if phi(end) < 1/exp(1)
        index = find(phi>0.01);
        coeff = polyfit(result(i).time(index), log(phi(index)),1);
        tauVT = -1/coeff(1);
    else
        tauVT = result(i).tauColl(end);
    end
    result(i).Zv = tauVT/mean(result(i).tauColl);
    Zv(i) = result(i).Zv;
end
Tinv3 = T_run.^(-1/3);
ZvT = 1.767075e+05*Tinv3.^4 -6.347191e+04*Tinv3.^3 + 8.422653e+03*Tinv3.^2 -4.179710e+02*Tinv3.^1 + 8.891025e+00;
ZvT = 10.^(ZvT);

% sort T_run
[T_run, I] = sort(T_run);
Zv = Zv(I);
ZvT = ZvT(I);
result = result(I);
Tinv3 = Tinv3(I);

semilogy(Tinv3, Zv, 'o')
hold on
semilogy(Tinv3, ZvT)


