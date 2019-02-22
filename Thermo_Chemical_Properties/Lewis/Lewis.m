function [ Cp,H,G,S ] = Lewis(sp,T,data)
%% evaluate thermodynmic properties for given temperature
%Cp J/K/mol
% H  kJ/mol
% S   J/K/mol
% G   kJ/mol
R = 8.314510; %J/mol/K
T = reshape(T,1,[]);
splist = {data.sp};
m = find(strcmp(splist,sp));
if isempty(m)
    error('No such specie');
end
nT = data(m).nT;
c = data(m).coeff(:,1:7);
Trange =  data(m).Trange;
Hdiff =  data(m).Hdiff;
b = data(m).coeff(:,8:9);
mass = data(m).m0;
%% group by temperature 
index = cell(1,nT);
findex  = [];
for i = 1:nT
    index{i} = find( T>Trange(i,1) & T<=Trange(i,2));
    findex = [findex, index{i} ];
end
noindex = setdiff(1:length(T),findex);
Cp = zeros(1,length(T));
H = Cp; S = Cp; G = Cp; 
%% calculate
for i = 1:nT
    cT = T(index{i});
    if ~isempty(cT)
        cTpow = repmat(cT,6,1);
        cTpow(2,:) = 1./cT;      %cT^(-1)
        cTpow(4,:) = cT.*cT;     %cT^2
        cTpow(5,:) = cTpow(4,:).*cT;  %cT^3
        cTpow(6,:) = cTpow(5,:).*cT;  %cT^4
        cTpow(1,:) = 1./cTpow(4,:);   %cT^-2
        
        
        Cp(index{i}) = R*(c(i,1)*cTpow(1,:) + c(i,2)*cTpow(2,:) + c(i,3) + c(i,4)*cT + ...
            c(i,5)*cTpow(4,:)+c(i,6)*cTpow(5,:)+c(i,7)*cTpow(6,:));
        H(index{i}) = R*cT.*(-c(i,1)*cTpow(1,:) + c(i,2)*log(cT)./cT + c(i,3) + c(i,4)*cT/2 + ...
            c(i,5)*cTpow(4,:)/3+c(i,6)*cTpow(5,:)/4+c(i,7)*cTpow(6,:)/5+b(i,1)./cT);
        S(index{i}) = R.*(-c(i,1)*cTpow(1,:)/2 - c(i,2)./cT + c(i,3)*log(cT) + c(i,4)*cT + ...
            c(i,5)*cTpow(4,:)/2+c(i,6)*cTpow(5,:)/3+c(i,7)*cTpow(6,:)/4+b(i,2));
        G(index{i}) = H(index{i})- cT.*S(index{i});
    end
end
H = H/1000;
G = G/1000;

end
