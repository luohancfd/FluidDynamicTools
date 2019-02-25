function [ Cp,H,G,S ] = Lewis7ParFun(T,c,Trange)
if nargin == 2
    Trange = [1000,5000; 300, 1000];
end
R = 8.314510; %J/mol/K

Cp = zeros(1,length(T));
H = Cp; S = Cp; G = Cp; 
nT = 2;
index = cell(1,nT);
findex  = [];
for i = 1:nT
    index{i} = find( T>Trange(i,1) & T<=Trange(i,2));
    findex = [findex, index{i} ];
end
noindex = setdiff(1:length(T),findex);
for i = 1:nT
    cT = T(index{i});
    if ~isempty(cT)
        cTpow = repmat(cT,6,1);
        cTpow(1,:) = 1./cT;      %cT^(-1)
        cTpow(3,:) = cT.*cT;     %cT^2
        cTpow(4,:) = cTpow(3,:).*cT;  %cT^3
        cTpow(5,:) = cTpow(4,:).*cT;  %cT^4
        % cTpow(2,:) % cT
        
        
        Cp(index{i}) = R*(c(i,1) + c(i,2)*cTpow(2,:) + c(i,3)*cTpow(3,:) + c(i,4)*cTpow(4,:) + ...
            c(i,5)*cTpow(5,:));
        H(index{i}) = R*cT.*(c(i,1) + c(i,2)*cTpow(2,:)/2 + c(i,3)*cTpow(3,:)/3 + c(i,4)*cTpow(4,:)/4 + ...
            c(i,5)*cTpow(5,:)/5 + c(i,6)*cTpow(1,:));
        S(index{i}) = R.*(c(i,1)*log(cTpow(2,:)) + c(i,2)*cTpow(2,:) + c(i,3)*cTpow(3,:)/2 + c(i,4)*cTpow(4,:)/3 + ...
            c(i,5)*cTpow(5,:)/4 + c(i,7));
        G(index{i}) = H(index{i})- cT.*S(index{i});
    end
end
H = H/1000;
G = G/1000;
end

