function [ E ] = SpecTableEnergy( Spec,diss )
% Spec in cm-1, E in eV
[n_elec,~ ] = size(Spec);
E = repmat(struct('Te',[],'ns',[],'Erv',[],'Jmax',[]),[n_elec,1]);
conver = [8065.479 349.755 83.5935 0.695028792]; % 1 eV | kcal/mol | kJ/mol | K |to cm-1

for i=1:n_elec
    par = Spec(i,3:end);
    % find max V, max J
    MaxV = 0;
    Ev = SpecEnergy(MaxV,0,par);
    prev = -0.1;
    while Ev < diss && Ev > prev
        prev = Ev;
        MaxV = MaxV +1;
        Ev = SpecEnergy(MaxV,0,par);      
    end
    MaxV = MaxV-1;
    
    MaxJ = 0;
    Ej = SpecEnergy(0,MaxJ,par);
    prev = -0.1;
    while Ej < diss && Ej > prev
        prev = Ej;
        MaxJ = MaxJ +1;
        Ej = SpecEnergy(0,MaxJ,par);       
    end
    MaxJ = MaxJ-1;
    
    Erv = zeros(MaxV+1,MaxJ+1);
    Jmax = ones(MaxV+1,1)*MaxJ;
    for v = 0:MaxV
        prev = 0;
        for j = 0:MaxJ
            Erv(v+1,j+1) = SpecEnergy(v,j,par);
            if Erv(v+1,j+1) > diss || Erv(v+1,j+1) < prev
                Erv(v+1,j+1) = 0;
                Jmax(v+1) = j-1;
                break;
            end
            prev =  Erv(v+1,j+1);
        end
    end
    E(i).Jmax = Jmax;
    E(i).Erv = Erv/conver(1);
    E(i).Te = Spec(i,2)/conver(4);
    E(i).ne = Spec(i,1);
    E(i).Vmax = MaxV;
end
end

