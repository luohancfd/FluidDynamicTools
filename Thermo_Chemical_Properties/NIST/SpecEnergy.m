function Erv = SpecEnergy(v,J,par)
% calculate rovibrational energy based on spectrosopy measurement data
% par is a vector contains: we wexe weye Be alphae gamme De betae
vph = v+0.5;
jp1 = J*(J+1);
Bv = par(4) - par(5)*vph + par(6)*vph*vph;
Dv = par(7) + par(8)*vph;
Fv = Bv*jp1 - Dv*jp1*jp1;
Gv = par(1)*vph-par(2)*vph*vph+par(3)*vph*vph*vph;


E0 = par(1)*0.5-par(2)*0.5*0.5+par(3)*0.25*0.5;

Erv = Gv+Fv-E0;
end

