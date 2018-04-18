function data = ReadNIST(filename)
data = repmat(struct('Particle',[],'data',zeros(100,10)),[100,1]); %buffer
fid = fopen(filename);
n_species = 0;
n_elec = 0;
cline = fgetl(fid);
while ischar(cline)
    cline = strtrim(cline);
    if ~isempty(cline)
        if ~strcmp(cline(1),'#')
            if strcmp(cline(1:8),'Molecule')
                if n_species > 0
                    data(n_species).data = data(n_species).data(1:n_elec,:);
                end
                n_species = n_species+1;
                data(n_species).Particle = strtrim(cline(10:end));
                n_elec = 0;
            elseif strcmp(cline(1:4),'Atom')
                if n_species > 0
                    data(n_species).data = data(n_species).data(1:n_elec,:);
                end
                n_species = n_species+1;
                data(n_species).Particle = strtrim(cline(6:end));
                n_elec = 0;
            else
                n_elec = n_elec + 1;
                A=sscanf(cline,'%d%s%f%f%f%f%f%f%f%f%f%f');
                data(n_species).data(n_elec,2:end)=A(end-8:end);
                
                TermSymbol = [char(A(2:end-10))]';
                
                Lam = -1; S = (A(1)-1)/2;
                if strcmp(TermSymbol,'Sigma')
                    Lam = 0;
                    data(n_species).data(n_elec,1)=1*A(1);
                elseif strcmp(TermSymbol,'Pi')
                    Lam = 1;
                    data(n_species).data(n_elec,1)=2*A(1);
                elseif strcmp(TermSymbol,'Delta')
                    Lam = 2;
                    data(n_species).data(n_elec,1)=2*A(1);
                else
                    data(n_species).data(n_elec,1)=A(end-9); % this is atom
                end
                
                if Lam ~= -1  %for molecule 
                    Jmin = abs(Lam -S); Jmax = Lam+S;
                    if Jmin ~= Jmax
                        nJ = Jmax - Jmin + 1;
                        data(n_species).data(n_elec,1) = data(n_species).data(n_elec,1)/nJ;
                    end
                end
            end
        end
    end
    cline = fgetl(fid);
end
data(n_species).data = data(n_species).data(1:n_elec,:);
data = data(1:n_species);
fclose(fid);

end

