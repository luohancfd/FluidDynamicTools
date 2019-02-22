function [prop] = read_ladder( file )

%file = strcat(ladder_dir,'/nn_ladder.txt');
fid = fopen(file,'r');

fprintf([fgetl(fid),'\n']);
sp = strtrim(fgetl(fid));
fprintf([sp,'\n']);

fprintf([fgetl(fid),'\n']);
Vmax = sscanf(fgetl(fid),'%d');
fprintf('%d\n', Vmax);

fprintf([fgetl(fid),'\n']);
Jmax = fscanf(fid,'%d',Vmax+1);
for i = 1:length(Jmax)
    fprintf('%d ', Jmax(i));
end
fprintf('\n');

l = strtrim(fgetl(fid));
while ~contains(l,'WellDepth')
    l = strtrim(fgetl(fid));
end
fprintf([l,'\n']);
De = sscanf(fgetl(fid),'%f');
fprintf('%f\n', De);

fprintf([fgetl(fid),'\n']);
Re = sscanf(fgetl(fid),'%f');
fprintf('%f\n', Re);

fprintf([fgetl(fid),'\n']);
Ev = zeros(Vmax+1,1);
for i = 1:Vmax+1
    Ev(i) = sscanf(fgetl(fid),'%f');
end

fprintf([fgetl(fid),'\n']);
Erv = zeros(Vmax+1,Jmax(1)+1);
for i = 1:Vmax+1
    l = fgetl(fid);
    Erv(i,1:Jmax(i)+1) = sscanf(l,'%f');
end

fprintf([fgetl(fid),'\n']);
Rin = zeros(Vmax+1,Jmax(1)+1);
for i = 1:Vmax+1
    l = fgetl(fid);
    Rin(i,1:Jmax(i)+1) = sscanf(l,'%f');
end

fprintf([fgetl(fid),'\n']);
Rout = zeros(Vmax+1,Jmax(1)+1);
for i = 1:Vmax+1
    l = fgetl(fid);
    Rout(i,1:Jmax(i)+1) = sscanf(l,'%f');
end

fclose(fid);


prop.Vmax = Vmax;
prop.Jmax = Jmax;
prop.Ev = Ev;
prop.Erv = Erv;
prop.De = De;
prop.Rin = Rin;
prop.Rout = Rout;
prop.Re = Re;
