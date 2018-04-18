clear
clc
MN2 = 28.016;
RgasN2 = 8314.4598/MN2;
thetaV = 3395;
%% following data from gupta
Omega_C =   [ -6.0614558000E-03  1.2689102000E-01 -1.0616948000E+00  8.0955466000E+02;...
           -7.6303990000E-03  1.6878089000E-01 -1.4004234000E+00  2.1427708000E+03;...
            0.0000000000E+00  0.0000000000E+00  0.0000000000E+00  1.1000000000E+00];
Omega11 = @(T) 1/pi*Omega_C(1,4)*T.^(Omega_C(1,1)*log(T).*log(T)+Omega_C(1,2).*log(T)+Omega_C(1,3));
Omega22 = @(T) 1/pi*Omega_C(2,4)*T.^(Omega_C(2,1)*log(T).*log(T)+Omega_C(2,2).*log(T)+Omega_C(2,3));
Bij     = @(T) Omega_C(3,4)*T.^(Omega_C(3,1)*log(T).*log(T)+Omega_C(3,2).*log(T)+Omega_C(3,3));

mu = @(T) 266.93e-7*sqrt(T*MN2)./Omega22(T)/1000*100;   %kg/m/s
kTtr = @(T) 1989.1e-7*sqrt(T/MN2)./Omega22(T)*4.184*100; %W/(m*K)
%kTtr = @(T) mu(T).*RgasN2*3/2*5/2; %W/(m*K)
kTrot = @(T)6.3605e-5*sqrt(T/MN2)./Omega11(T)*4.184*100;
%https://ntrs.nasa.gov/search.jsp?R=19900017748 2018-03-22T16:16:03+00:00Z
kTvib = @(T,s)6.3605e-5*sqrt(T/MN2)./Omega11(T).*s*4.184*100;

% Eucken's formula kTtr = mu *5/2*Cv_tr, Cv_tr = 3/2*RgasN2

grid_f = 'grid.h5';
flow_f = 'data.h5';
conn_f = 'conn.h5';
info = h5info(flow_f,'/info/grid');
ncell = info.Attributes(2).Value;  %number of actuall cell

data0 = h5read(flow_f,'/solution/run_1/interior');
data0 = data0(:,1:ncell);
svnames = h5read(flow_f,'/info/solver/svnames');
svnames = cellfun(@strtrim,svnames,'UniformOutput',false);
data0 = mat2cell(data0,ones(1,length(svnames)));
data = cell2struct(data0,svnames);
data.rho = data.rhosN2;
data.p = data.T.*data.rho.*RgasN2;
vel = sqrt(data.u.^2+data.v.^2+data.w.^2);
data.m = vel./sqrt(1.4*data.T*RgasN2);
data.mu = mu(data.T);
data.kt = kTtr(data.T) + kTrot(data.T);
w = thetaV./data.Tv;
cpv = w.^2.*exp(w)./(exp(w)-1).^2;
data.kv = kTvib(data.T,cpv);

xcn  =   h5read(grid_f,'/xcn');
conn = h5read(conn_f,'/ien');
conn = conn(2:end,:)';

names = fieldnames(data);
tdata=[];
tdata.varnames={'x','y','z', names{:}}; %#ok<CCAT>
tdata.FEvolumes(1).zonename='US3D_Data';
tdata.FEvolumes(1).x=xcn(1,:);
tdata.FEvolumes(1).y=xcn(2,:);
tdata.FEvolumes(1).z=xcn(3,:);
%temperature
for i = 1:length(names)
    tdata.FEvolumes(1).v(i,:)=data.(names{i});
end
tdata.FEvolumes(1).e2n=conn;
tdata.FEvolumes(1).varloc =1;
mat2tecplot(tdata,'fvs.plt')