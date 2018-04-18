function [ x,y,v,z ] = fem2ij( conn,datav,nvar,icell0,datax,datay,dataz,Imax,Jmax,direct)
%% Author: Han Luo
% Github repository: https://github.com/luohancfd/FluidDynamicTools/tree/master/Tecplot_Tools
%% convert FEM grid to IJ ordered
x = zeros(Jmax,Imax);
y = zeros(Jmax,Imax);
z = zeros(Jmax,Imax);
v = zeros(Jmax-1,Imax-1,nvar);
%% change order of conn to make sure always like the following
            %  ^ jp(Imax per +1)
            %  |          RHS order
            %  |--->ip
            %  p3-----p2
            %   |      |
            %   |   1  |
            %  p4-----p1
p1 = conn(icell0,:);
p2 = conn(icell0+1,:);
ix = [datax(p1(2))-datax(p1(1)),datay(p1(2))-datay(p1(1)),dataz(p1(2))-dataz(p1(1))];
iy = [datax(p1(4))-datax(p1(1)),datay(p1(4))-datay(p1(1)),dataz(p1(4))-dataz(p1(1))];
iz = cross(ix,iy);
if iz(3) < 0
    conn = fliplr(conn);
    p1 = conn(icell0,:);
    p2 = conn(icell0+1,:);
end
[~,ia,~]=intersect(p1,p2);
if ia(1) ~= 1
    neworder = ones(1,4);
    for i=1:4
        if i>= ia(1)
            neworder(i) = 1+i-ia(1);
        else
            neworder(i) = 5-ia(1)+i;
        end
    end
    ib = neworder;
    for i = 1:4
        ib(neworder(i)) = i;
    end
    conn = conn(:,ib);
end

% icell0: index of first cell
if direct == 1
    icell_matrix = zeros(Jmax-1,Imax-1);
    for j = 1:Jmax-1
        for i = 1:Imax-1
            ip = i;
            jp = j;
            i_cell = icell0-1+ip+(jp-1)*(Imax-1);
            icell_matrix(j,i) = i_cell;
            %  ^ jp(Imax per +1)
            %  |          RHS order
            %  |--->ip
            %  p3-----p2
            %   |      |
            %   |   1  |
            %  p4-----p1
            p1 = conn(i_cell,1);
            p2 = conn(i_cell,2);
            p3 = conn(i_cell,3);
            p4 = conn(i_cell,4);
            x(j,i) = datax(p4); y(j,i) = datay(p4); z(j,i) = dataz(p4);
            v(j,i,:)  = datav(:,i_cell);
            if i == Imax-1
                x(j,i+1) = datax(p1); y(j,i+1) = datay(p1); z(j,i+1) = dataz(p1);
            end
            if j == Jmax-1
                x(j+1,i) = datax(p3); y(j+1,i) = datay(p3); z(j+1,i) = dataz(p3);
                if i == Imax-1
                    x(j+1,i+1) = datax(p2); y(j+1,i+1) = datay(p2);
                    z(j+1,i+1) = dataz(p2);
                end
            end
        end
    end
elseif direct == 2
    icell_matrix = zeros(Jmax-1,Imax-1);
    for jp = 1:Imax-1
        for ip = 1:Jmax-1
            j = ip;
            i = Imax-jp;
            i_cell = icell0-1+ip+(jp-1)*(Jmax-1);
            icell_matrix(j,i) = i_cell;
            %                     ^ ip
            %                      |
            % jp (Jmax per +1) <---
            %  p2-----p1
            %   |      |
            %   |   1  |
            %  p3-----p4
            p1 = conn(i_cell,1);
            p2 = conn(i_cell,2);
            p3 = conn(i_cell,3);
            p4 = conn(i_cell,4);
            x(j,i) = datax(p3); y(j,i) = datay(p3); z(j,i) = dataz(p3);
            v(j,i,:)  = datav(:,i_cell);
            if i == Imax-1
                x(j,i+1) = datax(p4); y(j,i+1) = datay(p4); z(j,i+1) = dataz(p4);
            end
            if j == Jmax-1
                x(j+1,i) = datax(p2); y(j+1,i) = datay(p2); z(j+1,i) = dataz(p2);
                if i == Imax-1
                    x(j+1,i+1) = datax(p1); y(j+1,i+1) = datay(p1); z(j+1,i+1) = dataz(p1);
                end
            end
        end
    end
end
icell_matrix = flipud(icell_matrix);  %debug use
end

