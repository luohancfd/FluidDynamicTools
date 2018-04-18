function [zonetype,itype] = plt_read_zonetype(fid)
itype = fread(fid,1,'int32');
switch itype
    case     0
        zonetype = "ORDERED";
    case     1
        zonetype = "FELINESEG";
    case     2
        zonetype = "FETRIANGLE";
    case     3
        zonetype = "FEQUADRILATERAL";
    case     4
        zonetype = "FETETRAHEDRON";
    case     5
        zonetype = "FEBRICK";
    case     6
        zonetype = "FEPOLYGON";
    case     7
        zonetype = "FEPOLYHEDRON";
end
end