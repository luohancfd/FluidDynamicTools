function string= plt_read_string(fid,m,type)
%read a string
if nargin == 1
    string = fread(fid,1,'int32');
    while string(end) ~= 0
        string = [string, fread(fid,1,'int32')];
    end
    string = string(1:end-1);
else
    %type gives the size of one char
    string = fread(fid,m,type);
    string = char(string);
end
string = char(string);
end

