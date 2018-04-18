function [b,c] = plt_check_section( fid,a )
% check if section id is the one
% if yes, jump to the bit after the value and return true
% else,   jump back the bit before that and return false
floc = ftell(fid);
c = fread(fid,1,'float32');
if length(a) == 1
    if c == a
        b = true;
    else
        b = false;
        fseek(fid,floc,'bof');
    end
else
    b = false(size(a));
    for i = 1:length(a)
        if c == a(i)
            b(i) = true;
            break
        end
    end
    if sum(b) == 0
        fseek(fid,floc,'bof');
    end
end

