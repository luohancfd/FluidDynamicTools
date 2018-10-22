%% read capitelli data
%http://arc.aiaa.org/doi/10.2514/2.6517
clear
clc
atommass =struct('N',14.007,'O',15.999,'H',1.008,'C',12.011,'Ar',39.948);
% atommass =struct('N',14,'O',16,'H',1,'C',12,'Ar',40);

fid = fopen('DATA_Capitelli.txt','r');
data = repmat(struct('sp1',[],'sp2',[],'m1',[],'m2',[],'a',[]),100,1);
for i = 1:100
    if ~feof(fid)
        line = fgetl(fid);
        [temp,~]=regexp(line,'([A-Z]{1,2}\d?[+-]?)-([A-Z]{1,2}\d?[+-]?)\s*(.*)','tokens','once');
        sp1 = temp{1};
        sp2 = temp{2};
        % analysze mass
        [msp1,~]=regexp(sp1,'([A-Z]\d?)','tokens','once');
        ms1 = 0;
        for ii = 1:length(msp1)
            if length(msp1{ii}) == 2
                ms1 = atommass.(msp1{ii}(1))*str2double(msp1{ii}(2));
            else
                ms1 = atommass.(msp1{ii}(1));
            end
        end
        [msp2,~]=regexp(sp2,'([A-Z]\d?)','tokens','once');
        ms2 = 0;
        for ii = 1:length(msp1)
            if length(msp2{ii}) == 2
                ms2 = atommass.(msp2{ii}(1))*str2double(msp2{ii}(2));
            else
                ms2 = atommass.(msp2{ii}(1));
            end
        end            
        data(i).sp1 = sp1; data(i).sp2 = sp2;
        data(i).m1 = ms1; data(i).m2 = ms2;
        data(i).a = zeros(4,6);
        data(i).a(1,:) = sscanf(temp{3},'%e');
        for j = 2:4
            line = fgetl(fid);
            data(i).a(j,:)= sscanf(line,'%e');
        end
        fprintf('find %s - %s\n',temp{1},temp{2});
    else
        break
    end
end
data = data(1:i-1);
fclose(fid);
save('DATA_capitelli.mat','data')
% 
% fid = fopen('update.txt','w');
% for i = 1:length(data)
%     fprintf(fid,'%-8s',[data(i).sp1,'-',data(i).sp2]);
%     for j = 1:6
%         fprintf(fid,'%15.4E',data(i).a(1,j));
%     end
%     fprintf(fid,'\n');
%     for jj = 2:4
%         fprintf(fid,'        ');
%         for j = 1:6
%             fprintf(fid,'%15.4E',data(i).a(jj,j));
%         end
%         fprintf(fid,'\n');
%     end
% end
% fclose(fid);
