% This function reads the torsion angle restraint file in
% CYANA .aco format
%
% Babak Alipanahi
% University of Waterloo
% June 21, 2010

function [phi psi] = ang_reader(file_name,num)

if nargin ~= 2
    error('only two input arguments')
end

len_seq = length(num);

phi = ones(len_seq,2);
psi = ones(len_seq,2);

phi(:,1) = phi(:,1)*-180;
phi(:,2) = phi(:,2)* 180;

psi(:,1) = psi(:,1)*-180;
psi(:,2) = psi(:,2)* 180;

fid = fopen(file_name);
if fid < 0
    disp('warning! aco file could not be opend, returning [-180 180] values.')
else
    
%     raw_data = textscan(fid,'%d %s %s %f %f');
%     fclose(fid);
%     
    raw_data = textscan(fid,'%s','Delimiter','\r\n');
    raw_data = raw_data{1};
    fclose(fid);
    
    % removing empty lines
    max_i = numel(raw_data);
    index_empty = false(1,max_i);
    for i = 1:max_i
        if isempty(raw_data{i}) || ~isempty(regexp(raw_data{i},'\#', 'once'))  
            index_empty(i) = true;
        end
    end
    raw_data(index_empty) = [];
    max_i = numel(raw_data);
    
    % removing the residues with ambigous
    % torsion angle restraints
    bad_residues = nan(2,len_seq);
    bad_residues(1,:) = num;
    num_fields = nan(max_i,1);
    for i = 1:max_i
        temp_str = textscan(raw_data{i},'%s');
        temp_str = temp_str{1};
        num_fields(i) = numel(temp_str);
        if strcmp(temp_str{num_fields(i)},'OR')
           temp_res = str2double(temp_str{1});
           bad_residues(2,bad_residues(1,:) == temp_res) = 1; 
        end
    end
    
    
       
    for i = 1:max_i
        switch num_fields(i)
            case 5
                temp_data = textscan(raw_data{i},'%d %s %s %f %f');
            case 6
                temp_data = textscan(raw_data{i},'%d %s %s %f %f %s');
            case 7
                temp_data = textscan(raw_data{i},'%d %s %s %f %f %s %s');
            case 8
                temp_data = textscan(raw_data{i},'%d %s %s %f %f %s %s %s');

        end
        temp_res = temp_data{1};
        temp_res_name = char(temp_data{2});
        temp_ang_name = char(temp_data{3});
        temp_ang_low  = temp_data{4};
        temp_ang_hi   = temp_data{5};
        index_num = find(num == temp_res);
        if ~isempty(index_num) && Res2Num(Three2One(temp_res_name))
            if isnan(bad_residues(2,index_num))
                switch temp_ang_name
                    case 'PHI'
                        phi(index_num,1) = temp_ang_low;
                        phi(index_num,2) = temp_ang_hi;
                    case 'PSI'
                        psi(index_num,1) = temp_ang_low;
                        psi(index_num,2) = temp_ang_hi;
                end
            end
        end
    end
end


