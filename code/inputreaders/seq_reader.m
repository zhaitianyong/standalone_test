% This function reads the input sequence file in
% CYANA .seq format
% 
% Babak Alipanahi
% University of Waterloo
% June 21, 2010

function [seq num] = seq_reader(file_name)

if nargin ~= 1
    error('only one input argument')
end
    
fid = fopen(file_name);
if fid < 0
    error('seq file could not be opend')
end

raw_data = textscan(fid,'%s');
fclose(fid);

first_residue = 0;
residue_cnt = 0;
max_i = size(raw_data{1},1);
seq = nan(1,max_i);
num = nan(1,max_i);
i = 1;
while i <= max_i
    temp_name = upper(raw_data{1}{i});    
    if isnan(str2double(temp_name))
       residue_cnt = residue_cnt + 1;       
       seq(residue_cnt) = Res2Num(Three2One(temp_name));
       num(residue_cnt) = -1;
       if seq(residue_cnt) < 0
           error('unkown residue name');
       end
    elseif first_residue == 0
        first_residue = str2double(temp_name);
        num(residue_cnt) = first_residue;
    else
        num(residue_cnt) = str2double(temp_name);
    end
    i = i + 1;
end
seq(isnan(seq)) = [];
num(isnan(num)) = [];
index_positive = find(num > 0);
index_positive(end+1) = residue_cnt + 1;

for i = 1:length(index_positive)-1
    len = index_positive(i+1) - index_positive(i);
    num(index_positive(i):index_positive(i+1)-1) = num(index_positive(i)) + (0:len-1);
end    