% This function reads the input distance restraint file in
% CYANA .upl and .lol format
%
% Babak Alipanahi
% University of Waterloo
% June 21, 2010
% -edit: March 22, 2011

% generates constraints for full-H and omitted-H cases
% type= 0: full atom, = 1; H-omitted


function dist_con = dist_reader(file_name,A,type)

if nargin < 2
    error('Two inputs argument are required')
end

if nargin < 3
    type = 0;
end

fid = fopen(file_name);
if fid < 0
    error('dist file could not be opend')
end

%  format:
%  res_num res_name atom_name   res_num res_name atom_name distance #peak
%  peak_num

raw_data = textscan(fid,'%s','Delimiter','\r\n');
raw_data = raw_data{1};
fclose(fid);

% removing empty lines
XPLOR_format = 0;
max_i = numel(raw_data);
index_empty = false(1,max_i);
for i = 1:max_i
    if isempty(raw_data{i})
        index_empty(i) = true;
    else
        if ~isempty(regexp(raw_data{i},'HB1[^1-9]', 'once')) || ~isempty(regexp(raw_data{i},'HA1[^1-9]', 'once'))
            XPLOR_format = 1;
        end
    end
end
raw_data(index_empty) = [];
max_i = numel(raw_data);

% info = info_reader();

dist_con(max_i) = struct('sres',[],'stype',[],'satom','','tres',[],'ttype',[],'tatom','','dist',[],'peak',[]);

index_bad = false(max_i,1);

basic_format = '%d %s %s %d %s %s %f';

for i = 1:max_i
    temp_line = raw_data{i};
    dummy = textscan(temp_line,'%s');
    dummy = dummy{1};
    num_fields = numel(dummy);
    if num_fields == 7
        temp_line =  textscan(temp_line,basic_format);
    elseif num_fields == 8
        temp_line =  textscan(temp_line,[basic_format ' %f']);
    elseif num_fields == 9
        temp_line =  textscan(temp_line,[basic_format ' %s %f']);
    elseif num_fields > 9
        temp_format = [basic_format ' %s %f'];
        for j = num_fields - 9:num_fields
           temp_format = [temp_format ' %s'];  %#ok<*AGROW>
        end
    end
    %*********************************************************************
    % First handle the S part
    dist_con(i).sres  = temp_line{1};
    stype = char(temp_line{2});
    dist_con(i).stype = Res2Num(Three2One(stype));
    dist_con(i).satom = char(temp_line{3});
    dist_con(i).peak  = 0;
    
    delta_s = 0;
    if strcmp(dist_con(i).satom,'HN')
        dist_con(i).satom = 'H';
    end
    if XPLOR_format
        temp = atom_nom(stype,dist_con(i).satom,'xplor');
        if ~isempty(temp)
            dist_con(i).satom = temp;
        end
    end
    old_atom = dist_con(i).satom;
    if type > 0
        [temp flag] = atom_nom(stype,dist_con(i).satom,'homitted');
        if ~isempty(temp)
            dist_con(i).satom = temp;
        end
    else
        [temp flag] = atom_nom(stype, dist_con(i).satom,'full');
        if ~isempty(temp)
            dist_con(i).satom = temp;
        end
    end
    if flag > 0
        inds = index_finder({A(dist_con(i).stype).atom.name},{dist_con(i).satom,old_atom});
        delta_s = norm(A(dist_con(i).stype).X(:,inds(1)) - A(dist_con(i).stype).X(:,inds(2)));
    end
    %*********************************************************************
    % Then handle the T part
    dist_con(i).tres  = temp_line{4};
    ttype = char(temp_line{5});
    dist_con(i).ttype = Res2Num(Three2One(ttype));
    dist_con(i).tatom = char(temp_line{6});
        
    delta_t = 0;    
    if strcmp(dist_con(i).tatom,'HN')
        dist_con(i).tatom = 'H';
    end
    if XPLOR_format
        temp = atom_nom(ttype,dist_con(i).tatom,'xplor');
        if ~isempty(temp)
            dist_con(i).tatom = temp;
        end
    end
    
    old_atom = dist_con(i).tatom;
    if type > 0
        [temp flag] = atom_nom(ttype,dist_con(i).tatom,'homitted');
        if ~isempty(temp)
            dist_con(i).tatom = temp;
        end
    else
        [temp flag] = atom_nom(ttype,dist_con(i).tatom,'full');
        if ~isempty(temp)
            dist_con(i).tatom = temp;
        end
    end
    if flag > 0
        inds = index_finder({A(dist_con(i).ttype).atom.name},{dist_con(i).tatom,old_atom});
        delta_t = norm(A(dist_con(i).ttype).X(:,inds(1)) - A(dist_con(i).ttype).X(:,inds(2)));
    end         
    %*********************************************************************
    if temp_line{7} == 0
        index_bad(i) = 1;
    end
    dist_con(i).dist = temp_line{7} + delta_s + delta_t;
    
    if numel(dummy) >= 9
        if ~strcmp(temp_line{8},'#peak') 
            error('dist file format is wrong.')
        end
        dist_con(i).peak = temp_line{9};
    end
          
    if ~isempty(regexp(dist_con(i).satom,'O', 'once')) || ~isempty(regexp(dist_con(i).tatom,'O', 'once')) ...
            || ~isempty(regexp(dist_con(i).satom,'S', 'once')) || ~isempty(regexp(dist_con(i).tatom,'S', 'once'))
        dist_con(i).peak  = -1;  
    end
end
dist_con(index_bad) = [];
max_new = numel(dist_con);
if max_new ~= max_i
    fprintf('%d zero-dist constraints omitted.\n',max_i - max_new);
end
%**********************************************************************************
%**********************************************************************************
function ind = index_finder(list,query)

q = numel(query);
ind = -1*ones(1,q);
for j = 1:q
    for i = 1:numel(list)
        if strcmp(list{i},query{j})
            ind(j) = i;
            break
        end
    end
end
%**********************************************************************************

