% This function writes the output into a
% PDB file
%
% Babak Alipanahi
% University of Waterloo
% July 26, 2010
% -edit: Feb 15, 2011
% -edit: Aug 14, 2011
%
% Some ideas taken from pdbwrite.m (MathWorks)

function pdb_writer(file_name,X,Comp)

amino_3L = {'ALA', 'ARG', 'ASN','ASP','CYS', 'GLU','GLN', 'GLY', 'HIS', 'ILE',...
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL','CYS'};

atom_types = {'H', 'C', 'N', 'O', 'S', 'P','M'};


fid = fopen(file_name,'wt');
if fid < 0
    error('there was an error opening the file...')
end

cnt = 0;
max_i = size(X,2);

residue = Comp.residue;
info = Comp.info;
atom_names = Comp.atom_names;
residue_type = Comp.residue_type;

uniq_residues = unique(residue);
max_res = max(uniq_residues);
num_residues = numel(uniq_residues);
num_seq_res_lines = ceil(num_residues/13);




%chain_id = 'A';
chain_id = ' ';

try   
%     res_cnt = 1;
%     for i = 1:num_seq_res_lines
%        fprintf(fid,'SEQRES%4d %s%5d ',i,chain_id,max_res); 
%        for j = 1:13
%           if res_cnt <= num_residues
%               fprintf(fid,' %s',amino_3L{Comp.seq(res_cnt)});
%               res_cnt = res_cnt + 1;
%           end
%        end
%        fprintf(fid,'\n');
%     end
    
    for i = 1:max_i
        if info(5,i) > 0
            shift_left = 0;
            cnt = cnt + 1;            
            filler = '';
            tname = atom_names{i};
            
            if residue(i) == max_res && strcmp(tname,'O')
               tname = 'OXT'; 
            end
            
            if tname(1) == 'H' && isstrprop(tname(end),'digit') && length(tname) > 3
                shift_left = 1;
                left_space = '';
            end
            
            for j = length(tname)+1:4
                filler(end+1) = ' '; %#ok<AGROW>
            end
            if shift_left
                filler(end+1) = ' '; %#ok<AGROW>
            else
                left_space = ' ';
            end
            
            
            fprintf(fid,'ATOM  %5d%s %s%s%3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s\n',...
                cnt,left_space,tname,filler,amino_3L{residue_type(i)},chain_id,residue(i),X(1,i),X(2,i),X(3,i),atom_types{info(1,i)});
            TER_data.Res_name = amino_3L{residue_type(i)};
            TER_data.chain_id = chain_id;
            TER_data.Res_num  = residue(i);
            TER_data.atom_num = cnt + 1;
        end
    end         
    fprintf(fid,'TER   %5d      %3s %s%4d\n',TER_data.atom_num,TER_data.Res_name,TER_data.chain_id,TER_data.Res_num);
    fprintf(fid,'END');
    fclose(fid);
catch except
    fclose(fid);
    delete(file_name);
    rethrow(except);
end

