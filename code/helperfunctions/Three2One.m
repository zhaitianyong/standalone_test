function type = Three2One(input)

amino_3L = {'ALA', 'ARG', 'ASN','ASP','CYS', 'GLU', 'GLN','GLY', 'HIS', 'ILE',...
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL','CYSS'};

amino_1L = {'A','R','N','D','C','E','Q','G','H','I','L','K','M'...
    ,'F','P','S','T','W','Y','V','Z'};

type = -1;
for i = 1:21
    if strcmp(amino_3L{i},input)
        type = amino_1L{i};
        break;
    end
end
