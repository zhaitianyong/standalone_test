function type = Res2Num(input)

amino_1L = {'A','R','N','D','C','E','Q','G','H','I','L',...
    'K','M','F','P','S','T','W','Y','V','Z'};

if size(input,2) == 1
    type = -1;
    for i = 1:21
        if strcmp(amino_1L{i},input)
            type = i;
            break;
        end
    end
elseif size(input,2) > 1
    max_j = size(input,2);
    type = -1*ones(max_j,1);
    for j = 1:max_j
        for i = 1:21
            if strcmp(amino_1L{i},input(j))
                type(j) = i;
                break;
            end
        end
    end
end