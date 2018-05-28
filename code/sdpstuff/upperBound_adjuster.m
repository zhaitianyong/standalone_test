% This function adjusts distance constraints
%
% Babak Alipanahi
% University of Waterloo
% February 23, 2011

function oC = upperBound_adjuster(C,A,hydrogen_omission)

oC = C;

num = size(C,2);
for i = 1:num
    satom  = C(i).satom;
    tatom  = C(i).tatom;
    stype  = C(i).stype;
    ttype  = C(i).ttype;
    
    sswaps = A(stype).swaps;
    tswaps = A(ttype).swaps;
    max_s  = size(sswaps,2);
    max_t  = size(tswaps,2);
    
    % When length(sswaps) = 2, it means that with or
    % without hydrogen omission, this mapping should
    % be done    
    if hydrogen_omission
        max_s = size(sswaps,2);
        max_t = size(tswaps,2);
        for s = 1:max_s
            if strcmp(sswaps{s}{1},satom)
                if length(sswaps{s}) > 2
                    oC(i).satom = sswaps{s}{3};
                else
                    oC(i).satom = sswaps{s}{2};
                end
                break
            end
        end
        for t = 1:max_t
            if strcmp(tswaps{t}{1},tatom)
                if length(tswaps{t}) > 2
                    oC(i).tatom = tswaps{t}{3};
                else
                    oC(i).tatom = tswaps{t}{2};
                end
            end
        end
        % This is for the case that for example
        % HB1 has already been mapped to QB but
        % we want to remap that to CB
        for s = 1:max_s
            if strcmp(sswaps{s}{2},satom)
                if length(sswaps{s}) > 2
                    oC(i).satom = sswaps{s}{3};
                end
                break
            end
        end
        for t = 1:max_t
            if strcmp(tswaps{t}{2},tatom)
                if length(tswaps{t}) > 2
                    oC(i).tatom = tswaps{t}{3};                
                end
            end
        end
    else
        for s = 1:max_s
            if strcmp(sswaps{s}{1},satom)
                oC(i).satom = sswaps{s}{2};
                break
            end
        end
        for t = 1:max_t
            if strcmp(tswaps{t}{1},tatom)
                oC(i).tatom = tswaps{t}{2};
            end
        end
    end
end
