function Compile_mex_files()

clc

cd D:/standalone_test/code/refinement
source_files = {'bondcheckmex.c','distmex.c','gradfunmex.c',...
    'loboundmex.c','objfunmex.c','upboundmex.c'};


disp('Compiling C-Mex files')
disp('')

COMPUTER = computer;
IS64BIT  = strcmp(COMPUTER(end-1:end),'64');


for i= 1:numel(source_files)
    if IS64BIT
        str =  ['mex -largeArrayDims -Dchar16_t=UINT16_T ' source_files{i}];
    else
        str =  ['mex -Dchar16_t=UINT16_T ' source_files{i}];
    end    
    
    disp(str);
    eval(str);
end
cd ..
cd D:/standalone_test/matrix_completion/Auxiliary

% compute a subset entries of matrix product
mex partXY.c

% update a sparse matrix 
mex updateSval.c

cd ..
cd ..
disp('Done!')


