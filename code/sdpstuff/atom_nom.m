% mode: full, homitted, xplor
function [natom delta] = atom_nom(residuetype,atom,mode)

delta = -1;
natom = '';
switch mode
    %**********************************************************************
    case 'full'        
        switch residuetype
            case 'ALA'
                switch atom
                    case 'HB1',  natom = 'QB'; delta = 1;
                    case 'HB2',  natom = 'QB'; delta = 1;
                    case 'HB3',  natom = 'QB'; delta = 1;
                    case 'QB'
                end
            case 'ARG'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                    case 'HD2'
                    case 'HD3'
                    case 'QD'
                    case 'HH11'
                    case 'HH12'
                    case 'QH1'
                    case 'HH21'
                    case 'HH22'
                    case 'QH2'
                end
            case 'ASN'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HD21'
                    case 'HD22'
                    case 'QD2'
                end
            case 'ASP'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                end
            case 'CYS'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'HG'
                    case 'QB'
                end
            case 'CYSS'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                end
            case 'GLN'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                    case 'HE21'
                    case 'HE22'
                    case 'QE2'
                end
            case 'GLU'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                end
            case 'GLY'
                switch atom
                    case 'HA2'
                    case 'HA3'
                    case 'QA'
                end
            case 'HIS'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HD1'
                    case 'HD2'
                    case 'HE1'
                end
            case 'ILE'
                switch atom
                    case 'HB'
                    case 'HG21', natom = 'QG2'; delta = 1;
                    case 'HG22', natom = 'QG2'; delta = 1;
                    case 'HG23', natom = 'QG2'; delta = 1;
                    case 'QG2'
                    case 'HG12'
                    case 'HG13'
                    case 'QG1'
                    case 'HD11', natom = 'QD1'; delta = 1;
                    case 'HD12', natom = 'QD1'; delta = 1;
                    case 'HD13', natom = 'QD1'; delta = 1;
                    case 'QD1'
                end
            case 'LEU'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG'
                    case 'HD11', natom = 'QD1'; delta = 1;
                    case 'HD12', natom = 'QD1'; delta = 1;
                    case 'HD13', natom = 'QD1'; delta = 1;
                    case 'QD1'
                    case 'HD21', natom = 'QD2'; delta = 1;
                    case 'HD22', natom = 'QD2'; delta = 1;
                    case 'HD23', natom = 'QD2'; delta = 1;
                    case 'QD2'
                    case 'QQD'
                end
            case 'LYS'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                    case 'HD2'
                    case 'HD3'
                    case 'QD'
                    case 'HE2'
                    case 'HE3'
                    case 'QE'
                    case 'HZ1',  natom = 'QZ'; delta = 1;
                    case 'HZ2',  natom = 'QZ'; delta = 1;
                    case 'HZ3',  natom = 'QZ'; delta = 1;
                    case 'QZ'
                end
            case 'MET'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                    case 'HE1',  natom = 'QE'; delta = 1;
                    case 'HE2',  natom = 'QE'; delta = 1;
                    case 'HE3',  natom = 'QE'; delta = 1;
                    case 'QE'
                end
            case 'PHE'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HD1'
                    case 'HD2'
                    case 'QD'
                    case 'HE1'
                    case 'HE2'
                    case 'QE'
                    case 'QR'
                    case 'HZ'
                end
            case 'PRO'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                    case 'HD2'
                    case 'HD3'
                    case 'QD'
                end
            case 'SER'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG'
                end
            case 'THR'
                switch atom
                    case 'HB'
                    case 'HG21', natom = 'QG2'; delta = 1;
                    case 'HG22', natom = 'QG2'; delta = 1;
                    case 'HG23', natom = 'QG2'; delta = 1;
                    case 'QG2'
                    case 'HG1'
                end
            case 'TRP'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HD1'
                    case 'HE1'
                    case 'HE3'
                    case 'HZ2'
                    case 'HZ3'
                    case 'HH2'
                end
            case 'TYR'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HD1'
                    case 'HD2'
                    case 'QD'
                    case 'HE1'
                    case 'HE2'
                    case 'QE'
                    case 'QR'
                    case 'HH'
                end
            case 'VAL'
                switch atom
                    case 'HB'
                    case 'HG11', natom = 'QG1'; delta = 1;
                    case 'HG12', natom = 'QG1'; delta = 1;
                    case 'HG13', natom = 'QG1'; delta = 1;
                    case 'QG1'
                    case 'HG21', natom = 'QG2'; delta = 1;
                    case 'HG22', natom = 'QG2'; delta = 1;
                    case 'HG23', natom = 'QG2'; delta = 1;
                    case 'QG2'
                    case 'QQG'
                end
        end
        %**********************************************************************
        case 'homitted'        
        switch residuetype
            case 'ALA'
                switch atom
                    case 'HB1',  natom = 'CB'; delta = 1;
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB'
                end
            case 'ARG'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG2',  natom = 'CG'; delta = 1;
                    case 'HG3',  natom = 'CG'; delta = 1;
                    case 'QG',   natom = 'CG'; delta = 1;
                    case 'HD2',  natom = 'CD'; delta = 1;
                    case 'HD3',  natom = 'CD'; delta = 1;
                    case 'QD',   natom = 'CD'; delta = 1;
                    case 'HH11'
                    case 'HH12'
                    case 'QH1'
                    case 'HH21'
                    case 'HH22'
                    case 'QH2'
                end
            case 'ASN'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HD21'
                    case 'HD22'
                    case 'QD2'
                end
            case 'ASP'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                end
            case 'CYS'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'HG',   natom = 'SG'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                end
            case 'CYSS'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                end
            case 'GLN'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG2',  natom = 'CG'; delta = 1;
                    case 'HG3',  natom = 'CG'; delta = 1;
                    case 'QG',   natom = 'CG'; delta = 1;
                    case 'HE21'
                    case 'HE22'
                    case 'QE2'
                end
            case 'GLU'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG2',  natom = 'CG'; delta = 1;
                    case 'HG3',  natom = 'CG'; delta = 1;
                    case 'QG',   natom = 'CG'; delta = 1;
                end
            case 'GLY'
                switch atom
                    case 'HA2'
                    case 'HA3'
                    case 'QA'
                end
            case 'HIS'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HD1'
                    case 'HD2'
                    case 'HE1'
                end
            case 'ILE'
                switch atom
                    case 'HB'
                    case 'HG21', natom = 'QG2'; delta = 1;
                    case 'HG22', natom = 'QG2'; delta = 1;
                    case 'HG23', natom = 'QG2'; delta = 1;
                    case 'QG2'
                    case 'HG12', natom = 'CG1'; delta = 1;
                    case 'HG13', natom = 'CG1'; delta = 1;
                    case 'QG1',  natom = 'CG1'; delta = 1;
                    case 'HD11', natom = 'QD1'; delta = 1;
                    case 'HD12', natom = 'QD1'; delta = 1;
                    case 'HD13', natom = 'QD1'; delta = 1;
                    case 'QD1'
                end
            case 'LEU'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG'
                    case 'HD11', natom = 'QD1'; delta = 1;
                    case 'HD12', natom = 'QD1'; delta = 1;
                    case 'HD13', natom = 'QD1'; delta = 1;
                    case 'QD1'
                    case 'HD21', natom = 'QD2'; delta = 1;
                    case 'HD22', natom = 'QD2'; delta = 1;
                    case 'HD23', natom = 'QD2'; delta = 1;
                    case 'QD2'
                    case 'QQD'
                end
            case 'LYS'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG2',  natom = 'CG'; delta = 1;
                    case 'HG3',  natom = 'CG'; delta = 1;
                    case 'QG',   natom = 'CG'; delta = 1;
                    case 'HD2',  natom = 'CD'; delta = 1;
                    case 'HD3',  natom = 'CD'; delta = 1;
                    case 'QD',   natom = 'CD'; delta = 1;
                    case 'HE2',  natom = 'CE'; delta = 1;
                    case 'HE3',  natom = 'CE'; delta = 1;
                    case 'QE',   natom = 'CE'; delta = 1;
                    case 'HZ1',  natom = 'QZ'; delta = 1;
                    case 'HZ2',  natom = 'QZ'; delta = 1;
                    case 'HZ3',  natom = 'QZ'; delta = 1;
                    case 'QZ'
                end
            case 'MET'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG2',  natom = 'CG'; delta = 1;
                    case 'HG3',  natom = 'CG'; delta = 1;
                    case 'QG',   natom = 'CG'; delta = 1;
                    case 'HE1',  natom = 'QE'; delta = 1;
                    case 'HE2',  natom = 'QE'; delta = 1;
                    case 'HE3',  natom = 'QE'; delta = 1;
                    case 'QE'
                end
            case 'PHE'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HD1'
                    case 'HD2'
                    case 'QD'
                    case 'HE1'
                    case 'HE2'
                    case 'QE'
                    case 'QR'
                    case 'HZ'
                end
            case 'PRO'
                switch atom
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                    case 'HG2'
                    case 'HG3'
                    case 'QG'
                    case 'HD2'
                    case 'HD3'
                    case 'QD'
                end
            case 'SER'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HG',   natom = 'OG'; delta = 1;
                end
            case 'THR'
                switch atom
                    case 'HB'
                    case 'HG21', natom = 'QG2'; delta = 1;
                    case 'HG22', natom = 'QG2'; delta = 1;
                    case 'HG23', natom = 'QG2'; delta = 1;
                    case 'QG2'
                    case 'HG1',  natom = 'OG1'; delta = 1;
                end
            case 'TRP'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HD1'
                    case 'HE1'
                    case 'HE3'
                    case 'HZ2'
                    case 'HZ3'
                    case 'HH2'
                end
            case 'TYR'
                switch atom
                    case 'HB2',  natom = 'CB'; delta = 1;
                    case 'HB3',  natom = 'CB'; delta = 1;
                    case 'QB',   natom = 'CB'; delta = 1;
                    case 'HD1'
                    case 'HD2'
                    case 'QD'
                    case 'HE1'
                    case 'HE2'
                    case 'QE'
                    case 'QR'
                    case 'HH',   natom = 'OH'; delta = 1;
                end
            case 'VAL'
                switch atom
                    case 'HB'
                    case 'HG11', natom = 'QG1'; delta = 1;
                    case 'HG12', natom = 'QG1'; delta = 1;
                    case 'HG13', natom = 'QG1'; delta = 1;
                    case 'QG1'
                    case 'HG21', natom = 'QG2'; delta = 1;
                    case 'HG22', natom = 'QG2'; delta = 1;
                    case 'HG23', natom = 'QG2'; delta = 1;
                    case 'QG2'
                    case 'QQG'
                end
        end
        
        %**********************************************************************
    case 'xplor'        
        switch residuetype
            case 'ALA'
                switch atom
                    case 'HB1'
                    case 'HB2'
                    case 'HB3'
                    case 'QB'
                end
            case 'ARG'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG2'
                    case 'HG1', natom = 'HG3'; delta = -1;
                    case 'QG'
                    case 'HD2'
                    case 'HD1', natom = 'HD3'; delta = -1;
                    case 'QD'
                    case 'HH11'
                    case 'HH12'
                    case 'QH1'
                    case 'HH21'
                    case 'HH22'
                    case 'QH2'
                end
            case 'ASN'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HD21'
                    case 'HD22'
                    case 'QD2'
                end
            case 'ASP'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                end
            case 'CYS'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'HG'
                    case 'QB'
                end
            case 'CYSS'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                end
            case 'GLN'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG2'
                    case 'HG1', natom = 'HG3'; delta = -1;
                    case 'QG'
                    case 'HE21'
                    case 'HE22'
                    case 'QE2'
                end
            case 'GLU'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG2'
                    case 'HG1', natom = 'HG3'; delta = -1;
                    case 'QG'
                end
            case 'GLY'
                switch atom
                    case 'HA2'
                    case 'HA1', natom = 'HA3'; delta = -1;
                    case 'QA'
                end
            case 'HIS'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HD1'
                    case 'HD2'
                    case 'HE1'
                end
            case 'ILE'
                switch atom
                    case 'HB'
                    case 'HG21'
                    case 'HG22'
                    case 'HG23'
                    case 'QG2'
                    case 'HG12'
                    case 'HG11', natom = 'HG13'; delta = -1;
                    case 'QG1'
                    case 'HD11'
                    case 'HD12'
                    case 'HD13'
                    case 'QD1'
                end
            case 'LEU'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG'
                    case 'HD11'
                    case 'HD12'
                    case 'HD13'
                    case 'QD1'
                    case 'HD21'
                    case 'HD22'
                    case 'HD23'
                    case 'QD2'
                    case 'QQD'
                end
            case 'LYS'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG2'
                    case 'HG1', natom = 'HG3'; delta = -1;
                    case 'QG'
                    case 'HD2'
                    case 'HD1', natom = 'HD3'; delta = -1;
                    case 'QD'
                    case 'HE2'
                    case 'HE1', natom = 'HE3'; delta = -1;
                    case 'QE'
                    case 'HZ1'
                    case 'HZ2'
                    case 'HZ3'
                    case 'QZ'
                end
            case 'MET'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG2'
                    case 'HG1', natom = 'HG3'; delta = -1;
                    case 'QG'
                    case 'HE1'
                    case 'HE2'
                    case 'HE3'
                    case 'QE'
                end
            case 'PHE'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HD1'
                    case 'HD2'
                    case 'QD'
                    case 'HE1'
                    case 'HE2'
                    case 'QE'
                    case 'QR'
                    case 'HZ'
                end
            case 'PRO'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG2'
                    case 'HG1', natom = 'HG3'; delta = -1;
                    case 'QG'
                    case 'HD2'
                    case 'HD1', natom = 'HD3'; delta = -1;
                    case 'QD'
                end
            case 'SER'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HG'
                end
            case 'THR'
                switch atom
                    case 'HB'
                    case 'HG21'
                    case 'HG22'
                    case 'HG23'
                    case 'QG2'
                    case 'HG1'
                end
            case 'TRP'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HD1'
                    case 'HE1'
                    case 'HE3'
                    case 'HZ2'
                    case 'HZ3'
                    case 'HH2'
                end
            case 'TYR'
                switch atom
                    case 'HB2'
                    case 'HB1', natom = 'HB3'; delta = -1;
                    case 'QB'
                    case 'HD1'
                    case 'HD2'
                    case 'QD'
                    case 'HE1'
                    case 'HE2'
                    case 'QE'
                    case 'QR'
                    case 'HH'
                end
            case 'VAL'
                switch atom
                    case 'HB'
                    case 'HG11'
                    case 'HG12'
                    case 'HG13'
                    case 'QG1'
                    case 'HG21'
                    case 'HG22'
                    case 'HG23'
                    case 'QG2'
                    case 'QQG'
                end
        end
        
end

% switch residuetype
%     case 'ALA'
%         switch atom
%             case 'HB1'
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%         end
%     case 'ARG'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG2'
%             case 'HG3'
%             case 'QG'
%             case 'HD2'
%             case 'HD3'
%             case 'QD'
%             case 'HH11'
%             case 'HH12'
%             case 'QH1'
%             case 'HH21'
%             case 'HH22'
%             case 'QH2'
%         end
%     case 'ASN'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HD21'
%             case 'HD22'
%             case 'QD2'
%         end
%     case 'ASP'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%         end
%     case 'CYS'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'HG'
%             case 'QB'
%         end
%     case 'CYSS'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%         end
%     case 'GLN'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG2'
%             case 'HG3'
%             case 'QG'
%             case 'HE21'
%             case 'HE22'
%             case 'QE2'
%         end
%     case 'GLU'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG2'
%             case 'HG3'
%             case 'QG'
%         end
%     case 'GLY'
%         switch atom
%             case 'HA2'
%             case 'HA3'
%             case 'QA'
%         end
%     case 'HIS'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HD1'
%             case 'HD2'
%             case 'HE1'
%         end
%     case 'ILE'
%         switch atom
%             case 'HB'
%             case 'HG21'
%             case 'HG22'
%             case 'HG23'
%             case 'QG2'
%             case 'HG12'
%             case 'HG13'
%             case 'QG1'
%             case 'HD11'
%             case 'HD12'
%             case 'HD13'
%             case 'QD1'
%         end
%     case 'LEU'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG'
%             case 'HD11'
%             case 'HD12'
%             case 'HD13'
%             case 'QD1'
%             case 'HD21'
%             case 'HD22'
%             case 'HD23'
%             case 'QD2'
%             case 'QQD'
%         end
%     case 'LYS'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG2'
%             case 'HG3'
%             case 'QG'
%             case 'HD2'
%             case 'HD3'
%             case 'QD'
%             case 'HE2'
%             case 'HE3'
%             case 'QE'
%             case 'HZ1'
%             case 'HZ2'
%             case 'HZ3'
%             case 'QZ'
%         end
%     case 'MET'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG2'
%             case 'HG3'
%             case 'QG'
%             case 'HE1'
%             case 'HE2'
%             case 'HE3'
%             case 'QE'
%         end
%     case 'PHE'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HD1'
%             case 'HD2'
%             case 'QD'
%             case 'HE1'
%             case 'HE2'
%             case 'QE'
%             case 'QR'
%             case 'HZ'
%         end
%     case 'PRO'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG2'
%             case 'HG3'
%             case 'QG'
%             case 'HD2'
%             case 'HD3'
%             case 'QD'
%         end
%     case 'SER'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HG'
%         end
%     case 'THR'
%         switch atom
%             case 'HB'
%             case 'HG21'
%             case 'HG22'
%             case 'HG23'
%             case 'QG2'
%             case 'HG1'
%         end
%     case 'TRP'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HD1'
%             case 'HE1'
%             case 'HE3'
%             case 'HZ2'
%             case 'HZ3'
%             case 'HH2'
%         end
%     case 'TYR'
%         switch atom
%             case 'HB2'
%             case 'HB3'
%             case 'QB'
%             case 'HD1'
%             case 'HD2'
%             case 'QD'
%             case 'HE1'
%             case 'HE2'
%             case 'QE'
%             case 'QR'
%             case 'HH'
%         end
%     case 'VAL'
%         switch atom
%             case 'HB'
%             case 'HG11'
%             case 'HG12'
%             case 'HG13'
%             case 'QG1'
%             case 'HG21'
%             case 'HG22'
%             case 'HG23'
%             case 'QG2'
%             case 'QQG'
%         end
% end
%
