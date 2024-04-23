function viszkcso = subsasgn(viszkcso,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'L',   viszkcso.L = val; 
            case 'lambda',   viszkcso.lambda = val;
            case 'p', viszkcso.p= val; 
            case 'v', viszkcso.v = val;         
            case 'p1', viszkcso.p1= val; 
            case 'pN', viszkcso.pN = val;                
            case 'v1', viszkcso.v1= val; 
            case 'vN', viszkcso.vN = val;                
            case 't',  viszkcso.t = val; 
            case 'all',  viszkcso.all = val; 
            case 'dtki',  viszkcso.dtki = val; 
                
            otherwise
                viszkcso.tranziens_agelem_2csp=subsasgn(viszkcso.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end