function cso = subsasgn(cso,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'L'
                cso.L = val;
            otherwise
                cso.tranziens_agelem_2csp=subsasgn(cso.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end