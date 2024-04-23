function cso = subsasgn(cso,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'K0', cso.K0 = val;
            case 'K1', cso.K1 = val;
            otherwise
                cso.tranziens_agelem_2csp=subsasgn(cso.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end