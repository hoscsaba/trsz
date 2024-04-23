function ellenallas = subsasgn(ellenallas,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'K0', ellenallas.K0 = val;
            case 'K1', ellenallas.K1 = val;
            otherwise
                ellenallas.tranziens_agelem_2csp=subsasgn(ellenallas.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Baj van...');
end