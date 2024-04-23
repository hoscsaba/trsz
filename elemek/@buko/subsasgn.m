function buko = subsasgn(buko,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'Cd'
                buko.Cd = val;
            case 'is_working'
                buko.is_working=val;
            otherwise
                buko.tranziens_agelem_2csp=subsasgn(buko.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end