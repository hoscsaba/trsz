function sziv = subsasgn(sziv,index,val)

switch index.type
    case '.'
        switch index.subs
            otherwise
                sziv.tranziens_agelem_2csp=subsasgn(sziv.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end