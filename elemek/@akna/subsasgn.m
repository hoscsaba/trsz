function akna = subsasgn(akna,index,val)

switch index.type
    case '.'
        switch index.subs    
            case 'fignum'
                akna.tranziens_agelem_1csp.fignum=val;
            otherwise
                akna.tranziens_agelem_1csp = subsasgn(akna.tranziens_agelem_1csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end
