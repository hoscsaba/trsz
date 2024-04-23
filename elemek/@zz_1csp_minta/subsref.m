function out = subsref(med,index)

switch index.type
    case '.'
        switch index.subs
            case 'p',   out = med.p;
            otherwise
                out=subsref(med.tranziens_agelem_1csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end