function out = subsref(szip,index)

switch index.type
    case '.'
        switch index.subs
            case 'pnyit',   out = szip.pnyit;
            otherwise
                out=subsref(szip.tranziens_agelem_1csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end