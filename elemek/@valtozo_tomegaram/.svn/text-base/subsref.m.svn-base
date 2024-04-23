function out = subsref(tomegaram,index)

switch index.type
    case '.'
        switch index.subs
            case 'p',   out = tomegaram.p;
            otherwise
                out=subsref(tomegaram.tranziens_agelem_1csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end