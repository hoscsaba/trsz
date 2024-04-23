function out = subsref(ak,index)

switch index.type
    case '.'
        switch index.subs
            case 'p',     out = ak.p;
            case 'y',     out = ak.y;
            case 'A',     out = ak.A;
            case 'hmin',  out = ak.hmin;           
            case 'hmax',  out = ak.hmax;
            case 'ures',  out = ak.ures;
            case 'tele',  out = ak.tele;
            otherwise
                out = subsref(ak.tranziens_agelem_1csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end