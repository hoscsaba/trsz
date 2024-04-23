function out = subsref(cso,index)

switch index.type
    case '.'
        switch index.subs
            case 'Cd'
                out = cso.Cd;
            case 'init'
                out = cso.init;
            case 'h0'
                out = cso.h0;
            case 'kitevo'
                out = cso.kitevo;
                
            otherwise
                out=subsref(cso.tranziens_agelem_2csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end