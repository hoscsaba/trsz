function out = subsref(sziv,index)

switch index.type
    case '.'
        switch index.subs
            case 'n'
                out=sziv.n;
            case 'res'
                out=sziv.res;
            case 'cspe'
                out=sziv.cspe;
            case 'cspv'
                out=sziv.cspv;
            case 'tranziens'
                out=sziv.tranziens;
            case 'uzem'
                out=sziv.uzem;
            otherwise
                out=subsref(sziv.tranziens_agelem_2csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end