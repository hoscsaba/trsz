function out = subsref(trag2,index)

switch index.type
    case '.'
        switch index.subs
            case 'csp',  out = trag2.csp;
            case 'p',    out = trag2.p;
            case 'nev',  out = trag2.nev;
            case 'ro',   out = trag2.ro;
            case 'Q',    out = trag2.Q;
            case 'fignum',   out = trag2.fignum;
            otherwise, error(['Ismeretlen mezo:',index.subs])
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)');
end