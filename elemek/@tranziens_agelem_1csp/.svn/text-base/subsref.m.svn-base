function out = subsref(trag1,index)

switch index.type
    case '.'      
        switch index.subs         
            case 'csp', out = trag1.csp;
            case 'p',   out = trag1.p;
            case 'nev', out = trag1.nev;
            case 'ro',  out = trag1.ro;
            case 'Q',   out = trag1.Q;
            case 'fignum',   out = trag1.fignum;
            otherwise, error(['Ismeretlen mezo:',index.subs])
        end   
    otherwise
        error('HIBA!!! mezohivatkozas .-al (pl. elem.ro)')
end