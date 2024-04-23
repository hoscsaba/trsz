function out = subsref(elem,index)

switch index(1).type
    case '.'
        switch index(1).subs
            case 'nev'
                out=elem.nev;
            case 'x'
                out=elem.x;
            case 'y'
                out=elem.y;
            otherwise
                error('Ismeretlen mezo: ',index.type);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al!')
end