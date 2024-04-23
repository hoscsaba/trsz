function trag2 = subsasgn(trag2,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'csp', trag2.csp = val;            
            case 'p',   trag2.p  = val;
            case 'ro',  trag2.ro = val;
            case 'Q',   trag2.Q  = val;
            case 'Qr',  trag2.Qr = val;
            otherwise, error(['Ismeretlen mezo: ',index.subs]);
        end
end
