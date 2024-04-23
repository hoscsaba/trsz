function trag1 = subsasgn(trag1,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'csp', trag1.csp = val;            
            case 'p',   trag1.p  = val;
            case 'ro',  trag1.ro = val;
            case 'Q',   trag1.Q  = val;
            case 'Qr',  trag1.Qr = val;
            otherwise, error(['Ismeretlen mezo: ',index.subs]);
        end
end