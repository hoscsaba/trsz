function cso = subsasgn(cso,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'L',   cso.L = val;
            case 'lambda',   cso.lambda = val;
            case 'p', cso.p= val;
            case 'v', cso.v = val;
            case 'p1', cso.p1= val;
            case 'pN', cso.pN = val;
            case 'v1', cso.v1= val;
            case 'vN', cso.vN = val;
            case 't',  cso.t = val;
            case 'dtki',  cso.dtki = val;
            case 'pf_eleje_tipus', cso.pf_eleje_tipus = val;
            case 'pf_vege_tipus', cso.pf_vege_tipus = val;
            otherwise
                cso.tranziens_agelem_2csp=subsasgn(cso.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end