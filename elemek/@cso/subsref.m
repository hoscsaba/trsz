function out = subsref(cso,index)

switch index.type
    case '.'
        switch index.subs
            case 'p1', out = cso.p1;
            case 'pN', out = cso.pN;
            case 'v1', out = cso.v1;
            case 'vN', out = cso.vN;
            case 'v',  out = cso.v;
            case 'p',  out = cso.p;
            case 'x',  out = cso.x;
            case 'dt', out = cso.dt;
            case 'dtki',  out = cso.dtki;
            case 't',  out = cso.t;
            case 'A',  out = cso.A;
            case 'a',  out = cso.a;
            case 'lambda',  out = cso.lambda;
            case 'D',  out = cso.D;
            case 'L',  out = cso.L;
            case 'N',  out = cso.N;
            case 'pf_eleje_tipus', out = cso.pf_eleje_tipus;
            case 'pf_vege_tipus',  out = cso.pf_vege_tipus;
            otherwise
                out=subsref(cso.tranziens_agelem_2csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al!')
end