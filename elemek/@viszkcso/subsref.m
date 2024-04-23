function out = subsref(viszkcso,index)

switch index.type
    case '.'
        switch index.subs
            case 'p1', out = viszkcso.p1;
            case 'pN', out = viszkcso.pN;
            case 'v1', out = viszkcso.v1;
            case 'vN', out = viszkcso.vN;
            case 'v',  out = viszkcso.v;
            case 'p',  out = viszkcso.p;
            case 'x',  out = viszkcso.x;
            case 'dt', out = viszkcso.dt;
            case 'all', out = viszkcso.all;
            case 'dtki', out = viszkcso.dtki;
            case 'dtuj', out = viszkcso.dtuj;
            case 't',  out = viszkcso.t;
            case 'A',  out = viszkcso.A;
            case 'a',  out = viszkcso.a;
            case 'lambda',  out = viszkcso.lambda;
            case 'D',  out = viszkcso.D;
            case 'L',  out = viszkcso.L;
            case 'N',  out = viszkcso.N;
            case 'tipus',  out = viszkcso.tipus;
            case 'pf_eleje_tipus', out = viszkcso.pf_eleje_tipus;
            case 'pf_vege_tipus',  out = viszkcso.pf_vege_tipus;
            otherwise
                out=subsref(viszkcso.tranziens_agelem_2csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al!')
end