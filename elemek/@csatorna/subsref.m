function out = subsref(csatorna,index)

switch index.type
    case '.'
        switch index.subs
            case 'p1', out = csatorna.p1;
            case 'pN', out = csatorna.pN;
            case 'v1', out = csatorna.v1;
            case 'vN', out = csatorna.vN;
            case 'v',  out = csatorna.v;
            case 'p',  out = 1000*9.81*csatorna.y;
            case 'x',  out = csatorna.x;
            case 'dt', out = csatorna.dt;
            case 't',  out = csatorna.t;
            case 'tki',  out = csatorna.tki;
            case 'dtki',  out = csatorna.dtki;
            case 'A',  out = csatorna.A;
            case 'a',  out = csatorna.a;
            case 'lambda',  out = csatorna.lambda;
            case 'D',  out = csatorna.D;
            case 'L',  out = csatorna.L;
            case 'N',  out = csatorna.N;
            case 'y',  out = csatorna.y;
            case 'ze',  out = csatorna.ze;
            case 'zv',  out = csatorna.zv;                
            case 'dvB',  out = csatorna.dvB;
            case 'tipus',  out = csatorna.tipus;
            case 'pf_eleje_tipus', out = csatorna.pf_eleje_tipus;
            case 'pf_vege_tipus',  out = csatorna.pf_vege_tipus;
            otherwise
                out=subsref(csatorna.tranziens_agelem_2csp,index);
        end
    otherwise
        error('HIBA!!! mezohivatkozas .-al!')
end