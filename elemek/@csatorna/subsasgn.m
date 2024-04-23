function csatorna = subsasgn(csatorna,index,val)

switch index.type
    case '.'
        switch index.subs
            case 'L',   csatorna.L = val;
            case 'lambda',   csatorna.lambda = val;
            case 'p', csatorna.p= val;
            case 'v', csatorna.v = val;
            case 'p1', csatorna.p1= val;
            case 'pN', csatorna.pN = val;
            case 'v1', csatorna.v1= val;
            case 'vN', csatorna.vN = val;
            case 't',  csatorna.t = val;
            case 'dtki', csatorna.dtki = val;
            case 'A',  csatorna.A = val;
            case 'tki',  csatorna.tki = val;
            case 'pf_eleje_tipus', csatorna.pf_eleje_tipus = val;
            case 'pf_vege_tipus', csatorna.pf_vege_tipus = val;
            otherwise
                csatorna.tranziens_agelem_2csp=subsasgn(csatorna.tranziens_agelem_2csp,index,val);
        end
    otherwise
        error('Bajjj van...');
end