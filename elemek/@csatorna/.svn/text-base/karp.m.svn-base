function out = karp(csatorna,dt,dtmin)
ro = csatorna.tranziens_agelem_2csp.ro;
pv = 1000*9.81*(csatorna.y(end)+csatorna.zv) + csatorna.p0 ;
%AA = get_A(csatorna,csatorna.y(end));
AA = get_A(csatorna,csatorna.y(end-1));
%m = AA*csatorna.v(end)*ro;
m = AA*csatorna.v(end-1)*ro;
out = [pv, m, -1];
