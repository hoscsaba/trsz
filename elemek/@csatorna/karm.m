function out = karm(csatorna,dt,dtmin)

ro = csatorna.tranziens_agelem_2csp.ro;
p1 = 1000*9.81*(csatorna.y(1)+csatorna.ze) + csatorna.p0 ;
AA = get_A(csatorna,csatorna.y(1));
m = AA*csatorna.v(1)*ro;
out = [p1, m, 1];