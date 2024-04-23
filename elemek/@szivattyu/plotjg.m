function plotjg(sziv)

figure
plot(sziv.jgQ(2:length(sziv.jgQ)-1),sziv.jgH(2:length(sziv.jgQ)-1),'ro',sziv.jgQ(2:length(sziv.jgQ)-1),sziv.jgH(2:length(sziv.jgQ)-1));
xlabel('Q [m^3/s]'),ylabel('H [m]'), grid on
% axis([0 1.1*max(sziv.jgQ) 0 1.2*max(sziv.jgH)])
title(['A "',sziv.nev,'" szivattyu jelleggorbeje'])
