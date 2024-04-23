function fy = get_B(csatorna,y)
% Atfolyo keresztmetszet meghatarozasa
% eps=0.999;
x=csatorna.dvB;
% switch csatorna.tipus
%     case 'teglalap'
%         fy=x;
%     case 'kor'
%         if y<x*eps
%             r = x/2; theta=acos(1-y/r);
%             fy=2*r*sin(pi-theta);
%         else            
%             %y=x*eps;
%             %r = x/2; theta=acos(1-y/r);
%             %fy=2*r*sin(pi-theta);
%             fy=x/10;
%         end
%     otherwise
%         error(['Ismeretlen csatorna tipus:',csatorna.type]);
% end

if y<csatorna.dvB*0.999
    r = csatorna.dvB/2; theta=acos(1-y/r);
    fy=2*r*sin(pi-theta);
else
    %y=x*eps;
    %r = x/2; theta=acos(1-y/r);
    %fy=2*r*sin(pi-theta);
    fy=csatorna.dvB/10;
end

if ~isreal(fy)
    csatorna.tranziens_agelem_2csp.nev;
    x
    y
    fy
    error('get_B');
end