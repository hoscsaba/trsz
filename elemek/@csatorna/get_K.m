function fy = get_K(csatorna,y)
% Nedvesitett kerulet meghatarozosa
% eps=0.999;
x=csatorna.dvB;
% switch csatorna.tipus
%     case 'teglalap'
%         fy=2*y+x;
%     case 'kor'
%         if y<x*eps
%             r = x/2; theta=acos(1-y/r);
%             fy=2*r*theta;
%         else
%             %y=x*eps;
%             %r = x/2; theta=acos(1-y/r);
%             %fy=2*r*theta;
%             fy=2*x*pi+2*(y-x*eps);
%         end
%     otherwise
%         error(['Ismeretlen csatorna tipus:',csatorna.type]);
% end

if y<csatorna.dvB*0.999
    r = csatorna.dvB/2; theta=acos(1-y/r);
    fy=2*r*theta;
else
    %y=x*eps;
    %r = x/2; theta=acos(1-y/r);
    %fy=2*r*theta;
    fy=2*csatorna.dvB*pi+2*(y-csatorna.dvB*eps);
end

if ~isreal(fy)
    csatorna.tranziens_agelem_2csp.nev;
    x
    y
    fy
    error('get_K');
end