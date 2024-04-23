function fy = get_A(csatorna,y)
% Atfolyo keresztmetszet meghatarozasa
% eps=0.999;
x=csatorna.dvB;
% switch csatorna.tipus
%     case 'teglalap'
%         fy = y*x;
%     case 'kor'
%         if y < x*eps
%             r = x/2; 
%             theta = acos(1-y/r);
%             fy = r^2*(theta-sin(2*theta)/2);
%         else
%             %y=x*eps;
%             %r = x/2; theta = acos(1-y/r);
%             %fy = r^2*(theta-sin(2*theta)/2);
%             fy = x^2*pi/4 + (y-x)*x/10;
%         end
% end

if y < csatorna.dvB*0.999
    %r = csatorna.dvB/2;
    theta = acos(1-y/(csatorna.dvB/2));
    fy = (csatorna.dvB/2)^2*(theta-sin(2*theta)/2);
else
    %y=x*eps;
    %r = x/2; theta = acos(1-y/r);
    %fy = r^2*(theta-sin(2*theta)/2);
    fy = csatorna.dvB^2*pi/4 + (y-csatorna.dvB)*csatorna.dvB/10;
end

if ~isreal(fy)
    csatorna.tranziens_agelem_2csp.nev;
    x
    y
    fy
    error('get_A');
end
