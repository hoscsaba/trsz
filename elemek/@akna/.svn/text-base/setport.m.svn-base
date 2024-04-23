function [out,akna] = setport(flag,akna,varargin)

ro  = akna.tranziens_agelem_1csp.ro;
csp = akna.tranziens_agelem_1csp.csp;
g = 9.81;
p0 = 1e5;

switch flag
    case 1 % Newton-Raphson szamitas, nincs belso update
        xr = varargin{1};
        t  = varargin{2}{1}{1}(1);
        dt = varargin{2}{1}{1}(2);
        mar = varargin{3};
     
        out{1} = {0, 0, 0, p0 + ro*g*akna.y, [csp, -1e5]};

    case 3 % NR vege, belso update

        xr = varargin{1};
        t  = varargin{2}(1);
        dt = varargin{2}(2);
        mar = varargin{3};
        Q = xr;
        Qr = akna.Qr;
        Qrr = akna.Qrr;
        
        Q1 = Q;
        %%%
        %Q = (Q + Qr + Qrr)/3;
        
        Q = interp1([0, dt, 2*dt],[Qrr, Qr, Q], 2*dt+dt,'spline','extrap');
        dV = (xr + Q)*dt/2;
        Q2 = Q;

        %%%
        
        %yuj = akna.y + Q/akna.A*dt;
        yuj = akna.y +dV/akna.A;

        if yuj <= akna.hmin && akna.ures == 0
            yuj = akna.hmin;
            akna.ures = 1;
            fprintf('\nA(z) %s akna kiurult!\n',akna.tranziens_agelem_1csp.nev);
            Q = 0;
        elseif yuj >= akna.hmax && akna.tele == 0
            akna.tele = 1;
            fprintf('\nMegtelt a(z) %s akna!\n',akna.tranziens_agelem_1csp.nev);
        elseif yuj > akna.hmin && yuj < akna.hmax
            akna.ures = 0;
            akna.tele = 0;
        end
    
        % fprintf('\n%s akna: yregi=%8.3f, Qr=%5.3f, Q=%5.3f, dt=%5.3e, A=%5.3f, yuj=%8.3f\n',akna.tranziens_agelem_1csp.nev,akna.y,Qr,Q,dt,akna.A,yuj);

        akna.yr = akna.y;
        akna.Q = Q;
        akna.Qr = xr;
        akna.Qrr = Qr;
        akna.y = yuj;
        out = 1;
    
        %% RAJZ 
        if strcmp(akna.rajz,'rajz')
            figure(akna.tranziens_agelem_1csp.fignum)
            if akna.y > akna.hmax
                plot([t-dt t]/60,[akna.yr,akna.y],'r-')
            else
                plot([t-dt t]/60,[akna.yr,akna.y],'b-')
            end
            hold on, grid on
            title([akna.tranziens_agelem_1csp.nev,', h=',num2str(akna.y),...
                ' m, hfedlap=',num2str(akna.hmax),' m, Q=',num2str(akna.Q*3600),' m^3/h'])
            xlabel('t [min.]'), ylabel('vizszint, m'), hold on
        end
               
    otherwise
        error(['%s akna: flag=',num2str(flag)]);
end
