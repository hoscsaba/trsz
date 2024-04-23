function [out,cso]=setport(flag,cso,varargin)

ro  = cso.tranziens_agelem_2csp.ro;
csp = cso.tranziens_agelem_2csp.csp;
cspe = csp(1);
cspv = csp(2);
Qr  = cso.tranziens_agelem_2csp.Q;

L = cso.L; D = cso.D; A=cso.A; lam=cso.lambda;
ze = cso.ze; zv = cso.zv;
g = 9.81;
%p0 = 1e5;

switch flag
    case -1 % inicializalas
        out{1}=1;
        
    case 0 % stacioner szamitas
        out{1}={0, 0, lam*L/D/2/A^2/ro, 0, [csp(1),-1e5] [csp(2),1e5]};
        
    case 1 % tranziens szamitas
        % xr = varargin{1};
        % t  = varargin{2}{1}{1}(1);
        dt = varargin{2}{1}{1}(2);
        mar = varargin{3};
        
        if cso.init == 0
            k = 0;
            for i = 1:length(mar.elemek)
                if isa(mar.elemek{i},'akna')
                    csp_akna = mar.elemek{i}.csp;
                    if csp_akna(1) == cspe
                        he = mar.elemek{i}.y;
                        cso.aknae = i;
                    elseif csp_akna(1) == cspv
                        hv = mar.elemek{i}.y;
                        cso.aknav = i;
                    end
                elseif isa(mar.elemek{i},'szivattyu')
                    %h = mar.elemek{i};
                    csp_sziv = mar.elemek{i}.csp;
                    if csp_sziv(2) == csp(1)
                        k = k+1;
                        cso.sziv{k} = i;
                        for j = 1:length(mar.elemek)
                            if isa(mar.elemek{j},'akna')
                                csp_akna=mar.elemek{j}.csp;
                                if csp_akna(1) == csp_sziv(1)
                                    he = mar.elemek{j}.y;
                                    %mar.elemek{j}.nev
                                    cso.aknae = j;
                                end
                            end
                        end
                    end
                end
            end
            cso.init = 1;
        else
            he = mar.elemek{cso.aknae}.y;
            hv = mar.elemek{cso.aknav}.y;
        end
        
        uzem = 0;
        for j = 1:length(cso.sziv)
            uzem = max(uzem,mar.elemek{cso.sziv{j}}.uzem);
        end
        
        %if strcmp(cso.nev,'ny5')
        
        %    he
        %    hv
        %    uzem
        %end
        
        %         if he > ze && hv > zv
        %             dummy = 0;
        %         elseif he > ze && hv <= zv
        %             if ze > zv
        %                 dummy = ro*g*(zv - hv);
        %             else
        %                 dummy = ro*g*(zv - hv)- ro*g*(zv-ze);
        %                 if he < zv
        %                     dummy = ro*g*(he-hv);
        %                     if uzem == 0, Qr = 0; end
        %                 end
        %             end
        %         elseif he <= ze && hv > zv
        %             if ze > zv
        %                 if ze > hv
        %                     dummy = ro*g*(he-hv);
        %                 else
        %                     dummy = ro*g*(he-hv) - ro*g*(ze-hv);
        %                 end
        %             else
        %                 dummy = ro*g*(he - ze);
        %             end
        %         else
        %             dummy = ro*g*(he-hv);
        %             if uzem == 0, Qr = 0; end
        %         end
        
        
        
        if ze > zv
            if he > ze && hv > zv
                dummy = 0;
            elseif he > ze && hv <= zv
                dummy = - ro*g*hv;
            elseif he <= ze && hv > ze
                dummy = ro*g*he;
            elseif he <= ze && hv <= ze
                dummy = ro*g*(he - hv);
                if uzem == 0, Qr = 0.0; end
            end
        else
            if he > ze && hv > zv
                if he < hv
                    dummy = ro*g*(he - hv);
                    if uzem == 0, Qr = 0.0; end
                else
                    dummy = 0;
                end
                %                 if hv-zv < 0.01
                %                     dummy = 0.5*ro*g*(he - hv );
                %                 end
            elseif he > zv && hv <= zv
                dummy = - ro*g*(hv);
            elseif he <= zv && hv > zv
                dummy = ro*g*(he - hv);
                if uzem == 0, Qr = 0.0; end
            elseif he <= zv && hv <= zv
                dummy = ro*g*(he - hv);
                if uzem == 0, Qr = 0.0; end
            end
        end
        
        out{1}={L/A/dt, 0, lam*L/D/2/A^2/ro, dummy-L/A/dt*Qr*ro, [csp(1),-1e5] [csp(2),1e5]};
    case 3
        out{1}=1;
        
    otherwise
        error('Ilyen opcio nincs; vagy stac vagy instac!');
end

