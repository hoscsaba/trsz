function [out,cso]=setport(flag,cso,varargin)

ro  = cso.tranziens_agelem_2csp.ro;
csp = cso.tranziens_agelem_2csp.csp;
Qr  = cso.tranziens_agelem_2csp.Q;
L = cso.L; D = cso.D; A=cso.A; lam=cso.lambda;
ze = cso.ze; zv = cso.zv;
g = 9.81;
p0 = 1e5;

switch flag
    case -1 % inicializálás
        out{1}=1;
        
    case 0 % stacioner számítás
        out{1}={0, 0, lam*L/D/2/A^2/ro, 0, [csp(1),-1e5] [csp(2),1e5]};
        
    case 1 % tranziens számítás
%        xr = varargin{1};
	    t  = varargin{2}{1}{1}(1);
	    dt = varargin{2}{1}{1}(2);
        mar = varargin{3};
        
        
        if cso.init == 0
            k = 0;
            for i = 1:length(mar.elemek)
                if isa(mar.elemek{i},'akna')
                    if mar.elemek{i}.csp == cso.cspe
                        he = mar.elemek{i}.y;
                        %mar.elemek{i}.nev
                        cso.aknae = i;
                    elseif mar.elemek{i}.csp == cso.cspv
                        hv = mar.elemek{i}.y;
                        %mar.elemek{i}.nev
                        cso.aknav = i;
                    end
                elseif isa(mar.elemek{i},'szivattyu')
                    h = mar.elemek{i};
                    if mar.elemek{i}.cspv == cso.cspe
                        k = k+1;
                        cso.sziv{k} = i;
                        for j = 1:length(mar.elemek)
                            if isa(mar.elemek{j},'akna')
                                if mar.elemek{j}.csp == mar.elemek{i}.cspe
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
                      
        if he > ze && hv > zv
            dummy = 0;
        elseif he > ze && hv <= zv  
            if ze > zv
                dummy = ro*g*(zv - hv);
            else
                dummy = ro*g*(zv - hv)- ro*g*(zv-ze);
                if he < zv
                    dummy = ro*g*(he-hv);
                    if uzem == 0, Qr = 0; end
                end
            end
        elseif he <= ze && hv > zv
            if ze > zv
                if ze > hv
                    dummy = ro*g*(he-hv);
                else
                    dummy = ro*g*(he-hv) - ro*g*(ze-hv);
                end
            else
                dummy = ro*g*(he - ze);
            end
        else
            dummy = ro*g*(he-hv);
            if uzem == 0, Qr = 0; end
        end
        
        out{1}={L/A/dt, 0, lam*L/D/2/A^2/ro, dummy-L/A/dt*Qr*ro, [csp(1),-1e5] [csp(2),1e5]};
   case 3
        out{1}=1;
	
    otherwise
        error('Ilyen opcio nincs; vagy stac vagy instac!');
        end

