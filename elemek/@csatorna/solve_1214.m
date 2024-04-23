function csatorna = solve(csatorna,pf)
%% Peremfeltetelek atadasa
bct1 = pf{1}{1}; bcv1 = pf{1}{2};
bct2 = pf{2}{1}; bcv2 = pf{2}{2};

ro = csatorna.tranziens_agelem_2csp.ro;
nev = csatorna.tranziens_agelem_2csp.nev;
show_csat='csat22';

if strcmp(bct1,'p')
    y_akna_beg = (bcv1 - csatorna.p0)/(ro*9.81) - csatorna.ze;
else
    error('Rossz peremfeltetel az elejen!')
end

if strcmp(bct2,'p')
    y_akna_end = (bcv2 - csatorna.p0)/(ro*9.81) - csatorna.zv;
else
    error('Rossz peremfeltetel a vegen!')
end


%% Csatorna adatok, halogeneralas
type = csatorna.tipus;
%type 1-Vegtelen magas teglalap keresztmetszetu csatorna.
%     2-Kor keresztmetszetu csatorna. (1% Hasitek)
X = csatorna.dvB; % atmero v. szelesseg
He = csatorna.ze; % csatorna elje folyasi szint
Hv = csatorna.zv; % csatorna vege folyasi szint
n = csatorna.n;   % manning allando
L = csatorna.L;   % csatorna hosz
s = (He-Hv)/L;    % lejtes
g = 9.81;         % gravitacios gyorsulas
N = csatorna.N+1; % osztaspontok szama
Xmesh = linspace(0,L,N); dx = Xmesh(2)-Xmesh(1);

%% Kezdeti feltetelek, idolepes beallitas
dt = csatorna.dt;
Ytr = zeros(1,N); Vtr = zeros(1,N);
Ytr(1:N) = csatorna.y; Vtr(1:N) = csatorna.v;
Ytr_temp = Ytr; Ytr_temp( Ytr_temp>X*0.999 )=X*0.999;
% Keresztmetszet

theta = acos(1-Ytr_temp/(X/2));
AA = (X/2)^2*(theta-sin(2*theta)/2);
BB = X*sin(pi-theta);
KK = X*theta;

% AA = zeros(1,N);
% BB = zeros(1,N);
% for i = 1:length(Ytr)
%     if Ytr(i) < csatorna.dvB*0.999
%         %r = csatorna.dvB/2;
%         theta = acos(1-Ytr(i)/(csatorna.dvB/2));
%         AA(i) = (csatorna.dvB/2)^2*(theta-sin(2*theta)/2);
%         BB(i) = X*sin(pi-theta);
%         KK(i) = X*acos(1-Ytr(i)*2/X);
%     else
%         AA(i) = csatorna.dvB^2*pi/4 + (Ytr(i)-csatorna.dvB)*csatorna.dvB/10;
%         BB(i) = X/10;
%         KK(i) = 2*X*pi+2*(Ytr(i)-X*0.999);
%     end
% end

svs=g*AA./BB;
Rh=AA./KK;

%% Eleje perem

Aeq=zeros(2*N,2*N); Beq=zeros(2*N,1);

% Kontinuitas
Aeq(1,1)=1; Aeq(1,N+1)=0;
Aeq(1,2)=0; Aeq(1,N+2)=0;
Aeq(1,3)=0; Aeq(1,N+3)=0;

% Eleje vizszint leellenorzese kritikus kifolyasra es urulesre
% Q_beg=abs(Vtr(1)*get_A(csatorna,Ytr(1)));
Q_beg=abs(Vtr(1)*AA(1));
yc_beg=crit(csatorna,max([Q_beg 1e-6]));

% % kritikus kifolyas
% if (yc_beg > y_akna_beg) && (Vtr(1)<0)
%     y_value_beg = yc_beg;
%     
% else
%     y_value_beg=y_akna_beg;
% end
% % kiurules
% if y_value_beg < X*0.001;
%     y_value_beg = X*0.001;
% end

y_value_beg = max([y_akna_beg   Ytr(1)-g*dt^2/2   -Vtr(1)*yc_beg/abs(Vtr(1))   y_akna_beg   X*0.001]);
Beq(1) = y_value_beg;

% Mozgasegyenlet
% Aeq(N+1,1)=(-3*g/dx/2);     Aeq(N+1,N+1)=(g*n^2*abs(Vtr(1))/get_Rh(csatorna,Ytr(1))^(4/3) - 3*Vtr(1)/dx/2 + 1/dt)*1; % ??? *1
Aeq(N+1,1)=(-3*g/dx/2);     Aeq(N+1,N+1)=(g*n^2*abs(Vtr(1))/Rh(1)^(4/3) - 3*Vtr(1)/dx/2 + 1/dt)*1; % ??? *1

Aeq(N+1,2)=2*g/dx;          Aeq(N+1,N+2)=2*Vtr(1)/dx;
Aeq(N+1,3)=-g/dx/2;         Aeq(N+1,N+3)=-Vtr(1)/dx/2;
Beq(N+1)=g*s+Vtr(1)/dt;

%% Kozbenso pontok
% svs = sonic velocity squared
% svs2=zeros(1,N);
% for k=1:N
%     %svs(k)=min(10^2,g*get_A(csatorna,Ytr(k))/get_B(csatorna,Ytr(k)));   
%     svs2(k)=g*get_A(csatorna,Ytr(k))/get_B(csatorna,Ytr(k));  
% end



for k=2:N-1
    %% Kontinuitas
    Aeq(k,k-1)=-Vtr(k)/dx/2;    Aeq(k,N+k-1)=-svs(k)/dx/g/2;
    Aeq(k,k)  =1/dt;            Aeq(k,N+k)  =0;
    Aeq(k,k+1)=Vtr(k)/dx/2;     Aeq(k,N+k+1)=svs(k)/dx/g/2;
    Beq(k)=Ytr(k)/dt;
    
    %% Mozgas egyenlet
    Aeq(N+k,k-1)=-g/dx/2;   Aeq(N+k,N+k-1)=-Vtr(k)/dx/2;
    Aeq(N+k,k)  =0;         Aeq(N+k,N+k)  =1/dt+g*n^2*abs(Vtr(k))/Rh(k)^(4/3);
    Aeq(N+k,k+1)=g/dx/2;    Aeq(N+k,N+k+1)=Vtr(k)/dx/2;
    Beq(N+k)=g*s+Vtr(k)/dt;   
end

%% Vege Perem

% Vege vizszint leellenorzese kritikus kifolyasra
%Q_end=abs(Vtr(end-1)*get_A(csatorna,Ytr(end-1)))
Q_end=abs(Vtr(end-1)*AA(end-1));

yc_end=crit(csatorna,max([Q_end 1e-6]));
yn_end=normal(csatorna,max([Q_end 1e-6]),s,n);

% Q_end
% 
% yc_end

% normal szint meghatarozas
if isempty(yn_end)
    yn_end=2*yc_end;
else
    yn_end=yn_end(1);
end
% kritikus kifolyas

% %Vtr(end-10:end)
% if ((yc_end>y_akna_end) && (Vtr(end)>0)) || ...
%         (y_akna_end < X*0.001)
%     if yn_end>yc_end
%         y_value_end=1.1*yc_end;
%     else
%         if yn_end>y_akna_end
%             y_value_end=Ytr(end-1);
%         else
%             y_value_end=y_akna_end;
%         end
%     end
% else
%     y_value_end=y_akna_end;
% end
% 
% % kiurules
% if y_value_end<X*0.001;
%     y_value_end=X*0.001;
% end

y_value_end=max([y_akna_end   Ytr(end)-g*dt^2/2   min([yc_end yn_end])   X*0.001]);

vege_matyak = 'konti_off';

switch vege_matyak
    case 'konti_off'
        % PF.
        Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
        Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
        Aeq(N,N)  =1;   Aeq(N,N+N)  =0;
        Beq(N) = y_value_end;
        
        % Mozgas egyenlet
        Aeq(2*N,N-2) = g/dx/2;      Aeq(2*N,N+N-2) = Vtr(end)/dx/2;
        Aeq(2*N,N-1) = -2*g/dx;     Aeq(2*N,N+N-1) = -2*Vtr(end)/dx;
        Aeq(2*N,N)   = 3*g/dx/2;    Aeq(2*N,N+N)   = g*n^2*abs(Vtr(N))/get_Rh(csatorna,Ytr(end))^(4/3) + 3*Vtr(end)/dx/2 + 1/dt;
%        Aeq(2*N,N)   = 3*g/dx/2;    Aeq(2*N,N+N)   = g*n^2*abs(Vtr(N))/Rh(end)^(4/3) + 3*Vtr(end)/dx/2 + 1/dt;
        Beq(2*N)=g*s+Vtr(end)/dt;
        
    case 'kar'
        
        %sv_end=sqrt(g*get_A(csatorna,csatorna.y(end))/get_B(csatorna,csatorna.y(end)));     
        sv_end=sqrt(g*AA(end)/BB(end));      
        
        if Vtr(end)<-sv_end
            % Hangsebesseg feletti BEaramlas
            % y_akna = y(end)+v(end)^2/2g
            % y(end) = ykrit
            if strcmp(nev,show_csat)
                fprintf('\n\n %s end: Hangsebesseg feletti bearamlas, t=%gs, dt=%gs',nev,csatorna.t,dt);
            end
            
            % y(end) = ykrit
            Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
            Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
            Aeq(N,N)  =1;   Aeq(N,N+N)  =0;
            Beq(N) = y_value_end;
            
            % y_akna = y(end)+v(end)^2/2g
            Aeq(2*N,N-2) = 0;  Aeq(2*N,N+N-2) = 0;
            Aeq(2*N,N-1) = 0;  Aeq(2*N,N+N-1) = 0;
            Aeq(2*N,N)   = 0;  Aeq(2*N,N+N)   = 1;
            Beq(2*N) = -sqrt(max(0,(y_akna_end-yc_end))*2*g);            
            
            if strcmp(nev,show_csat)
                fprintf('\n\t y_akna=%gm, y_krit=%gm',y_akna_end,yc_end);
            end
            
        elseif Vtr(end)>sv_end
            % Hangsebesseg feletti KIaramlas
            % C+ es C- karakterisztika menten a szelre interpolalunk
            
            if strcmp(nev,show_csat)
            fprintf('\n\n %s end: Hangsebesseg feletti kiaramlas, t=%gs, dt=%gs',nev,csatorna.t,dt);
            end
            
            % Karakterisztika
            for i=1:length(Ytr)
                %sv(i)=sqrt(g*get_A(csatorna,csatorna.y(i))/get_B(csatorna,csatorna.y(i)));
                sv(i)=sqrt(g*AA(i)/BB(i));
                dxp(i)=(Vtr(i)+sv(i))*dt;
                dxm(i)=(Vtr(i)-sv(i))*dt;
            end
            xtp=interp1(dxp-(L-Xmesh),Xmesh,0);
            xtm=interp1(dxm-(L-Xmesh),Xmesh,0);
            vtp=interp1(Xmesh,Vtr,xtp);
            vtm=interp1(Xmesh,Vtr,xtm);
            ytp=interp1(Xmesh,Ytr,xtp);
            ytm=interp1(Xmesh,Ytr,xtm);
            svtp=interp1(Xmesh,sv,xtp);
            svtm=interp1(Xmesh,sv,xtm);
            
            alfa=(ytp+svtp/g*vtp)+dt*svtp*(s-vtp*abs(vtp)*n^2/get_Rh(csatorna,ytp)^(4/3));
            beta=(ytm-svtm/g*vtm)-dt*svtm*(s-vtm*abs(vtm)*n^2/get_Rh(csatorna,ytm)^(4/3));
            y_end=(alfa+beta)/2;
            sv_end=sqrt(g*get_A(csatorna,y_end)/get_B(csatorna,y_end));
            v_end=(alfa-y_end)*g/sv_end;           
            
            if strcmp(nev,show_csat)
%                 figure(1)
%                 for i=1:length(Ytr)
%                     plot([Xmesh(i),Xmesh(i)+dxp(i)],[0 dt],'r-x'), hold on
%                     plot([Xmesh(i),Xmesh(i)+dxm(i)],[0 dt],'b--o'), hold on
%                     plot([xtp,xtp+dt*(vtp+svtp)],[0 dt],'k--'), hold on
%                     plot([xtm,xtm+dt*(vtm-svtm)],[0 dt],'k-'), hold on
%                 end
%                 hold off, grid on, xlabel('x'), ylabel('t')
                fprintf('\n\t C+: x=%6.1fm, y=%6.3fm, v=%6.3fm/s, a=%6.3fm/s',xtp,ytp,vtp,svtp);
                fprintf('\n\t C-: x=%6.1fm, y=%6.3fm, v=%6.3fm/s, a=%6.3fm/s',xtm,ytm,vtm,svtm);
                %fprintf('\n\t perem:         y=%6.3fm, v=%6.3fm/s',yend,vend);               
                %pause
            end            
            if y_akna_end<0 && v_end<0
                %v_end = 0.001;
                %y_end = 0.001*X;
            end           
            % PF.
            Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
            Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
            Aeq(N,N)  =0;   Aeq(N,N+N)  =1;
            Beq(N) = v_end;

            
            Aeq(2*N,N-2) = 0;  Aeq(2*N,N+N-2) = 0;
            Aeq(2*N,N-1) = 0;  Aeq(2*N,N+N-1) = 0;
            Aeq(2*N,N)   = 1;  Aeq(2*N,N+N)   = 0;
            Beq(2*N) =  y_end;
        else
            % Hangsebesseg alatti KI/BEaramlas
            % C+ menten a szelre interpolalunk                        
           
            % Karakterisztika
            for i=1:length(Ytr)
                %sv(i)=sqrt(g*get_A(csatorna,csatorna.y(i))/get_B(csatorna,csatorna.y(i)));
                sv(i)=sqrt(g*AA(i)/BB(i));
                dxp(i)=(Vtr(i)+sv(i))*dt;                
            end
            xtp=interp1(Xmesh+dxp,Xmesh,L);
            vtp=interp1(Xmesh,Vtr,xtp);
            ytp=interp1(Xmesh,Ytr,xtp);
            svtp=interp1(Xmesh,sv,xtp);
            
            alfa_end=(ytp+svtp/g*vtp)+dt*svtp*(s-vtp*abs(vtp)*n^2/get_Rh(csatorna,ytp)^(4/3));            
            
            % y(end)+a/g*v(end)=alfa
            Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
            Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
            Aeq(N,N)  =0;   Aeq(N,N+N)  =1;
            sv_end=sqrt(g*get_A(csatorna,y_value_end)/get_B(csatorna,y_value_end));
            v_end=(alfa_end-y_value_end)*g/sv_end;

            %y_value_end
            %v_end
            %Vtr(end)
            %v_end = (v_end + Vtr(end))/2
            %pause
            
            if y_akna_end<0 && v_end<0
                %v_end = 0.001;
                %y_value_end = alfa_end - v_end*sv_end/g;
            end
            
            Beq(N) = v_end;
            
            % y(end)=y_akna
            Aeq(2*N,N-2) = 0;  Aeq(2*N,N+N-2) = 0;
            Aeq(2*N,N-1) = 0;  Aeq(2*N,N+N-1) = 0;
            Aeq(2*N,N)   = 1;  Aeq(2*N,N+N)   = 0;                                    
            Beq(2*N) =  y_value_end;
            
            if strcmp(nev,show_csat)
                fprintf('\n\n %s end: Hangsebesseg alatti ki/bearamlas, t=%gs, dt=%gs, y_akna=%5.3fm',nev,csatorna.t,dt,y_akna_end);                
                fprintf('\n\t xtp: x=%6.1fm, y=%6.3fm, v=%6.3fm/s, a=%6.3fm/s, Q=%5.3fm3/h, A=%5.3fm2',xtp,ytp,vtp,svtp,vtp*get_A(csatorna,ytp)*3600,get_A(csatorna,ytp));
                fprintf('\n\t end: x=%6.1fm, y=%6.3fm, v=%6.3fm/s, a=%6.3fm/s, Q=%5.3fm3/h, A=%5.3fm2',L,y_value_end,v_end,sv_end,v_end*get_A(csatorna,y_value_end)*3600,get_A(csatorna,y_value_end));
            end
        end
        
    otherwise
        % PF.
        Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
        Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
        Aeq(N,N)  =1;   Aeq(N,N+N)  =0;
        Beq(N) = y_value_end;
        
        % Konti
        %svs = sonic velocity squared
        svs=g*get_A(csatorna,Ytr(end))/get_B(csatorna,Ytr(end));
        Aeq(2*N,N-2) =  Vtr(end)/2/dx;        Aeq(2*N,N+N-2) = svs/g/2/dx;
        Aeq(2*N,N-1) = -2*Vtr(end)/dx;        Aeq(2*N,N+N-1) = -2*svs/g/dx;
        Aeq(2*N,N)   = 1/dt+3*Vtr(end)/2/dx;  Aeq(2*N,N+N)   = svs/g*3/2/dx;
        Beq(2*N) = Ytr(end)/dt;
        
end
%% Megoldas
[mo flag]=bicgstab(Aeq,Beq,1e-3,[],[],[],[Ytr(1,:), Vtr(1,:)]');
Ytr=mo(1:N)'; Vtr=mo(N+1:2*N)';

if strcmp(nev,show_csat)
    fprintf('\n \t y(end)=%g, v(end)=%g',Ytr(end),Vtr(end));
%     sparse(imag(Aeq))
%     sparse(imag(Beq))
end

%% Simitas
%---------------------------------------------------
%% Mesterseges viszkozitas
% av_c1=1; av_c2=av_c1/10;
% dv=zeros(size(Vtr)); hull_seb=dv; av_dy=dv;
% for i=2:length(Vtr)-1
%     if Vtr(i)>0
%         dv(i)=Vtr(i)-Vtr(i-1);
%     else
%         dv(i)=Vtr(i+1)-Vtr(i);
%     end
%     if dv(i)<0
%         AA = get_A(csatorna,Ytr(i));
%         BB = get_B(csatorna,Ytr(i));
%         hull_seb(i) = sqrt(9.81*AA/BB);
%         av_dy(i) = (av_c1*hull_seb(i)*ro*abs(dv(i)) + av_c2*dv(i)^2)/g/ro;
%         %av_dy(i) = av_c2*dv(i)^2/g/ro;
%         av_dy(i) = min(0.1,av_dy(i));
%     end
% end
% Ytr=Ytr+av_dy;

%% Harom pont atlaga
Ytr=max(Ytr,1e-3);
Q=zeros(size(csatorna.y)); A=Q;
% for i=1:length(csatorna.y)
%     A(i)=get_A(csatorna,csatorna.y(i));
%     Q(i)=Vtr(i)*A(i);
% end

Q = Vtr.*AA;

Dyj=Ytr(3:end)-Ytr(2:end-1); 
Dyb=Ytr(2:end-1)-Ytr(1:end-2);
Y_mean=(Ytr(3:end)+Ytr(2:end-1)+Ytr(1:end-2))/3; 

%DYtemp=Y_mean-Ytr(2:end-1);
DY=Y_mean-Ytr(2:end-1);
Dytot=Dyj.*Dyb;
mask=find(Dytot>0);
DY(mask) = 0;

Ytr(2:end-1)=Ytr(2:end-1)+DY;

% for i=1:length(csatorna.y)
%     Vtr(i)=Q(i)/get_A(csatorna,Ytr(i));
% end

Dvj=Vtr(3:end)-Vtr(2:end-1); 
Dvb=Vtr(2:end-1)-Vtr(1:end-2);
V_mean=(Vtr(3:end)+Vtr(2:end-1)+Vtr(1:end-2))/3; 

DV=V_mean-Vtr(2:end-1);
Dvtot=Dvj.*Dvb; mask=find(Dvtot>0); 
DV(mask)=0;
Vtr(2:end-1)=Vtr(2:end-1)+DV;



%% Spline
% dV=diff(Vtr(j,:))./diff(Xmesh); dVtreshold=0;
% if max(dV.^2) > dVtreshold
%     p=ones(1,N); pef=[inf inf];
%     Vtr(j,:) = simi_spline(Xmesh',Vtr(j,:)',p,pef,Xmesh');
%     Ytr(j,:) = simi_spline(Xmesh',Ytr(j,:)',p,pef,Xmesh');
% end
%     if Vtr(j,1)>Vtr(j,2)
%         Vtr(j,1)=Vtr(j,2);
%     end
%% adatok visszairasa es rajzolas

csatorna.y = Ytr; csatorna.v = Vtr;
csatorna.t = csatorna.t + dt;

if strcmp(csatorna.rajz,'rajz')
    
    V=trapz(Xmesh,A);
    
    figure(csatorna.tranziens_agelem_2csp.fignum)
    
    subplot(3,1,1)
    plot([0 L],[He Hv],'k',[0 L],[He+X Hv+X],'k',Xmesh,(Ytr+He)-s*Xmesh,'b-+');
    %plot(Xmesh,Ytr,'b-+');
    title([csatorna.tranziens_agelem_2csp.nev,...
        ', t = ',num2str(round(10*csatorna.t/60)/10),...
        ' min., Q_e = ',num2str(round(10*Q(1)*3600)/10),...
        ' m^3/h, Q_v = ',num2str(round(10*Q(end)*3600)/10),...
        ' m^3/h, V = ',num2str(round(10*V)/10),...
        ' m3, he = ',num2str(round(10*(Ytr(1)+He))/10),...
        ' m, hv = ',num2str(round(10*(Ytr(end)+Hv))/10),' m']);
    axis tight
    
    subplot(3,1,2)
    plot(Xmesh,Q*3600,'b+-'), ylabel('Q [m3/h]')
    
    subplot(3,1,3)
    plot(Xmesh,Vtr,'b+-'), ylabel('v [m/s]')
    
    %subplot(4,1,4)
    %plot(Xmesh,A,'b+-'), ylabel('A [m2]')
end

%% dt update
csatorna.dt=update_dt(csatorna);

%% ------------------------------------------------------------------------
%--------------------------------------------------------------------------
% functions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% spline
function yy = simi_spline(x,y,p,pf,xx)
%
% yy = simi_spline(x,y,p,pf,xx)
% keozelito spline
% az (x,y) pontokra illeszt a p k�zel�t� param�terrel
% pf - peremfeltetel
%      ha veges, akkor pf(1) = g'(x_1) �s pf(2) = g'(x_N)
%      ha inf, akkor g''(x_1) = 0 illetve g''(x_N) = 0
% visszateresi erteke yy = g(xx)
%
% Vaik Istvan, 2005.
%
N = length(x);
if length(y) ~= N || length(p) ~= N
    exit;
end;

% x szerint novekvo sorba rendezes
[x, s] = sort(x);
y = y(s);
p = p(s);


dx = x(2:N) - x(1:N-1);

A=zeros(N); B=zeros(N);
sparse(A);  sparse(B);

A = diag( dx/6, -1) + diag( [ dx(1)/3 ; ( dx(1:N-2) + dx(2:N-1) )/3; dx(N-1)/3 ] ) + diag( dx/6,1);
B = -diag( 1./dx, -1 ) + diag([ 1/dx(1); 1./dx(1:N-2) + 1./dx(2:N-1); 1/dx(N-1) ] ) - diag( 1./dx, 1);

b = [-pf(1); zeros( N-2,1); -pf(2)];

if pf(1) == inf
    A(1,1) = 1;
    A(1,2) = 0;
    B(1,1) = 0;
    B(1,2) = 0;
    b(1) = 0;
end;

if pf(2) == inf
    A(N,N) = 1;
    A(N,N-1) = 0;
    B(N,N) = 0;
    B(N,N-1) = 0;
    b(N) = 0;
end;

D=zeros(N);
sparse(D);
D = diag(p);

%%
sparse(b);
gdd = (A+B/D*B)\(b-B*y);
g = y + D\B*gdd;
%%
M = length(xx);
yy = zeros(M,1);
for k=1:M
    i = max( find( x < xx(k) ) );
    if isempty(i)
        if xx(k) == x(1)
            yy(k) = g(1);
        else
            continue;
        end;
    else
        if xx(k) > x(N)
            continue;
        else
            yy(k) = g(i)*( dx(i) - ( xx(k)-x(i) ) ) / dx(i);
            yy(k) = yy(k) + g(i+1)*( xx(k)-x(i) )/ dx(i);
            yy(k) = yy(k) + gdd(i)*( -(xx(k)-x(i))^3 + 3*(xx(k)-x(i))^2*dx(i) - 2*(xx(k)-x(i))*dx(i)^2 )/ 6 / dx(i);
            yy(k) = yy(k) + gdd(i+1)*( (xx(k)-x(i))^3 - (xx(k)-x(i))*dx(i)^2 )/6/dx(i);
        end;
    end;
end;


[xx, s] = sort(xx);
yy = yy(s);

%% normal
function fy=normal(csatorna,Q,s,n)
type=csatorna.tipus;
x=csatorna.dvB;
relax=0.8; eps=1e-6;

switch type
    case 'teglalap'
        type_old=1;
    case 'kor'
        type_old=2;
end

minusz=x*0.5;
while (Fnorm(minusz,x,type,Q,s,n)>0) || (dFnorm(minusz,x,type,Q,s,n)<0)
    minusz=minusz/2;
end

old=minusz; hiba=1;
while abs(hiba)>eps
    new=old-relax*Fnorm(old,x,type,Q,s,n)/dFnorm(old,x,type,Q,s,n);
    hiba=(new-old); hiba=Fnorm(old,x,type,Q,s,n);
    old=new;
    if ((old>x) || (dFnorm(old,x,type,Q,s,n)<0)) && (type_old==2)
        fy=[];
        break
    end
    fy(1)=new;
end

if ((~isempty(fy)) && (Fnorm(x,x,type,Q,s,n)<0)) && (type_old==2)
    plusz=new; minusz=x;
    while abs(plusz-minusz)>eps
        half=(plusz+minusz)/2;
        if Fnorm(half,x,type,Q,s,n)>0
            plusz=half;
        else
            minusz=half;
        end
    end
    fy(2)=(plusz+minusz)/2;
end
%% normal function
%function fy = Fnorm(csatorna,y,x,type,Q,s,n)
function fy = Fnorm(y,x,type,Q,s,n)

theta = acos(1-y/(x/2));
AA = (x/2)^2*(theta-sin(2*theta)/2);
KK = x*theta;
Rh = AA./KK;
% C  =get_Rh(csatorna,y)^(1/6)/n;
C  =Rh^(1/6)/n;
% fy = s - Q^2/get_A(csatorna,y)^2/C^2/get_Rh(csatorna,y);
fy = s - Q^2/AA^2/C^2/Rh;
%% normal function der
%function fy = dFnorm(csatorna,y,x,type,Q,s,n)
function fy = dFnorm(y,x,type,Q,s,n)

switch type
    case 'teglalap'
        dA=x; dK=2;
    case 'kor'
        r = x/2; theta=acos(1-y/r);
        dtheta=1/( r*(sqrt(1-(1-y/r)^2)) );
        dA=r^2*(dtheta-cos(2*theta)*dtheta);
        dK=2*r*dtheta;
end

theta = acos(1-y/(x/2));
AA = (x/2)^2*(theta-sin(2*theta)/2);
KK = x*theta;
Rh = AA/KK;
% 
% dRh=dA/get_K(csatorna,y) - get_A(csatorna,y)*dK/get_K(csatorna,y)^2;
% dC=dRh/(get_Rh(csatorna,y)^(5/6))/6/n;
% 
% C=get_Rh(csatorna,y)^(1/6)/n;
% 
% fy = 2*Q^2*dA/get_A(csatorna,y)^3/C^2/get_Rh(csatorna,y) + ...
%     2*Q^2*dC/get_A(csatorna,y)^2/C^3/get_Rh(csatorna,y) + ...
%     Q^2*dRh/get_A(csatorna,y)^2/C^2/get_Rh(csatorna,y)^2;

dRh=dA/KK - AA*dK/KK^2;
dC=dRh/(Rh^(5/6))/6/n;

C=Rh^(1/6)/n;

fy = 2*Q^2*dA/AA^3/C^2/Rh + ...
    2*Q^2*dC/AA^2/C^3/Rh + ...
    Q^2*dRh/AA^2/C^2/Rh^2;

%% kritikus
function fy=crit(csatorna,Q)
type=csatorna.tipus;
x=csatorna.dvB;
switch type
    case 'teglalap'
        type_old=1;
    case 'kor'
        type_old=2;
end

switch type_old
    case 1
        g=9.81;
        fy=(Q^2/g/x^2)^(1/3);
    case 2
        % Vezerlo parameterek
        eps=1e-6;
        
        % kiindulasi pozitiv ertek beallitasa
        plusz=x;
        t=1; p=Fkrit(x,x,type,Q);
        % kiindulasi negativ ertek meghatarozasa
        while p>0
            t=t+1;
            p=Fkrit(x/t,x,type,Q);
        end
        minusz=x/t;
        
        % iteracio
        while abs(plusz-minusz)>eps
            half=(plusz+minusz)/2;
            if Fkrit(half,x,type,Q)>0
                plusz=half;
            else
                minusz=half;
            end
        end
        fy=(plusz+minusz)/2;
end
%% kritikus function
%function fy = Fkrit(csatorna,y,x,type,Q)
function fy = Fkrit(y,x,type,Q)
g=9.81;

%fy = 1-Q^2*get_B(csatorna,y)/get_A(csatorna,y)^3/g;

theta = acos(1-y/(x/2));
AA = (x/2)^2*(theta-sin(2*theta)/2);
BB =x*sin(pi-theta);
fy = 1-Q^2*BB/AA^3/g;
%% ymin es ymax beallitasa rajzolashoz
function range = set_y_minmax(y)
if max(y)>0
    ymax=1.01*max(y);
else
    ymax=0.99*max(y);
end
if min(y)<0
    ymin=1.01*min(y);
else
    ymin=0.99*min(y);
end
range = [ymin ymax];
