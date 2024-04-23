function single_channel_forward_v03
clc; clear all;
%% Csatorna adatok megadasa;
type=2; % 1-Vegtelen magas teglalap keresztmetszetu csatorna.
        % 2-Kor keresztmetszetu csatorna. (1% Hasitek)
X=0.5; He=216.80; Hv=214.13; n=0.013; L=440; s=(He-Hv)/L; g=9.81;

N=400; Xmesh=linspace(0,L,N); dx=Xmesh(2)-Xmesh(1);
%% Tranziens futtatas
tic

Total_Time=200; dt=0.1; NT=Total_Time/dt; T=zeros(NT,1);
Yinit=0.25*ones(1,N); Vinit=0*ones(1,N);
Ytr=zeros(NT,N); Ytr(1,:)=Yinit;
Vtr=zeros(NT,N); Vtr(1,:)=Vinit;

percentage=1;
for j=2:NT;
    if j>percentage*NT/100
        fprintf('%d\n', percentage);
        percentage=percentage+1;
    end
    T(j)=T(j-1)+dt;
    Aeq=zeros(2*N,2*N); Beq=zeros(2*N,1); sparse(Aeq);
%% Eleje Perem
        % Kontinuitás
        Aeq(1,1)=1; Aeq(1,N+1)=0;
        Aeq(1,2)=0; Aeq(1,N+2)=0;
        Aeq(1,3)=0; Aeq(1,N+3)=0;
            % Eleje vízszint leellenõrzése kritikus kifolyásra
            Q_beg=abs(Vtr(j-1,1)*A(Ytr(j-1,1),X,type)); yc_beg=crit(type,X,max([Q_beg 1e-6]));
            y_akna_beg=TS6(T(j),Yinit(1),0,2.5*Yinit(1),1)+0*TS6(T(j),0,20,-3.0*Yinit(1),50);
            %--------------------------------------------------------------
                y_value_beg=max([y_akna_beg   Ytr(j-1,1)-g*dt^2/2   -Vtr(j-1,1)*yc_beg/abs(Vtr(j-1,1))   y_akna_beg   X*0.001]);
                %----------------------------------------------------------
                Beq(1)=y_value_beg;
        % Mozgás egyenlet
        Aeq(N+1,1)=(-3*g/dx/2);     Aeq(N+1,N+1)=(g*n^2*abs(Vtr(j-1,1))/Rh(Ytr(j-1,1),X,type)^(4/3) - 3*Vtr(j-1,1)/dx/2 + 1/dt)*1;
        Aeq(N+1,2)=2*g/dx;          Aeq(N+1,N+2)=2*Vtr(j-1,1)/dx;
        Aeq(N+1,3)=-g/dx/2;         Aeq(N+1,N+3)=-Vtr(j-1,1)/dx/2;
            Beq(N+1)=g*s+Vtr(j-1,1)/dt;
%% Kozbenso pontok
                for k=2:N-1
                        % Kontinuitás
                            sound_speed=g*A(Ytr(j-1,k),X,type)/B(Ytr(j-1,k),X,type);
                        Aeq(k,k-1)=-Vtr(j-1,k)/dx/2;    Aeq(k,N+k-1)=-sound_speed/dx/g/2;
                        Aeq(k,k)  =1/dt;                Aeq(k,N+k)  =0;
                        Aeq(k,k+1)=Vtr(j-1,k)/dx/2;     Aeq(k,N+k+1)=sound_speed/dx/g/2;
                            Beq(k)=Ytr(j-1,k)/dt;
                        % Mozgásegyenlet
                        Aeq(N+k,k-1)=-g/dx/2;   Aeq(N+k,N+k-1)=-Vtr(j-1,k)/dx/2;
                        Aeq(N+k,k)  =0;         Aeq(N+k,N+k)  =1/dt+g*n^2*abs(Vtr(j-1,k))/Rh(Ytr(j-1,k),X,type)^(4/3);
                        Aeq(N+k,k+1)=g/dx/2;    Aeq(N+k,N+k+1)=Vtr(j-1,k)/dx/2;
                            Beq(N+k)=g*s+Vtr(j-1,k)/dt;
                end
%% Vege Perem
        % Kontinuitás
        Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
        Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
        Aeq(N,N)  =1;   Aeq(N,N+N)  =0;
            % Vége vízszint leellenõrzése kritikus kifolyásra
            Q_end=abs(Vtr(j-1,end)*A(Ytr(j-1,end),X,type));
            yc_end=crit(type,X,max([Q_end 1e-6]));
            yn_end=normal(type,X,max([Q_end 1e-6]),s,n); if length(yn_end)==0, yn_end=2*yc_end; else yn_end=yn_end(1); end
            y_akna_end=0*TS6(T(j),0*Yinit(N),0,0.0*Yinit(N),50);
            %--------------------------------------------------------------
                y_value_end=max([y_akna_end   Ytr(j-1,end)-g*dt^2/2   min([yc_end yn_end])   y_akna_end   X*0.001]);
                %----------------------------------------------------------
                Beq(N)=y_value_end;
        % Mozgás egyenlet
        Aeq(2*N,N-2)=g/dx/2;        Aeq(2*N,N+N-2)=Vtr(j-1,N)/dx/2;
        Aeq(2*N,N-1)=-2*g/dx;       Aeq(2*N,N+N-1)=-2*Vtr(j-1,N)/dx;
        Aeq(2*N,N)  =(3*g/dx/2);    Aeq(2*N,N+N)  =g*n^2*abs(Vtr(j-1,N))/Rh(Ytr(j-1,N),X,type)^(4/3) + 3*Vtr(j-1,N)/dx/2 + 1/dt;
            Beq(2*N)=g*s+Vtr(j-1,N)/dt;
%% Megoldas
        [mo flag]=bicgstab(Aeq,Beq,1e-6,[],[],[],[Ytr(j-1,:), Vtr(j-1,:)]');
        Ytr(j,:)=mo(1:N)'; Vtr(j,:)=mo(N+1:2*N)';
%% Simitas
Dyj=Ytr(j,3:end)-Ytr(j,2:end-1); Dyb=Ytr(j,2:end-1)-Ytr(j,1:end-2);
Y_mean=(Ytr(j,3:end)+Ytr(j,2:end-1)+Ytr(j,1:end-2))/3; DY=Y_mean-Ytr(j,2:end-1);
Dytot=Dyj.*Dyb; mask=find(Dytot>0); DY(mask)=0;
Ytr(j,2:end-1)=Ytr(j,2:end-1)+DY;

Dvj=Vtr(j,3:end)-Vtr(j,2:end-1); Dvb=Vtr(j,2:end-1)-Vtr(j,1:end-2);
V_mean=(Vtr(j,3:end)+Vtr(j,2:end-1)+Vtr(j,1:end-2))/3; DV=V_mean-Vtr(j,2:end-1);
Dvtot=Dvj.*Dvb; mask=find(Dvtot>0); DV(mask)=0;
Vtr(j,2:end-1)=Vtr(j,2:end-1)+DV;
end
toc
%% Rajzolas
pause
for k=1:NT
    figure(5)

    Q_end=abs(Vtr(k,end)*A(Ytr(k,end),X,type));
%     Q_end=abs(Vtr(k,1)*A(Ytr(k,1),X,type));
    % normál szint számítás
    Yn=normal(type,X,max([Q_end 1e-6]),s,n);
            if length(Yn)==0
                Yn=2*yc_end;
            else
                Yn=Yn(1);
            end
    % kritikus szint számítás
    yc_end=crit(type,X,max([Q_end 1e-6]));
    
    subplot(2,1,1)
    plot([0 L],[He Hv],'k',[0 L],[He+X Hv+X],'k',[0 L],[He+Yn Hv+Yn],'g',[0 L],[He+yc_end Hv+yc_end],'r',Xmesh,(Ytr(k,:)+He)-s*Xmesh);
    xlim([0 L]); ylim([min((Ytr(k,:)+He)-s*Xmesh)-max(Ytr(:,end)), 1.001*max((Ytr(k,:)+He)-s*Xmesh)]);
    
    subplot(2,1,2)
    plot(Xmesh,Vtr(k,:));
    xlim([0 L]); ylim([1.1*min(min(Vtr)), 1.1*max(max(Vtr))]);
    
            tit = title( strcat(num2str(ceil(T(k))),' s') );
            set(tit,'Fontsize',26,'Fontweight','demi');

%     pause
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Függvények
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------






function dy = Ly(x,y,type,X,Q,s,n,Yc);
% A differenciál egyenletet leíró függvény. Ha a szinte eléri a D szintet
% akkor a jobboldal konstans lesz, ami a teltszelvényû áramlást jellemzi.

g=9.81;
dy = (s - Q^2/A(y,X,type).^2/(Rh(y,X,type)^(1/6)/n)^2/Rh(y,X,type)) / (1-Q^2*B(y,X,type)./A(y,X,type).^3/g);

%--------------------------------------------------------------------------
% Timeseries
%--------------------------------------------------------------------------
function fy=TS1(t,xo)
fy=xo;

function fy=TS2(t,xo)
fy=xo*exp(-t/100);

function fy=TS3(t,xo,x1)
fy=(xo-x1)*exp(-t/100)+x1;

function fy=TS4(t,xo,x1,T)
fy=(x1-xo)/T*t+xo;

function fy=TS5(t,xo,to,A,T)
if (t<to) || (t>to+T)
    fy=xo;
else
    fy=xo+A*(1-cos(2*pi*(t-to)/T));
end

function fy=TS6(t,xo,to,x1,T)
if t<to
    fy=xo;
else
    fy=(xo-x1)*exp(-(t-to)/T)+x1;
end