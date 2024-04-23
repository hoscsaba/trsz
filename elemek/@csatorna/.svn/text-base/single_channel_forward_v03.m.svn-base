function single_channel_forward_v03
clc; clear all;
%--------------------------------------------------------------------------
% Csatorna adatok megadasa;
%--------------------------------------------------------------------------
type=2; % 1-Vegtelen magas teglalap keresztmetszetu csatorna.
        % 2-Kor keresztmetszetu csatorna. (1% Hasitek)
X=1; He=4; Hv=3.5; n=0.013; L=100; s=(He-Hv)/L;
Q=3; g=9.81;

N=41;
Xmesh=linspace(0,L,N); dx=Xmesh(2)-Xmesh(1);
%--------------------------------------------------------------------------
% Tranziens futtatások
%--------------------------------------------------------------------------
Total_Time=200; dt=0.2;
NT=Total_Time/dt;
T=zeros(NT,1);
Yinit=0.5*ones(1,N); Vinit=0*ones(1,N);
Ytr=zeros(NT,N); Ytr(1,:)=Yinit;
Vtr=zeros(NT,N); Vtr(1,:)=Vinit;

percentage=1;
for j=2:NT;
    if j>percentage*NT/100
        fprintf('%d\n', percentage);
        percentage=percentage+1;
    end
    T(j)=T(j-1)+dt;
    Aeq=zeros(2*N,2*N); Beq=zeros(2*N,1);
        % Eleje Perem------------------------------------------------------
        % Kontinuitás
        Aeq(1,1)=1; Aeq(1,N+1)=0;
        Aeq(1,2)=0; Aeq(1,N+2)=0;
        Aeq(1,3)=0; Aeq(1,N+3)=0;
            % Eleje vízszint leellenõrzése kritikus kifolyásra
            yc_beg=Vtr(j-1,1)^2/g; y_akna_beg=TS6(T(j),1*Yinit(1),0,0*Yinit(1),50);
            if (yc_beg>y_akna_beg) && (Vtr(j-1,1)<0)
                y_value_beg=yc_beg;
            else
                y_value_beg=y_akna_beg;
            end
            if y_value_beg<X*0.01;
                y_value_beg=X*0.01;
            end
            Beq(1)=y_value_beg;
        % Mozgás egyenlet
        Aeq(N+1,1)=(-3*g/dx/2);     Aeq(N+1,N+1)=(g*n^2*abs(Vtr(j-1,1))/Rh(Ytr(j-1,1),X,type)^(4/3) - 3*Vtr(j-1,1)/dx/2 + 1/dt)*1;
        Aeq(N+1,2)=2*g/dx;          Aeq(N+1,N+2)=2*Vtr(j-1,1)/dx;
        Aeq(N+1,3)=-g/dx/2;         Aeq(N+1,N+3)=-Vtr(j-1,1)/dx/2;
            Beq(N+1)=g*s+Vtr(j-1,1)/dt;
                % Közbensõ pontok------------------------------------------
                for k=2:N-1
                    % Kontinuitás
                        if B(Ytr(j-1,k),X,type)==0
                            sound_speed=g*A(Ytr(j-1,k),X,type)/1e-10;
                        else
                            sound_speed=g*A(Ytr(j-1,k),X,type)/B(Ytr(j-1,k),X,type);
                        end;
                        sound_bound=600;
                        if sound_speed>sound_bound
                            sound_speed=sound_bound;
                        end
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
        % Vége
        % Perem%-----------------------------------------------------------
        % Kontinuitás
        Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
        Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
        Aeq(N,N)  =1;   Aeq(N,N+N)  =0;
            % Vége vízszint leellenõrzése kritikus kifolyásra
            yc_end=Vtr(j-1,end)^2/g; y_akna_end=TS6(T(j),1*Yinit(N),5,4*Yinit(N),50);
            if (yc_end>y_akna_end) && (Vtr(j-1,end)>0)
                y_value_end=yc_end;
            else
                y_value_end=y_akna_end;
            end
            if y_value_end<X*0.01;
                y_value_end=X*0.01;
            end
            Beq(N)=y_value_end;
        % Mozgás egyenlet
        Aeq(2*N,N-2)=g/dx/2;        Aeq(2*N,N+N-2)=Vtr(j-1,N)/dx/2;
        Aeq(2*N,N-1)=-2*g/dx;       Aeq(2*N,N+N-1)=-2*Vtr(j-1,N)/dx;
        Aeq(2*N,N)  =(3*g/dx/2);    Aeq(2*N,N+N)  =g*n^2*abs(Vtr(j-1,N))/Rh(Ytr(j-1,N),X,type)^(4/3) + 3*Vtr(j-1,N)/dx/2 + 1/dt;
            Beq(2*N)=g*s+Vtr(j-1,N)/dt;
        
        % Megoldás
        Aeq;
        mo=Aeq\Beq; Ytr(j,:)=mo(1:N)'; Vtr(j,:)=mo(N+1:2*N)';
        
        % Simítás
                dV=diff(Vtr(j,:))./diff(Xmesh); dVtreshold=0;
                if max(dV.^2) > dVtreshold
                    p=ones(1,N); pf=[inf inf];
                    Vtr(j,:) = simi_spline(Xmesh',Vtr(j,:)',p,pf,Xmesh');
                    Ytr(j,:) = simi_spline(Xmesh',Ytr(j,:)',p,pf,Xmesh');
                end
end

for k=1:NT
    figure(4)
    

    subplot(2,1,1)
    plot([0 L],[He Hv],[0 L],[He+X Hv+X],Xmesh,(Ytr(k,:)+He)-s*Xmesh);
    xlim([0 L]); ylim([min((Ytr(k,:)+He)-s*Xmesh)-max(Ytr(:,end)), 1.1*max((Ytr(k,:)+He)-s*Xmesh)]);
    
    subplot(2,1,2)
    plot(Xmesh,Vtr(k,:));
    xlim([0 L]); ylim([1.1*min(min(Vtr)), 1.1*max(max(Vtr))]);
    
            tit = title( strcat(num2str(ceil(T(k))),' s') );
            set(tit,'Fontsize',26,'Fontweight','demi');

    pause
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