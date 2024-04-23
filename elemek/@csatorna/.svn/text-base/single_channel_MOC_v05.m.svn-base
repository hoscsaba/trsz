function single_channel_MOC_v05
clc; clear all;
%% Csatorna adatok megadasa;
type=2; % 1-Vegtelen magas teglalap keresztmetszetu csatorna.
        % 2-Kor keresztmetszetu csatorna. (1% Hasitek)
X=0.5; He=217.00; Hv=215.00; n=0.013; L=800; s=(He-Hv)/L; g=9.81;

N=800; Xmesh=linspace(0,L,N);
%% Tranziens futtatas
tic

NT=4000; T=zeros(NT,1); h_min=0.015;
Yinit=0.25*ones(1,N); Vinit=0.0*ones(1,N);
Yinit=X*h_min*ones(1,N);

Ytr=zeros(NT,N); Ytr(1,:)=Yinit;
Vtr=zeros(NT,N); Vtr(1,:)=Vinit;

percentage=1;
for j=2:NT;
    if j>percentage*NT/100
        fprintf('%d\n', percentage);
        percentage=percentage+1;
    end
    
    % keresztmetszet jellemzoinek kiszamolasa
    Y_temp=Ytr(j-1,:); Y_temp( Y_temp>X )=X;
    r=X/2; theta=acos(1-Y_temp/r);
        AA=r^2*(theta-sin(2*theta)/2);
        BB=2*r*sin(pi-theta); BB( BB<0.01 )=0.01;
        KK=2*r*theta; Rhh=AA./KK;
        sp=sqrt(g*AA./BB);
%% Inner points
%     figure(2);
%     subplot(2,2,1); hold on;
    tP=zeros(1,N); xP=zeros(1,N);
    yP=zeros(1,N); vP=zeros(1,N);
    % Vectorized solver----------------------------------------------------
    yR=Ytr(j-1,3:N); vR=Vtr(j-1,3:N); yL=Ytr(j-1,1:N-2); vL=Vtr(j-1,1:N-2);
    aR=sp(3:N);      XR=Xmesh(3:N);   aL=sp(1:N-2);      XL=Xmesh(1:N-2);
    JR=n^2*abs(vR).*vR./Rhh(3:N).^(4/3); JL=n^2*abs(vL).*vL./Rhh(1:N-2).^(4/3);
    mR=1./(vR-aR); mL=1./(vL+aL);
    
    tP(2:N-1)=(-mR*T(j-1)+ mL*T(j-1)+mL.*mR.*(XL-XR))./(mL-mR);
    dt=tP(2:N-1)-T(j-1);
    xP(2:N-1)=(mL.*XL-mR.*XR)./(mL-mR);
    
    yP(2:N-1)=(aR*g.*yL + aL.*(aR.*(dt*g.*(-JL+JR)+vL-vR)+g*yR))./((aL+aR)*g);
    vP(2:N-1)=(aL.*(dt*g.*(-JL+s)+vL)+aR.*(dt*g.*(-JR+s)+vR)+g*(yL-yR))./(aL+aR);
    %----------------------------------------------------------------------
%     pause
%     for k=2:N-1
%         % datas of the left and right sides
%         yR=Ytr(j-1,k+1); vR=Vtr(j-1,k+1); yL=Ytr(j-1,k-1); vL=Vtr(j-1,k-1);
%         aR=sp(k+1);      XR=Xmesh(k+1);   aL=sp(k-1);      XL=Xmesh(k-1);
%         JR=n^2*abs(vR)*vR/Rhh(k+1)^(4/3); JL=n^2*abs(vL)*vL/Rhh(k-1)^(4/3);
%         mR=1/(vR-aR); mL=1/(vL+aL);
%         % location of point P
%         tP(k)=(-mR*T(j-1)+mL*(T(j-1)+mR*(XL-XR)))/(mL-mR); dt=tP(k)-T(j-1);
%         xP(k)=(mL*XL-mR*XR)/(mL-mR);
%         % properties of point P
%         yP(k)=(aR*g*yL+aL*(aR*(dt*g*(-JL+JR)+vL-vR)+g*yR))/((aL+aR)*g);
%         vP(k)=(aL*(dt*g*(-JL+s)+vL)+aR*(dt*g*(-JR+s)+vR)+g*(yL-yR))/(aL+aR);
%         % visualization of the characteristics points
% %         plot([XL xP(k)],[T(j-1),tP(k)],[XR xP(k)],[T(j-1) tP(k)]);
%     end
%     hold off;
%% Interpolation
% Vectorized interpolator--------------------------------------------------
% t_int=zeros(1,N);  yPint=zeros(1,N);  vPint=zeros(1,N);
% t_int1=zeros(1,N); yPint1=zeros(1,N); vPint1=zeros(1,N);
% t_int2=zeros(1,N); yPint2=zeros(1,N); vPint2=zeros(1,N);
% 
% yR=Ytr(j-1,3:N); vR=Vtr(j-1,3:N); yL=Ytr(j-1,1:N-2); vL=Vtr(j-1,1:N-2);
% 
% t_int1(2:N-1)=( (Xmesh(2:N-1)-Xmesh(1:N-2)).*tP(2:N-1)-(Xmesh(2:N-1)-Xmesh(1:N-2))*T(j-1) )./(xP(2:N-1)-Xmesh(1:N-2))+T(j-1);
% yPint1(2:N-1)=( (yP(2:N-1)-yL).*t_int1(2:N-1)- (yP(2:N-1)-yL)*T(j-1) )./(tP(2:N-1)-T(j-1))+yL;
% vPint1(2:N-1)=( (vP(2:N-1)-vL).*t_int1(2:N-1)- (vP(2:N-1)-vL)*T(j-1) )./(tP(2:N-1)-T(j-1))+vL;
% 
% t_int2(2:N-1)=( (Xmesh(2:N-1)-Xmesh(3:N)).*tP(2:N-1)-(Xmesh(2:N-1)-Xmesh(3:N))*T(j-1) )./(xP(2:N-1)-Xmesh(3:N))+T(j-1);
% yPint2(2:N-1)=( (yP(2:N-1)-yR).*t_int2(2:N-1)- (yP(2:N-1)-yR)*T(j-1) )./(tP(2:N-1)-T(j-1))+yR;
% vPint2(2:N-1)=( (vP(2:N-1)-vR).*t_int2(2:N-1)- (vP(2:N-1)-vR)*T(j-1) )./(tP(2:N-1)-T(j-1))+vR;
% 
% mask=find( t_int1<=t_int2 );
%     yPint(mask)=yPint1(mask); vPint(mask)=vPint1(mask); t_int(mask)=t_int1(mask);
% mask=find( t_int2<t_int1 );
%     yPint(mask)=yPint2(mask); vPint(mask)=vPint2(mask); t_int(mask)=t_int2(mask);
% 
% tmin=min(t_int(2:end-1)); T(j)=tmin; dt=T(j)-T(j-1);
%     
% Ytr(j,2:N-1)=(yPint(2:N-1)-Ytr(j-1,2:N-1))*(tmin-T(j-1))./(t_int(2:N-1)-T(j-1))+Ytr(j-1,2:N-1);
% Vtr(j,2:N-1)=(vPint(2:N-1)-Vtr(j-1,2:N-1))*(tmin-T(j-1))./(t_int(2:N-1)-T(j-1))+Vtr(j-1,2:N-1);
%--------------------------------------------------------------------------
t_int=zeros(1,N); yPint=zeros(1,N); vPint=zeros(1,N);
    for k=2:N-1
        yR=Ytr(j-1,k+1); vR=Vtr(j-1,k+1); yL=Ytr(j-1,k-1); vL=Vtr(j-1,k-1);
        if Xmesh(k)<xP(k)
            t_int(k)=(Xmesh(k)-Xmesh(k-1))*(tP(k)-T(j-1))/(xP(k)-Xmesh(k-1))+T(j-1);
            yPint(k)=(yP(k)-yL)*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+yL;
            vPint(k)=(vP(k)-vL)*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+vL;
        else
            t_int(k)=(Xmesh(k)-Xmesh(k+1))*(tP(k)-T(j-1))/(xP(k)-Xmesh(k+1))+T(j-1);
            yPint(k)=(yP(k)-yR)*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+yR;
            vPint(k)=(vP(k)-vR)*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+vR;
        end
    end
tmin=min(t_int(2:end-1)); T(j)=tmin; dt=T(j)-T(j-1);
    for k=2:N-1
        Ytr(j,k)=(yPint(k)-Ytr(j-1,k))*(tmin-T(j-1))/(t_int(k)-T(j-1))+Ytr(j-1,k);
        Vtr(j,k)=(vPint(k)-Vtr(j-1,k))*(tmin-T(j-1))/(t_int(k)-T(j-1))+Vtr(j-1,k);
    end
    % minimal water level
    Ytr( Ytr<X*h_min )=X*h_min;
    % visualization of the intarpolant
%     subplot(2,2,2); plot(xrend,yPrend,xrend,vPrend);
%% First boundary condition
    Q_beg=abs(Vtr(j-1,1)*AA(1)); yc_beg=crit(type,X,max([Q_beg 1e-6])); Qbegin(j)=Q_beg;
    y_akna_beg=1*TS6(T(j),1*Yinit(1),0,20*Yinit(1),1)-...
               1*TS6(T(j),0*Yinit(1),20,19*Yinit(1),1);
    Ytr(j,1)=max([y_akna_beg   Ytr(j-1,1)-g*dt^2/2   -Vtr(j-1,1)*yc_beg/abs(Vtr(j-1,1))   X*h_min]);
    Vtr(j,1)=Vtr(j,2);
%% End boundary condition
    Q_end=abs(Vtr(j-1,end)*AA(end)); Qend(j)=Q_end;
    yc_end=crit(type,X,max([Q_end 1e-6]));
    yn_end=normal(type,X,max([Q_end 1e-6]),s,n); if length(yn_end)==0, yn_end=2*yc_end; else yn_end=yn_end(1); end
    y_akna_end=1*TS6(T(j),Yinit(N),0,1.0*Yinit(N),10)+...
               0*TS1(T(j),3*Yinit(N));
    Ytr(j,end)=max([y_akna_end   Ytr(j-1,end)-g*dt^2/2   min([yc_end yn_end])   X*h_min]);
    Vtr(j,end)=Vtr(j,end-1);
%% Visualization
% subplot(2,2,3);
% plot([0 L],[He Hv],'k',[0 L],[He+X Hv+X],'k',Xmesh,(Ytr(j,:)+He)-s*Xmesh);
% subplot(2,2,4);
% plot(Xmesh,Vtr(j,:),Xmesh,sp);

% pause
% subplot(2,2,1);
% plot([0 0],[T(j) T(j)]);

end
toc

figure(3)
plot(T,Qbegin,T,Qend);
V1=trapz(T,Qbegin)
V2=trapz(T,Qend)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Függvények
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------






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