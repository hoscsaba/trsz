function single_channel_MOC_v06
clc; clear all;
%% Csatorna adatok megadasa;
type=2; % 1-Vegtelen magas teglalap keresztmetszetu csatorna.
        % 2-Kor keresztmetszetu csatorna. (1% Hasitek)
X=0.2; He=217.00; Hv=215.00; n=0.013; L=100; s=(He-Hv)/L; g=9.81;

N=500; Xmesh=linspace(0,L,N);
%% Tranziens futtatas
tic

NT=1000; T=zeros(NT,1); h_min=0.015;
Yinit=0.1*ones(1,N); Vinit=0.0*ones(1,N);
Yinit=1*X*h_min*ones(1,N);

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
%     for k=2:N-1
%         plot([Xmesh(k-1) xP(k)],[T(j-1) tP(k)],[Xmesh(k+1) xP(k)],[T(j-1) tP(k)]);
%     end
%     hold off;
%% Interpolation
t_int=zeros(1,N); yPint=zeros(1,N); vPint=zeros(1,N);
    for k=2:N-1
        if Xmesh(k)<xP(k)
            t_int(k)=(Xmesh(k)-Xmesh(k-1))*(tP(k)-T(j-1))/(xP(k)-Xmesh(k-1))+T(j-1);
            yPint(k)=(yP(k)-yL(k-1))*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+yL(k-1);
            vPint(k)=(vP(k)-vL(k-1))*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+vL(k-1);
        else
            t_int(k)=(Xmesh(k)-Xmesh(k+1))*(tP(k)-T(j-1))/(xP(k)-Xmesh(k+1))+T(j-1);
            yPint(k)=(yP(k)-yR(k-1))*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+yR(k-1);
            vPint(k)=(vP(k)-vR(k-1))*(t_int(k)-T(j-1))/(tP(k)-T(j-1))+vR(k-1);
        end
    end
    yPint(N)=yP(N-1); vPint(N)=vP(N-1); t_int(N)=tP(N-1);
    yPint(1)=yP(2);   vPint(1)=vP(2);   t_int(1)=tP(2);
%--------------------------------------------------------------------------
tmin=min(t_int(1:end)); T(j)=tmin; dt=T(j)-T(j-1);
    for k=1:N
        Ytr(j,k)=(yPint(k)-Ytr(j-1,k))*(tmin-T(j-1))/(t_int(k)-T(j-1))+Ytr(j-1,k);
        if Ytr(j,k)<X*h_min, Ytr(j,k)=X*h_min; end
        Vtr(j,k)=(vPint(k)-Vtr(j-1,k))*(tmin-T(j-1))/(t_int(k)-T(j-1))+Vtr(j-1,k);
    end
%% First boundary condition
yc_beg=crit_new(X,Vtr(j-1,1));
    y_akna_beg=1*TS6(T(j),1*Yinit(1),0,20*Yinit(1),1)-...
               1*TS6(T(j),0*Yinit(1),20,19*Yinit(1),1);
    Ytr(j,1)=max([y_akna_beg   Ytr(j-1,1)-g*dt^2/2   -Vtr(j-1,1)*yc_beg/abs(Vtr(j-1,1))   X*h_min]);
%% Recheck the end boundary water level
yc_end=crit_new(X,Vtr(j-1,end)); yn_end=norm_new(X,Vtr(j-1,end),n,s);
    y_akna_end=1*TS6(T(j),Yinit(N),0,1*Yinit(N),10)+...
               0*TS1(T(j),3*Yinit(N));
    Ytr(j,end)=max([y_akna_end   Ytr(j-1,end)-g*dt^2/2   min([yc_end yn_end])   X*h_min]);
%% Visualization
% subplot(2,2,3);
% plot([0 L],[He Hv],'k',[0 L],[He+X Hv+X],'k',Xmesh,(Ytr(j,:)+He)-s*Xmesh,[0 L],[He+yc_end Hv+yc_end],'r',[0 L],[He+yn_end Hv+yn_end],'g');
% subplot(2,2,4);
% plot(Xmesh,Vtr(j,:),Xmesh,sp);

% pause
% subplot(2,2,1);
% plot([0 0],[T(j) T(j)]);

end
toc

% figure(3)
% plot(T,Qbegin,T,Qend);
% V1=trapz(T,Qbegin)
% V2=trapz(T,Qend)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Függvények
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function fy=norm_new(X,v,n,s)

r=X/2; eps=0.01;
Rh_r_norm_end=(( v^2*n^2/s )^(3/4))/r;

theta=eps:eps:2.247;
Rh_r=(theta-sin(2*theta)/2)./theta/2; theta(1)=0; Rh_r(1)=0; Rh_r_max=0.6086;
if Rh_r_norm_end>=Rh_r_max
    fy=2*X;
else
    m1=find( Rh_r<=Rh_r_norm_end ); m2=find( Rh_r>Rh_r_norm_end );
    theta_norm = (theta(m2(1))-theta(m1(end)))*(Rh_r_norm_end-Rh_r(m1(end)))/(Rh_r(m2(1))-Rh_r(m1(end)))+theta(m1(end));
    
%     theta_norm = interp1(Rh_r,theta,Rh_r_norm_end);
    fy=r*(1-cos(theta_norm));
end

function fy=crit_new(X,v)

r=X/2; eps=0.01;
AK_end=v^2/r/9.81;

theta=0:eps:pi-eps;
AK=(theta-sin(2*theta)/2)./sin(pi-theta)/2;
theta(end)=pi; AK(end)=1000;

m1=find( AK<=AK_end ); m2=find( AK>AK_end );
theta_crit = (theta(m2(1))-theta(m1(end)))*(AK_end-AK(m1(end)))/(AK(m2(1))-AK(m1(end)))+theta(m1(end));

% theta_crit = interp1(AK,theta,AK_end);
fy=r*(1-cos(theta_crit));

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