function single_channel_MOC_v04
clc; clear all;
%% Csatorna adatok megadasa;
type=2; % 1-Vegtelen magas teglalap keresztmetszetu csatorna.
        % 2-Kor keresztmetszetu csatorna. (1% Hasitek)
X=0.5; He=217.00; Hv=215.00; n=0.013; L=800; s=(He-Hv)/L; g=9.81;

N=600; Xmesh=linspace(0,L,N);
%% Tranziens futtatas
tic

NT=3000; T=zeros(NT,1); h_min=0.01;
Yinit=X*h_min*ones(1,N); Vinit=0.0*ones(1,N);

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
    figure(2);
%     subplot(2,2,1); hold on;
    tP=zeros(1,N); xP=zeros(1,N);
    yP=zeros(1,N); vP=zeros(1,N);
    for k=2:N-1
        % datas of the left and right sides
        yR=Ytr(j-1,k+1); vR=Vtr(j-1,k+1); yL=Ytr(j-1,k-1); vL=Vtr(j-1,k-1);
        aR=sp(k+1);      XR=Xmesh(k+1);   aL=sp(k-1);      XL=Xmesh(k-1);
        JR=n^2*abs(vR)*vR/Rhh(k+1)^(4/3); JL=n^2*abs(vL)*vL/Rhh(k-1)^(4/3);
        mR=1/(vR-aR); mL=1/(vL+aL);
        % location of point P
        tP(k)=(-mR*T(j-1)+mL*(T(j-1)+mR*(XL-XR)))/(mL-mR); dt=tP(k)-T(j-1);
        xP(k)=(mL*XL-mR*XR)/(mL-mR);
        % properties of point P
        yP(k)=(aR*g*yL+aL*(aR*(dt*g*(-JL+JR)+vL-vR)+g*yR))/((aL+aR)*g);
        vP(k)=(aL*(dt*g*(-JL+s)+vL)+aR*(dt*g*(-JR+s)+vR)+g*(yL-yR))/(aL+aR);
        % visualization of the characteristics points
%         plot([XL xP(k)],[T(j-1),tP(k)],[XR xP(k)],[T(j-1) tP(k)]);
    end
%     hold off;
%% Interpolation
tmin=min(tP(2:end-1)); T(j)=tmin; count=1; dt=T(j)-T(j-1);
xint=zeros(1,2*N); yPint=zeros(1,2*N); vPint=zeros(1,2*N);
xint(1)=Xmesh(1);
    for k=2:N-1
        yR=Ytr(j-1,k+1); vR=Vtr(j-1,k+1); yL=Ytr(j-1,k-1); vL=Vtr(j-1,k-1);
        xint_temp =(xP(k)-Xmesh(k-1))*(tmin-T(j-1))/(tP(k)-T(j-1))+Xmesh(k-1);
            if xint_temp>xint(count)
                count=count+1;
                xint(count)=xint_temp;
                yPint(count)=(yP(k)-yL)*(tmin-T(j-1))/(tP(k)-T(j-1))+yL;
                vPint(count)=(vP(k)-vL)*(tmin-T(j-1))/(tP(k)-T(j-1))+vL;
            end
        xint_temp =(xP(k)-Xmesh(k+1))*(tmin-T(j-1))/(tP(k)-T(j-1))+Xmesh(k+1);
            if xint_temp>xint(count)
                count=count+1;
                xint(count)=xint_temp;
                yPint(count)=(yP(k)-yR)*(tmin-T(j-1))/(tP(k)-T(j-1))+yR;
                vPint(count)=(vP(k)-vR)*(tmin-T(j-1))/(tP(k)-T(j-1))+vR;
            end
    end
    yPint(1)=yPint(2); vPint(1)=vPint(2);
        count=count+1;
        xint(count)=Xmesh(end); yPint(count)=yPint(count-1); vPint(count)=vPint(count-1);
    xint(count+1:end)=[]; yPint(count+1:end)=[]; vPint(count+1:end)=[];
    [xint index dummy]=unique(xint);
    yPint=yPint(index);
    vPint=vPint(index);

    [xrend,index]=sort(xint);
    yPrend=yPint(index); vPrend=vPint(index);
    
    yPnew = interp1(xrend,yPrend,Xmesh,'linear','extrap');
    vPnew = interp1(xrend,vPrend,Xmesh,'linear','extrap');
    
    Ytr(j,2:end-1)=yPnew(2:end-1); Vtr(j,2:end-1)=vPnew(2:end-1);
    % minimal water level
    Ytr( Ytr<X*h_min )=X*h_min;
    % visualization of the intarpolant
    subplot(2,2,2); plot(xrend,yPrend,xrend,vPrend);
%% First boundary condition
    Q_beg=abs(Vtr(j-1,1)*AA(1)); yc_beg=crit(type,X,max([Q_beg 1e-6])); %Qbegin(j)=Q_beg;
    y_akna_beg=1*TS6(T(j),Yinit(1),0,20*Yinit(1),1)-...
               0*TS6(T(j),0*Yinit(1),20,0*Yinit(1),1);
    Ytr(j,1)=max([y_akna_beg   Ytr(j-1,1)-g*dt^2/2   -Vtr(j-1,1)*yc_beg/abs(Vtr(j-1,1))   X*h_min]);
    Vtr(j,1)=Vtr(j,2);
%% End boundary condition
    Q_end=abs(Vtr(j-1,end)*AA(end)); %Qend(j)=Q_end;
    yc_end=crit(type,X,max([Q_end 1e-6]));
    yn_end=normal(type,X,max([Q_end 1e-6]),s,n); if length(yn_end)==0, yn_end=2*yc_end; else yn_end=yn_end(1); end
    y_akna_end=1*TS6(T(j),Yinit(N),0,1.0*Yinit(N),10);
    Ytr(j,end)=max([y_akna_end   Ytr(j-1,end)-g*dt^2/2   min([yc_end yn_end])   X*h_min]);
    Vtr(j,end)=Vtr(j,end-1);
%% Visualization
subplot(2,2,3);
plot([0 L],[He Hv],'k',[0 L],[He+X Hv+X],'k',Xmesh,(Ytr(j,:)+He)-s*Xmesh);
subplot(2,2,4);
plot(Xmesh,Vtr(j,:),Xmesh,sp);

% pause
subplot(2,2,1);
plot([0 0],[T(j) T(j)]);

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