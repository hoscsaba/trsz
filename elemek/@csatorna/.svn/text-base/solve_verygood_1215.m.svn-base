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

%% Csatorna adatok megadasa;
%type=2; % 1-Vegtelen magas teglalap keresztmetszetu csatorna.
% 2-Kor keresztmetszetu csatorna. (1% Hasitek)
%X=0.5; He=217.00; Hv=215.00; n=1*0.013; L=250; s=(He-Hv)/L; g=9.81;

%N=800;

Xmesh=linspace(0,L,N);
%% Tranziens futtatas

% NT=1000;
% T=zeros(NT,1);
% Yinit=0.01*ones(1,N); Vinit=0.2*ones(1,N);
% Ytr=zeros(NT,N);
% Ytr=Yinit;
% Vtr=zeros(NT,N);
% Vtr=Vinit;

% dt = csatorna.dt;
t = csatorna.t;
Ytr = zeros(1,N); Vtr = zeros(1,N);
Ytr(1:N) = csatorna.y; Vtr(1:N) = csatorna.v;

%percentage=1;
%for j=2:NT;
%    if j>percentage*NT/100
%        fprintf('%d\n', percentage);
%        percentage=percentage+1;
%    end

% keresztmetszet jellemzoinek kiszamolasa
Y_temp=Ytr(:); Y_temp( Y_temp>X )=X;
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
for k=2:N-1
    % data of the left and right sides
    yR=Ytr(k+1); vR=Vtr(k+1); yL=Ytr(k-1); vL=Vtr(k-1);
    aR=sp(k+1);      XR=Xmesh(k+1);   aL=sp(k-1);      XL=Xmesh(k-1);
    JR=n^2*abs(vR)*vR/Rhh(k+1)^(4/3); JL=n^2*abs(vL)*vL/Rhh(k-1)^(4/3);
    mR=1/(vR-aR); mL=1/(vL+aL);
    % location of point P
    tP(k)=(-mR*t+mL*(t+mR*(XL-XR)))/(mL-mR);
    dt=tP(k)-t;
    xP(k)=(mL*XL-mR*XR)/(mL-mR);
    % properties of point P
    yP(k)=(aR*g*yL+aL*(aR*(dt*g*(-JL+JR)+vL-vR)+g*yR))/((aL+aR)*g);
    vP(k)=(aL*(dt*g*(-JL+s)+vL)+aR*(dt*g*(-JR+s)+vR)+g*(yL-yR))/(aL+aR);
    if dt < 0
        error('dt<0, a karakterisztikak a multban metszik egymast')
        tP(k)
        tP(k) = tP(k-1);
        tP(k)
        xP(k) = xP(k-1);
        yP(k) = yP(k-1);
        vP(k) = vP(k-1);
    end
end


%% Interpolation
tmin=min(tP(2:end-1)); T=tmin; count=1; dt=T-t;

if (T-t) < 0
   tP-t
   error('hiba') 
end

xint=zeros(1,2*N); yPint=zeros(1,2*N); vPint=zeros(1,2*N);
for k=2:N-1
    yR=Ytr(k+1); vR=Vtr(k+1); yL=Ytr(k-1); vL=Vtr(k-1);
    count=count+1;
    xint(count) =(xP(k)-Xmesh(k-1))*(tmin-t)/(tP(k)-t)+Xmesh(k-1);
    yPint(count)=(yP(k)-yL)*(tmin-t)/(tP(k)-t)+yL;
    vPint(count)=(vP(k)-vL)*(tmin-t)/(tP(k)-t)+vL;
    if tP(k)-t>1.001*dt
        count=count+1;
        xint(count) =(xP(k)-Xmesh(k+1))*(tmin-t)/(tP(k)-t)+Xmesh(k+1);
        yPint(count)=(yP(k)-yR)*(tmin-t)/(tP(k)-t)+yR;
        vPint(count)=(vP(k)-vR)*(tmin-t)/(tP(k)-t)+vR;
    end
end

xint(1)=Xmesh(1); yPint(1)=yPint(2); vPint(1)=vPint(2);

xint(count+1:end)=[]; yPint(count+1:end)=[]; vPint(count+1:end)=[];
[xint index dummy]=unique(xint);
yPint=yPint(index);
vPint=vPint(index);


[xrend,index]=sort(xint);
yPrend = yPint(index);
vPrend = vPint(index);

yPnew = interp1(xrend,yPrend,Xmesh,'linear','extrap');
vPnew = interp1(xrend,vPrend,Xmesh,'linear','extrap');

Ytr(2:end-1) = yPnew(2:end-1);
Vtr(2:end-1) = vPnew(2:end-1);
% minimal water level
Ytr( Ytr<X*0.001 )=X*0.001;
% visualization of the intarpolant
%     subplot(2,2,2); plot(xrend,yPrend,xrend,vPrend);
%% First boundary condition
Q_beg=abs(Vtr(1)*AA(1)); yc_beg=crit(csatorna,type,X,max([Q_beg 1e-6])); %Qbegin(j)=Q_beg;
%y_akna_beg=1*TS6(T,Yinit(1),0,20.0*Yinit(1),10)-...
%    1*TS6(T,0*Yinit(1),20,19.0*Yinit(1),10);

Ytr(1)=max([y_akna_beg   Ytr(1)-g*dt^2/2   -Vtr(1)*yc_beg/abs(Vtr(1))   X*0.001]);
Vtr(1)=Vtr(2);
%% End boundary condition
Q_end=abs(Vtr(end)*AA(end)); %Qend(j)=Q_end;
yc_end=crit(csatorna,type,X,max([Q_end 1e-6]));
yn_end=normal(csatorna,max([Q_end 1e-6]),s,n); if length(yn_end)==0, yn_end=2*yc_end; else yn_end=yn_end(1); end
%y_akna_end=1*TS6(T,Yinit(N),0,1.0*Yinit(N),10);
Ytr(end)=max([y_akna_end   Ytr(end)-g*dt^2/2   min([yc_end yn_end])   X*0.001]);
Vtr(end)=Vtr(end-1);

if T < csatorna.t
    error('negativ idolepes!')
end

csatorna.t = T;
csatorna.y = Ytr;
csatorna.v = Vtr;

%% Visualization
% subplot(2,2,3);
% plot([0 L],[He Hv],'k',[0 L],[He+X Hv+X],'k',Xmesh,(Ytr(j,:)+He)-s*Xmesh);
% subplot(2,2,4);
% plot(Xmesh,Vtr(j,:),Xmesh,sp);

% pause
% subplot(2,2,1);
% plot([0 0],[T T]);

%end

% figure(3)
% plot(T,Qbegin,T,Qend);
% V1=trapz(T,Qbegin)
% V2=trapz(T,Qend)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Függvények
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% kritikus
function fy=crit(csatorna,type,X,Q)

x = X;

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
        t=1; p=Fkrit(csatorna,x,x,type,Q);
        % kiindulasi negativ ertek meghatarozasa
        while p>0
            t=t+1;
            p=Fkrit(csatorna,x/t,x,type,Q);
        end
        minusz=x/t;
        
        % iteracio
        while abs(plusz-minusz)>eps
            half=(plusz+minusz)/2;
            if Fkrit(csatorna,half,x,type,Q)>0
                plusz=half;
            else
                minusz=half;
            end
        end
        fy=(plusz+minusz)/2;
end

%% kritikus function
function fy = Fkrit(csatorna,y,x,type,Q)
g=9.81;
fy = 1-Q^2*get_B(csatorna,y)/get_A(csatorna,y)^3/g;

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
while (Fnorm(csatorna,minusz,x,type,Q,s,n)>0) || (dFnorm(csatorna,minusz,x,type,Q,s,n)<0)
    minusz=minusz/2;
end

old=minusz; hiba=1;
while abs(hiba)>eps
    new=old-relax*Fnorm(csatorna,old,x,type,Q,s,n)/dFnorm(csatorna,old,x,type,Q,s,n);
    hiba=(new-old); hiba=Fnorm(csatorna,old,x,type,Q,s,n);
    old=new;
    if ((old>x) || (dFnorm(csatorna,old,x,type,Q,s,n)<0)) && (type_old==2)
        fy=[];
        break
    end
    fy(1)=new;
end

if ((~isempty(fy)) && (Fnorm(csatorna,x,x,type,Q,s,n)<0)) && (type_old==2)
    plusz=new; minusz=x;
    while abs(plusz-minusz)>eps
        half=(plusz+minusz)/2;
        if Fnorm(csatorna,half,x,type,Q,s,n)>0
            plusz=half;
        else
            minusz=half;
        end
    end
    fy(2)=(plusz+minusz)/2;
end
%% normal function
function fy = Fnorm(csatorna,y,x,type,Q,s,n)

C  =get_Rh(csatorna,y)^(1/6)/n;
fy = s - Q^2/get_A(csatorna,y)^2/C^2/get_Rh(csatorna,y);

%% normal function der
function fy = dFnorm(csatorna,y,x,type,Q,s,n)

switch type
    case 'teglalap'
        dA=x; dK=2;
    case 'kor'
        r = x/2; theta=acos(1-y/r);
        dtheta=1/( r*(sqrt(1-(1-y/r)^2)) );
        dA=r^2*(dtheta-cos(2*theta)*dtheta);
        dK=2*r*dtheta;
end

dRh=dA/get_K(csatorna,y) - get_A(csatorna,y)*dK/get_K(csatorna,y)^2;
dC=dRh/(get_Rh(csatorna,y)^(5/6))/6/n;

C=get_Rh(csatorna,y)^(1/6)/n;

fy = 2*Q^2*dA/get_A(csatorna,y)^3/C^2/get_Rh(csatorna,y) + ...
    2*Q^2*dC/get_A(csatorna,y)^2/C^3/get_Rh(csatorna,y) + ...
    Q^2*dRh/get_A(csatorna,y)^2/C^2/get_Rh(csatorna,y)^2;

%%

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