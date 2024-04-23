function csatorna = solve(csatorna,pf)

bct1=pf{1}{1};
bcv1=pf{1}{2};
bct2=pf{2}{1};
bcv2=pf{2}{2};

fprintf('\n\n Csatorna lepes:  PFE tipus:  %s  ertek: %g    PFV tipus:  %s  ertek: %g',bct1,bcv1,bct2,bcv2);

ro = csatorna.tranziens_agelem_2csp.ro;
dt = csatorna.dt;
N  = length(csatorna.v);

dh=(csatorna.ze+bcv1/1000/9.81)-(csatorna.zv+bcv2/1000/9.81);
dzeta=0.1;
A=1;
mp=sqrt(abs(dh)/dzeta*2*9.81*A^2*ro^2);
if dh<0
    mp=-mp;
end

AA=get_A(csatorna,mean(csatorna.y),csatorna.dvB,csatorna.tipus);
csatorna.t=csatorna.t+csatorna.dt;
mp;
%pause
% Aeq=zeros(2*N,2*N); 
% Beq=zeros(2*N,1);
% 
% Vtr(1,:)=csatorna.y
% pause
% 
% % Eleje Perem
% % Kontinuitas
% Aeq(1,1)=1; Aeq(1,N+1)=0;
%         Aeq(1,2)=0; Aeq(1,N+2)=0;
%         Aeq(1,3)=0; Aeq(1,N+3)=0;
%             % Eleje vizszint leellenorzese kritikus kifolyasra
%             yc_beg=Vtr(j-1,1)^2/g; y_akna_beg=TS1(T(j),1*Yinit(1));
%             if (yc_beg>y_akna_beg) && (Vtr(j-1,1)<0)
%                 y_value_beg=yc_beg;
%             else
%                 y_value_beg=y_akna_beg;
%             end
%             if y_value_beg<X*0.001;
%                 y_value_beg=X*0.001;
%             end
%             Beq(1)=y_value_beg;
%         % Mozg�s egyenlet
%         Aeq(N+1,1)=(-3*g/dx/2);   Aeq(N+1,N+1)=(g*n^2*abs(Vtr(j-1,1))/Rh(Ytr(j-1,1),X,type)^(4/3) - 3*Vtr(j-1,1)/dx/2 + 1/dt)*1;
%         Aeq(N+1,2)=2*g/dx;          Aeq(N+1,N+2)=2*Vtr(j-1,1)/dx;
%         Aeq(N+1,3)=-g/dx/2;         Aeq(N+1,N+3)=-Vtr(j-1,1)/dx/2;
%             Beq(N+1)=g*s+Vtr(j-1,1)/dt;
%                 for k=2:N-1
%                     % Kontinuit�s
%                     Aeq(k,k-1)=-Vtr(j-1,k)/dx/2;    Aeq(k,N+k-1)=-A(Ytr(j-1,k),X,type)/B(Ytr(j-1,k),X,type)/dx/2;
%                     Aeq(k,k)  =1/dt;                Aeq(k,N+k)  =0;
%                     Aeq(k,k+1)=Vtr(j-1,k)/dx/2;     Aeq(k,N+k+1)=A(Ytr(j-1,k),X,type)/B(Ytr(j-1,k),X,type)/dx/2;
%                         Beq(k)=Ytr(j-1,k)/dt;
%                     % Mozg�segyenlet
%                     Aeq(N+k,k-1)=-g/dx/2;   Aeq(N+k,N+k-1)=-Vtr(j-1,k)/dx/2;
%                     Aeq(N+k,k)  =0;         Aeq(N+k,N+k)  =1/dt+g*n^2*abs(Vtr(j-1,k))/Rh(Ytr(j-1,k),X,type)^(4/3);
%                     Aeq(N+k,k+1)=g/dx/2;    Aeq(N+k,N+k+1)=Vtr(j-1,k)/dx/2;
%                         Beq(N+k)=g*s+Vtr(j-1,k)/dt;
%                 end
%         % V�ge Perem
%         % Kontinuit�s
%         Aeq(N,N-2)=0;   Aeq(N,N+N-2)=0;
%         Aeq(N,N-1)=0;   Aeq(N,N+N-1)=0;
%         Aeq(N,N)  =1;   Aeq(N,N+N)  =0;
%             % V�ge v�zszint leellen�rz�se kritikus kifoly�sra
%             yc_end=Vtr(j-1,end)^2/g; y_akna_end=TS3(T(j),1*Yinit(N),0);
%             if (yc_end>y_akna_end) && (Vtr(j-1,end)>0)
%                 y_value_end=yc_end;
%             else
%                 y_value_end=y_akna_end;
%             end
%             if y_value_end<X*0.001;
%                 y_value_end=X*0.001;
%             end
%             Beq(N)=y_value_end;
%         % Mozg�s egyenlet
%         Aeq(2*N,N-2)=g/dx/2;        Aeq(2*N,N+N-2)=Vtr(j-1,N)/dx/2;
%         Aeq(2*N,N-1)=-2*g/dx;       Aeq(2*N,N+N-1)=-2*Vtr(j-1,N)/dx;
%         Aeq(2*N,N)  =(3*g/dx/2);    Aeq(2*N,N+N)  =g*n^2*abs(Vtr(j-1,N))/Rh(Ytr(j-1,N),X,type)^(4/3) + 3*Vtr(j-1,N)/dx/2 + 1/dt;
%             Beq(2*N)=g*s+Vtr(j-1,N)/dt;
%         
%         % Megold�s
%         Aeq;
%         mo=Aeq\Beq; Ytr(j,:)=mo(1:N)'; Vtr(j,:)=mo(N+1:2*N)';
%         
%         % Sim�t�s
%         mo_Y=mo(1:N)'; N_half = ceil((N+1)/2); N_cut = ceil(N_half/1.5);
%         signal_Y=mo_Y-(mo_Y(end)-mo_Y(1))/L*Xmesh-mo_Y(1);
%             coefs_Y = fft(signal_Y);
%                 if (mod(N,2)==0)
%                     coefs_Y(N_half-N_cut:N_half+N_cut)=0;
%                 else
%                     coefs_Y(N_half-N_cut:N_half+N_cut+1)=0;
%                 end
%         signal_Y_cut = ifft(coefs_Y); mo_Y_cut=signal_Y_cut+(mo_Y(end)-mo_Y(1))/L*Xmesh+mo_Y(1);
%         Ytr(j,:)=mo_Y_cut;
%         
%         mo_V=mo(N+1:2*N)';
%         signal_V=mo_V-(mo_V(end)-mo_V(1))/L*Xmesh-mo_V(1);
%             coefs_V = fft(signal_V);
%                 if (mod(N,2)==0)
%                     coefs_V(N_half-N_cut:N_half+N_cut)=0;
%                 else
%                     coefs_V(N_half-N_cut:N_half+N_cut+1)=0;
%                 end
%         signal_V_cut = ifft(coefs_V); mo_V_cut=signal_V_cut+(mo_V(end)-mo_V(1))/L*Xmesh+mo_V(1);
%         Vtr(j,:)=mo_V_cut;
% 
% end
%     
