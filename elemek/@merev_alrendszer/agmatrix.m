function [A,B,C,D,mar]=agmatrix(mar,pf,flag,xr,varargin)

% pf: peremfelt�telek
% pf{:}{1}: csomopont szama

elem = mar.elemek;
csp  = mar.csp;
n_csp =length(csp);
n_elem=length(elem);

sor=1;
%fprintf('\n--------------');
for i=1:length(elem)
%    fprintf('\n %d',i);
    [temp1,elem{i}] = setport(flag,elem{i},xr(i),varargin,mar);  
    %fprintf('flag,elem: %8.3f % 8.3f\n',flag,elem{i}.y)
    %pause
    for k=1:length(temp1)
        A(sor,sor)=temp1{k}{1}; % A 
        B(sor,sor)=temp1{k}{2}; % B
        C(sor,sor)=temp1{k}{3}; % C
        D(sor,1)  =temp1{k}{4}; % D    
        
        ll=length(temp1{k})-4;
        for j=1:ll
            temp3=temp1{k}{4+j};
            for i=1:length(csp)
                if temp3(1)==i, P(sor,i)=temp3(2); end
            end
        end    
        sor=sor+1;
    end
end
mar.elemek=elem;

ll=1;
for i=1:length(elem)
    temp1=elem{i}.csp;
    switch length(temp1)
        case 1
            % Ez a sor megkeresi a csp m�trixban a csp sz�m�t
            ize=1; while ~(temp1(1)==ize), ize=ize+1; end
            K(ize,ll)=-1;
        case 2
            ize1=1; while ~(temp1(1)==ize1), ize1=ize1+1; end
            ize2=1; while ~(temp1(2)==ize2), ize2=ize2+1; end
            if isa(elem{i},'cso')	
                K(ize1,ll) =-1; 
                ll=ll+1; 
                K(ize2,ll) = 1;
            else
                K(ize1,ll)=-1; K(ize2,ll)=1;
            end
        end
        ll=ll+1;
    end

% csom�pontok hidrosztatikus magass�ga:
for i=1:length(elem)
    temp1=elem{i}.csp;
    if length(temp1)==2
        ize1=1; while ~(temp1(1)==ize1), ize1=ize1+1; end
        ize2=1; while ~(temp1(2)==ize2), ize2=ize2+1; end
        dh=mar.csp{ize1}{2}-mar.csp{ize2}{2};
        D(i,1) = D(i,1) - dh*elem{i}.ro*9.81;
    end
end

A=[A P; K zeros(length(csp),length(csp))];
B=[B zeros(length(elem),length(csp)); zeros(length(csp),length(csp)+length(elem))];
C=[C zeros(length(elem),length(csp)); zeros(length(csp),length(csp)+length(elem))];
D=[D; zeros(length(csp),1)];

% Kontinuit�si egyenletek

for i=1:length(csp)
%  fprintf('\n csp{%d}{3} hossza: %d',i,length(csp{i}{3})); disp(csp{i}{3});
%  fprintf('\n elem hossza:       %d\n',length(elem));
    for j=1:length(csp{i}{3})
        %temp = elem{abs(csp{i}{3}(j))};
        A(n_elem+i,abs(csp{i}{3}(j))) = sign(csp{i}{3}(j)); %*temp.ro;
    end
    D(n_elem+i) = -csp{i}{5};
end

% Peremfelt�telek
k=0; % Ez a hozz�adott t�rfogat�ramokat jelenti, rugalmas cs�vek miatt
for i=1:length(pf)

    user_cspszam=pf{i}{1};
    j=1; while ~(user_cspszam==mar.csp{j}{1}), j=j+1; end
    %   fprintf('Peremfelt�tel �t�ll�t�s...  %d --> %d\n\n',mar.csp{j}{1},j);  
    switch pf{i}{2}
        case 'p'
            fprintf('');
            % A csomponti konti helyett p=eloirt tip. egyenlet van.
            no_of_eq=n_elem + j;
            A(no_of_eq,:) = 0;
            A(no_of_eq,n_elem+j) = 1;
            B(no_of_eq,:) = 0;
            C(no_of_eq,:) = 0;
            D(no_of_eq)   = -pf{i}{3};
            
        case {'cso','viszkcso'}
            % Ujabb agaegyenlet hozzaadasa...
            k = k+1;
            no_of_eq = n_elem + n_csp+k;
            % Matrixok bovitese ujabb sorral es oszloppal
            A(no_of_eq,:) = 0; A(:,no_of_eq) = 0;
            B(no_of_eq,:) = 0; B(:,no_of_eq) = 0;
            C(:,no_of_eq) = 0; C(no_of_eq,:) = 0;
            % Kitoltes
            A(no_of_eq,n_elem+n_csp+k) = pf{i}{3}(2);            
            A(no_of_eq,n_elem+j) = pf{i}{3}(1);            
            D(no_of_eq)   = -pf{i}{3}(3);
	    
            % es a csomoponti kontiba ujabb tag.
            A(n_elem + j,n_elem+n_csp+k) = pf{i}{3}(4);
					      
        case 'vizszint_&_konti'
            % Ujabb agaegyenlet - p=eloirt - hozzaadasa...
            k=k+1; %% Itt nyultam bele!
            no_of_eq = n_elem + n_csp+k;
            % Matrixok bovitese ujabb sorral es oszloppal
            A(no_of_eq,:) = 0; A(:,no_of_eq) = 0;
            B(no_of_eq,:) = 0; B(:,no_of_eq) = 0;
            C(:,no_of_eq) = 0; C(no_of_eq,:) = 0;
            A(no_of_eq,no_of_eq) = 1;               % !!!!!!!
            D(no_of_eq) = pf{i}{3}(2);         % !!!!!!!
      

            % a csomponti kontiba ujabb tag

            A(n_elem + j,no_of_eq) = pf{i}{3}(3);
		
        
        otherwise
            error(['ismeretlen PF tipus:',pf{i}{2}]);
    end
end

%if strcmp(mar.nev,'sziv')
% mar.nev
% elem
%  A
%  B
%  C
%  D(1)
% pause
%end

