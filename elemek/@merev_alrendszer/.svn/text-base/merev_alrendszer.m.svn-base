function mar = merev_alrendszer(be_nev,be_elemek,be_cspok,be_wdir)



mar.plot_iter =0;
mar.save_level=0;
mar.wdir = be_wdir;
mar.dtki = 0.0;
mar.dtkiorig = 0.0;

fnev=be_nev;
mar.fp=fopen(strcat(fnev,'.dat'),'r');
mar.fnev=strcat(fnev,'.out');
mar.nev=fnev;
mar.n_elem=length(be_elemek);
mar.n_csp =length(be_cspok);
for i=1:mar.n_elem, mar.elemek{i}=be_elemek{i}; end
for i=1:mar.n_csp, mar.csp{i}=be_cspok{i}; end

% A csomopont struktura {3} feltoltese: erkezo (+) es tavozo (-) agak
% felsorolasa
for i=1:length(mar.elemek)
    elem = mar.elemek{i};
    cspv = elem.csp;
    mar.csp{cspv(1)}{3} = [mar.csp{cspv(1)}{3}, -i];
    
    if length(cspv)==2
        mar.csp{cspv(2)}{3} = [mar.csp{cspv(2)}{3}, i];
    end
end

mar.warnings=0; % 0->off, 1->on

if mar.save_level>0 mar.fout=fopen(strcat(fnev,'.out'),'w'); end

%_____________________________________________________________________
% Meg a csomopontok atszamozasa elott elmentjuk a bemenetet           |
if mar.save_level>0 
    workdir = mar.wdir;
    fnev=fullfile(workdir,strcat(mar.nev,'.out'));
    %mar.fout=fopen(strcat(fnev,'.out'),'w'); 
    mar.fout=fopen(fnev,'w'); 
end
if mar.save_level>0
    save_merev(fnev,mar.elemek,mar.csp); fclose(mar.fout);
end

if mar.save_level>1    
    for i=1:mar.n_elem, info(mar.elemek{i},2,fnev); end    
end

% %% Rendszermatrix felepitese
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ll=1;
% for j=1:mar.n_elem
%     for i=1:length(mar.elemek{j}.csp)
%         temp=mar.elemek{j}.csp; 
% 	k=temp(i);
% 	cspszam=1;
% 	while ~(k==cspszam), cspszam=cspszam+1; end
%         RM(cspszam,j)=i; ll=ll+1;
%     end
% end
% 
% %fprintf('\n A rendszerm�trix: (sor -> csom�pont, oszlop -> �g)\n'); disp(RM);
% flag
% pause
% if flag==0
%   for i=1:length(mar.csp)
%     for j=1:length(mar.elemek)
%         if RM(i,j)==1, mar.csp{i}{3}=[mar.csp{i}{3}, -j]; end
%         if RM(i,j)==2, mar.csp{i}{3}=[mar.csp{i}{3}, j]; end
%     end
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mar.x=[];
mar.plot_it=[0 0];
mar.lepes=1; mar.iter=1;
mar.t=-1e-10; % Ez nagyon fontos!
mar.gorbek=[];

mar = class(mar,'merev_alrendszer');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_merev(fnev,elem,csp)

%ncsp=length(csp); nelem=length(elem);

fout=fopen(fnev,'a');
fprintf(fout,'\nMerev alrendszer adatai\n');
clockk=clock;
fprintf(fout,'\nletrehozas datuma: %g/%g/%g  %g:%g:%g',clockk(1),clockk(2),clockk(3),clockk(4),clockk(5),round(clockk(6)));
fprintf(fout,'\nnev:               %s',fnev);
fprintf(fout,'\ncsomopontok szama: %d',length(csp));
fprintf(fout,'\nelemek szama:      %d',length(elem));

 fprintf(fout,'\n\n\nElemek felsorol�sa:');
 fprintf(fout,'\n\n ssz.  nev           csp1          csp2 ');
 fprintf(fout,'\n-----------------------------------------------');
 for i=1:length(elem)
     fprintf(fout,'\n%3d.  %s',i,elem{i}.nev);
     temp=elem{i}.csp;
     switch length(temp)
         case 1
             fprintf(fout,'\t %5s (%2d)',csp{temp(1)}{6},temp(1));
         case 2
            fprintf(fout,'\t %5s (%2d)  -> %5s (%2d)',csp{temp(1)}{6},temp(1),csp{temp(2)}{6},temp(2));
     end
end

fprintf(fout,'\n\n\nCsomopontok felsorolasa:');
fprintf(fout,  '\n\nssz.    nev         h[m]      m_be[kg/s]      befuto(+) es tavozo(-) agak');
fprintf(fout,    '\n-----------------------------------------------------------------------------');
for i=1:length(csp)
    fprintf(fout,'\n%3d. %8s    %+6.4e    %+5.2e       ',i,csp{i}{6},csp{i}{2},csp{i}{5});
    temp=csp{i}{3};
    for j=1:length(temp), fprintf(fout,' %+d',temp(j));  end    
end

fprintf(fout,'\n\n*****************************************************************************');
fprintf(fout,'\n*                     Elemek reszletes leirasa                              *');
fprintf(fout,'\n*****************************************************************************\n');
