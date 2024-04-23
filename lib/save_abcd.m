function save_ABCD(fnev,elem,csp,A,B,C,D,flag,varargin)

ncsp=length(csp);
nQ=length(elem);

switch flag
    case 'x'
        fp=fopen(fnev,'a');
        fprintf(fp,'\n\nlepes: %d',varargin{1});
        ikszel(fp,A,'A',nQ,ncsp);
        ikszel(fp,B,'B',nQ,ncsp);
        ikszel(fp,C,'C',nQ,ncsp);
        ikszel_b(fp,D,'D',nQ,ncsp);
                fprintf(fp,'\n*******************************************************************************************');
        
    case 'val'
        xr=varargin{1}; xu=varargin{2};
        fp=fopen(fnev,'a+');
        fprintf(fp,'\n\nlepes: %d  iteracio: %d',varargin{3},varargin{4});
        kitolt(fp,A,'A',nQ,ncsp);        
        kitolt(fp,B,'B',nQ,ncsp);        
        kitolt(fp,C,'C',nQ,ncsp);        
        kitolt_b(fp,D,'D',nQ,ncsp);
        kitolt_x(fp,xr,xu,nQ,ncsp);
        fprintf(fp,'\n\n*******************************************************************************************');
        
    otherwise
        error('Bajjj van...');    
end

fclose(fp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ikszel(fp,M,neve,nQ,ncsp)
n1=length(M(:,1));
n2=length(M(:,1))-ncsp-nQ;

fprintf(fp,'\n\n\nRendszermátrixok: %s\n\n',neve); 
fprintf(fp,'    |'); for i=1:nQ,   fprintf(fp,' Q%2d |',i); end 
for i=1:ncsp, fprintf(fp,' p%2d |',i); end
for i=1:n2, fprintf(fp,' y%2d |',i); end

fprintf(fp,'\n----+'); for i=1:nQ, fprintf(fp,'-----+'); end
for i=1:ncsp, fprintf(fp,'-----+'); end
for i=1:n2, fprintf(fp,'-----+'); end

for i=1:n1
    if i<=nQ,       fprintf(fp,'\nQ%2d |',i);
    elseif i<=nQ+ncsp, fprintf(fp,'\np%2d |',i-nQ);
    else,           fprintf(fp,'\ny%2d |',i-nQ-ncsp); end

    for j=1:n1
        if M(i,j)==0, fprintf(fp,'     ');
        else,         fprintf(fp,'  x  ');  end
        
        if j==nQ, fprintf(fp,'|');
        else, if j==nQ+ncsp, fprintf(fp,'|');
            else, if j==n1, fprintf(fp,'|');
                else fprintf(fp,' ');
        end,end,end

    end
    if i==nQ | i==nQ+ncsp | i==n1
        fprintf(fp,'\n----+'); for i=1:n1, fprintf(fp,'-----+'); end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kitolt(fp,M,neve,nQ,ncsp)
n1=length(M(:,1));
n2=length(M(:,1))-ncsp-nQ;
fprintf(fp,'\n\n\nRendszermátrixok: %s\n\n',neve);
fprintf(fp,'    ||'); for i=1:nQ,   fprintf(fp,'     Q%2d     |',i); end
fprintf(fp,'|'); for i=1:ncsp, fprintf(fp,'     p%2d     |',i); end
fprintf(fp,'|'); for i=1:n2, fprintf(fp,'     y%2d     |',i); end

fprintf(fp,'\n----++'); for i=1:nQ, fprintf(fp,'-------------+'); end
fprintf(fp,'+'); for i=1:ncsp, fprintf(fp,'-------------+'); end
fprintf(fp,'+'); for i=1:n2, fprintf(fp,'-------------+'); end

for i=1:n1
    if i<=nQ, fprintf(fp,'\nQ%2d ||',i);
    elseif i<=nQ+ncsp,     fprintf(fp,'\np%2d ||',i-nQ); 
    else, fprintf(fp,'\ny%2d ||',i-nQ-ncsp); end
    
    for j=1:n1
        if M(i,j)==0
            fprintf(fp,'             ');
        else
            fprintf(fp,' %+5.3e ',M(i,j));
        end
        fprintf(fp,'|');
        if j==nQ, fprintf(fp,'|'); end
        if j==nQ+ncsp, fprintf(fp,'|'); end        
    end
    if i==nQ | i==nQ+ncsp | i==n1
        fprintf(fp,'\n----++'); for i=1:nQ, fprintf(fp,'-------------+'); end
        fprintf(fp,'+'); for i=1:ncsp, fprintf(fp,'-------------+'); end    
        fprintf(fp,'+'); for i=1:n2, fprintf(fp,'-------------+'); end    
    end
end

% fprintf(fp,'\n----++'); for i=1:nQ, fprintf(fp,'------------+'); end
% fprintf(fp,'+'); for i=1:ncsp, fprintf(fp,'------------+'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ikszel_b(fp,D,neve,nQ,ncsp)

n1=length(D); n2=length(D)-ncsp-nQ;

fprintf(fp,'\n\n\nInhomogén vektor: D\n\n'); 
fprintf(fp,  '    | D |');
fprintf(fp,'\n----+---+');

for i=1:n1
    if i<=nQ
        fprintf(fp,'\nQ%2d |',i);
        if D(i)==0, fprintf(fp,'   |'); else fprintf(fp,' x |'); end
    elseif i<=nQ+ncsp
        fprintf(fp,'\np%2d |',i-nQ);
        if D(i)==0, fprintf(fp,'   |'); else fprintf(fp,' x |'); end
    else
        fprintf(fp,'\ny%2d |',i-nQ-ncsp);
        if D(i)==0, fprintf(fp,'   |'); else fprintf(fp,' x |'); end
    end
    
    if i==nQ | i==nQ+ncsp | i==n1
        fprintf(fp,'\n----+---+');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kitolt_b(fp,D,neve,nQ,ncsp)

n1=length(D); n2=length(D)-ncsp-nQ;

fprintf(fp,'\n\n\nInhomogén vektor: D\n\n'); 
fprintf(fp,  '    ||      D      |');
fprintf(fp,'\n----++-------------+');

for i=1:n1
    if i<=nQ, fprintf(fp,'\nQ%2d || %+5.3e |',i,D(i));
    elseif i<=nQ+ncsp, fprintf(fp,'\np%2d || %+5.3e |',i-nQ,D(i));
    else, fprintf(fp,'\ny%2d || %+5.3e |',i-nQ-ncsp,D(i)); end
    
    if i==nQ | i==nQ+ncsp | i==n1, fprintf(fp,'\n----++-------------+');  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kitolt_x(fp,xr,xu,nQ,ncsp)

n1=length(xr); n2=length(xr)-ncsp-nQ;

fprintf(fp,'\n\n\nAz egyenletrendszer régi és javított megoldása:\n\n'); 
fprintf(fp,  '    ||    xrégi    |     xúj     |');
fprintf(fp,'\n----++-------------+- -----------+');

for i=1:n1
  if i<=nQ, fprintf(fp,'\nQ%2d || %+5.3e | %+5.3e |',i,xr(i),xu(i));
  elseif i<=nQ+ncsp, fprintf(fp,'\np%2d || %+5.3e | %+5.3e |',i-nQ,xr(i),xu(i));
  else fprintf(fp,'\ny%2d || %+5.3e | %+5.3e |',i-nQ-ncsp,xr(i),xu(i)); end

  if i==nQ | i==nQ+ncsp | i==n1, fprintf(fp,'\n----++-------------+-------------+'); end
end
