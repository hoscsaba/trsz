function elem=setfluid(elem,foly)

ok=0; fp=fopen('fluids.dat','r');
for i=1:4, temp=fgetl(fp); end

while ~feof(fp) & ok==0
    temp=fgetl(fp);
    [nev,mar]=strtok(temp,',');
    elem.folynev=nev;
%    fprintf('\n %s   %s   %d',nev,foly,strncmp(nev,foly,length(foly)));
    if strncmp(nev,foly,length(foly))
        [ro,mar]=strtok(mar,','); elem.ro=str2num(ro);
        [nu,mar]=strtok(mar,','); elem.nu=str2num(nu);
        [B,mar] =strtok(mar,','); elem.B =str2num(B);
        ok=1;
    end
end
fclose(fp);

if ok==0
    error(['Ismeretlen folyad�k: ',foly]); 
else
    mu=str2num(nu)*str2num(ro);
    elem.mu=mu;
end
 %ize
