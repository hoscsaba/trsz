function save_res(mar)

fp = fopen(mar.fnev,'a');
nQ = length(mar.elemek);
csp = mar.csp;
xu = mar.x;

%fprintf(fp,'\n\n*****************************************************************************\n');
fprintf(fp,'\nEREDMÉNYEK:\n');
fprintf(fp,'\ntérfogatáramok és tömegáramok:');
    nQ = length(mar.elemek);
    elemek = mar.elemek;
    for i=1:length(elemek)
        fprintf(fp,'\nQ%2d = %+8.5f [m^3/s]   m%2d = %+5.2f [kg/s]',i,elemek{i}.Q,i,elemek{i}.Q*elemek{i}.ro);
    end
    
    fprintf(fp,'\n\nnyomások:');
    for i=1:length(mar.csp)
        fprintf(fp,'\np%2d = %+5.2f [bar] = %+6.1f [m]',mar.csp{i}{1},mar.csp{i}{4}/1e5,mar.csp{i}{4}/1e3/9.81);
    end
    fprintf(fp,'\n');

fclose(fp);
   