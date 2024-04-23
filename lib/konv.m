% Adatfile konverter

function konv(fname)

wdir = what;

fname_in = fullfile(wdir.path,fname);

[data text raw] = xlsread(fname_in,1);

[maxsor maxoszlop] = size(raw);

fname_out = strcat(fname(1:end-4),'.tpr');
f01 = fopen(fname_out,'w');

for i = 1:maxsor
    tmp = raw(i,:);
    for j = 1:length(tmp)
        if isnan(tmp{j})
            break
        elseif j == 1
            fprintf(f01,'%s',num2str(tmp{j}));
        else
            fprintf(f01,',%s',num2str(tmp{j}));
        end
    end
    fprintf(f01,'\n');
end

fclose(f01);

