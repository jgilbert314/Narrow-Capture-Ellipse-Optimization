%%what LGO does

clear all;
close all
clc;

%%
lgo__license = regexprep(which('matlablgo'), 'matlablgo.mexw32','lgo_license.dat');

lgo__fid = fopen(lgo__license,'rt');
lgo__in = fscanf(lgo__fid,'%s')
fclose(lgo__fid);

%%
[lgo__s,lgo__w] = system('ipconfig /all | findstr [0-9]-[0-9]')

lgo__x = findstr(lgo__w,':')


%%

lgo__out = [];
for lgo__indexer=1:min(3,length(lgo__x));
    lgo__out = [lgo__out; regexprep(lgo__w(lgo__x(lgo__indexer)+2 : lgo__x(lgo__indexer)+18), '-', '')];
end;
lgo__out


lgo__out = dec2base(hex2dec(lgo__out)*7,29)

lgo__out = sprintf('%s',lgo__out)

lgo__vld = strcmp(lgo__in,lgo__in)

