%%
clear;

numOpt = 311;
engArr = zeros(1, numOpt);

bestEng = 0;
bestInd = 0;

for itr = 1:numOpt
    filename = ['C:\Projects\Narrow Escape\Elliptical Domain\Ellipse\LogFiles\Log_eps-diff500-ecc0p125\log', ...
        num2str(itr), '.csv' ];
    
    dataArr = readLogFile(filename);
    thisEng = dataArr{1}(end-1);
    engArr(itr) = thisEng;
    if ( (thisEng < bestEng) && (thisEng > -0.2) )
       bestEng = thisEng;
       bestInd = itr;
    end
        
    
end

%%
plot(engArr);
ylim([-0.16, -0.14]);