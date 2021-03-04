%%
clear;

% good_file = 'GoodLog0.csv';
% comp_file = 'log1.csv';S
% new_file = 'GoodLog1.csv';
% 
% 
% goodData = readLogFile(good_file); % Load most recent archive
% compData = readLogFile(comp_file);
% 
% % Write current optimums to log file
% [ writeData ] = updateDataSet(compData, goodData);
% writeLogFile(writeData, new_file);

%%
clear;
filename_list = {'GoodLog0.csv', 'GoodLog9.csv'};

% Domain setup
ecc = 1;                  % Eccentricity of ellipse
epsilon = 0;                  % Trap size
InputStruct = buildInputStruct(ecc, epsilon); % Initialize arguments to pass to optimizer


hold on
for itr = 1:length(filename_list)
    filename = filename_list{itr};
    dataArr = readLogFile(filename);
    
%     nData = length(dataArr);
%     p = zeros(1, nData);
%     N_list = zeros(1, nData);
%     for itrN = 1:nData
%         N = dataArr{itrN}(1);
%         xOpt = dataArr{itrN}(2:2*N+1);
%         [A_0, nu_0] = calcCartesian2Elliptical(xOpt(1:N), xOpt(N+1:end), InputStruct.mu);
%         x = [A_0, nu_0];
%         
%         InputStruct.N = N;
%         p(itrN) = meritFuncEllipse(x, InputStruct);
%         const = -1/2/N*log(epsilon); % Constant term used to calculate AMFPT
%         AMFPT = const + 2*pi/N*p(itr);
%         if (AMFPT <= 0)
%             disp(['Unphysical Result:', ' ', num2str(N), ' ', num2str(AMFPT)]);
%         end
%         N_list(itrN) = N;
%         p(itrN) = AMFPT;
%     end
    
    p = getIntEng(dataArr);
    N_list = getN(dataArr);
    AMFPT = -1/2*log(epsilon)./N_list + 2*pi*p./N_list;
    plot(N_list, AMFPT, '- .');
end
hold off
legend('Old', 'New');