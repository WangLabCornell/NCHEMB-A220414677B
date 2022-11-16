% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

function [Times, jindex, ForcepN, fmax, fmean, jend] = Show_each_trace(filename, filename_t, f)
%% As of 07/15/2015 - Tung Le

    %% Open the file with experimental data
    fileID = fopen(filename,'r');
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    Times = dataArray{:, 1};
    AOMSentV = dataArray{:, 13};
    ForcepN = dataArray{:, 23};
    jindex = dataArray{:, 24};
    clearvars delimiter startRow formatSpec fileID dataArray ans;
    
        %% Open the theoretical curve
    %filename_t = 'force_kkb_1p6kbafterA20_100mMNaCl_4mMMg';
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%[^\n\r]';
    fileID_t = fopen(filename_t,'r');
    dataArray = textscan(fileID_t, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID_t);

    extensionnm = dataArray{:, 1};
    forcepN_t = dataArray{:, 2};
    jindex_t = dataArray{:, 3};
    
    %clearvars filename_t delimiter startRow formatSpec fileID_t dataArray ans;
    clearvars  delimiter startRow formatSpec fileID_t dataArray ans;
    
    
    
    %% Filter data based on the AOMSent Magnitude and the length of the sequence
    index1 = find(AOMSentV > 0.1 & ForcepN > 10) ; 
    %% Discard the last np data points
    
    jindex = jindex(index1); ForcepN = ForcepN(index1);Times = Times(index1); 
    jindex = jindex - median(jindex(ForcepN < 12 & ForcepN > 10));
    
    index2 =  ForcepN > 10;
    jindex = jindex(index2); Times = Times(index2); ForcepN = ForcepN(index2);
    
    jindex_s = movingmean(jindex,5); ForcepN_s =   movingmean(ForcepN,5);
    
    %np = 10;
    %Times_n(end-np:end) = []; jindex_n(end-np:end) = []; ForcepN_n(end-np:end) = [];
    
    
    fmax = max(ForcepN_s);
    fmean = mean(ForcepN_s);
    jend= jindex(end-4);
    
%% Ploting
figure(f);
            

plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5,'LineStyle','-');   

plot(jindex_s,ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

axis([-100 4200 5 65]);

textp1 = { ['Max force: ' , num2str(fmax) ' pN'] ,['Mean force: ' num2str(fmean) ' pN'] ,['j-index end @: ' num2str(jend) ' bp']};   
text(0, 65, textp1, 'HorizontalAlignment', 'left','VerticalAlignment', 'top','FontSize',11, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','k');

set(gca,'XGrid','on');
title(filename);

%set(gca, 'XTick', -100:200:4600);
xlabel('j-index (bp)', 'FontSize',13); 
ylabel('Force (pN)', 'FontSize',13);
title(['Force vs. j-index ' filename], 'FontSize',13);



end
