% Copyright 2022 Wang Lab
% 
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

%load('data_scII_50uMETO_neg_all_pause_analysis.mat');

load('dat_temp_box92.mat');

i = 58;% 224; %180; % 231; %287;
z_level = {};

f1 = figure;
subplot(2,1,1)
plot( trace{i}.time_holding_new, trace{i}.z_holding_new, ':' ,'Color',[0.5, 0.5, 0.5]);
hold all

for j = 1 : size(trace{i}.index_pause_t,2)
    meanz(j) = mean(trace{i}.z_holding_new(trace{i}.index_pause_t{j}));
end

for j = 1 : size(trace{i}.index_pause_t,2)
    
    z_level{j} = ones(1, length(trace{i}.index_pause_t{j})) * meanz(j);
    
end


plot(trace{i}.time_holding_new(cell2mat(trace{i}.index_pause_t)), cell2mat(z_level),'-k');


subplot(2, 1,2)

for j = 1 : size(trace{i}.index_pause_t,2)
    
    turn_level{j} = ones(1, length(trace{i}.index_pause_t{j})) * trace{i}.turn_mean_u(j);
    
end

plot( trace{i}.time_holding_new, trace{i}.turn_relaxed2_new, ':' ,'Color',[0.5, 0.5, 0.5]);
hold all
plot(trace{i}.time_holding_new(cell2mat(trace{i}.index_pause_t)), cell2mat(turn_level),'-k');
