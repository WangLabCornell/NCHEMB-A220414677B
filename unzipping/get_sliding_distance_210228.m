% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

filename = uigetfile('*.*','Select the DNA file','multiselect','on');
%if isequal(filename,0)
%    disp('User selected Cancel');
%else
%    if iscell(filename)
%        K = length(filename);
%        tracenuc_a = [];
%        for i = 1 : K
%            load(filename{i});
%            trace_a = [trace_a trace];
%        end
%    else
        load(filename);
%        trace_a = trace;
%    end
%end

%load('data1.mat');

P = length(trace);

for j = 1 : P
    
    if strcmp(trace{j}.selected ,'yes')
        if ~strcmp(trace{j}.type, 'DNA')
            %% find sliding distance
            name = trace{j}.name;
            Times = trace{j}.Times ;
            jindex =             trace{j}.jindex ;
            ForcepN =             trace{j}.ForcepN ;
            jindex_s =             trace{j}.jindex_s ;
            ForcepN_s =             trace{j}.ForcepN_s ;
            jindex_s_n = trace{j}.jindex_s_n;

            jindex_t = trace{j}.jindex_t ;
            forcepN_t =  trace{j}.forcepN_t ;

            dt = mean(diff(Times));
            N = round(1/40/dt);     % resampling to 40 Hz
            jindex_s2 = movingmean(jindex,N);
            ForcepN_s2 = movingmean(ForcepN,N);

            [xData_theo, yData_theo] = prepareCurveData( jindex_t, forcepN_t );
            yData_theo = yData_theo; 
            fitresult_theo = splinefit(xData_theo, yData_theo,round(max(xData_theo)/5));

            [xData_exp, yData_exp] = prepareCurveData( jindex_s2, ForcepN_s2 );
            fitresult_exp = splinefit( xData_exp, yData_exp, round(max(xData_exp)/3));

            j_above = jindex_s_n(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+2 & jindex_s_n >= 0);
            f_above = ForcepN_s(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+2 & jindex_s_n >= 0);

            %% Total sliding distance
            djr = 20;
            jr = 0 : djr : 4000;
            jrn = jr + djr/2;
            njr = histc(j_above, jr);

            n_thres = 1;

            index_sliding = find(njr > n_thres);
            jrn_sliding = jrn(index_sliding);

            index1 = jrn_sliding(diff(jrn_sliding) > djr);
            index2 = jrn_sliding(find(diff(jrn_sliding) > djr)+1);
            
            if any(index1) && any(index2)

                j_slidingr_sf =  [jrn_sliding(1) index1(1) ; index2(1:end-1)' index1(2:end)'; index2(end) jrn_sliding(end)];

                j_sliding_total = sum(j_slidingr_sf(:,2) - j_slidingr_sf(:,1)+1);
            else
                j_sliding_total = 0;
            end
            
            
            trace{j}.j_sliding_total = j_sliding_total;
        else
            trace{j}.j_sliding_total = 0;
            
        end
        
    end


end

save(filename, 'trace');

clear
