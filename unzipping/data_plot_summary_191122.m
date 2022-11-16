% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.


%filename = uigetfile('*.*','Select the DNA file','multiselect','on');

        %load(filename);


%P = length(trace);

%for j = 1 : P
    
%    if strcmp(trace{j}.selected ,'yes')
%        if ~strcmp(trace{j}.type, 'DNA')
%            %% find sliding distance
%            name = trace{j}.name;
%            Times = trace{j}.Times ;
%            jindex =             trace{j}.jindex ;
%            ForcepN =             trace{j}.ForcepN ;
%            jindex_s =             trace{j}.jindex_s ;
%            ForcepN_s =             trace{j}.ForcepN_s ;
%            jindex_s_n = trace{j}.jindex_s_n;

%            jindex_t = trace{j}.jindex_t ;
%            forcepN_t =  trace{j}.forcepN_t ;

%            dt = mean(diff(Times));
%            N = round(1/40/dt);     % resampling to 40 Hz
%            jindex_s2 = movingmean(jindex,N);
%            ForcepN_s2 = movingmean(ForcepN,N);

%            [xData_theo, yData_theo] = prepareCurveData( jindex_t, forcepN_t );
%            yData_theo = yData_theo; 
%            fitresult_theo = splinefit(xData_theo, yData_theo,round(max(xData_theo)/5));

%            [xData_exp, yData_exp] = prepareCurveData( jindex_s2, ForcepN_s2 );
%            fitresult_exp = splinefit( xData_exp, yData_exp, round(max(xData_exp)/3));

%            j_above = jindex_s_n(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+2 & jindex_s_n >= 0);
%            f_above = ForcepN_s(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+2 & jindex_s_n >= 0);

            %% Total sliding distance
 %           djr = 20;
%            jr = 0 : djr : 4000;
%            jrn = jr + djr/2;
%            njr = histc(j_above, jr);

%            n_thres = 1;

%            index_sliding = find(njr > n_thres);
%            jrn_sliding = jrn(index_sliding);

%            index1 = jrn_sliding(diff(jrn_sliding) > djr);
%            index2 = jrn_sliding(find(diff(jrn_sliding) > djr)+1);
            
%            if any(index1) && any(index2)

%                j_slidingr_sf =  [jrn_sliding(1) index1(1) ; index2(1:end-1)' index1(2:end)'; index2(end) jrn_sliding(end)];

%                j_sliding_total = sum(j_slidingr_sf(:,2) - j_slidingr_sf(:,1)+1);
%            else
%                j_sliding_total = 0;
%            end
            
            
%            trace{j}.j_sliding_total = j_sliding_total;
%        else
%            trace{j}.j_sliding_total = 0;
            
%        end
        
%    end


%end

%save(filename, 'trace');

%clear

filename = uigetfile('*.*','Select the DNA file','multiselect','on');
load(filename);
P = length(trace);

j= 0;

j_sliding_total_a = [];
for i = 1 : P
    if strcmp(trace{i}.selected ,'yes')
        
           j = j+1;
           
           name_a{j} = trace{i}.name;
           
           j_sliding_total_a(j) = trace{i}.j_sliding_total;

           if ~strcmp(trace{i}.type, 'DNA')
               j_slidingr_sf{j} = trace{i}.j_slidingr_sf;
           else
               j_slidingr_sf{j} = [];
           end
           
           Times{j} = trace{i}.Times;
           jindex_s_n{j} = trace{i}.jindex_s_n;
           ForcepN_s{j} = trace{i}.ForcepN_s;
           jindex_t{j} = trace{i}.jindex_t;
           forcepN_t{j} = trace{i}.forcepN_t;
           
           
           j_first(j) = trace{i}.j_first;
           f_first(j) = trace{i}.f_first;
           fmax(j) = trace{i}.fmax;
           fmean(j) = trace{i}.fmean;
           jend(j) = trace{i}.jend;
           j_fmax(j) = trace{i}.j_fmax;

               type{j} = trace{i}.type;

    end
end

M = j;

index_dna = [];
index_protein = [];
index_notallign = [];
j = 0;
for i = 1 : M
    
    if strcmp(type{i}, 'DNA')
        index_dna = [index_dna i];
    elseif strcmp(type{i}, 'protein')
        index_protein = [index_protein i];
    else
        index_notallign = [index_notallign i];
    end
end

%% Histogram of the final J-index of naked DNA only
f1 = figure;
hist(jend(index_dna),40);
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('Final j-index (bp)','FontSize',14);
ylabel('Count','FontSize',14);
title(['Final j-index of Naked DNA, N = ' num2str(length(index_dna))], 'FontSize',14, 'FontName', 'Calibri','FontWeight','bold');

%% Plot all naked DNA traces
f2 = figure;
hold all
for i = 1 : length(index_dna)
    plot(jindex_s_n{index_dna(i)}, ForcepN_s{index_dna(i)}, '-', 'Color', rand(1,3), 'LineWidth', 0.5);
end
plot(jindex_t{1},  forcepN_t{1}, 'Color',[0.5, 0.5, 0.5], 'LineWidth', 2);

set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('J-index (bp)','FontSize',14);
ylabel('Force (pN)','FontSize',14);
title(['Naked DNA, N = ' num2str(length(index_dna))], 'FontSize',14, 'FontName', 'Calibri','FontWeight','bold');
xlim([-100 4500]);
ylim([9 20]);
grid on

disp('# of DNA traces'); length(index_dna)

jmax_theory = jindex_t{1}(end);

index_notdna = [index_protein index_notallign];

j_first_notdna = j_first(index_notdna);
f_first_notdna = f_first(index_notdna);
fmax_notdna = fmax(index_notdna);
fmean_notdna = fmean(index_notdna);
jend_notdna = jend(index_notdna);
j_fmax_notdna = j_fmax(index_notdna);

j_end_theory = 4000;
p_threshold = 0.8;
%% Breakage probability of traces with bound protein
jend_notdna_within = jend_notdna(j_first_notdna >= 50);
disp('# of traces with DNA-bound protein excluding fork binding that breaks');n_breakage_notdna = sum(jend_notdna_within <= p_threshold * j_end_theory)
disp('# of traces with DNA-bound protein excluding fork binding');length(jend_notdna_within)
disp('Tether breakage fraction');p_breakage_notdna = n_breakage_notdna / length(jend_notdna_within)
dp_breakage_notdna = sqrt( p_breakage_notdna * (1- p_breakage_notdna) / length(jend_notdna_within))

%% Breakage probability of traces with bound protein and total
jend_all_within = [jend_notdna(j_first_notdna >= 50) jend(index_dna)];
n_breakage_all = sum(jend_all_within <= p_threshold * j_end_theory);
p_breakage_all = n_breakage_all / length(jend_all_within);
dp_breakage_all = sqrt( p_breakage_all * (1- p_breakage_all) / length(jend_all_within));

%% Breakage probability of traces with DNA only
jend_dna = [jend(index_dna)];
n_breakage_dna = sum(jend_dna <= p_threshold * j_end_theory);
p_breakage_dna = n_breakage_dna / length(jend_dna);
dp_breakage_dna = sqrt( p_breakage_dna * (1- p_breakage_dna) / length(jend_dna));

j_sliding_total_a_nodna = j_sliding_total_a(index_notdna);
%% Sliding distance for traces with boudn protein
j_sliding_total_a_nodna_within = j_sliding_total_a_nodna(j_first_notdna >= 50);
j_sliding_total_a_nodna_within = j_sliding_total_a_nodna_within(isfinite(j_sliding_total_a_nodna_within));
fs= figure;
hist(j_sliding_total_a_nodna_within);
xlabel('Total sliding distance (bp)','FontSize',14);
ylabel('Count','FontSize',14);
xlim([0 4000]);
title(['For J >= 50bp, N = ' num2str(length(j_sliding_total_a_nodna_within )) '; mean: ' num2str(mean(j_sliding_total_a_nodna_within)) ' \pm ' num2str(std(j_sliding_total_a_nodna_within)) ', sem: ' num2str(std(j_sliding_total_a_nodna_within)/sqrt( length(j_sliding_total_a_nodna_within) ))]);


%% First binding location vs. max force
f3 = figure;
plot(j_first_notdna, fmax_notdna, 'ow', 'MarkerSize',6, 'MarkerFaceColor','k', 'LineWidth',0.5);
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('J-index of first binding (bp)','FontSize',14);
ylabel('Max Force of the trace (pN)','FontSize',14);
xlim([-100 4000]);
title(['All, N = ' num2str(length(j_first_notdna)) '; J >= 50bp, N = ' num2str(length(j_first_notdna(j_first_notdna>=50)))]);

%% unzipping distance from first binding
dj = jend_notdna - j_first_notdna;
jr = 0:100:4000;
jrn = jr+50;
ndj = histc(dj,jr);
f4 = figure;
bar(jrn, ndj, 'BarWidth',1);
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('Unzipping distance from first binding (bp)','FontSize',14);
ylabel('Count','FontSize',14);
xlim([-100 4000]);
title(['All, N = ' num2str(length(dj))]);

%% unzipping distance from first binding for binding over 50 bp
dj_withindna = jend_notdna((j_first_notdna >=50)) - j_first_notdna((j_first_notdna >=50));
jr = 0:100:4000;
jrn = jr+50;
ndj_withindna = histc(dj_withindna,jr);
f4_2 = figure;
bar(jrn, ndj_withindna, 'BarWidth',1);
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('Unzipping distance from first binding (bp)','FontSize',14);
ylabel('Count','FontSize',14);
xlim([-100 4000]);
title(['For binding within DNA, N = ' num2str(length(dj_withindna))]);

%% Histogram of Location of first binding
f5 = figure;
n_jfirst = histc(j_first_notdna,jr);
bar(jrn, n_jfirst, 'BarWidth',1);
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('Location of first binding (bp)','FontSize',14);
ylabel('Count','FontSize',14);
xlim([-100 4000]);

%% Max force histogram
fmax_notdna_withinDNA = fmax_notdna(j_first_notdna >=50);
fmax_notdna_withinDNA = fmax_notdna_withinDNA(isfinite(fmax_notdna_withinDNA));
fr = 0:4:120;
frn = fr+2;
n_fmax = histc(fmax_notdna, fr);
n_fmax_withinDNA = histc(fmax_notdna_withinDNA, fr);
f6 = figure;
hold all
b1 = bar(frn, n_fmax, 'BarWidth',1, 'FaceColor', [0.5, 0.5, 0.5]);
b2 = bar(frn, n_fmax_withinDNA, 'BarWidth',1, 'FaceColor', 'r');
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('Max Force (pN)','FontSize',14);
ylabel('Count','FontSize',14);
title({['All: N = ' num2str(length(fmax_notdna)) ', With J-index >= 50 bp: N = ' num2str(length(fmax_notdna_withinDNA))], [num2str(mean(fmax_notdna_withinDNA)) ', sd: ' num2str(std(fmax_notdna_withinDNA)) ', sem: ' num2str(std(fmax_notdna_withinDNA)/sqrt(length(fmax_notdna_withinDNA)) )]});

%% All traces with protein binding
f7 = figure;
hold all
for i = 1 : length(index_protein)
    plot(jindex_s_n{index_protein(i)}, ForcepN_s{index_protein(i)}, '-', 'Color', rand(1,3), 'LineWidth', 0.5);
end
for i = 1 : length(index_notallign)
    plot(jindex_s_n{index_notallign(i)}, ForcepN_s{index_notallign(i)}, '-', 'Color', rand(1,3), 'LineWidth', 0.5);
end
plot(jindex_t{1},  forcepN_t{1}, 'Color',[0.5, 0.5, 0.5], 'LineWidth', 2);

set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('J-index (bp)','FontSize',14);
ylabel('Force (pN)','FontSize',14);
title(['With bound protein, N = ' num2str(length(index_protein)+ length(index_notallign))], 'FontSize',14, 'FontName', 'Calibri','FontWeight','bold');
xlim([-100 4500]);
ylim([9 65]);
grid on

%% Number of sliding region for nearly full length tether
length_max = 0;%4000*0.9; % for eto only
j_slidingr_sf_protein  = j_slidingr_sf(index_notdna);
j_slidingr_sf_protein_full = j_slidingr_sf_protein(j_first_notdna >= 50 & jend_notdna >= length_max);
name_a_notdna = name_a(index_notdna);
name_a_notdna_full = name_a_notdna(j_first_notdna >= 50 & jend_notdna >= length_max);

L = length(j_slidingr_sf_protein_full);

for i = 1 : L
    num_region(i) = sum( j_slidingr_sf_protein_full{i}(:,2) - j_slidingr_sf_protein_full{i}(:,1) > 0);     
end

f9 = figure;
num_r = 0: 20;
n_num = histc(num_region, num_r);
bar(num_r,n_num, 'BarWidth',1);
xlabel('Number of sliding regions','FontSize',14);
ylabel('Count','FontSize',14);
xlim([0 15]);
title(['Number of sliding regions, full-length (N = ' num2str(length(num_region )) '); mean: ' num2str(mean(num_region)) ' \pm ' num2str(std(num_region)) ', sem: ' num2str(std(num_region)/sqrt( length(num_region) ))]);

disp('# of sliding regions'); mean(num_region)
%% Sliding distance for traces with boudn protein, full length
j_sliding_total_a_nodna_within_full = j_sliding_total_a_nodna(j_first_notdna >= 50 & jend_notdna >= length_max);
j_sliding_total_a_nodna_within_full = j_sliding_total_a_nodna_within_full(isfinite(j_sliding_total_a_nodna_within_full));
f8= figure;
hist(j_sliding_total_a_nodna_within_full);
xlabel('Total sliding distance, full-length tether (bp)','FontSize',14);
ylabel('Count','FontSize',14);
xlim([0 4000]);
title(['For J >= 50bp, full-length, N = ' num2str(length(j_sliding_total_a_nodna_within_full )) '; mean: ' num2str(mean(j_sliding_total_a_nodna_within_full)) ' \pm ' num2str(std(j_sliding_total_a_nodna_within)) ', sem: ' num2str(std(j_sliding_total_a_nodna_within)/sqrt( length(j_sliding_total_a_nodna_within) ))]);
