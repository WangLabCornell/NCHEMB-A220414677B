% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

filename = uigetfile('*.*','Select the DNA file','multiselect','on');

if isequal(filename,0)
    disp('User selected Cancel');
else
    if iscell(filename)
        K = length(filename);
        trace_a = [];
        for i = 1 : K
            load(filename{i});
            trace_a = [trace_a trace];
        end
    else
       load(filename);
        trace_a = trace;
    end
end

trace = trace_a;
%load(filename);
P = length(trace);

j= 0;

lbp_thres = 50;

j_sliding_total_a = [];
for i = 1 : P
    if strcmp(trace{i}.selected ,'yes')
        
           j = j+1;
           
           name_a{j} = trace{i}.name;
           
           j_sliding_total_a(j) = trace{i}.j_sliding_total;
           
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

jmax_theory = jindex_t{1}(end);

index_notdna = [index_protein index_notallign];

j_first_notdna = j_first(index_notdna);
f_first_notdna = f_first(index_notdna);
fmax_notdna = fmax(index_notdna);
fmean_notdna = fmean(index_notdna);
jend_notdna = jend(index_notdna);
j_fmax_notdna = j_fmax(index_notdna);

j_sliding_total_a_nodna = j_sliding_total_a(index_notdna);
%% Sliding distance for traces with boudn protein
j_sliding_total_a_nodna_within = j_sliding_total_a_nodna(j_first_notdna >= lbp_thres);
fs= figure;
dxs = 100;
x_slide = 0:dxs : 4000;
n_slide = histc(j_sliding_total_a_nodna_within, x_slide);
x_sliden = x_slide + dxs/2;
bar(x_sliden, n_slide/sum(n_slide));
xlabel('Total sliding distance (bp)','FontSize',14);
ylabel('Fraction','FontSize',14);
xlim([0 4000]);
title(['For J >= ' num2str(lbp_thres) 'bp, N = ' num2str(length(j_sliding_total_a_nodna_within )) '; mean: ' num2str(mean(j_sliding_total_a_nodna_within)) ' \pm ' num2str(std(j_sliding_total_a_nodna_within)) ', sem:' num2str(std(j_sliding_total_a_nodna_within)/sqrt(length(j_sliding_total_a_nodna_within)))]);


%% First binding location vs. max force
f3 = figure;
plot(j_first_notdna, fmax_notdna, 'ow', 'MarkerSize',6, 'MarkerFaceColor','k', 'LineWidth',0.5);
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('J-index of first binding (bp)','FontSize',14);
ylabel('Max Force of the trace (pN)','FontSize',14);
xlim([-100 4000]);
title(['All, N = ' num2str(length(j_first_notdna)) '; J >= ' num2str(lbp_thres) 'bp, N = ' num2str(length(j_first_notdna(j_first_notdna>=lbp_thres)))]);

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

%% unzipping distance from first binding for binding over lbp_thres bp
dj_withindna = jend_notdna((j_first_notdna >=lbp_thres)) - j_first_notdna((j_first_notdna >=lbp_thres));
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
fmax_notdna_withinDNA = fmax_notdna(j_first_notdna >=lbp_thres);
fr = 0:4:120;
frn = fr+2;
n_fmax = histc(fmax_notdna, fr);
n_fmax_withinDNA = histc(fmax_notdna_withinDNA, fr);
f6 = figure;
hold all
b1 = bar(frn, n_fmax/sum(n_fmax), 'BarWidth',1, 'FaceColor', [0.5, 0.5, 0.5]);
b2 = bar(frn, n_fmax_withinDNA/sum(n_fmax_withinDNA), 'BarWidth',1, 'FaceColor', 'r');
set(gca,'FontSize',14,'LineWidth',1.5, 'Box','off');
xlabel('Max Force (pN)','FontSize',14);
ylabel('Fracion','FontSize',14);
title({['All: N = ' num2str(length(fmax_notdna))], [ 'With J-index >= ' num2str(lbp_thres) ' bp: N = ' num2str(length(fmax_notdna_withinDNA)) '; mean:' num2str(mean(fmax_notdna_withinDNA)) ', sd:' num2str(std(fmax_notdna_withinDNA)) ', sem:' num2str(std(fmax_notdna_withinDNA)/sqrt(length(fmax_notdna)))]});

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
