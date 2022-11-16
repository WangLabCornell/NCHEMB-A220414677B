% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

function [peakforce_all, n_loop, name, dL_peak_mean, l0p5_a, fmax, ft, lt, pair_peak_width_all, peak_width_trace_all, location_break] = data_statistic_stretching_all            
                    
%filename_t='force_kkb_1p6kbafterA20_100mMNaCl_4mMMg';
filename = uigetfile('*.*','Select the DNA file','multiselect','on');
if isequal(filename,0)
    disp('User selected Cancel');
else
    %if iscell(filename)
    %    K = length(filename);
    %    tracenuc_a = [];
    %    for i = 1 : K
    %        load(filename{i});
    %        tracenuc_a = [tracenuc_a tracenuc];
    %    end
    %else
        load(filename);
        %tracenuc_a = tracenuc;
    %end
    
    %load(filename);
    
    L = length(tracenuc_a);
    
    rise = 0.338;  % nm
    load('x_F_WLC');  % L_r: relative length, Force_WLC: WLC force
    Loo = 12667;
    Lo = Loo * 0.338;
    
    
    name = {};
    LDNAnm_s = {};
    ForcepN_s = {};
    Times_s = {};

    peak_force= {};
    peak_location= {};
    peak_width= {};
    peak_width_trace_all = {};
    L_bp_peak= {};
    t_on = [];
    
    dL_bp_peak_all = [];
    peakforce_all = [];
    
    t_hold = [];
    location_break = [];
    force_break = [];
    
    %f_thr = 65;
    
    fmax = [];
    l0p5_a = [];
    dL_peak_mean = [];
    pair_peak_width_all = [];
    
    
    j=0;
    for i = 1 : L 
        if strcmp(tracenuc_a{i}.converted , 'yes')%isfield(tracenuc_a{i}, 'comments') %&& any(strcmp(tracenuc{i}.name , goodtrace)) 
            if 1%strcmp(tracenuc_a{i}.comments , 'yes')
                if 1%(-tracenuc_a{i}.L_bp_fulllength + Loo )/Loo <= 0.15 && ...
                    %abs(tracenuc_a{i}.num_outerturn - tracenuc_a{i}.no_innerturn_extension)/tracenuc_a{i}.no_innerturn_extension <= 0.2 && ...
                    %(tracenuc_a{i}.no_innerturn_extension + tracenuc_a{i}.num_outerturn)/2 >0
                    
                    j = j+1;
                    t_on(j) = tracenuc_a{i}.t_on_pick;
                   
                    name{j} = tracenuc_a{i}.name;
                    LDNAnm_s{j} = tracenuc_a{i}.LDNAnm_s;
                    ForcepN_s{j} = tracenuc_a{i}.ForcepN_s;
                    Times_s{j} = tracenuc_a{i}.Times_s;
                    
                    
                    l0p5_a(j) = tracenuc_a{i}.l0p5;
                    
                    ft = tracenuc_a{i}.ft;
                    lt(j,:) = tracenuc_a{i}.lt;
                    
                    peak_force{j} = tracenuc_a{i}.peak_force;
                    peak_location{j} = tracenuc_a{i}.peak_location;
                    
                    L_bp_peak{j} = tracenuc_a{i}.L_bp_peak;
                    
                    peak_width{j} = diff(tracenuc_a{i}.L_bp_peak);
                    
%                    location_break(j) = tracenuc_a{i}.location_break;
                    
                    fmax(j) = tracenuc_a{i}.f_max;

                    
                    %% ignore high force than 58 pN
                    f_thres = 55;
                    lbp_temp = L_bp_peak{j};
                    peak_force_temp = peak_force{j};
                    
                    lbp_temp = lbp_temp(peak_force{j}<= f_thres);
                    peak_force_temp = peak_force_temp(peak_force{j}<= f_thres);

                    peak_width_trace = diff(lbp_temp);

                    n_loop(j) = length(lbp_temp);

                    peak_width_trace_all{j} = peak_width_trace;

                    if length(peak_width_trace) >= 2
                        pair_peak_width = [peak_width_trace(1 : end-1), peak_width_trace(2 : end)];
                    else
                        pair_peak_width = [0/0 0/0];
                    end

                    pair_peak_width_all = [pair_peak_width_all; pair_peak_width];
                    
                    dL_bp_peak_all = [dL_bp_peak_all; peak_width_trace];

                    dlbp_temp = diff(lbp_temp);

                    if mean(dlbp_temp) >= 0 
                        dL_peak_mean = [dL_peak_mean; mean(dlbp_temp)];
                    end

                    peakforce_all = [peakforce_all; peak_force_temp];
                  
                    %%% only include significant forcepeak (with width = 1kb)
                    %index_major_peak = find(peak_width{j} >= 1000);
                    %n_major(j) = length(index_major_peak);
                    %F_peak_major{j} = peak_force{j}(index_major_peak);
                    %Lbp_peak_major{j} = L_bp_peak{j}(index_major_peak);
                    %peak_width_major{j} = peak_width{j}(index_major_peak);
                else
                    %tracenuc_a{i}.anchoringgeometry = 'bad';
                end
            end
        end
    end
    
    
    dL_peak_mean = dL_peak_mean(isfinite(dL_peak_mean));
    
    M = j;




    %save('data.mat','tracenuc_a');
        if M >= 1
            
            % major peaks
            %f10 = figure;
            %subplot(1,2,1)
            %hist(n_major);
            %xlabel('Number of major peaks');
            %ylabel('Count');
            %set(gca,'FontSize',15,'LineWidth',1.5, 'Box','off');
            
            %subplot(1,2,2)
            %hist(cell2mat(F_peak_major(:)));
            %xlabel('Peak force of major peaks (pN)');
            %ylabel('Count');
            %xlim([0 65]);
            %set(gca,'FontSize',15,'LineWidth',1.5, 'Box','off');
            
            % number of loops
            dn_bin = 1;
            n_bin = 0:dn_bin:100;
            n_loop_bin = histc(n_loop , n_bin);
            fl = figure;
            hold all
            bar(n_bin + dn_bin/2, n_loop_bin);
            title(['Loop number, mean: ' num2str(mean(n_loop)) ', median: ' num2str(median(n_loop))  ', SD: ' num2str(std(n_loop)) ', SEM: ' num2str(std(n_loop)/sqrt(length(n_loop)))]);

            
            %% plot all traces
            f7 = figure;
            %pos = [10 70 700 900];
            %set(f7, 'Pos', pos);
            %subplot(3,1,1:2);
            hold all
            plot(L_r * Lo  , Force_WLC,'--k');
            
            for j = 1:M
                plot(LDNAnm_s{j}(Times_s{j} >= t_on(j)), ForcepN_s{j}(Times_s{j} >= t_on(j)), '.','MarkerFaceColor',rand(1,3), 'MarkerSize',1);
            end
            axis([500 5000 0 75]);
            ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
            xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
            set(gca,'FontSize',13,'FontName','Calibri');   
            title(['Summary, N = ' num2str(M)], 'interpreter','none');
            
            savefig(f7,'All_traces.fig');
            print('All_traces','-dpng','-r0');

            %% plot innerturn release
            dxw = 25;
            xw = 0: dxw : Loo;
            xwn = xw+dxw/2;
            nw = histc(dL_bp_peak_all, xw);
            
            pw = nw / sum(nw);
            dpw = sqrt(pw .* (1- pw) ./ sqrt(sum(nw)));
            
            f8 = figure;
            subplot(2,1,1)
            hold on
            errorbar(xwn, pw, dpw, 'k', 'LineStyle', 'none') ;
            h1 = stairs(xw, pw);
            
            axis([0 2000 0 max(pw)+1.2*max(dpw)]);
            set(h1, 'Color','r', 'LineWidth',1.5);
            ylabel('Count', 'FontSize',11, 'FontName', 'Calibri');
            xlabel('Adjacent disruption size (bp)', 'FontSize',11, 'FontName', 'Calibri'); 
            set(gca,'FontSize',11,'FontName','Calibri');   
            %title(['Mean: ' num2str(mean(dL_bp_peak_all)) ' \pm ' num2str(std(dL_bp_peak_all)) ' bp']);
            title(['N = ' num2str(sum(nw)) ]);
            
            subplot(2,1,2)
            nw_mean = histc(dL_peak_mean, xw);
            
            pw_mean = nw_mean / sum(nw_mean);
            dpw_mean = sqrt(pw_mean .* (1- pw_mean) ./ sqrt(sum(nw_mean)));
            
            hold on
            h12 = stairs(xw, pw_mean);
            xlim([0 2000]);
            set(h12, 'Color','r', 'LineWidth',1.5);
            ylabel('Count', 'FontSize',11, 'FontName', 'Calibri');
            xlabel('Mean adjacent disruption size per trace (bp)', 'FontSize',11, 'FontName', 'Calibri'); 
            set(gca,'FontSize',11,'FontName','Calibri');   
            title(['N = ' num2str(sum(nw_mean)) ', Mean: ' num2str(mean(dL_peak_mean)) ' \pm ' num2str(std(dL_peak_mean)) ' bp']);
            
                        
            savefig(f8,'DNA_loop_size.fig');
            print('DNA_loop_size','-dpng','-r0');
            
            %% plot innerturn release
           % peakforce_all = peakforce_all(peakforce_all < 55);
            dfw = 4;
            fw = 0: dfw : 70;
            fwn = fw+dfw/2;
            nf = histc(peakforce_all, fw);
            f9 = figure;
            hold on
            h2 = bar(fwn, nf);
            hold all
            lthr = line([f_thres f_thres], [0 max(nf)+2]);
            set(lthr, 'LineStyle','--', 'Color','k');
            axis([0 70 0 max(nf)+2]);
            set(h2, 'FaceColor','r', 'BarWidth',1);
            ylabel('Count', 'FontSize',11, 'FontName', 'Calibri');
            xlabel('Disruption force (pN)', 'FontSize',11, 'FontName', 'Calibri'); 
            set(gca,'FontSize',11,'FontName','Calibri');   
            title(['Disruption force: ' num2str(mean(peakforce_all)) ' \pm ' num2str(std(peakforce_all)) ', sem: ' num2str(std(peakforce_all)/sqrt(length(peakforce_all))) ' pN']);
                       
            savefig(f9,'DNA_loop_disrupt_force.fig');
            print('DNA_loop_disrupt_force','-dpng','-r0');

            %subplot(3, 2, 5)
            %hold on
            %plot(t_hold, force_break,'o', 'MarkerSize', 4, 'Color','k', 'MarkerFaceColor','r');
            %xlim([0 10.5]);
            %ylim([0 75]);
            %xlabel('Hold time (s)', 'FontSize',11, 'FontName', 'Calibri');
            %ylabel('Disruption force (pN)', 'FontSize',11, 'FontName', 'Calibri'); 
            %set(gca,'FontSize',11,'FontName','Calibri');   
            %title(['Disruption force: ' num2str(mean(peakforce_all)) ' \pm ' num2str(std(peakforce_all)) ' pN']);
                       
            %force_break_s = force_break(t_hold <= 10);
            %subplot(3, 2, 6)
            %hist(force_break_s,20);
            %xlim([0 70]);
            %xlabel('Force at DSB (pN)', 'FontSize',11, 'FontName', 'Calibri');
            %ylabel('Count', 'FontSize',11, 'FontName', 'Calibri'); 
            %set(gca,'FontSize',11,'FontName','Calibri');   
            %title(['DSB force: ' num2str(mean(force_break_s)) ' \pm ' num2str(std(force_break_s)) ' pN']);
            
            
            
            %no_peak = cellfun(@numel,peak_force);
            %no_peak_s = no_peak(t_hold <= 10);
            
            %m = 1:25;
            %f_dsb_a = zeros(length(m), 2);
            %for k = 1 : length(m)
            %    f_dsb_a(m(k), :) = [mean(force_break_s(no_peak_s == m(k))) std(force_break_s(no_peak_s == m(k)))];
            %end
            
            %figure;
            %plot(no_peak_s, force_break_s,'x');
            %hold all
            %errorbar(m, f_dsb_a(:, 1), f_dsb_a(:, 2) ,'s-');
            %set(gca,'FontSize',13,'FontName','Calibri');
            %xlabel('Number of peaks identified', 'FontSize',14, 'FontName', 'Calibri');
            %ylabel('Peak force (pN)', 'FontSize',14, 'FontName', 'Calibri');
            %ylim([0 75]);
            
            
            f9 = figure;
            dfr = 2;
            fr = 0:dfr:150;
            frn = fr+dfr/2;
            nfm = histc(fmax, fr);
            pfm = nfm/sum(nfm);
            line([65, 65], [0, max(nfm)+2]);
            hold all

            bar(frn, pfm, 'BarWidth',1);
            pl65 = sum(fmax <= 65)/sum(fmax >= 0 );
            dpl65 = sqrt(pl65 * (1 - pl65)/sum(fmax >= 0 ));
            title({['< 65-pN fraction: ' num2str(pl65) ' \pm ' num2str(dpl65)], ['mean force: ' num2str(mean(fmax)) , ', sd: ' num2str(std(fmax)), ', sem: ' num2str(std(fmax)/sqrt(length(fmax))) ' pN']});
            xlim([0 130]);
            %ylim([0 max(nfm)+2]);
            ylim([0 0.5]);
            xlabel('Breakage force (pN)', 'FontSize',12, 'FontName', 'Calibri');
            ylabel('Fraction (binsize: 2 pN)', 'FontSize',12, 'FontName', 'Calibri'); 
            set(gca,'FontSize',12,'FontName','Calibri');   
            savefig(f9,'Breakage_force_l65pNfraction.fig');

            
            %l_0p5_plot = lt(:,1);
            l_0p5_plot = l0p5_a;
            f10 = figure;
            dext = 100;
            extr = 0:dext:5000;
            extn = extr+dext/2;
            next = histc(l_0p5_plot, extr);
            
            line([3000, 3000], [0, max(next)+2]);
            hold all
            bar(extn, next, 'BarWidth',1);
            
            pe3000 = sum(l_0p5_plot >= 3000)/sum(l_0p5_plot >= 0 );
            dpe3000 = sqrt(pe3000 * (1 - pe3000)/sum(l_0p5_plot >= 0 ));
            title(['> 3um 0.5-pN extension: ' num2str(pe3000) ' \pm ' num2str(dpe3000)]);
            xlim([0 4500]);
            ylim([0 max(next)+2]);
            xlabel('Extension at 0.5 pN (nm)', 'FontSize',12, 'FontName', 'Calibri');
            ylabel('Count', 'FontSize',12, 'FontName', 'Calibri'); 
            set(gca,'FontSize',12,'FontName','Calibri');   
            savefig(f10,'extension0p5_m3umfraction.fig');
        end
        
       
end
end

