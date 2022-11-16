filename = uigetfile('*.*','Select the DNA file','multiselect','on');

if isequal(filename,0)
    disp('User selected Cancel');
else
    if iscell(filename)
        K = length(filename);
        tracenuc_a = [];
        for i = 1 : K
            load(filename{i});
            tracenuc_a = [tracenuc_a tracenuc];
        end
    else
        load(filename);
        tracenuc_a = tracenuc;
    end
    
    L = length(tracenuc_a);
    
    rise = 0.338;  % nm
    load('x_F_WLC');  % L_r: relative length, Force_WLC: WLC force
    Loo = 12688;
    Lo = Loo * 0.338;

    
     j=0;
    for i = 1 : L 
        if isfield(tracenuc_a{i}, 'converted') %&& any(strcmp(tracenuc{i}.name , goodtrace)) 
            if strcmp(tracenuc_a{i}.converted , 'yes')
                
                j = j+1;
                name{j} = tracenuc_a{i}.name;
                Times1{j} = tracenuc_a{i}.Times1;
                ForcepN1{j} = tracenuc_a{i}.ForcepN1;
                LDNAnm1{j} = tracenuc_a{i}.LDNAnm1;
                fit_Xcentering{j} = tracenuc_a{i}.fit_Xcentering;
                l0p5(j) = tracenuc_a{i}.l0p5;
                
                if isfield(tracenuc_a{i},'jbp_0p5')
                    jbp_0p5(j) = tracenuc_a{i}.jbp_0p5;
                else
                    jbp_0p5(j) = l0p5(j)./(interp1(Force_WLC, L_r, 0.5) * 0.338);
                end
                                
                t_on(j) = tracenuc_a{i}.t_on_pick;
                                  
                LDNAnm_s{j} = tracenuc_a{i}.LDNAnm_s;
                ForcepN_s{j} = tracenuc_a{i}.ForcepN_s;
                Times_s{j} = tracenuc_a{i}.Times_s;
                
                if isfield(tracenuc_a{i},'j_bp')
                    j_bp{j} = tracenuc_a{i}.j_bp;
                else
                    j_bp{j} = LDNAnm_s{j}./(interp1(Force_WLC, L_r, ForcepN_s{j}) * 0.338);
                end
                         
                ft = tracenuc_a{i}.ft;
                lt(j,:) = tracenuc_a{i}.lt;
                f_max(j) = tracenuc_a{i}.f_max;
                t_max(j) = tracenuc_a{i}.t_max;
                l_max(j) = tracenuc_a{i}.l_max;
                %peak_force{j} =  tracenuc_a{i}.peak_force;
                %peak_location{j} = tracenuc_a{i}.peak_location;
                %peak_width{j} = tracenuc_a{i}.peak_width;
                %L_bp_peak{j} = tracenuc_a{i}.L_bp_peak;
                
                %%find peak again to locate in time
                index_start = find(Times_s{j} >= t_on(j));
                index_start = index_start(1);
                tracenuc_a{i}.index_start = index_start;
                
                index_max = tracenuc_a{i}.index_max;
                window_size2 = 20;   % 30ms
                ForcepN_ss = movingmean(ForcepN_s{j}(index_start:index_max), window_size2);
                LDNAnm_ss = movingmean(LDNAnm_s{j}(index_start:index_max), window_size2);
                Times_ss = Times_s{j}(index_start:index_max);
                
                tracenuc_a{i}.ForcepN_ss = ForcepN_ss;
                tracenuc_a{i}.LDNAnm_ss = LDNAnm_ss;
                tracenuc_a{i}.Times_ss = Times_ss;
                
                indexp = peakfinder(ForcepN_ss,0.5,1,1,false,false);
                peak_force = ForcepN_ss(indexp);
                peak_location = LDNAnm_ss(indexp);
                peak_time = Times_ss(indexp);
                L_bp_peak = peak_location./(0.338 * interp1(Force_WLC, L_r, peak_force));
                
                tracenuc_a{i}.indexp = indexp;
                tracenuc_a{i}.peak_force = peak_force;
                tracenuc_a{i}.peak_location = peak_location;
                tracenuc_a{i}.peak_time = peak_time;
                tracenuc_a{i}.L_bp_peak = L_bp_peak;
                
                f_peak_max = 65;
                peak_location = peak_location(peak_force <= f_peak_max);
                peak_time = peak_time(peak_force <= f_peak_max);
                L_bp_peak = L_bp_peak(peak_force <= f_peak_max);
                peak_force = peak_force(peak_force <= f_peak_max);
                                
       
                f1 = figure;
                pos = [10 70 1100 800];
                set(f1, 'Pos', pos);  
                
                subplot(2,2,1)
                hold all
                plot(LDNAnm1{j}, ForcepN1{j},'o-');
                %plot(fit_Xcentering{j});
                plot(l0p5(j), 0.5, 'or', 'MarkerFaceColor','r');
                ylim([-0.1, 1])
                xlim([0 4000]);
                set(gca,'FontSize',13,'FontName','Calibri'); legend off;
                xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
                ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
                title({name{j}, ['Extension @ 0.5 pN: ' num2str(l0p5(j)) ' nm'] },'Interpreter','None','FontSize',12);  

                subplot(2,2,3)
                hold all
                plot(L_r * Loo *0.338, Force_WLC,':', 'Color','k');
                plot(LDNAnm_s{j}, ForcepN_s{j},'o', 'Color','r', 'MarkerFaceColor','r', 'MarkerSize',1);
                plot(LDNAnm_ss, ForcepN_ss, '-k');
                
                plot(peak_location, peak_force,'Marker','o','MarkerSize',4,'MarkerFaceColor','green','Color','k','LineStyle','none');
                axis([400 5500 0 125 ]);
                set(gca,'FontSize',13,'FontName','Calibri');
                xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
                ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
                title(['Extension @ 1,2, 3 pN: ' num2str(lt(j, 1)) ', ' num2str(lt(j, 2)) ', ' num2str(lt(j, 3)) ' nm'],'Interpreter','None','FontSize',10);        

                subplot(2,2,2)
                hold all
                plot(Times_s{j}, ForcepN_s{j},'--', 'Color',[0.5,0.5,0.5]);
                plot(Times_s{j}, ForcepN_s{j},'o', 'Color','k', 'MarkerFaceColor','k', 'MarkerSize',1);
                plot(t_max(j), f_max(j),'or','MarkerFaceColor','r');
                ylim([0 125]);
                xlim([t_on(j)+3.5 Times_s{j}(end)]);
                title(['F_max: ' num2str(f_max(j)) ' pN'],'Interpreter','None','FontSize',10);        
                set(gca,'FontSize',13,'FontName','Calibri');
                xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
                ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');


                subplot(2,2,4)
                hold all
                lm = line([min(Times_s{j}) max(Times_s{j})], [Loo Loo]);
                set(lm,'LineStyle','--','Color','b');
                lminn = line([min(Times_s{j}) max(Times_s{j})], [jbp_0p5(j) jbp_0p5(j)]);
                set(lminn,'LineStyle','--','Color','b');

                plot(Times_s{j}, j_bp{j},':.', 'Color','k', 'MarkerSize', 3, 'MarkerFaceColor','k');

                plot(peak_time, L_bp_peak,'Marker','o','MarkerSize',4,'MarkerFaceColor','green','Color','k','LineStyle','none');

                ylim([ min([0.9*jbp_0p5(j) 1.2*Loo]) max([0.9*jbp_0p5(j) 1.2*Loo]) ]);
                xlim([t_on(j)+3.5 Times_s{j}(end)]);

                set(gca,'FontSize',13,'FontName','Calibri');
                xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
                ylabel('Converted bp no. (bp)', 'FontSize',13, 'FontName', 'Calibri');


                savefig(f1,[name{j}(1:end-4) '_sum,.fig']);
                export_fig(f1,[name{j}(1:end-4) '_sum.png']);
                close(f1);
            end
        end
    end
    
    save('data_all.mat', 'tracenuc_a');
    
end