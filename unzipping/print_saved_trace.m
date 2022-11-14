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
P = length(trace);

for j = 1 : P
    
    if strcmp(trace{j}.selected ,'yes')
            name = trace{j}.name;
            Times = trace{j}.Times ;
            jindex =             trace{j}.jindex ;
            ForcepN =             trace{j}.ForcepN ;
            jindex_s =             trace{j}.jindex_s ;
            ForcepN_s =             trace{j}.ForcepN_s ;
            jindex_s_n = trace{j}.jindex_s_n;

            jindex_t = trace{j}.jindex_t ;
            forcepN_t =  trace{j}.forcepN_t ;

            
        if ~strcmp(trace{j}.type, 'DNA')
            
            j_first = trace{j}.j_first;
            
            dt = mean(diff(Times));
            N = round(1/40/dt);     % resampling to 40 Hz
            jindex_s2 = movingmean(jindex,N);
            ForcepN_s2 = movingmean(ForcepN,N);

            [xData_theo, yData_theo] = prepareCurveData( jindex_t, forcepN_t );
            yData_theo = yData_theo; 
            fitresult_theo = splinefit(xData_theo, yData_theo,round(max(xData_theo)/5));

            [xData_exp, yData_exp] = prepareCurveData( jindex_s2, ForcepN_s2 );
            fitresult_exp = splinefit( xData_exp, yData_exp, round(max(xData_exp)/3));

            j_above = jindex_s_n(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+2 & jindex_s_n >= 0 & ForcepN_s <= 55);
            f_above = ForcepN_s(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+2 & jindex_s_n >= 0 & ForcepN_s <= 55);

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
                
                j_slidingr_sf = [jrn_sliding(1) jrn_sliding(end)];
                j_sliding_total =  j_slidingr_sf(2) - j_slidingr_sf(1)+1;
                
                %j_sliding_total = 0;
            end
            
            trace{j}.j_sliding_total = j_sliding_total;
            trace{j}.j_above = j_above;
            trace{j}.f_above = f_above;
            trace{j}.j_slidingr_sf = j_slidingr_sf;
            
            comment = ['Protein; ' name];
            
            f0 = figure;
            pos = [100 100 1550 300];
            set(f0, 'Pos', pos);
            hold on;
            if any(j_slidingr_sf)
                for i = 1 : size(j_slidingr_sf,1)
                    p1 = patch([j_slidingr_sf(i,1) j_slidingr_sf(i,2) j_slidingr_sf(i,2) j_slidingr_sf(i,1) ], [0 0 80 80], [0.3 0.75 0.93]);
                    set(p1, 'LineStyle','none');
                end
            end
            plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5,'LineStyle','-');   
            plot(jindex_s_n,ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            
            plot(j_above, f_above, '.b');
           
            
            lfirst = line([j_first, j_first],[0 80]);
            set(lfirst, 'LineStyle','--');
            axis([-190 4500 10 65]);

            
            set(gca,'XGrid','on');
            xlabel('j-index (bp)', 'FontSize',13); 
            ylabel('Force (pN)', 'FontSize',13);
            title([comment ', sliding distance: ' num2str(j_sliding_total) ' bp'], 'FontSize',13, 'Interpreter', 'none');
            j_slidingr_sf = [];
            j_sliding_total = [];
        else
            comment = ['DNA; ' name];
            
            f0 = figure;
            pos = [100 100 1550 300];
            set(f0, 'Pos', pos);
            hold on;
            plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5,'LineStyle','-');   
            plot(jindex_s_n,ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

            axis([-190 4500 10 65]);

            
            set(gca,'XGrid','on');
            xlabel('j-index (bp)', 'FontSize',13); 
            ylabel('Force (pN)', 'FontSize',13);
            title(comment , 'FontSize',13, 'Interpreter', 'none');
            
            
        end
        
        % save figure
        savefig(f0,[name(1:end-4) '_plot.fig']);
        export_fig(f0,[name(1:end-4) '_plot.png']);
    end
end

save('data_all.mat', 'trace');
