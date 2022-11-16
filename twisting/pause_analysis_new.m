filename = uigetfile('*.*','Select the DNA file','multiselect','on');
load(filename);
load('x_f_WLC.mat');

P = length(trace);

 turn_step_a = [];

for i =  1 : P
    if strcmp(trace{i}.comments ,'good')
        if trace{i}.height_c_before >= 2.4 && trace{i}.F_from_var_length_p < 1
            name = trace{i}.name
            
            time_holding = trace{i}.time_holding1;
            turn_holding = trace{i}.turn_holding1;
            z_holding = trace{i}.z_holding1;

            sign_winding = mean(diff(turn_holding))>0;
            d_0 = turn_holding(1);

            turn_shift = trace{i}.turn_shift;
           
            height_0 = trace{i}.height_0;

            d_wound = turn_holding(end);
            dd_wound = 0.05;
            index_hold = abs(turn_holding-d_wound)<dd_wound;  

            % data when holding only
            time_holding = time_holding(index_hold);
            turn_holding = turn_holding(index_hold);
            z_holding = z_holding(index_hold);

            rate = 1/30;  
            dt = mean(diff(time_holding));
            window = round(1/rate / dt);
            order = 2;
            if mod(window, 2) == 0 % more smooth, used to find pausing
                            window = window+1;
            end
            z_holding_s = sgolayfilt(z_holding,order,window);  % smoothed data
            
            rate2 = 1/5;   % less smoothed
            window2 = round(1/rate2 / dt);
            if mod(window2, 2) == 0
                            window2 = window2+1;
            end
            z_holding_s2 = sgolayfilt(z_holding,order,window2);  % smoothed data
            
            
            height_0 = trace{i}.height_0;
            Turn2_mean_1way = trace{i}.turn_m_before_positive;
            Z_beadc2_mean_1way =  trace{i}.Z_m_before_positive;
            
            Turn4_mean_1way = trace{i}.turn_m_after_positive;
            Z_beadc4_mean_1way =  trace{i}.Z_m_after_positive;
            
            F_from_var_length_p = trace{i}.F_from_var_length_p;
            
            %if ~isfield(trace{i}, 'paraf')
                % Fit the hat if needed
                [paraf, height_0, turn_0, ttempt, ztempt,fitresult_height_turn_right, fitresult_height_turn_left] = fit_naked_dna(Turn2_mean_1way,Z_beadc2_mean_1way);
                trace{i}.paraf = paraf;
                trace{i}.turn_0 = turn_0;
                trace{i}.fitresult_height_turn_right = fitresult_height_turn_right;
                trace{i}.fitresult_height_turn_left = fitresult_height_turn_left;
                trace{i}.ttempt = ttempt;
                trace{i}.ztempt = ztempt;
            %else
            
            %    paraf = trace{i}.paraf;
            %    fitresult_height_turn_left = trace{i}.fitresult_height_turn_left;
            %    fitresult_height_turn_right = trace{i}.fitresult_height_turn_right;
                
            %end
            

            if sign_winding
                        turn_convert = fitresult_height_turn_right(paraf, z_holding_s);
                        turn_relaxed = (turn_holding - d_0) - turn_convert; % the number of turn relaxed by topo, always positive          
                        
                        turn_convert2 = fitresult_height_turn_right(paraf, z_holding_s2);
                        turn_relaxed2 = (turn_holding - d_0) - turn_convert2;
            else
                        turn_convert = fitresult_height_turn_left(paraf, z_holding_s); 
                        turn_relaxed = -(turn_holding - d_0) + turn_convert; % the number of turn relaxed by topo, always positive
                        
                        turn_convert2 = fitresult_height_turn_left(paraf, z_holding_s2); 
                        turn_relaxed2 = -(turn_holding - d_0) + turn_convert2;
            end


            % Region with activity
            index_notactive = find(z_holding_s > 0.9 * height_0); % below 90% of maximum height
            if any(index_notactive)
                index_notactive = index_notactive(1);
            else
                index_notactive = length(turn_relaxed);
            end

            time_holding = time_holding(1:index_notactive);
            turn_holding = turn_holding(1:index_notactive);
            z_holding = z_holding(1:index_notactive);
            z_holding_s = z_holding_s(1:index_notactive);
            turn_convert = turn_convert(1:index_notactive);
            turn_convert2 = turn_convert2(1:index_notactive);
            turn_relaxed = turn_relaxed(1:index_notactive);
            turn_relaxed2 = turn_relaxed2(1:index_notactive);
            
            trace{i}.time_holding_new = time_holding;
            trace{i}.turn_holding_new = turn_holding;
            trace{i}.z_holding_new = z_holding;
            trace{i}.z_holding_s_new = z_holding_s;
            trace{i}.turn_convert_new = turn_convert;
            trace{i}.turn_convert2_new = turn_convert2;
            trace{i}.turn_relaxed_new = turn_relaxed;
            trace{i}.turn_relaxed2_new = turn_relaxed2;

            % wait-time analysis, pause analysis
            t_thres = 2;
            dturn = 0.2;
            turn_r = min(turn_relaxed) : dturn : max(turn_relaxed); turn_r = turn_r(:);
            turn_rn = turn_r + dt/2;
            n_r = histc(turn_relaxed, turn_r);

            %%% Caveat, finding the n_thres that balance long and short dwell???
            n_thres = t_thres/dt; % 1.5 s
            %n_thres = 120;

            % finding cluster
            index_bin = find(n_r >= n_thres);  % location of high enough bin in 'turn_rn'

            index_pause_t = {}; 
            index_cluster = {}; turn_rn_cluster = {}; n_r_cluster = {}; t_pause = []; turn_mean = []; dturn_mean = []; t_pause_u = [];
            index_pause_t_mm= []; turn_mean_u = [];
            
            if any(index_bin)
                index_edge = find(diff(index_bin) >= 2);
                i_end = [index_bin(index_edge); index_bin(end)];  % end point of a cluster
                i_start = [index_bin(1); index_bin(index_edge+1)];  % start point of a cluster

                m1 = length(i_end); 
                m2 = length(i_start); 

                if m1 == m2
                    for j = 1 : m1
                        index_cluster{j} = i_start(j):i_end(j);
                        turn_rn_cluster{j} = turn_rn(index_cluster{j});
                        n_r_cluster{j} = n_r(index_cluster{j});

                        % waittime at pause, % mean turn at pause with sd
                        t_pause(j) = sum(n_r_cluster{j})*dt;
                        turn_mean(j) = sum(turn_rn_cluster{j} .* n_r_cluster{j}) / sum(n_r_cluster{j});
                        dturn_mean(j) = sqrt(sum(turn_rn_cluster{j}.^2 .* n_r_cluster{j}) / sum(n_r_cluster{j}) - turn_mean(j)^2);

                    end

                    index_pause_long = t_pause >= 2;%t_thres;
                    turn_mean = turn_mean(index_pause_long);
                    dturn_mean = dturn_mean(index_pause_long);
                    t_pause = t_pause(index_pause_long);
                    
                    %%% joining regions with turn relaxed less than 2 turns
                    %turn_mean_temp = turn_mean;
                    %index_low_turn = turn_mean < 2;
                    %turn_mean_low = mean(turn_mean_temp(index_low_turn));
                    %turn_mean_temp(index_low_turn) = [];
                    %turn_mean = [turn_mean_low turn_mean_temp];

                    % map back to time space
                    M = length(turn_mean);
                    for j = 1 : M
                        %index_pause_t{j} =  find(turn_relaxed >= turn_mean(j)-dturn_mean(j) & turn_relaxed < turn_mean(j)+dturn_mean(j));
                        index_pause_t{j} =  find( abs(turn_relaxed - turn_mean(j)) <= 0.4);
                        index_pause_t_mm(j,:) = [min(index_pause_t{j})  max(index_pause_t{j})];
                    end
                    [~,index_sort ] = sort(index_pause_t_mm(:,1),1);
                    index_pause_t_mm = index_pause_t_mm(index_sort, :);
                    turn_mean = turn_mean(index_sort); dturn_mean = dturn_mean(index_sort); t_pause = t_pause(index_sort);
                                       
                    %%Join too close region
                    %turn_mean
                    %% Join overlapping region
                    dm = -1;
                    while dm < 0
                        M = size(index_pause_t_mm,1);
                        if M>=2
                            index_pause_t_mm_new = [];
                            j = 1;
                            while j <= M
                                index = find(index_pause_t_mm(:,1) < index_pause_t_mm(j,2) & index_pause_t_mm(:,1) > index_pause_t_mm(j,1));
                                if ~isempty(index)
                                    index_pause_t_mm_new = [index_pause_t_mm_new; [index_pause_t_mm(j,1) max(index_pause_t_mm(index,2))]];
                                    j = index(end)+1;
                                else
                                    index_pause_t_mm_new = [index_pause_t_mm_new; [index_pause_t_mm(j,1) index_pause_t_mm(j,2)]];
                                    j = j+1;
                                end
                            end
                        else
                            index_pause_t_mm_new = index_pause_t_mm;
                        end

                        dm = size(index_pause_t_mm_new,1) - M;
                        index_pause_t_mm = index_pause_t_mm_new;    
                    end
                    

                    Mp = size(index_pause_t_mm_new,1);
                    
                    for j = 1 : Mp
                        index_pause_t{j} = index_pause_t_mm_new(j,1): index_pause_t_mm_new(j,2);
                        t_pause_u(j) = time_holding(index_pause_t{j}(end)) - time_holding(index_pause_t{j}(1));
                        turn_mean_u(j) = mean( turn_relaxed(index_pause_t{j}));
                    end

                    turn_step_a = [turn_step_a diff(turn_mean_u)];
                else
                    disp('Error check!');
                end
            else
                Mp = 0;
            end
            
            %% assigning first pausing
            dT = 10; dM = round(dT/dt);
            [xData, yData] = prepareCurveData( time_holding(1:dM), turn_relaxed2(1:dM) );
            ft = fittype( 'poly1' );
            [fitresult, ~] = fit( xData, yData, ft );
            turn_relaxed0 = fitresult(time_holding(1));
            
            d_turnshift_threshold = 1.5;
            
            if Mp >= 1               
                if abs(turn_mean_u(1) -  turn_relaxed0) <= 3 || turn_mean_u(1)<= 3
                    waittime = time_holding(index_pause_t{1}(end)) - time_holding(1);

                    if Mp>= 2
                        index_burst = index_pause_t{1}(end)+1 : index_pause_t{2}-1;
                        pause_first  = t_pause_u(2);
                        turn_first = turn_mean_u(2) - turn_mean_u(1);
                    else
                        
                        pause_first  = 0/0;
                        if abs(turn_shift) >= d_turnshift_threshold
                            index_burst = index_pause_t{1}(end)+1 : length(turn_relaxed);
                            %turn_first = turn_relaxed(end) - turn_mean_u(1); % or 
                            turn_first = 30;
                        else
                            index_burst = [];
                            turn_first = 0/0;
                        end
            
                    end   
                else
                    waittime = 0;
                    index_burst = 1 : index_pause_t{1}-1;
                    pause_first  = t_pause_u(1);
                    turn_first = turn_mean_u(1) - turn_relaxed(1);
                end
            else
                 waittime = 0;
                 index_burst = 1 : length(turn_relaxed);
                 pause_first  = 0/0;
                 %turn_first = turn_relaxed(end) - turn_relaxed(1); % or 
                 turn_first = 30;
            end
            
            %% find the burst speed
            window_slope = window2;
            if any(index_burst) && abs(turn_shift) >= d_turnshift_threshold
                x = time_holding(index_burst);
                y = turn_relaxed2(index_burst);

                if length(y) > window_slope
                    v = movingslope(y, window_slope,1,dt); 
                else
                    v = movingslope(y, 2,1,dt);
                end
                
                dtu = 0.5;
                tu = min(y):dtu:max(y);
                tun = tu + dtu/2; tun(end) = [];
                vu = zeros(1, length(tun));
                for j = 1 : length(tu)-1
                   vu(j) = mean(v( y >= tu(j) & y < tu(j+1)));
                end
                if ~isempty(vu(isfinite(vu)))
                    v_burst = median(vu(isfinite(vu)));
                else
                    v_burst = mean(v);
                end
            
            else
                dtu = 0.5;
                v = 0/0;
                x = 0/0;
                y = 0/0;
                vu = 0/0;
                v_burst = 0/0;
            end
            
            trace{i}.index_pause_t = index_pause_t;
            trace{i}.t_pause_u = t_pause_u;
            trace{i}.turn_mean_u = turn_mean_u;
            trace{i}.d_turnshift_threshold = d_turnshift_threshold;
            trace{i}.waittime = waittime;
            trace{i}.index_burst = index_burst;
            trace{i}.pause_first = pause_first;
            trace{i}.turn_first = turn_first;
            trace{i}.turn_first = turn_first;
            trace{i}.v_burst = v_burst;
            trace{i}.v_burst_instanteneous = v;

                f1 = figure;
                pos = [10 70 1200 500];
                set(f1, 'Pos', pos);
                
                subplot(2,5,[1 2 6 7]);
                hold all
                line([0 0],[ -2 15],'LineStyle','--','LineWidth',1, 'Color','k');
                line([-55 150],[ height_0*0.9 height_0*0.9],'LineStyle','--','LineWidth',0.5, 'Color','k');
                line([d_wound d_wound],[ -2 15],'LineStyle','--','LineWidth',1, 'Color',[0.87, 0.49, 0]);
                plot( ttempt, ztempt,'-c','LineWidth',2);
                plot(Turn4_mean_1way, Z_beadc4_mean_1way, '-o','MarkerSize',2,'LineWidth',1,'MarkerFaceColor',[0.5, 0.5, 0.5],'Color',[0.5, 0.5, 0.5]);
                plot(turn_holding, movingmean(z_holding,10), '.','Color',[0.87, 0.49, 0]);
                plot(Turn2_mean_1way, Z_beadc2_mean_1way, '-ok','MarkerSize',4,'LineWidth',1.5,'MarkerFaceColor','k');

                force_topo = num2str(F_from_var_length_p);
                turn_shift_topo = num2str(turn_shift);
                height_0_topo = num2str(height_0);
                ylabel('Z (\mum)','FontSize',12);
                xlabel('Turn added','FontSize',12);
                set(gca,'FontSize',12,'FontName','Calibri');
                hmm = 1.1*max(Z_beadc2_mean_1way);
                xlim([-50 Turn4_mean_1way(end)+10]); ylim([0 hmm+0.2]);
                title([name ', F = '  force_topo 'pN, Turn shift: ' turn_shift_topo ' turns, h_max: ' height_0_topo ' um'],'Interpreter','None','FontSize',10,'FontWeight','bold');
                
                subplot(2,5,5)
                plot(x,y,'-k', 'LineWidth',1.5);
                xlabel('Time (s)'); ylabel('# of turn relaxed');
                
                subplot(2, 5,10)
                hist(vu,20);
                xlabel(['Resampled speed (tps), bin: ' num2str(dtu) ' turns']); ylabel('Count');
                title(['Burst speed: ' num2str(v_burst) ' tps']);
                
                subplot(2,5,3:4)
                hold all
                plot(time_holding, turn_relaxed,'LineWidth',1.5, 'Color',[0.5, 0.5, 0.5]);
                if ~isempty(index_pause_t)
                    for j = 1 : M
                        plot(time_holding(index_pause_t{j}), turn_relaxed(index_pause_t{j}), '.b-','LineWidth',2.5,'MarkerSize',2.5);
                    end
                end
                xlabel('Time (s)');  ylabel('Turn relaxed'); 
                title({name, ['Turn @ pausing: ' num2str(turn_mean_u)] , ['First duration: ' num2str(pause_first) ' s, # turn to pausing: ' num2str(turn_first)]});

                subplot(2,5,8:9)
                hold all
                plot(time_holding, z_holding,'.','MarkerSize',1.5, 'Color',[0.5, 0.5, 0.5]);
                if ~isempty(index_pause_t)

                for j = 1 : M
                    plot(time_holding(index_pause_t{j}), z_holding_s(index_pause_t{j}), '.b-','LineWidth',2.5,'MarkerSize',2.5);
                end
                end
                title( ['Pause duration (s): ' num2str(t_pause_u) ]);
                xlabel('Time (s)');  ylabel('Z (um)'); 
                
                savefig(f1,[name '_' num2str(i) '_pause_analysis.fig']);
                export_fig(f1,[name '_' num2str(i) '_pause_analysis.png']);
                

        end
    end
end

%save(filename,'trace');
save('dat_temp.mat','trace');

f2 = figure;
dt = 0.2; t = -0.5:dt: 20.5;
n = histc(turn_step_a, t);
tn = t+dt/2;
bar(tn,n);
xlabel('Turn between pausing (turns)');
ylabel('Count');
savefig(f2,'Step_between_pauses.fig');









