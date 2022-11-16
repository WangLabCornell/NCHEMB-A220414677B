% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

function analyze_topo_dna_with_pause_analysis
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Read data and input parameter
R_bead = 0.5;
kbt = 4.08; % pN.nm
rise = 0.338; % nm
Lo13p7 = 13700;    
f_e = 0.5; %pN



%load('data_IIA_1pM_noETO_pos.mat');
%j = 0;
%trace1 = trace;
%for i = 1 : length(trace1)
%    if strcmp(trace1{i}.comments, 'good')
%        j = j+1;
%        trace_selected{j} = trace1{i}.name;
%    end
%end

%trace = [];

load('x_f_WLC.mat');
[Times1, Turn1, Box_name1, X_bead1, Y_bead1, Z_bead1, L1, Z_piezo1] = get_mtdata4('MT data with MT DAQ data_forcecalibration.txt');
[Timesg, Turng, Box_nameg, X_beadg, Y_beadg, Z_beadg, Lg, Z_piezog] = get_mtdata4('MT data with MT DAQ data_geometrycorrection.txt');
[Times2, Turn2, Box_name2, X_bead2, Y_bead2, Z_bead2, L2, Z_piezo2] = get_mtdata4('MT data with MT DAQ data_initialHAT.txt');

[Times4, Turn4, Box_name4, X_bead4, Y_bead4, Z_bead4, L4, Z_piezo4] = get_mtdata4('MT data with MT DAQ data_finalHAT.txt');

% first winding
[Times3, Turn3, Box_name3, X_bead3, Y_bead3, Z_bead3, L3, Z_piezo3] = get_mtdata4('MT data with MT DAQ data_topoholding.txt');

% drift correction
% beadstuck name: 253

j = 0;
for i = 1 : L1
    if strcmp(Box_name1{i}, 'Box1.1 X (px)') || strcmp(Box_name1{i}, 'Box2.1 X (px)')  || strcmp(Box_name1{i}, 'Box4.1 X (px)')  || strcmp(Box_name1{i}, 'Box114.1 X (px)')  %1, 2, 4, 114
        j = j+1;
        index_stuck(j) = i;
    end
end
index_stuck

Z_beadg = Z_beadg - repmat(mean(Z_beadg(:, index_stuck),2),1,Lg);
Z_bead1 = Z_bead1 - repmat(mean(Z_bead1(:, index_stuck),2),1,L1);
Z_bead2 = Z_bead2 - repmat(mean(Z_bead2(:, index_stuck),2),1,L2);
Z_bead3 = Z_bead3 - repmat(mean(Z_bead3(:, index_stuck),2),1,L3);
Z_bead4 = Z_bead4 - repmat(mean(Z_bead4(:, index_stuck),2),1,L4);

%% Part 2: Select TC tethers, surface finding, geometry correction, find B-field

%% calculate z_variance during second twisting to identify non-tc tether
z_variance = mean(Z_bead2.^2,1) - mean(Z_bead2,1).^2;       
%% find tether with smaller variance and assign them as tether with chromatin
index_chromatin =  z_variance>0.02 ;% & z_variance<0.0850*(0.8*12644*0.338/1000)^2 ;
L = length(find(index_chromatin));
Box_namec1 = Box_name1(index_chromatin); Box_namec2 = Box_name2(index_chromatin); Box_namec3 = Box_name3(index_chromatin); Box_namec4 = Box_name4(index_chromatin);

X_beadc1 = X_bead1(:,index_chromatin); X_beadc2 = X_bead2(:,index_chromatin); X_beadc3 = X_bead3(:,index_chromatin); X_beadc4 = X_bead4(:,index_chromatin);
Y_beadc1 = Y_bead1(:,index_chromatin); Y_beadc2 = Y_bead2(:,index_chromatin); Y_beadc3 = Y_bead3(:,index_chromatin); Y_beadc4 = Y_bead4(:,index_chromatin);
Z_beadc1 = Z_bead1(:,index_chromatin); Z_beadc2 = Z_bead2(:,index_chromatin); Z_beadc3 = Z_bead3(:,index_chromatin); Z_beadc4 = Z_bead4(:,index_chromatin);


%% find the surface, applied geometry correction from the stuck chromatin
surface = mean(Z_beadc4(Times4>Times4(end)-5-5 & Times4<Times4(end)-5,:),1); % choose turn where the chromatin were clearly stuck
%surface = mean(Z_beadc4(Turn4>190 & Turn4 < 195),1); % choose turn where the chromatin were clearly stuck
Z_beadc1 = Z_beadc1 - repmat(surface,size(Z_beadc1,1),1);
Z_beadc2 = Z_beadc2 - repmat(surface,size(Z_beadc2,1),1);
Z_beadc3 = Z_beadc3 - repmat(surface,size(Z_beadc3,1),1);
Z_beadc4 = Z_beadc4 - repmat(surface,size(Z_beadc4,1),1);


Z_beadc2 = z_correct(Times3, Z_beadc2, L);
Z_beadc3 = z_correct(Times3, Z_beadc3, L);
Z_beadc4 = z_correct(Times3, Z_beadc4, L);

%% find anchoring pos
pos_tether = zeros(L,2); 
%% dealing with drift
for j = 1 : L
        [xData, yData] = prepareCurveData( Times1, X_beadc1(:,j));
        ft1 = fittype( 'poly1' );
        [fitresult1, ~] = fit( xData, yData, ft1 );
        [xData, yData] = prepareCurveData( Times1, Y_beadc1(:,j));
        ft2 = fittype( 'poly1' );
        [fitresult2, ~] = fit( xData, yData, ft2 );
        pos_tether(j, :) =  [fitresult1(Times1(end)), fitresult2(Times1(end))];
end
               
%% find relative anchoring point
% partition the data according to turn number for the twisting in step 2 (Geometry correction)
k=1;
time_patg={}; turn_patg={};index_patg={};
time_patg{1} = Timesg(1);
turn_patg{1} = Turng(1);
index_patg{1} = 1;
for i = 1 : length(Turng)-1
        if abs(Turng(i+1)- turn_patg{k}(end))<0.01
            time_patg{k} = [time_patg{k} Timesg(i+1)];   
            turn_patg{k} = [turn_patg{k} Turng(i+1)];   
            index_patg{k} = [index_patg{k} i+1];   
        else
            if Turng(i+1)-turn_patg{k}(end) >= 0.12
                k = k+1;
                time_patg{k} = Timesg(i+1);
                turn_patg{k} = Turng(i+1);
                index_patg{k} = i+1;
            end
        end
end
    
K = length(index_patg);
colorK = rand(K,3);    
n1 = length(time_patg);
Turng_mean  = zeros(1,n1);
X_beadcg_mean  = [];
Y_beadcg_mean  = []; 
X_beadcg = X_beadg(:,index_chromatin);
Y_beadcg = Y_beadg(:,index_chromatin);
for i = 1 : n1
        Turng_mean(i) = mean(turn_patg{i}); 
        X_beadcg_mean(i,:) = mean(X_beadcg(index_patg{i},:),1);
        Y_beadcg_mean(i,:) = mean(Y_beadcg(index_patg{i},:),1);
end

%% fit circle to find center_rotation usign Kasa algorithm
% source: https://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit
% more accurate algorithm: Pratt and Taubin circle fit
center_rotation = zeros(size(X_beadcg_mean,2), 2);
Rfit = zeros(size(X_beadcg_mean,2),1);
for i = 1 : size(X_beadcg_mean,2)
       a=[X_beadcg_mean(:,i) Y_beadcg_mean(:,i) ones(size(X_beadcg_mean(:,i)))]\(-(X_beadcg_mean(:,i).^2+Y_beadcg_mean(:,i).^2));
       center_rotation(i,:) = [-.5*a(1), -.5*a(2)];
       Rfit(i)  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end
pos_anchoring_relative_center = center_rotation - pos_tether;        % vector pointing from bead center to the anchoring point
radius = sqrt(pos_anchoring_relative_center(:,1).^2 + pos_anchoring_relative_center(:,2).^2);
height_correction = R_bead -  heaviside(R_bead-radius) .* sqrt(R_bead^2 - radius.^2);
%% correction to height 
Z_beadc1 = Z_beadc1+ repmat(height_correction',size(Z_beadc1,1),1);
Z_beadc2 = Z_beadc2+ repmat(height_correction',size(Z_beadc2,1),1);
Z_beadc3 = Z_beadc3+ repmat(height_correction',size(Z_beadc3,1),1);
Z_beadc4 = Z_beadc4+ repmat(height_correction',size(Z_beadc4,1),1);

%% finding B field
x_anchor = pos_anchoring_relative_center(:,1);
y_anchor = pos_anchoring_relative_center(:,2);

index_anchor = abs(x_anchor) <=1 & abs(y_anchor) <=1;
x_anchor = x_anchor(index_anchor);
y_anchor = y_anchor(index_anchor);

[xData, yData] = prepareCurveData( x_anchor, y_anchor );
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
[fitresult_Bfield, ~] = fit( xData, yData, ft, opts );
uB = [1, fitresult_Bfield.p1]/norm([1, fitresult_Bfield.p1]); % direction of B field
  
f1 = figure;
plot(0, 0,'xk','MarkerSize',10, 'LineWidth',5);
hold all
xas = -0.8:0.01:0.8;
yas = fitresult_Bfield(xas);
plot(xas, yas, '-k','LineWidth',1.5);
hold all
for i = 1 : L
        plot([0 pos_anchoring_relative_center(i,1)], [0, pos_anchoring_relative_center(i,2)],'LineWidth',1.5, 'Color', 'r', 'Marker','v', 'MarkerFaceColor','w');
end
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
axis equal
xlabel('\DeltaX(\mum)');
ylabel('\DeltaY(\mum)');
set(gca,'FontSize',12,'FontName','Calibri');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4: Average before and after data during the pulse
%% average data during the pulse
%% partition the before topo data

% before hat
[Turn2_mean_1way, Z_beadc2_mean_1way, Turn2_mean_1way_negative, Z_beadc2_mean_1way_negative] = average_at_pause(Times2, Turn2, Z_beadc2);
% final hat
[Turn4_mean_1way, Z_beadc4_mean_1way, Turn4_mean_1way_negative, Z_beadc4_mean_1way_negative] = average_at_pause(Times4, Turn4, Z_beadc4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5: Run through each trace, manually select good trace, and identify HAT shape and save statistics
%% run through each trace and save data
uB_var1 = zeros(1,L);
uB_var2 = zeros(1,L);
tv = 60;
d_wound1 = Turn3(end);

trace = cell(1,L);
for j = 1:L
    
    %if 1 %any(strcmp(trace_selected, Box_namec1{j})) %1%
    if 1%strcmp(Box_namec1{j}, 'Box228.0 X (px)')  %|| strcmp(Box_namec1{j}, 'Box78.0 X (px)')  || strcmp(Box_namec1{j}, 'Box179.0 X (px)') 
        trace{j}.name = Box_namec1{j};
        ftemp = figure;
        hold all   
        
        %indexf = abs(Times3 - 524) < 5; % |  abs(Times3 - 126) < 4 ; 
        %Z_beadc3(indexf, j) =  interp1(Times3(~indexf), Z_beadc3(~indexf, j), Times3(indexf));
        
        line([0 0],[ -2 15],'LineStyle','--','LineWidth',2, 'Color','k');
        line([d_wound1 d_wound1],[ -2 15],'LineStyle','--','LineWidth',2, 'Color',[0.87, 0.49, 0]);
                
        plot(movingmean(Turn3,10), movingmean(Z_beadc3(:,j),10), '-','LineWidth',2,'Color',[0.87, 0.49, 0]);    
                
        plot(Turn2_mean_1way, Z_beadc2_mean_1way(:,j), '-ok','MarkerSize',5,'LineWidth',2,'MarkerFaceColor','k');
        
        plot(Turn4_mean_1way, Z_beadc4_mean_1way(:,j), '-or','MarkerSize',5,'LineWidth',1.5,'Color',[1, 0.6, 0.78]);
       
        title(Box_namec1{j} ,'Interpreter','None','FontSize',12);
        ylabel('Z (\mum)','FontSize',12);
        xlabel('Turn','FontSize',12);
        set(gca,'FontSize',12,'FontName','Calibri');
        hmm = 2*max(Z_beadc2_mean_1way(:,j));
        xlim([-50 140]); ylim([min([-0.5 hmm+0.2]), max([-0.5 hmm+0.2]) ]);
        grid on;

        qstring = 'Does the trace look good?';
        choice = questdlg(qstring,'Question?','Yes','No','Yes');

        if strcmp(choice, 'Yes') 
            close(ftemp);

            %% Force calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            [height_c_before, F_from_var_length, F_from_var_length_p] = force_calculation(Times1, X_beadc1(:, j), Y_beadc1(:,j), Z_beadc1(:,j), uB, height_correction(j), tv);
            %trace{j}.height_c_before = height_c_before;

            % initial HAT
            trace{j}.turn_m_before_positive = Turn2_mean_1way;
            trace{j}.Z_m_before_positive = Z_beadc2_mean_1way(:,j);
            trace{j}.turn_m_before_negative = Turn2_mean_1way_negative;
            trace{j}.Z_m_before_negative = Z_beadc2_mean_1way_negative(:,j);
            trace{j}.height_c_before = max(Z_beadc2_mean_1way(:,j)) ;

            % after hat
            trace{j}.turn_m_after_positive = Turn4_mean_1way;
            trace{j}.Z_m_after_positive = Z_beadc4_mean_1way(:,j);
            trace{j}.turn_m_after_negative = Turn4_mean_1way_negative;
            trace{j}.Z_m_after_negative = Z_beadc4_mean_1way_negative(:,j);
            trace{j}.height_c_after = max(Z_beadc4_mean_1way(:,j)) ;

            % force 
            trace{j}.radius = radius(j);
            trace{j}.F_from_var_length = F_from_var_length;
            trace{j}.F_from_var_length_p = F_from_var_length_p;

            % topo activity, first winding
            trace{j}.time_holding1 = Times3; time_holding1 = Times3;
            trace{j}.turn_holding1 = Turn3; turn_holding1 = Turn3;
            trace{j}.z_holding1 = Z_beadc3(:,j); z_holding1 = Z_beadc3(:,j);
            
            %% Fit the hat, deduce height-turn state conversion on the right and the left sides
            [paraf, height_0, turn_0, ttempt, ztempt,fitresult_height_turn_right, fitresult_height_turn_left] = fit_naked_dna(Turn2_mean_1way,Z_beadc2_mean_1way(:,j));

            % Height at (+) and (-) buckling 
            slope_left = abs(paraf(7) ); % um/turn
            slope_right = abs(paraf(6) ); % um/turn
            trace{j}.slope_left = slope_left;
            trace{j}.slope_right = slope_right;
            trace{j}.height_0 = height_0;
            
            trace{j}.fitresult_height_turn_right = fitresult_height_turn_right;
            trace{j}.fitresult_height_turn_left = fitresult_height_turn_left;
            
            trace{j}.paraf = paraf;

            %% Analyze topo activity, deduce wait time and average rate
            
           % [v_waitinclude_mean, Z_waitinclude, time_waitinclude, Z_topo_pause_mean_1, v_mean_1, Z_pausefree_1, v_pausefree_1, dwelltime_pause_1, sign_winding_1, t_on_topo_1, Z_holding_s1_1, 
           %     Z_holding_s2_1, v_1, time_1, Z_1, turn_1, index_pause_time_1, pause_density_1, distance_between_pauses_1, v_pausefree_mean_1, index_hold_1, index_no_topo_hold_1] = ...
           %     analayze_topo_activity(time_holding1,turn_holding1,z_holding1, trace{j}.height_c_before);
            
            [v_waitinclude_mean, v_waitinclude_tps_mean, Z_waitinclude, time_waitinclude, turn_re_waitinclude, turn_waitinclude, Z_topo_pause_mean_1, ...
    Z_pausefree_1, v_pausefree_1, v_pausefree_tps_1, dwelltime_pause_1, sign_winding_1, t_on_topo_1, ...
    Z_holding_s1_1, Z_holding_s2_1, turn_relaxed_1, v_1, v_mean_1, v_tps_1, v_tps_mean_1, time_1, Z_1, turn_1, turn_re_1, index_pause_time_1, pause_density_1, ...
    distance_between_pauses_1, distance_between_pauses_turn_1, v_pausefree_mean_1, v_pausefree_tps_mean_1, turn_re_pausefree_1, index_hold_1, index_no_topo_hold_1] = ...
    analayze_topo_activity(time_holding1,turn_holding1, z_holding1, trace{j}.height_c_before, paraf, fitresult_height_turn_right, fitresult_height_turn_left);

           
            
            trace{j}.t_on_topo = t_on_topo_1 ;
            
            trace{j}.sign_winding =  sign_winding_1 ;
            
            trace{j}.dwelltime_pause =  dwelltime_pause_1 ;
            
            trace{j}.index_hold_1 = index_hold_1;
            trace{j}.index_no_topo_hold_1 = index_no_topo_hold_1;
            
            trace{j}.v = v_1;
            trace{j}.v_tps = v_tps_1;
            
            trace{j}.v_mean = v_mean_1;
            
            trace{j}.v_mean_tps = v_tps_mean_1 ;
            
            trace{j}.time = time_1; trace{j}.Z = Z_1; trace{j}.turn = turn_1; trace{j}.turn_re = turn_re_1;
            
            
            trace{j}.Z_holding_s1 = Z_holding_s1_1;
            trace{j}.turn_relaxed = turn_relaxed_1;
            
            
            trace{j}.Z_pausefree = Z_pausefree_1;
            trace{j}.Z_topo_pause_mean = Z_topo_pause_mean_1;
            
            t_pause_threshold = 10;
            num_pause_buckling = sum(Z_topo_pause_mean_1 <=  trace{j}.height_c_before * 0.9  & dwelltime_pause_1 >= t_pause_threshold);
            trace{j}.num_pause_buckling = num_pause_buckling;
            
            trace{j}.index_pause_time =  index_pause_time_1 ;
            trace{j}.distance_between_pauses =  distance_between_pauses_1 ;
            
            %distance_between_pauses_turn = distance_between_pauses_1 / slope_right;
            trace{j}.distance_between_pauses_turn =  distance_between_pauses_turn_1;
            
            trace{j}.v_pausefree =  v_pausefree_1 ;
            trace{j}.v_pausefree_tps =  v_pausefree_tps_1 ;
            
            trace{j}.v_pausefree_mean =  v_pausefree_mean_1 ;
            
            %v_pausefree_mean_tps = v_pausefree_mean_1/slope_right;
            trace{j}.v_pausefree_mean_tps =  v_pausefree_tps_mean_1;%v_pausefree_mean_tps;
            
            
            
            trace{j}.v_pausefree =  v_pausefree_1;
            
            %v_pausefree_tps = v_pausefree_1/slope_right;
            
            trace{j}.v_pausefree_tps =  v_pausefree_tps_1;
     
            trace{j}.pause_density =  pause_density_1;
            
            trace{j}.v_waitinclude_mean = v_waitinclude_mean; 
            %trace{j}.v_waitinclude_mean_tps = v_waitinclude_mean/slope_right;
            trace{j}.v_waitinclude_mean_tps = v_waitinclude_tps_mean;
            
            trace{j}.Z_waitinclude = Z_waitinclude; 
            
            trace{j}.time_waitinclude = time_waitinclude; 
            trace{j}.turn_waitinclude = turn_waitinclude;
            trace{j}.turn_re_waitinclude = turn_re_waitinclude;
            
             %% find topological shift %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tb = Turn2_mean_1way; zb = Z_beadc2_mean_1way(:,j); 
            ta = Turn4_mean_1way; za = Z_beadc4_mean_1way(:,j);

            ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.SmoothingParam = 0.5;
            %ia = za(:)>=0.2 & ta(:) < 200;
            ia =  ta(:) < 110;
            ta = ta(ia);za = za(ia);
            [xData, yData] = prepareCurveData( ta(:), za(:) ); [fitr_a, ~] = fit( xData, yData, ft, opts );
            %tb = tb(zb>=0.2);zb = zb(zb>=0.2);
            [xData, yData] = prepareCurveData( tb(:), zb(:) ); [fitr_b, ~] = fit( xData, yData, ft, opts );

            tcl = 0;

            tcr = tb(end);
            tc = linspace(tcl, tcr,100);
            
            if sign_winding_1
                dc = -30:0.5:80;
            else
                dc = -80:0.5:30;
            end
            
            deltad = 3;
            sd = zeros(1,length(dc));
            for i = 1 : length(dc)
                sd(i) = sum((fitr_b(tc) - fitr_a(tc+dc(i))).^2);
            end
            dmin = dc(sd == min(sd));
            x = dc(dc>=dmin-deltad & dc< dmin+deltad);
            y = sd(dc>=dmin-deltad & dc< dmin+deltad);
            [xData, yData] = prepareCurveData( x, y );
            ft = fittype( 'poly2' );
            [fitresult, ~] = fit( xData, yData, ft );
            turn_shift = -fitresult.p2/fitresult.p1/2;
            trace{j}.turn_shift = turn_shift;
            
            %% plot data
            f2 = figure;
            pos = [10 70 900 950];
            set(f2, 'Pos', pos);
            subplot(3,1,1:2);
            hold all
            
            line([0 0],[ -2 15],'LineStyle','--','LineWidth',1, 'Color','k');
            
            line([-55 150],[ height_0*0.9 height_0*0.9],'LineStyle','--','LineWidth',0.5, 'Color','k');
            
            line([d_wound1 d_wound1],[ -2 15],'LineStyle','--','LineWidth',1, 'Color',[0.87, 0.49, 0]);
            plot( ttempt, ztempt,'-c','LineWidth',2);
            
            plot(Turn4_mean_1way, Z_beadc4_mean_1way(:,j), '-o','MarkerSize',2,'LineWidth',1,'MarkerFaceColor',[0.5, 0.5, 0.5],'Color',[0.5, 0.5, 0.5]);
          
            plot(movingmean(Turn3,10), movingmean(Z_beadc3(:,j),10), '.','Color',[0.87, 0.49, 0]);
            
            
            plot(Turn2_mean_1way, Z_beadc2_mean_1way(:,j), '-ok','MarkerSize',4,'LineWidth',1.5,'MarkerFaceColor','k');
            

            force_topo = num2str(F_from_var_length_p);
            title([Box_namec1{j} ', F = '  force_topo(1:4) 'pN, Turn shift: ' num2str(turn_shift) ' turns, h_max: ' num2str(height_0) ' um'],'Interpreter','None','FontSize',12,'FontWeight','bold');
            ylabel('Z (\mum)','FontSize',12);
            xlabel('Turn added','FontSize',12);
            set(gca,'FontSize',12,'FontName','Calibri');
            hmm = 1.1*max(Z_beadc2_mean_1way(:,j));

            xlim([-50 Turn4_mean_1way(end)+10]); ylim([0 hmm+0.2]);
            textp1 = { ['Max height: ' , num2str(height_0) ' \mum'] };   
            
            t_no_topo = time_holding1(index_no_topo_hold_1);
            
            
            subplot(3,1,3);
            hold all
            line([0 time_holding1(end)+5],[height_0 height_0],'LineStyle','--','LineWidth',1, 'Color',[0.5, 0.5, 0.5]);
            line([0 time_holding1(end)+5],0.9 * [height_0 height_0],'LineStyle','--','LineWidth',1, 'Color',[0.5, 0.5, 0.5]);
            plot(time_holding1, Z_holding_s1_1, '--','Color',[0, 0.75, 0.75], 'LineWidth',1);
            
            plot(time_holding1(index_hold_1), Z_holding_s1_1(index_hold_1), '-','Color',[0, 0.75, 0.75], 'LineWidth',1);
            plot(time_waitinclude, Z_waitinclude, '-','Color','magenta', 'LineWidth',1);
         
                       
            plot(time_1(cell2mat(index_pause_time_1(:))), Z_1(cell2mat(index_pause_time_1(:))),'.k', 'MarkerSize',2);
            plot(time_holding1(index_no_topo_hold_1), Z_holding_s1_1(index_no_topo_hold_1),':','LineWidth',1.5, 'Color',[0.5,0.5,0.5]);
            xlabel('Time (s)', 'FontSize',12); 
            ylabel('Height (\mum)', 'FontSize',12);
            xlim([0 915]); ylim([min(Z_holding_s1_1)-0.1 hmm+0.1]);
            %title('Positive side','FontSize',13);
            set(gca,'FontSize',12,'FontName','Calibri');     
            
            title(['On-time: ' num2str(t_on_topo_1) ' s, mean rate: ' num2str(trace{j}.v_mean_tps) ' tps, wait included:' num2str(trace{j}.v_waitinclude_mean_tps) ' tps'],'Interpreter','None','FontSize',12,'FontWeight','bold');

            prompt = {'Comment?'};
            topheader = 'Add comment';
            num_line = 1;
            defaultvals = {'good'};
            comments = inputdlg(prompt, topheader,num_line,defaultvals);
            trace{j}.comments = comments;        

            savefig(f2,[Box_namec1{j} '.fig']);
            %print(Box_namec1{j},'-dpng','-r0');
            export_fig(f2,[Box_namec1{j} '.png']);


            
            data1 = [trace{j}.F_from_var_length_p height_0  sign_winding_1 t_on_topo_1     ];
            data1 =  [trace{j}.name trace{j}.comments num2cell(data1)];
            if ~isempty(data1)
                        xlswrite('summary.xls', data1, 1, ['A' num2str(j)]);
            end

            close(f2);
        else
            prompt = {'Comment?'};
            topheader = 'Add comment';
            num_line = 1;
            defaultvals = {'bad_hysteresis'};
            comments = inputdlg(prompt, topheader,num_line,defaultvals);
            trace{j}.comments = comments;
            data1 = [trace{j}.name comments];
            if ~isempty(data1)
                        xlswrite('summary.xls', data1, 2, ['A' num2str(j)]);
            end
            close(ftemp);
        end
    else
        trace{j}.name = Box_namec1{j};
        trace{j}.comments = {'bad_hysteresis'};
    end
    
    save('data.mat','trace');
end
    
%save('data.mat','trace');

end

function [Turn2_mean_1way, Z_beadc2_mean_1way, Turn2_mean_1way_negative, Z_beadc2_mean_1way_negative] = average_at_pause(Times2, Turn2, Z_beadc2)
    k=1;
    time_pat{1} = Times2(1);
    turn_pat{1} = Turn2(1);
    index_pat{1} = 1;
    for i = 1 : length(Turn2)-1
            if abs(Turn2(i+1)-turn_pat{k}(end))<0.01
                time_pat{k} = [time_pat{k} Times2(i+1)];   
                turn_pat{k} = [turn_pat{k} Turn2(i+1)];   
                index_pat{k} = [index_pat{k} i+1];   
            else
                if abs(Turn2(i+1) - turn_pat{k}(end)) == 1
                    k = k+1;
                    time_pat{k} = Times2(i+1);
                    turn_pat{k} = Turn2(i+1);
                    index_pat{k} = i+1;
                end
            end
    end    
    n1 = length(time_pat);
    Turn2_mean  = zeros(1,n1);
    L = size(Z_beadc2,2);
    Z_beadc2_mean  = zeros(n1,L);
    for i = 1 : n1
            Turn2_mean(i) = mean(turn_pat{i}); 
            Z_beadc2_mean(i,:) = mean(Z_beadc2(index_pat{i},:),1);
    end
    index_pos_step2 = diff(Turn2_mean)>0; index_pos_step2_negative = diff(Turn2_mean)<0;
    Turn2_mean_1way = Turn2_mean(index_pos_step2); Z_beadc2_mean_1way = Z_beadc2_mean(index_pos_step2,:);
    Turn2_mean_1way_negative = Turn2_mean(index_pos_step2_negative); Z_beadc2_mean_1way_negative = Z_beadc2_mean(index_pos_step2_negative,:);
end


function [height_c_before, F_from_var_length, F_from_var_length_p] = force_calculation(Times1, X_beadc1, Y_beadc1, Z_beadc1, uB, height_correction, tv)
        
        kbt = 4.08;
        R_bead = 0.5;
        height_c_before = mean(Z_beadc1);
               
        index = find(Times1<=tv);
        [xData, yData] = prepareCurveData( Times1(index), X_beadc1(index));
        
        ft1 = fittype( 'poly1' );
        [fitresult1, ~] = fit( xData, yData, ft1 );
        xx = X_beadc1(index) - fitresult1(Times1(index));
        [xData, yData] = prepareCurveData( Times1(index), Y_beadc1(index));
        ft2 = fittype( 'poly1' );
        [fitresult2, ~] = fit( xData, yData, ft2 );
        yy = Y_beadc1(index) - fitresult2(Times1(index));
        uBp = [-uB(2), uB(1)]; uBp = uBp/norm(uBp);
        rB = [xx yy] * uB'; % projection on the B field
        rBp = [xx yy] * uBp'; % projection on the perpendicular direction to B field
        uB_var1 = mean(rB.^2) - mean(rB)^2;    % variance along along B   
        uB_var2 = mean(rBp.^2) - mean(rBp)^2;

        F_from_var_length = kbt * height_c_before*1e3/(uB_var1*1e6);                        % parallel force
        F_from_var_length_p = kbt * (height_c_before + R_bead - height_correction)*1e3/(uB_var2*1e6);    % perpendicular force 
end

function [paraf, height_0, turn_0, ttempt, ztempt,fitresult_height_turn_right, fitresult_height_turn_left] = fit_naked_dna(turn_c,z_m_c)
        % use 3-piece function to fit the data for dna % get the residue of fitting
        range_f = [-40 40];
        turn_c = turn_c(:);
        z_m_c = z_m_c(:);
        turn_c_f = turn_c(turn_c >= range_f(1) & turn_c <= range_f(2));
        z_m_c_f = z_m_c(turn_c >= range_f(1) & turn_c <= range_f(2));
        
        index_max = find(z_m_c_f >= max(z_m_c_f)*0.9);
        xf = turn_c_f( index_max);
        yf = z_m_c_f( index_max);
        [xData, yData] = prepareCurveData( xf, yf );
        ft = fittype( 'poly2' );
        [fitparabol, ~] = fit( xData, yData, ft );
        turn_guess = -fitparabol.p2 / (2 * fitparabol.p1);
        h_max_guess = fitparabol.p3 - fitparabol.p2^2/(4*fitparabol.p1);

        % for DNA, piecewise
        [xData, yData] = prepareCurveData( turn_c_f, z_m_c_f );       
        F_positive = @(para, x) (para(1) + para(2) * (x - para(3)).^2) .* (x < para(4) & x > para(5)) + (para(1) + para(2)* (para(4) - para(3)).^2 +para(6).*(x-para(4))) .* (x > para(4)) + ...
                            + (para(1) + para(2)* (para(5) - para(3)).^2 +para(7).*(x - para(5))) .* (x < para(5)) ;
        para0 = [h_max_guess, -0.001 , turn_guess, turn_guess+10, turn_guess-10, -70/1000, 70/1000];
        [paraf,~,~,~,~, ~, ~] = lsqcurvefit(F_positive,para0,xData, yData);
        
        %height_0 = paraf(1);
        
        height_0 = mean(z_m_c_f(abs(turn_c_f - turn_c_f(z_m_c_f == max(z_m_c_f))) <= 2));
        turn_0 = paraf(3);
        
        ttempt = -65:0.1:65;
        ztempt = F_positive(paraf,ttempt);
        
        % positive turns
        fitresult_height_turn_right = @(para, h) (para(4) + 1/para(6) * (h - para(1) - para(2) * (para(4)- para(3))^2) ) .* (h < para(1)+para(2)*(para(4) - para(3))^2) + ...
                    (para(3) + sqrt( (h - para(1))/para(2) ) ) .* (h >= para(1)+para(2)*(para(4)-para(3))^2 & h < para(1)) + para(3)*( h >= para(1));
        % negative turns
        fitresult_height_turn_left = @(para, h) (para(5) + 1/para(7) * (h - para(1) - para(2) * (para(5)- para(3))^2) ) .* (h < para(1)+para(2)*(para(5) - para(3))^2) + ...
                    (para(3) - sqrt( (h - para(1))/para(2) ) ) .* (h >= para(1)+para(2)*(para(5)-para(3))^2 & h < para(1)) + para(3)*( h >= para(1));
end

                
             
function [v_waitinclude_mean, v_waitinclude_tps_mean, Z_waitinclude, time_waitinclude, turn_re_waitinclude, turn_waitinclude, Z_topo_pause_mean, ...
    Z_pausefree, v_pausefree, v_pausefree_tps, dwelltime_pause, sign_winding, t_on_topo, ...
    Z_holding_s1, Z_holding_s2, turn_relaxed, v, v_mean, v_tps, v_tps_mean, time, Z, turn, turn_re, index_pause_time, pause_density, ...
    distance_between_pauses, distance_between_pauses_turn, v_pausefree_mean, v_pausefree_tps_mean, turn_re_pausefree, index_hold, index_no_topo_hold] = ...
    analayze_topo_activity(time_holding,turn_holding, Z_holding, h_max, paraf, fitresult_height_turn_right, fitresult_height_turn_left)  
        
                    
        d_wound = turn_holding(end);
        d_0 = turn_holding(1);
        sign_winding = mean(diff(turn_holding))>0; % positive if winding positively
        
        rate = 1/2.5;  
        dt = mean(diff(time_holding));
        window = round(1/rate / dt);
        %order = 2;
        if mod(window, 2) == 0
                window = window+1;
        end
        %Z_holding_s1 = sgolayfilt(Z_holding,order,window);
        Z_holding_s1 = movingmean(Z_holding,window);
        
        if sign_winding
            turn_convert = fitresult_height_turn_right(paraf, Z_holding_s1);
            turn_relaxed = (turn_holding - d_0) - turn_convert; % the number of turn relaxed by topo, always positive          
        else
            turn_convert = fitresult_height_turn_left(paraf, Z_holding_s1); 
            turn_relaxed = -(turn_holding - d_0) + turn_convert; % the number of turn relaxed by topo, always positive
        end
        
         
        window2 = round(2.5/ dt);
        if mod(window2, 2) == 0
                window2 = window2+1;
        end
        %Z_holding_s2 = sgolayfilt(Z_holding,order,window2);
        Z_holding_s2 = movingmean(Z_holding,window2);
          
        %% data when holding    
        dd_wound = 0.05;
        index_hold = abs(turn_holding-d_wound)<dd_wound;   % holding
        
        %% v_magnet
        turn_magnet = turn_holding(~index_hold & abs(turn_holding)>=1);
        time_magnet = time_holding(~index_hold & abs(turn_holding)>=1);
        v_magnet = mean(diff(turn_magnet)./diff(time_magnet));
        
        %% wait time
        time_holding_magnetstop = time_holding(index_hold);
        Z_holding_magnetstop = Z_holding_s1(index_hold);
        turn_holding_magnetstop = turn_holding(index_hold);
        turn_holding_relaxed = turn_relaxed(index_hold);
         
        index_start_find = find( Z_holding_s1 == min(Z_holding_s1(abs(time_holding - time_holding_magnetstop(1)) <= 1)));        
        index_hold = time_holding>= time_holding(index_start_find);
        
        time_holding_magnetstop = time_holding(index_hold);
        Z_holding_magnetstop = Z_holding_s1(index_hold);
        turn_holding_magnetstop = turn_holding(index_hold);
        turn_holding_relaxed = turn_relaxed(index_hold);
        
        time_holding_magnetwind = time_holding(~index_hold);
        Z_holding_magnetwind = Z_holding_s1(~index_hold);
        turn_holding_magnetwind = turn_holding(~index_hold);
        turn_holding_relaxed_magnetwind = turn_relaxed(~index_hold);
        
        index_temp = find(Z_holding_s1-Z_holding_magnetwind(end) <= 75/1000);
        
        
        %index_no_topo_hold = abs(turn_holding-d_wound)<dd_wound & time_holding >= time_holding(index_temp(1)) & time_holding <= time_holding(index_temp(end)) ;    
        index_no_topo_hold = time_holding>= time_holding(index_start_find)  &  time_holding <= time_holding(index_temp(end)) ;   
        
        index_no_topo_hold_indices = find(index_no_topo_hold);
            
            if isempty(index_no_topo_hold_indices)
                t_on_topo = 0;
                time = time_holding_magnetstop;
                Z = Z_holding_magnetstop;
                turn = turn_holding_magnetstop;
                turn_re = turn_holding_relaxed;
                
                time_waitinclude = time_holding_magnetstop;
                Z_waitinclude = Z_holding_magnetstop;
                turn_waitinclude = turn_holding_magnetstop;
                turn_re_waitinclude = turn_holding_relaxed;
            else
                t_on_topo = time_holding(index_no_topo_hold_indices(end)) - time_holding_magnetstop(1);
                time = time_holding(~index_no_topo_hold & index_hold);
                Z = Z_holding_s1(~index_no_topo_hold & index_hold);
                turn = turn_holding(~index_no_topo_hold & index_hold);
                turn_re = turn_relaxed(~index_no_topo_hold & index_hold);
                
                time_waitinclude = time_holding_magnetstop;
                Z_waitinclude = Z_holding_magnetstop;
                turn_waitinclude = turn_holding_magnetstop;
                turn_re_waitinclude = turn_holding_relaxed;
            end 

        %% Relaxation speed
        window_slope = window2; v_pausefree_mean = []; v_pausefree = []; Z_pausefree = [];  
        time_pause = {}; index_pause_time = {}; t_pause_update = {}; Z_topo_pause = {}; dwelltime_pause = []; pause_density = []; distance_between_pauses = []; 
        v = [];   v_tps = [];  turn_topo_pause = {}; distance_between_pauses_turn = [];
        Z_topo_pause_mean = []; v_mean = [];   v_waitinclude = []; v_pausefree_tps = [];  v_pausefree_tps_mean = []; v_tps_mean = [];
        turn_re_pausefree = [];
        
        if length(Z) >= 2
                %index_end = find(Z == max(Z));
                condition = Z >= h_max * 0.9;
                                
                if sum(condition) == length(condition)  %% too fast
                    v_tps = 0/0; v = 0/0; v_pausefree_mean = 0/0; v_pausefree = 0/0; index_pause_time{1} = 0/0; t_pause_update{1} = 0/0; time_pause{1} = 0/0;
                    Z_topo_pause{1} = 0/0; dwelltime_pause = 0/0; pause_density =0/0; distance_between_pauses = 0/0; Z_pausefree = 0/0;
                    turn_topo_pause{1} = 0/0; distance_between_pauses_turn = 0/0; Z_topo_pause_mean = 0/0; v_mean = 0/0; v_waitinclude = 0/0;  v_waitinclude_tps = 0/0;
                    v_pausefree_tps = 0/0; turn_re_pausefree = 0/0; v_waitinclude_mean = 0/0;
                    
                    v_pausefree_tps_mean = v_magnet;
                    v_tps_mean = v_magnet;
                    v_waitinclude_tps_mean = v_magnet;
                    
                else
                    
                    index_end = find(condition);

                    if isempty(index_end)
                        index_end = length(Z);
                    end

                    tend = time(index_end(1));
                    Z = Z(time <= tend);
                    turn = turn(time <= tend);
                    time = time(time <= tend);
                    turn_re = turn_re(time <= tend);

                    turn_waitinclude = turn_waitinclude(time_waitinclude <= tend);
                    Z_waitinclude = Z_waitinclude(time_waitinclude <= tend);
                    time_waitinclude = time_waitinclude(time_waitinclude <= tend);
                    turn_re_waitinclude = turn_re_waitinclude(time_waitinclude <= tend);

                    window_slope = min([window2, length(Z)]);
                    %if length(Z) > window2

                        v_waitinclude = movingslope(Z_waitinclude, window_slope,1,dt);   
                        v_waitinclude_tps = movingslope(turn_re_waitinclude, window_slope,1,dt);  
                        v_waitinclude_mean = mean(v_waitinclude);
                        v_waitinclude_tps_mean = mean(v_waitinclude_tps);

                        v = movingslope(Z, window_slope,1,dt);
                        v_mean = mean(v);
                        v_tps = movingslope(turn_re, window_slope,1,dt); 
                        v_tps_mean = mean(v_tps);

                        %% Dwell-time analysis
                        binsize = 0.04;
                        range = min(Z) : binsize : max(Z);
                        P = length(range)-1;
                        dwell_time = zeros(1, P);
                        for i = 1 : P
                                index_jtemp = Z > range(i) & Z<range(i+1);
                                if ~any(index_jtemp)
                                    dwell_time(i) = 0;
                                else
                                    timetemp = time(index_jtemp); %timetemp = sort(timetemp);
                                    difftime = diff(timetemp); 
                                    index2 = find(difftime >= 2*dt); % at least separated by 2 steps to be considered two seperate regions
                                    if isempty(index2)
                                        dwell_time(i) = timetemp(end) - timetemp(1); % in ms
                                    else
                                        index2n = [0; index2(:)];
                                        index3 = [index2n(2:end); length(timetemp)];
                                       dwell_time(i) = sum(timetemp(index3) -  timetemp(index2n+1));
                                    end
                                end
                        end
                        range_new = range + binsize/2; range_new(end) = [];

                        dwelltime_threshold = 1.5;
                        index_dwell = dwell_time >= dwelltime_threshold;

                        if any(index_dwell)
                            count = 1;
                            dwelltemp = [];
                            Z_pause = {};
                            for i = 1 : length(index_dwell)-1
                                    if index_dwell(i) == 1
                                        dwelltemp = [dwelltemp range_new(i)];
                                    end
                                    Z_pause{count} = dwelltemp;   % contain all j of paused states
                                    if index_dwell(i)-index_dwell(i+1) == -1
                                        dwelltemp = [];
                                        count = count +1;
                                    end
                            end
                            if index_dwell(end) == 1
                                dwelltemp = [dwelltemp range_new(length(index_dwell))];
                                Z_pause{count} = dwelltemp;
                            end

                            Z_pause = Z_pause(~cellfun(@isempty, Z_pause));


                            for i = 1 : length(Z_pause)
                                    index_p = Z>=min(Z_pause{i})-binsize/4 & Z<=max(Z_pause{i})+binsize/4;
                                    time_pause{i} = time(index_p);
                                    index_pause_time{i} = find(time>=min(time_pause{i}) & time<=max(time_pause{i}));    % save
                                    t_pause_update{i} = time(index_pause_time{i});  % save
                                    Z_topo_pause{i} = Z(index_pause_time{i}); % save  
                                    turn_topo_pause{i} = turn_re(index_pause_time{i}); % save  

                                    dwelltime_pause(i) = t_pause_update{i}(end) - t_pause_update{i}(1);
                            end

                            Z_topo_pause_mean = cellfun(@mean, Z_topo_pause);
                            turn_topo_pause_mean = cellfun(@mean, turn_topo_pause);

                            pause_density = (length(Z_topo_pause_mean)-1)/(Z(end)-Z(1));
                            if length(Z_topo_pause_mean) >=2
                                distance_between_pauses = diff(Z_topo_pause_mean);

                                distance_between_pauses_turn = diff(turn_topo_pause_mean);


                            else
                                distance_between_pauses = Z(end)-Z(1);

                                distance_between_pauses_turn = turn_re(end)-turn_re(1);
                            end

                            v_pausefree = v;
                            v_pausefree_tps = v_tps;

                            time_pausefree = time;
                            Z_pausefree = Z;
                            turn_re_pausefree = turn_re;

                            v_pausefree(cell2mat(index_pause_time(:))) = [];
                            v_pausefree_tps(cell2mat(index_pause_time(:))) = [];
                            Z_pausefree(cell2mat(index_pause_time(:))) = [];
                            turn_re_pausefree(cell2mat(index_pause_time(:))) = [];
                            time_pausefree(cell2mat(index_pause_time(:))) = [];


                            %% mean pause free speed, 

                            v_pausefree_mean = mean(v_pausefree);
                            v_pausefree_tps_mean = mean(v_pausefree_tps);
                            
                        else
                            
                            v_pausefree = 0/0;
                            v_pausefree_mean = 0/0;
                            v_pausefree_tps = 0/0;
                            v_pausefree_tps_mean = 0/0;
                            
                        end

                            


                    %else
                    %        v = 0/0; 
                    %        v_tps = 0/0; 
                    %        v_pausefree = 0/0;
                    %        v_pausefree_mean = 0/0;
                    %        v_pausefree_tps = 0/0;%

                    %        v_pausefree_tps_mean = 0/0;
                    %        v_mean = 0/0;  
                    %        v_tps_mean = 0/0;
                            
                    %        v_waitinclude_mean = 0/0;
                    %        v_waitinclude_tps_mean = 0/0;
                    %end
                    
                end
   
        else
            
                v_waitinclude = movingslope(Z_waitinclude, window_slope,1,dt);   
                
                v_waitinclude_tps = movingslope(turn_re_waitinclude, window_slope,1,dt);   
                
                if any(v_waitinclude)
                    v_waitinclude_mean = mean(v_waitinclude);
                    
                    v_waitinclude_tps_mean = mean(v_waitinclude_tps);
                else
                    v_waitinclude_mean = 0/0;
                    
                    v_waitinclude_tps_mean = 0/0;
                end
        end
                
        
            

end

function Z_beadc3 = z_correct(Times3, Z_beadc3, L)

zmin = -1; zmax = 5;
for i = 1 :L
    ztemp = Z_beadc3(:,i);
    index_large = ztemp<zmin | ztemp>zmax;
    ttemp = Times3(~index_large);
    ztemptt = ztemp(~index_large);
    
    if sum(index_large) <= size(ttemp,1)
        if mean(ztemp)> zmin && mean(ztemp)< zmax
            if sum(index_large) >=1
                ztemp(index_large) = interp1(ttemp(:), ztemptt(:), Times3(index_large));
                Z_beadc3(:,i) = ztemp;
            end
        end
    end
end
end
