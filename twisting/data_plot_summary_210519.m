% Copyright 2022 Wang Lab
% 
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

filename = uigetfile('*.*','Select the DNA file','multiselect','on');
load(filename);
%load('data_IIA_0p5pM_noETO_pos_2.mat');
load('x_f_WLC.mat');

P = length(trace);


turn_m_before_positive = {};
Z_m_before_positive = {};

num = [];
for i= 1 : P
    
    num(i) = sscanf(trace{i}.name,'Box%d');
end

num_pos = find(diff(num) < 0);

P1 = num_pos;

j= 0;
for i = 1 : P1
    if strcmp(trace{i}.comments ,'good')
        if trace{i}.height_c_before >= 2.4 && trace{i}.F_from_var_length_p < 1
           j = j+1;
           name_a{j} = trace{i}.name;
           %force(j) = trace{i}.F_from_var_length_p;

           %height_c_before(j) = trace{i}.height_c_before;
           
           %height_c_after(j) = trace{i}.height_c_after;
           
           %turn_m_before_positive{j} = trace{i}.turn_m_before_positive;
           %Z_m_before_positive{j} = trace{i}.Z_m_before_positive;
                     
           %turn_m_after_positive{j} = trace{i}.turn_m_after_positive;
           %Z_m_after_positive{j} = trace{i}.Z_m_after_positive;
           
           %force(j) = trace{i}.F_from_var_length_p;
           
           %turn_shift(j) = trace{i}.turn_shift;
           
           %time_holding{j} = trace{i}.time_holding1;
           
           %z_holding{j} = trace{i}.z_holding1;

           %t_on_topo(j) =  trace{i}.t_on_topo;

           %waittime(j) =   trace{i}.waittime;
            
            pause_first(j) =   trace{i}.pause_first;
            
            turn_first(j) =   trace{i}.turn_first;
            
            v_burst(j) =   trace{i}.v_burst;

            if ~isempty(trace{i}.v_mean_tps)
                v_mean(j) =   trace{i}.v_mean_tps;
            else
                v_mean(j) =   NaN;
            end

            
            %turn_mean_u{j} =   trace{i}.turn_mean_u;
            
            %turn_step{j} = diff(trace{i}.turn_mean_u);

        end
        
    end
    
end

M = j;

v_burst = v_burst(turn_first >= 4);

turn_first = turn_first(isfinite(turn_first));

pause_first = pause_first(isfinite(pause_first));

v_burst = v_burst(isfinite(v_burst));

v_mean = v_mean(isfinite(v_mean));

dr = 0 : 0.5 : 29;
tr = 0 : 0.5 : 100;

d_cu = []; t_cu = [];
for i = 1 : length(dr)
    d_cu(i) = sum(turn_first < dr(i)) / length(turn_first);
    d_cu(i) = 1 - d_cu(i);
end

for i = 1 : length(tr)
    t_cu(i) = sum(pause_first < tr(i)) / length(pause_first);
    t_cu(i) = 1 - t_cu(i);
end

ft = fittype( 'exp(-a * x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = 0;
opts.StartPoint = 0.5;

[x1, y1] = prepareCurveData( dr, d_cu );
[fit1, ~] = fit( x1, y1, ft, opts );
fit1r = confint(fit1, 0.95); 
rate_first_fit = fit1.a;
drate_first_fit = diff(fit1r)/2;

ft = fittype( 'exp(-x/a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = 0;
opts.StartPoint = 1;

[x2, y2] = prepareCurveData( tr, t_cu );
[fit2, ~] = fit( x2, y2, ft, opts );
fit2r = confint(fit2, 0.95); 
dwelltime_first = fit2.a;
ddwelltime_first = diff(fit2r)/2;

% ------------- second experiment

j= 0;
for i = P1 + 1 : P
    if strcmp(trace{i}.comments ,'good')
        if trace{i}.height_c_before >= 2.4 && trace{i}.F_from_var_length_p < 1
           j = j+1;
           %name_a{j} = trace{i}.name;
           %force(j) = trace{i}.F_from_var_length_p;

           %height_c_before(j) = trace{i}.height_c_before;
           
           %height_c_after(j) = trace{i}.height_c_after;
           
           %turn_m_before_positive{j} = trace{i}.turn_m_before_positive;
           %Z_m_before_positive{j} = trace{i}.Z_m_before_positive;
                     
           %turn_m_after_positive{j} = trace{i}.turn_m_after_positive;
           %Z_m_after_positive{j} = trace{i}.Z_m_after_positive;
           
           %force(j) = trace{i}.F_from_var_length_p;
           
           %turn_shift(j) = trace{i}.turn_shift;
           
           %time_holding{j} = trace{i}.time_holding1;
           
           %z_holding{j} = trace{i}.z_holding1;

           %t_on_topo(j) =  trace{i}.t_on_topo;

           %waittime(j) =   trace{i}.waittime;
            
            pause_first2(j) =   trace{i}.pause_first;
            
            turn_first2(j) =   trace{i}.turn_first;
            
            v_burst2(j) =   trace{i}.v_burst;

            if ~isempty(trace{i}.v_mean_tps)
                v_mean2(j) =   trace{i}.v_mean_tps;
            else
                v_mean2(j) =   NaN;
            end
            
            %turn_mean_u{j} =   trace{i}.turn_mean_u;
            
            %turn_step{j} = diff(trace{i}.turn_mean_u);

        end
        
    end
    
end

M2 = j;

v_burst2 = v_burst2(turn_first2 >= 4);

turn_first2 = turn_first2(isfinite(turn_first2));

pause_first2 = pause_first2(isfinite(pause_first2));

v_burst2 = v_burst2(isfinite(v_burst2));

v_mean2 = v_mean2(isfinite(v_mean2));

dr = 0 : 0.5 : 29;
tr = 0 : 0.5 : 100;

d_cu2 = []; t_cu2 = [];
for i = 1 : length(dr)
    d_cu2(i) = sum(turn_first2 < dr(i)) / length(turn_first2);
    d_cu2(i) = 1 - d_cu2(i);
end

for i = 1 : length(tr)
    t_cu2(i) = sum(pause_first2 < tr(i)) / length(pause_first2);
    t_cu2(i) = 1 - t_cu2(i);
end

ft = fittype( 'exp(-a * x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = 0;
opts.StartPoint = 0.5;

[x1, y1] = prepareCurveData( dr, d_cu2 );
[fit12, ~] = fit( x1, y1, ft, opts );
fit12r = confint(fit12, 0.95); 
rate_first_fit2 = fit12.a;
drate_first_fit2 = diff(fit12r)/2;

ft = fittype( 'exp(-x/a)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = 0;
opts.StartPoint = 1;

[x2, y2] = prepareCurveData( tr, t_cu2 );
[fit22, ~] = fit( x2, y2, ft, opts );
fit22r = confint(fit22, 0.95); 
dwelltime_first2 = fit22.a;
ddwelltime_first2 = diff(fit22r)/2;


%% plot data

fm1 = figure;
subplot(2,1,1)
hold all
plot(dr, d_cu, 'xk-');
plot(dr, d_cu2, 'xr-');
hold all
plot(fit1);
plot(fit12);
legend off
lt = line([30 30], [0 1]); set(lt, 'LineStyle', '--');
xlabel('# of turns relaxed prior to the first pause');
ylabel('Survivor probability');
set(gca,'FontSize',15,'FontName','Calibri');
ylim([0, 1]);
xlim([0, 40]);

rate_first_fit_mean = 1/(M + M2) * (M * rate_first_fit + M2 * rate_first_fit2);
drate_first_fit_mean = sqrt(( M * drate_first_fit^2 + M2 * drate_first_fit2^2 +  M * (rate_first_fit_mean - rate_first_fit)^2 + M2 * (rate_first_fit_mean - rate_first_fit2)^2)/(M+M2) );

title({[num2str(rate_first_fit) ' \pm ' num2str(drate_first_fit) ' 1/turn'],[num2str(rate_first_fit2) ' \pm ' num2str(drate_first_fit2) ' 1/turn'], [num2str(rate_first_fit_mean) ' \pm ' num2str(drate_first_fit_mean) ' 1/turn']});


subplot(2,1,2)
hold all
plot(tr, t_cu, 'x-k');
plot(tr, t_cu2, 'x-r');
hold all
plot(fit2);
plot(fit22);
legend off
lb = line([100 100], [0 1]); set(lb, 'LineStyle', '--');
xlabel('Dwell-time to exit a pause (s)');
ylabel('Survivor probability');
set(gca,'FontSize',15,'FontName','Calibri');
ylim([0, 1]);
xlim([0, 320]);

dwelltime_first_mean = 1/(M + M2) * (M * dwelltime_first + M2 * dwelltime_first2);
ddwelltime_first_mean = sqrt(( M * ddwelltime_first^2 + M2 * ddwelltime_first2^2 +  M * (dwelltime_first_mean - dwelltime_first)^2 + M2 * (dwelltime_first_mean - dwelltime_first2)^2)/(M+M2) );

title({[num2str(dwelltime_first) ' \pm ' num2str(ddwelltime_first) ' s'],[num2str(dwelltime_first2) ' \pm ' num2str(ddwelltime_first2) ' s'], [num2str(dwelltime_first_mean) ' \pm ' num2str(ddwelltime_first_mean) ' s']});


pos = [10 70 500 700];
set(fm1, 'Pos', pos);

savefig(fm1,'turn_enter_pause__time_exit_pause2.fig');
print('turn_enter_pause__time_exit_pause2','-dpng','-r0');


%% Pause-free speed
dx = 0.1;
x=0:dx:20;
xn = x+dx/2;
n =histc(v_burst,x);

n2 =histc(v_burst2,x);

%[xData, yData] = prepareCurveData( xn, n );
%ft = fittype( 'gauss1' );
%opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Display = 'Off';
%opts.Lower = [-Inf -Inf 0];
%opts.StartPoint = [10 2 0.1];
%[fitresult, gof] = fit( xData, yData, ft, opts );

f2 = figure;
hold all
h = bar(xn, n,'Facecolor','g','LineWidth',1.5, 'BarWidth',1);
h2 = bar(xn, n2,'Facecolor','r','LineWidth',1.5, 'BarWidth',1);
%hold all
%plot(fitresult);
xlim([0,8]);
xlabel('Speed of first burst (turns/s)');
ylabel('Count');
set(gca,'FontSize',14,'LineWidth',1.5);

v_burst_mean = 1/(M + M2) * (M * mean(v_burst) + M2 * mean(v_burst2));
dv_burst_mean =   sqrt((M * (mean(v_burst) - v_burst_mean)^2 + M2 * (mean(v_burst2) - v_burst_mean)^2)/(M+M2) );

title({[num2str(mean(v_burst)) ' tps'], [num2str(mean(v_burst2)) ' tps'], [num2str(v_burst_mean) ' \pm ' num2str(dv_burst_mean) ' turn/s']});


%title( {['Mean +- SD: ' num2str(mean(v_burst)) ' \pm ' num2str(std(v_burst)) ', sem: ' num2str(std(v_burst)/sqrt(length(v_burst))) ' tps'], ['Fitted: ' num2str(fitresult.b1) ' \pm ' num2str(fitresult.c1/sqrt(2)) ', sem: ' num2str(std(v_burst)/sqrt(length(v_burst))) ' tps'] },'FontSize',14);

savefig(f2,'Pause_free_speed2.fig');
print('Pause_free_speed2','-dpng','-r0');

dx = 0.1;
x=0:dx:20;
xn = x+dx/2;
n =histc(v_mean,x);
n2 =histc(v_mean2,x);
f3 = figure;
hold all
h = bar(xn, n,'Facecolor','g','LineWidth',1.5, 'BarWidth',1);
h2 = bar(xn, n2,'Facecolor','r','LineWidth',1.5, 'BarWidth',1);
%hold all
%plot(fitresult);
xlim([0,8]);
xlabel('Mean speed (turns/s)');
ylabel('Count');
set(gca,'FontSize',14,'LineWidth',1.5);

v_mean_mean = 1/(M + M2) * (M * mean(v_mean) + M2 * mean(v_mean2));
dv_mean_mean =   sqrt((M * (mean(v_mean) - v_mean_mean)^2 + M2 * (mean(v_mean2) - v_mean_mean)^2)/(M+M2) );

title({[num2str(mean(v_mean)) ' tps'], [num2str(mean(v_mean2)) ' tps'], [num2str(v_mean_mean) ' \pm ' num2str(dv_mean_mean) ' turn/s']});


%title( {['Mean +- SD: ' num2str(mean(v_burst)) ' \pm ' num2str(std(v_burst)) ', sem: ' num2str(std(v_burst)/sqrt(length(v_burst))) ' tps'], ['Fitted: ' num2str(fitresult.b1) ' \pm ' num2str(fitresult.c1/sqrt(2)) ', sem: ' num2str(std(v_burst)/sqrt(length(v_burst))) ' tps'] },'FontSize',14);

savefig(f3,'Mean_speed2.fig');
print('Mean_speed2','-dpng','-r0');

