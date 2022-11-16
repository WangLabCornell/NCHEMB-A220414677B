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
        tracenuc_a = [];
        for i = 1 : K
            load(filename{i});
            length(tracenuc)
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
    name = {};
    time_hold = [];
    l0p5 = [];

    
    j=0;
    for i = 1 : L
        j = j+1;
        name{j} = tracenuc_a{i}.name;
        l0p5(j) = tracenuc_a{i}.l0p5;
        time_hold(j) = tracenuc_a{i}.time_hold;
    end
    
    time_hold = time_hold(isfinite(time_hold));
    tr = linspace(0,120,200);
    M = length(tr);
    p = zeros(1, M);
    dp = zeros(1, M);
    ns = zeros(1, M);
    for i = 1 : M
        ns(i) = sum(time_hold <= tr(i));
        p(i) = 1 - ns(i)/length(time_hold);
        dp(i) = sqrt(p(i)* (1- p(i)) /length(time_hold));
    end
    p_u = p + dp/2;
    p_d = p - dp/2;
    
    p = p(tr < 120);
    p_u = p_u(tr < 120);
    p_d = p_d(tr < 120);
    dp = dp(tr < 120);
    tr = tr(tr < 120);
    
    [xData, yData] = prepareCurveData( tr, p );
    ft = fittype( 'a*exp(-b*x)+(1-a)*exp(-c*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.2 0.1 2];
    [fitresult, ~] = fit( xData, yData, ft, opts );

    a1 = fitresult.a;
    k1 = fitresult.b;
    a2 = 1-fitresult.a;
    k2 = fitresult.c;

    if k1> k2
        a_fast = a1; k_fast = k1;
        a_slow = a2; k_slow = k2;
    else
        a_fast = a2; k_fast = k2;
        a_slow = a1; k_slow = k1;
    end

    % fraction breaks fast
    % error bar from 2 confidence intervals, 2 sigma
    a = confint(fitresult,0.95); 
    
    
    a = a(:,1);
    a_mean = mean(a); 
    
    fraction_sd = 1/2 * abs(a(1) - a(2));
    if fitresult.c > fitresult.b
        fraction_fast = 1 - a_mean;
    else
        fraction_fast = a_mean;
    end
    
    
    t_1o2 = fsolve(@(x)fitresult(x)-0.5,0.2);
    
    f1 = figure;
    hold all
    plot(tr, p, '.-r');
    plot(tr, p_u,'--g');
    plot(tr, p_d,'--g');
    ff = plot(fitresult);
    set(ff,'LineWidth',2.5);
    xlabel('Hold time (s)');
    ylabel('Survival fraction');
    ylim([0, 1]); xlim([0, 120]);
    title(['T_{1/2}: ' num2str(t_1o2) ' s; Fraction breaks fast: ' num2str(fraction_fast) ' \pm ' num2str(fraction_sd)]);
    savefig(f1,'Summary_2exp.fig');


end
