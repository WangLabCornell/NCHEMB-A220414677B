% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

function [l0p5_a, time_hold_a] = showtracestretch_hold

filename = uigetfile('*.*','Select the DNA file','multiselect','on');
%pathname = uigetfile_n_dir;
if isequal(filename,0)
    disp('User selected Cancel');
else
    if ~iscellstr(filename)
        L = 1;
    else
        L = length(filename);
    end
    
    button = 'Yes';
    j=0;
    
    
        load('x_F_WLC');  % L_r: relative length, Force_WLC: WLC force
        Loo = 12688;
        Lo = Loo * 0.338;

        
        ft = [0.5 1 2 3];
        lt = zeros(L,3);
        
        f0 = figure;
        f2 = figure;
        
        
        l0p5_a = [];
        time_hold_a = [];
        
    while j < L && ~strcmp(button,'Cancel')
        %close all;
        f1 = figure;
        pos = [10 70 600 800];
        set(f1, 'Pos', pos);  
        hold on;
        j = j+1;
        if ~iscellstr(filename)
            pn = filename;
        else
            pn = filename{j};
        end

        [Times, ForcepN, LDNAnm, AOMSentV, xpolym] = getdata(pn);
        
        if max(LDNAnm)-min(LDNAnm)<1
            tracenuc{j}.converted = 'no';
        else
            tracenuc{j}.converted = 'yes';
            
           

            Times1 = Times(AOMSentV < 0.1);
            ForcepN1 = ForcepN(AOMSentV < 0.1);
            LDNAnm1 = LDNAnm(AOMSentV < 0.1);
            AOMSentV1 = AOMSentV(AOMSentV < 0.1);

            % finding extension at intiail X centering
            index11 = ForcepN1 > 0 & LDNAnm1 > 0 & ForcepN1 <= 1;
            LDNAnm1 = LDNAnm1(index11);
            ForcepN1 = ForcepN1(index11);
            
            [xData, yData] = prepareCurveData( LDNAnm1, ForcepN1 );
            
            if any(xData)
                
                % Set up fittype and options.

                ftX = fittype( '4.1/43*(1/(4*(1-x/b)^2) -1/4 + x/b)', 'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.Lower = 0;
                opts.StartPoint = 5000;

                % Fit model to data.
                [fit_Xcentering, ~] = fit( xData, yData, ftX, opts );

                l0p5 = fsolve(@(x)fit_Xcentering(x)-0.5, max(xData));

                l0p5_a(j) = l0p5;
                
                jbp_0p5 = l0p5./(interp1(Force_WLC, L_r, 0.5) * 0.338);
            else
                l0p5 = 0/0;
                l0p5_a(j) = l0p5;
                fit_Xcentering = 0/0;
                
                jbp_0p5 = 0;
            end
            
            

            %% Final stretching
            Times = Times(AOMSentV >= 0.1);
            ForcepN = ForcepN(AOMSentV >= 0.1);
            LDNAnm = LDNAnm(AOMSentV >= 0.1);
            xpolym  = xpolym(AOMSentV >= 0.1);

            dt = mean(diff(Times));
             
             
            % no smoothing
            window_size = 1;%round(0.01/dt);   % 30ms
            ForcepN_s = movingmean(ForcepN, window_size);
            LDNAnm_s = movingmean(LDNAnm, window_size);
            Times_s = Times;


            %% Start
                    t_on_pick = Times_s(1)+0.03;
                    index_start = find(Times_s >= t_on_pick); 
                    index_start = index_start(1);
            
                    ForcepN_s2 = ForcepN_s(index_start:end);
                    LDNAnm_s2 = LDNAnm_s(index_start:end);
                    Times_s2 = Times_s(index_start:end);

                  
            %% Hold time 
            f_hold = 59; df_hold = 4;
            index_hold = abs(ForcepN_s2 - f_hold) <= df_hold;
            Times_s2_hold =  Times_s2(index_hold);
            LDNAnm_s2_hold = LDNAnm_s2(index_hold);
            ForcepN_s2_hold = ForcepN_s2(index_hold);
            
            time_hold = length(Times_s2_hold) * dt;
            
            time_hold_a(j) = time_hold;
 
            figure(f0);
            hold all
            plot(LDNAnm_s, ForcepN_s,'.:', 'Color',rand(1,3), 'MarkerSize',2);

            figure(f1);

            subplot(3,1,1)
            hold all
            plot(LDNAnm1, ForcepN1,'o-');
            plot(fit_Xcentering);
            plot(l0p5, 0.5, 'or', 'MarkerFaceColor','r');
            ylim([-0.1, 1])
            xlim([0 4000]);
            set(gca,'FontSize',13,'FontName','Calibri');
            legend off

            xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
            title({pn, ['Extension @ 0.5 pN: ' num2str(l0p5) ' nm'] },'Interpreter','None','FontSize',12);  

            subplot(3,1,2)
            hold all
            plot(L_r * Loo *0.338, Force_WLC,':', 'Color','k');

            plot(LDNAnm_s, ForcepN_s,'o', 'Color','k', 'MarkerFaceColor','k', 'MarkerSize',2);
            plot(LDNAnm_s2_hold, ForcepN_s2_hold,'o', 'Color','r', 'MarkerFaceColor','r', 'MarkerSize',4);
            
            axis([400 5500 0 100 ]);

            set(gca,'FontSize',13,'FontName','Calibri');
            xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
            title('Extension vs. Force: ','Interpreter','None','FontSize',10);        

            subplot(3,1,3)
            hold all
            plot(Times_s, ForcepN_s,'--', 'Color',[0.5,0.5,0.5]);
            plot(Times_s, ForcepN_s,'o', 'Color','k', 'MarkerFaceColor','k', 'MarkerSize',1);
            plot(Times_s2_hold, ForcepN_s2_hold,'o', 'Color','r', 'MarkerFaceColor','r', 'MarkerSize',4);
                        
           
            ylim([0 100]);
            xlim([0 120]);
            title(['Hold time: ' num2str(time_hold) ' s'],'Interpreter','None','FontSize',10);        

            set(gca,'FontSize',13,'FontName','Calibri');
            xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
            
            figure(f2)
            plot(Times_s,LDNAnm_s,'--');                     
           
            axis([0 120 400 5500])
            title(['Time vs Extantion'],'Interpreter','None','FontSize',10);        

            set(gca,'FontSize',13,'FontName','Calibri');
            xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Extension (pN)', 'FontSize',13, 'FontName', 'Calibri');


            savefig(f1,[pn '_check.fig']);
            export_fig(f1,[pn '_check.png']);

            tracenuc{j}.name = pn;
            tracenuc{j}.Times1 = Times1;
            tracenuc{j}.ForcepN1 = ForcepN1;
            tracenuc{j}.LDNAnm1 = LDNAnm1;
            tracenuc{j}.fit_Xcentering = fit_Xcentering;
            tracenuc{j}.l0p5 = l0p5;
            tracenuc{j}.jbp_0p5 = jbp_0p5;


            tracenuc{j}.Times = Times;
            tracenuc{j}.ForcepN = ForcepN;
            tracenuc{j}.LDNAnm = LDNAnm;
            tracenuc{j}.Times_s = Times_s;
            tracenuc{j}.ForcepN_s = ForcepN_s;
            tracenuc{j}.LDNAnm_s = LDNAnm_s;
            
            tracenuc{j}.ForcepN_s2 = ForcepN_s2;
            tracenuc{j}.Times_s2 = Times_s2;
            tracenuc{j}.LDNAnm_s2 = LDNAnm_s2;
            
            tracenuc{j}.Times_s2_hold = Times_s2_hold;
            tracenuc{j}.LDNAnm_s2_hold = LDNAnm_s2_hold ;
            tracenuc{j}.ForcepN_s2_hold = ForcepN_s2_hold;

            tracenuc{j}.t_on_pick = t_on_pick;
            
            tracenuc{j}.time_hold = time_hold;
            


        end
        
        save('data.mat','tracenuc');
        
    end
    
    figure(f0);
    hold all
    plot(L_r * Loo *0.338, Force_WLC,':', 'Color','k');
        axis([500 8000 0 125 ]);
        
        set(gca,'FontSize',13,'FontName','Calibri');
        xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
        ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
    figure(f2)
    hold all
    plot(Times_s,LDNAnm_s,'--');                     
           
    axis([0 120 400 5500])
    title(['Time vs Extantion'],'Interpreter','None','FontSize',10);        

    set(gca,'FontSize',13,'FontName','Calibri');
    xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
    ylabel('Extension (pN)', 'FontSize',13, 'FontName', 'Calibri');
    savefig(f2,[pn '_check2.fig']);
    export_fig(f2,[pn '_check2.png']);
        
end
end

function [Times, ForcepN,LDNAnm, AOMSentV, xpolym] = getdata(filename)
    fileID = fopen(filename,'r');
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    Times = dataArray{:, 1};
    AOMSentV = dataArray{:, 13};
    ForcepN = dataArray{:, 23};
    %jindex = dataArray{:, 23};
    LDNAnm = dataArray{:, 22};
    xpolym = dataArray{:, 17};
    clearvars delimiter startRow formatSpec fileID dataArray ans;
    
    %time_t = 5;
    %ForcepN = ForcepN(Times>time_t);
    %LDNAnm = LDNAnm(Times>time_t);
    %Times = Times(Times>time_t);
    
    %index0 = find(ForcepN == max(ForcepN));
    
    %Times = Times(1: index0(1)-1);
    %ForcepN = ForcepN(1: index0(1)-1);
    %LDNAnm = LDNAnm(1: index0(1)-1);
    
end
