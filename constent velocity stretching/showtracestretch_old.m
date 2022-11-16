% Copyright 2022 Wang Lab
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published % by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

function [lt, ft, f_max_a, l0p5_a, f_mean_a] = showtracestretch_old

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
        
        f_max_a = [];
        l0p5_a = [];
        
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
            
            dt = mean(diff(Times));

            Times1 = Times(AOMSentV < 0.1);
            ForcepN1 = ForcepN(AOMSentV < 0.1);
            LDNAnm1 = LDNAnm(AOMSentV < 0.1);
            AOMSentV1 = AOMSentV(AOMSentV < 0.1);

            % finding extension at intiail X centering
            index11 = ForcepN1 > 0 & LDNAnm1 > 0;
            LDNAnm1 = LDNAnm1(index11);
            ForcepN1 = ForcepN1(index11);
            
            [xData, yData] = prepareCurveData( LDNAnm1, ForcepN1 );
            
            if any(xData)
                
                % Set up fittype and options.
                %ftX = fittype( 'exp2' );
                %opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                %opts.Display = 'Off';
                %opts.Lower = [-Inf 0 -Inf 0];
                %opts.StartPoint = [0.05 1e-06 0.0005 0.001];

                %ftX = fittype( 'exp(b*x) - exp(c*x)', 'independent', 'x', 'dependent', 'y' );
                %opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                %opts.Display = 'Off';
                %opts.Lower = [0 0];
                %opts.StartPoint = [1e-06 0.001];

                ftX = fittype( '4.1/43*(1/(4*(1-x/b)^2) -1/4 + x/b)', 'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.Lower = 0;
                opts.StartPoint = 5000;

                % Fit model to data.
                [fit_Xcentering, ~] = fit( xData, yData, ftX, opts );

                l0p5 = fsolve(@(x)fit_Xcentering(x)-0.5, max(xData));

                l0p5_a(j) = l0p5;
            else
                l0p5 = 0/0;
                l0p5_a(j) = l0p5;
                fit_Xcentering = 0/0;
            end

            %% Final stretching
            Times = Times(AOMSentV >= 0.1);
            ForcepN = ForcepN(AOMSentV >= 0.1);
            LDNAnm = LDNAnm(AOMSentV >= 0.1);
            xpolym  = xpolym(AOMSentV >= 0.1);

            % no smoothing
            window_size = 1;%round(0.01/dt);   % 30ms
            ForcepN_s = movingmean(ForcepN, window_size);
            LDNAnm_s = movingmean(LDNAnm, window_size);
            Times_s = Times;

            if ForcepN_s(end) >= 5
                index_high = length(Times_s);
            else
                %dFdt = diff(ForcepN_s)/mean(diff(Times_s));
                indices_zero = find(ForcepN_s <= 5);
                dindices_zero = diff(indices_zero);
                indices_zero_end = find(dindices_zero>=2);
                indices_zero_end = indices_zero_end(end)+1;

                index_high = indices_zero(indices_zero_end);
            end
            %end
            t_e = Times_s(index_high); l_e = LDNAnm_s(index_high); f_e = ForcepN_s(index_high); 


            %% Final peak force
            [f_max, index_max] = max(ForcepN_s(index_high-100:index_high));
            index_max = index_high-100 + index_max - 1;
            t_max = Times_s(index_max);
            l_max = LDNAnm_s(index_max);

            f_max_a(j) = f_max;
            
            %% Find time-averaged force from transition to 60 pN: 
            vpolym = movingslope(xpolym, 100 ,1, dt);
            index_mean = find(vpolym <= mean(vpolym) & Times_s <= t_max  & ForcepN_s <= 60 ) ;
            f_mean = mean(ForcepN_s(index_mean));
            f_mean_a(j) = f_mean;

                    %f_test = figure;
                    %plot(Times_s, ForcepN_s, '.:');
                    %xlabel('Time (s)'); ylabel('Force (pN)');
                    %prompt = {'Pick t_on:'};
                    %tag = 'Pick start point:';
                    %num_line = 1;
                    %defaultvals = {num2str(Times_s(1)+1)};
                    %answer = inputdlg(prompt, tag,num_line,defaultvals);
                    %t_on_pick = cellfun(@str2num, answer)';
                    %close(f_test);
                    
                    t_on_pick = Times_s(1)+1;

                    index_start = find(Times_s >= t_on_pick); 


            %% Extension at 0.5, 1, 2, 3 pN
            window_size2 = 100;   % 30ms
            ForcepN_ss = movingmean(ForcepN(index_start:index_max), window_size2);
            LDNAnm_ss = movingmean(LDNAnm(index_start:index_max), window_size2);
            Times_ss = Times(index_start:index_max);


            for k = 1 : length(ft)
                indexc = abs(ForcepN_ss - ft(k)) <= 0.1;
                lr = LDNAnm_ss(indexc);
                fr = ForcepN_ss(indexc);


                xr = min(lr) : 50 : max(lr);
                xrn = xr(:)+25;
                nr = histc(lr, xr);
                indexse = find(nr> 0)%mean(nr(nr ~= 0))/2);
                indexgap = find(diff(indexse) > 1)

                if any(indexse)
                    if any(indexgap)
                        index_final = indexse(1): indexse(1)+indexgap(1);
                    else
                        index_final = indexse(1): indexse(end);
                    end
                    lt(j, k) = sum(nr(index_final) .* xrn(index_final)) / sum( nr(index_final));   
                else
                    lt(j, k) = 0/0;
                end

                       

            end


            j_bp = LDNAnm_s./(interp1(Force_WLC, L_r, ForcepN_s) * 0.338);
            % plot and save the stretching curve with multiple stretching lines

            figure(f0);
            hold all
            plot(LDNAnm_s, ForcepN_s,'.:', 'Color',rand(1,3), 'MarkerSize',2);

            figure(f1);

            subplot(4,1,1)
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

            subplot(4,1,2)
            hold all
            plot(L_r * Loo *0.338, Force_WLC,':', 'Color','k');

            plot(LDNAnm_s, ForcepN_s,'o', 'Color','r', 'MarkerFaceColor','r', 'MarkerSize',2);
            plot(LDNAnm_ss, ForcepN_ss, '-k');
            axis([400 5500 0 125 ]);

            set(gca,'FontSize',13,'FontName','Calibri');
            xlabel('Extension (nm)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');
            title(['Extension @ 1,2, 3 pN: ' num2str(lt(j, 1)) ', ' num2str(lt(j, 2)) ', ' num2str(lt(j, 3)) ' nm'],'Interpreter','None','FontSize',10);        

            subplot(4,1,3)
            hold all
            plot(Times, ForcepN_s,'--', 'Color',[0.5,0.5,0.5]);
            plot(Times, ForcepN_s,'o', 'Color','k', 'MarkerFaceColor','k', 'MarkerSize',1);
            plot(t_max, f_max,'or','MarkerFaceColor','r');
                        
            line([Times(index_mean(1)) Times(index_mean(1))], [0 120]);
            line([Times(index_mean(end)) Times(index_mean(end))], [0 120]);
            
            ylim([0 125]);
            title(['F_max: ' num2str(f_max) ' pN; F_mean: ' num2str(f_mean) ' pN'],'Interpreter','None','FontSize',10);        

            set(gca,'FontSize',13,'FontName','Calibri');
            xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Force (pN)', 'FontSize',13, 'FontName', 'Calibri');

            
            force_rate = movingslope(ForcepN_s, 100, 1,dt);
            subplot(4,1,4)
            hold all
            %l3 = line([Times(1) Times(end)], [10 10]); set(l3, 'Color',[0.5, 0.5, 0.5]);
            plot(Times, force_rate,'--.', 'Color','k', 'MarkerSize', 3, 'MarkerFaceColor','k');

            ylim([0 100]);

            set(gca,'FontSize',13,'FontName','Calibri');
            xlabel('Time (s)', 'FontSize',13, 'FontName', 'Calibri');
            ylabel('Force changing rate (pN/s)', 'FontSize',13, 'FontName', 'Calibri');


            savefig(f1,[pn '_check.fig']);

            tracenuc{j}.name = pn;
            tracenuc{j}.Times1 = Times1;
            tracenuc{j}.ForcepN1 = ForcepN1;
            tracenuc{j}.LDNAnm1 = LDNAnm1;
            tracenuc{j}.fit_Xcentering = fit_Xcentering;
            tracenuc{j}.l0p5 = l0p5;


            tracenuc{j}.Times = Times;
            tracenuc{j}.ForcepN = ForcepN;
            tracenuc{j}.LDNAnm = LDNAnm;
            tracenuc{j}.Times_s = Times_s;
            tracenuc{j}.ForcepN_s = ForcepN_s;
            tracenuc{j}.LDNAnm_s = LDNAnm_s;
            tracenuc{j}.Times_ss = Times_ss;
            tracenuc{j}.ForcepN_ss = ForcepN_ss;
            tracenuc{j}.LDNAnm_ss = LDNAnm_ss;
            tracenuc{j}.t_on_pick = t_on_pick;
            tracenuc{j}.index_high = index_high;

            tracenuc{j}.f_max = f_max;
            tracenuc{j}.t_max = t_max;
            tracenuc{j}.l_max = l_max;
            tracenuc{j}.index_max = index_max;
            
            tracenuc{j}.index_mean = index_mean;
            tracenuc{j}.f_mean = f_mean;
            
            tracenuc{j}.j_bp = j_bp;

            tracenuc{j}.ft = ft;
            tracenuc{j}.lt = lt(j, :);



            save('data.mat','tracenuc');
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
    
    time_t = 5;
    ForcepN = ForcepN(Times>time_t);
    LDNAnm = LDNAnm(Times>time_t);
    Times = Times(Times>time_t);
    
    %index0 = find(ForcepN == max(ForcepN));
    
    %Times = Times(1: index0(1)-1);
    %ForcepN = ForcepN(1: index0(1)-1);
    %LDNAnm = LDNAnm(1: index0(1)-1);
    
end
