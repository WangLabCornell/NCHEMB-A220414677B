function [jend_a, fmean_a, fmax_a, j_fmax_a, j_first_a, f_first_a, index_dna, index_protein, index_notallign, name_a] = alignment_stalldetection5(filename_t)
filename = uigetfile('*.*','Select the DNA file','multiselect','on');
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

    jend_a = zeros(1, L);
    fmean_a = zeros(1, L);
    fmax_a = zeros(1, L);
    j_fmax_a = zeros(1, L);
    name_a = {};
    j_first_a = zeros(1, L);
    f_first_a = zeros(1, L);
    index_dna = zeros(1, L);
    index_protein = zeros(1, L);
    index_notallign = zeros(1, L);

    trace = {};
    while j < L && ~strcmp(button,'Cancel')
        close all;
        j=j+1;
        
        if ~iscellstr(filename)
            fid = filename;
        else
            fid = filename{j};
        end
        
        trace{j}.name = fid;
            
        name_a{j} = fid;
        %% Time of the trace
        %namefile = name(1:11);
        %namefile_n = ['..\..\..\' namefile '\configuration.ini'];
        %delimiter = ' ';
        %startRow = 51;
        %endRow = 51;
        %formatSpec = '%*s%*s%*s%s%*s%*s%[^\n\r]';
        %fileID = fopen(namefile_n,'r');
        %textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
        %dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
        %fclose(fileID);   
        %VarName2 = dataArray{:, 1}; time = VarName2{1};time([1, end]) =[];
        
        %trace{j}.time = time;
        
        %% Show the trace     
             [Times, jindex, ForcepN, jindex_s, ForcepN_s,  forcepN_t, jindex_t] = Show_each_trace(fid, filename_t);
        
            trace{j}.Times = Times;
            trace{j}.jindex = jindex;
            trace{j}.ForcepN = ForcepN;
            trace{j}.jindex_s = jindex_s;
            trace{j}.ForcepN_s = ForcepN_s;
            
            trace{j}.jindex_t = jindex_t;
            trace{j}.forcepN_t = forcepN_t;
            

                        %% Fitting the theoretical curve using a spline
                        [xData_theo, yData_theo] = prepareCurveData( jindex_t, forcepN_t );
                        yData_theo = yData_theo; 
                        fitresult_theo = splinefit(xData_theo, yData_theo,round(max(xData_theo)/5));
                        
                        dt = mean(diff(Times));
                        N = round(1/40/dt);     % resampling to 40 Hz
                        jindex_s2 = movingmean(jindex,N);
                        ForcepN_s2 = movingmean(ForcepN,N);

                        % Fit exp. data using spline function
                        [xData_exp, yData_exp] = prepareCurveData( jindex_s2, ForcepN_s2 );
                        fitresult_exp = splinefit( xData_exp, yData_exp, round(max(xData_exp)/3));
                        
            %% Ploting
            f0 = figure;
            pos = [100 100 1550 300];
            set(f0, 'Pos', pos);
            hold on;
            plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5,'LineStyle','-');   
            plot(jindex_s,ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

            axis([-190 4500 10 65]);

            
            set(gca,'XGrid','on');
            title(fid);
            xlabel('j-index (bp)', 'FontSize',13); 
            ylabel('Force (pN)', 'FontSize',13);
            title(['Force vs. j-index ' fid], 'FontSize',13);
        
            qstring = 'Does the trace look good?';
            choice = questdlg(qstring,'Question?','Yes','No','Yes');  

        
            if strcmp(choice, 'Yes')
                figure(f0);
                title([fid ', good' ]);
                savefig(f0,[fid '.fig']);
                export_fig(f0,[fid '.png']);
                trace{j}.selected = 'yes';
                
                qstring2 = 'Can the trace be alligned?';
                choice2 = questdlg(qstring2,'Question?','Yes','No','Yes'); 

                if strcmp(choice2, 'Yes')
                    
                    choice3 = 'No';
                    while strcmp(choice3, 'No')
                        %% Choose a range
                        prompt = {'j1min:','j1max:','j2min:','j2max:'};
                        titlej = 'j_index range';
                        num_line = 1;
                        defaultvals = {'100', '3500', '100','3500'};
                        answer = inputdlg(prompt, titlej,num_line,defaultvals);
                        range = cellfun(@str2num, answer)';

                        trace{j}.range = range;

                        %% Search for the alignment
                        
                        a = linspace(0.5, 1.7,35);
                        j0 = linspace(-50, 50, 35);
                        N1 = length(a); N2 = length(j0);
                        R = zeros(N1, N2);
                        for i = 1 : N1
                            for k= 1 : N2
                                    R(i, k ) = crosscorrelationR([a(i), j0(k)], range, fitresult_theo, fitresult_exp);
                            end
                        end
                        R2 = R;
                        R2(R2<=0.5) = 0.5;
                        R2m = imregionalmax(R2,8);

                        [i_at, i_j0t] = ind2sub(size(R2m), find(R2m == 1));
                        at = a(i_at); j0t = j0(i_j0t);
                        if length(at) > 10
                            at(11:end) = []; j0t(11:end) = [];
                        end
                        coefi = zeros(length(at),2);
                        R3 = zeros(length(at),1);
                        for i = 1 : length(at)
                            coefi(i,:) = fminsearch(@(x)(1-crosscorrelationR(x, range, fitresult_theo, fitresult_exp)),[at(i), j0t(i)]);
                            R3(i) = crosscorrelationR(coefi(i,:),  range, fitresult_theo, fitresult_exp);
                        end
                        atm = coefi(:,1); j0tm = coefi(:,2);
                        a_f = atm(R3 == max(R3)); j0_f = j0tm(R3 == max(R3));
                        a_f = a_f(1); j0_f = j0_f(1);
                        
                        %% Moving the experimental data
                        jindex_n = (jindex - j0_f) ./ a_f;
                        jindex_s2_n = (jindex_s2 - j0_f) ./ a_f;
                        jindex_s_n = (jindex_s - j0_f) ./ a_f;

                        f1 = figure;
                        pos = [100 100 1550 300];
                        set(f1, 'Pos', pos);
                        hold on;
                        title(['Force vs. j-index ' fid], 'FontSize',13);
                        plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5);     % Plotting the fitted theoretical curve
                        plot(jindex_s_n, ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                         
                        axis([-190 4500 10 65]);
                        xlabel('j-index (bp)', 'FontSize',13); 
                        ylabel('Force (pN)', 'FontSize',13);
                        set(gca,'FontSize',13,'LineWidth',1.5,'YGrid','on');

                        qstring3 = 'Does the allignment look good?';
                        choice3 = questdlg(qstring3,'Question?','Yes','No','No');
                    end

                    trace{j}.jindex_s2 = jindex_s2;
                    trace{j}.ForcepN_s2 = ForcepN_s2;
                    
                    trace{j}.jindex_n = jindex_n;
                    trace{j}.jindex_s_n = jindex_s_n;
                    trace{j}.jindex_s2_n = jindex_s2_n;
    
                    close(f0);

                    qstring4 = 'Free DNA or not?';
                    choice4 = questdlg(qstring4,'Question?','Yes','No','No');
                    
                        fmax = max(ForcepN_s);
                        fmean = mean(ForcepN_s);
                        j_fmax = jindex_n(ForcepN_s == max(ForcepN_s));
                        jend = jindex_n(end-4);
                   

                    if strcmp(choice4, 'No')
                        trace{j}.type = 'protein';
                        index_protein(j) = 1;
                        
                        j_above = jindex_s_n(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+3 & jindex_s_n >= 0);
                        
                        if any(j_above)
                        
                            if length(j_above) >=50
                                j_first = median(j_above(1:50));
                            else
                                j_first = median(j_above(1:end));
                            end
                            f_first = median( ForcepN_s( abs(jindex_s_n-j_first )<= 25 ) ) ;
                        else
                            j_first = j_fmax;
                            f_first = fmax;
                        end
                        
                        trace{j}.j_first = j_first;
                        trace{j}.f_first = f_first;
                        
                        j_first_a(j) = j_first;
                        f_first_a(j) = f_first;
                        
                        figure(f1);
                        hold all;
                        plot(j_first, f_first, '+b', 'MarkerSize',10, 'LineWidth', 1.5);
                        
                    else
                        trace{j}.type = 'DNA';              
                        index_dna(j) = 1;
                        
                        
                        j_first = 0;
                        f_first = 15 ;
                        
                        trace{j}.j_first = j_first;
                        trace{j}.f_first = f_first;
                        
                        j_first_a(j) = j_first;
                        f_first_a(j) = f_first;
                        
                        figure(f1);
                        hold all;
                        plot(j_first, f_first, '+b', 'MarkerSize',10, 'LineWidth', 1.5);
                        
                        
                    end
                    
                        fmax = max(ForcepN_s);
                        fmean = mean(ForcepN_s);
                        j_fmax = jindex_n(ForcepN_s == max(ForcepN_s));
                        jend = jindex_n(end-4);
                            
                        trace{j}.fmax = fmax;
                        trace{j}.fmean = fmean;
                        trace{j}.jend = jend;
                        trace{j}.j_fmax = j_fmax;

                        figure(f1);
                        hold all;
                        plot(j_fmax, fmax, '+k', 'MarkerSize',10, 'LineWidth', 1.5);
                        textp1 = { ['Max force: ' , num2str(fmax) ' pN'] ,['Mean force: ' num2str(fmean) ' pN'] ,['j-index end @: ' num2str(jend) ' bp']};   
                        text(4000, 65, textp1, 'HorizontalAlignment', 'left','VerticalAlignment', 'top','FontSize',11, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','k');
                        savefig(f1,[fid '_aligned.fig']);
                        export_fig(f1,[fid '_aligned.png']);
                        
                        %close(f1);
                
                else
                        trace{j}.type = 'notallign'; 
                        index_notallign(j) = 1;
                        jindex_n = jindex;
                        jindex_s_n = jindex_s;
                        trace{j}.jindex_n = jindex_n;
                        trace{j}.jindex_s_n = jindex_s;
                        
                    
                        fmax = max(ForcepN_s);
                        fmean = mean(ForcepN_s);
                        j_fmax = jindex_n(ForcepN_s == max(ForcepN_s));
                        jend = jindex_n(end-4);
                            
                        trace{j}.fmax = fmax;
                        trace{j}.fmean = fmean;
                        trace{j}.jend = jend;
                        trace{j}.j_fmax = j_fmax;
                        close(f0);
                        
                        j_above = jindex_s_n(ForcepN_s > ppval(fitresult_theo, jindex_s_n)+3 & jindex_s_n >= 0);
                        j_first = median(j_above(1:50));
                        f_first = median( ForcepN_s( abs(jindex_s_n-j_first )<= 25 ) ) ;
                        
                        trace{j}.j_first = j_first;
                        trace{j}.f_first = f_first;
                        
                        j_first_a(j) = j_first;
                        f_first_a(j) = f_first;
                        
                        
                        f1 = figure;
                        pos = [100 100 1550 300];
                        set(f1, 'Pos', pos);
                        hold on;
                        plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5);     % Plotting the fitted theoretical curve
                        plot(jindex_s_n, ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                        plot(j_fmax, fmax, '+k', 'MarkerSize',10 , 'LineWidth', 1.5);
                        plot(j_first, f_first, '+b', 'MarkerSize',10, 'LineWidth', 1.5);
                         
                        axis([-190 4500 10 65]);
                        xlabel('j-index (bp)', 'FontSize',13); 
                        ylabel('Force (pN)', 'FontSize',13);
                        title(['Force vs. j-index ' fid], 'FontSize',13);
                        set(gca,'FontSize',13,'LineWidth',1.5,'YGrid','on');
                        
                        textp1 = { ['Max force: ' , num2str(fmax) ' pN'] ,['Mean force: ' num2str(fmean) ' pN'] ,['j-index end @: ' num2str(jend) ' bp']};   
                        text(4000, 65, textp1, 'HorizontalAlignment', 'left','VerticalAlignment', 'top','FontSize',11, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','k');
                        savefig(f1,[fid '_aligned.fig']);
                        export_fig(f1,[fid '_aligned.png']);
                        
                        %close(f1);
                end
                

                
                        
                    jend_a(j) = jend;
                    fmean_a(j) = fmean;
                    fmax_a(j) = fmax;
                    j_fmax_a(j) = j_fmax;
            
            else
                            
                figure(f0);
                axis([-190 4500 10 65]);
                title([fid ', bad' ]);
                savefig(f0,[fid '.fig']);
                export_fig(f0,[fid '.png']);
                trace{j}.selected = 'no';
                close(f0);
            end
        
        save('data','trace');
    end
    
    index_dna = logical(index_dna);
    index_protein = logical(index_protein);
    index_notallign = logical(index_notallign);
    
    %save('data1','trace');
end
end



function [Times, jindex, ForcepN, jindex_s, ForcepN_s,   forcepN_t, jindex_t] = Show_each_trace(filename, filename_t)
%% As of 07/15/2015 - Tung Le

    %% Open the file with experimental data
    fileID = fopen(filename,'r');
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    
    Times = dataArray{:, 1};
    AOMSentV = dataArray{:, 13};
    ForcepN = dataArray{:, 23};
    jindex = dataArray{:, 24};
    clearvars delimiter startRow formatSpec fileID dataArray ans;
    
        %% Open the theoretical curve
    %filename_t = 'force_kkb_1p6kbafterA20_100mMNaCl_4mMMg';
    delimiter = '\t';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%[^\n\r]';
    fileID_t = fopen(filename_t,'r');
    dataArray = textscan(fileID_t, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID_t);

    extensionnm = dataArray{:, 1};
    forcepN_t = dataArray{:, 2};
    jindex_t = dataArray{:, 3};
   
    clearvars  delimiter startRow formatSpec fileID_t dataArray ans;
 
    %% Filter data based on the AOMSent Magnitude and the length of the sequence
    index1 = find(AOMSentV > 0.1 & ForcepN > 9) ; 
    %% Discard the last np data points
    
    jindex = jindex(index1); ForcepN = ForcepN(index1);Times = Times(index1); 
    jindex = jindex - median(jindex(ForcepN < 11 & ForcepN > 9));
    
    index2 =  ForcepN > 9;
    jindex = jindex(index2); Times = Times(index2); ForcepN = ForcepN(index2);
    
    jindex_s = movingmean(jindex,5); ForcepN_s =   movingmean(ForcepN,5);
    

    
end

function R = crosscorrelationR(x, range0, fitresult_theo, fitresult_exp)
% Choose the range of the data
    a = x(1); j0 = x(2);
    bin = 1;
    N = length(range0);
    
    if N == 2 
            range1 = range0(1): bin : range0(2);
            ftheo_mean = mean(ppval(fitresult_theo, range1));
            fexp_mean = mean(ppval(fitresult_exp, a.*range1 + j0));
        R = (trapz(range1, (ppval( fitresult_theo, range1) - ftheo_mean) .* (ppval(fitresult_exp, a .* range1 + j0) - fexp_mean)) )...
            / sqrt( trapz(range1, (ppval(fitresult_theo, range1) - ftheo_mean).^2 ) * trapz(range1, (ppval(fitresult_exp, a .* range1 + j0) - fexp_mean).^2 ) );
    else
        range1 = range0(1): bin : range0(2);
        range2 = range0(3) : bin : range0(4);

        ftheo_mean = mean(ppval(fitresult_theo, range1));
        ftheo2_mean = mean(ppval(fitresult_theo, range2));
        fexp_mean = mean(ppval(fitresult_exp, a.*range1 + j0));
        fexp2_mean = mean(ppval(fitresult_exp, a.*range2 + j0));
        R = (trapz(range1, (ppval( fitresult_theo, range1) - ftheo_mean) .* (ppval(fitresult_exp, a .* range1 + j0) - fexp_mean)) + trapz(range2, (ppval( fitresult_theo, range2) - ftheo2_mean) .* (ppval(fitresult_exp, a .* range2 + j0) - fexp2_mean))  )...
            / sqrt( ( trapz(range1, (ppval(fitresult_theo, range1) - ftheo_mean).^2 ) + trapz(range2, (ppval(fitresult_theo, range2) - ftheo2_mean).^2 )) * ( trapz(range1, (ppval(fitresult_exp, a .* range1 + j0) - fexp_mean).^2 ) + trapz(range2, (ppval(fitresult_exp, a .* range2 + j0) - fexp2_mean).^2 )) );
    end
end