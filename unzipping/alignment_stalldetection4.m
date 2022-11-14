function [jend_a, fmean_a, fmax_a, name_a] = alignment_stalldetection4(filename_t)
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

    count = 0;
    countdna = 0;
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
        f0 = figure;
        pos = [100 100 1550 300];
        set(f0, 'Pos', pos);
        hold on;
        
        [Times, jindex, ForcepN, fmax, fmean, jend] = Show_each_trace(fid, filename_t, f0);
        
            trace{j}.Times = Times;
            trace{j}.jindex = jindex;
            trace{j}.ForcepN = ForcepN;
            trace{j}.fmax = fmax;
            trace{j}.fmean = fmean;
            trace{j}.jend = jend;
            
            jend_a(j) = jend;
            fmean_a(j) = fmean;
            fmax_a(j) = fmax;
        
        qstring = 'Does the trace look good?';
        choice = questdlg(qstring,'Question?','Yes','No','Yes');
        
            figure(f0)
            savefig(f0,[fid '.fig']);
            export_fig(f0,[fid '.png']);
        
        if strcmp(choice, 'Yes')
            trace{j}.selected = 'yes';
            
            choice2 = 'No';
            
            while strcmp(choice2, 'No')
                %% Choose a range
                prompt = {'j1min:','j1max:','j2min:','j2max:'};
                titlej = 'j_index range';
                num_line = 1;
                defaultvals = {'100', '3000', '100','3000'};
                answer = inputdlg(prompt, titlej,num_line,defaultvals);
                range = cellfun(@str2num, answer)';
                
                trace{j}.range = range;

                %% Perform the allignment and save data to new file
                f1 = figure;
                set(f1, 'Pos', pos);
                
                hold all;    
                
                [Times_n, jindex_n_old, jindex_n, ForcepN_n] = Unzipping_trace_nuc(fid,filename_t, f1, range);
                
                figure(f1);
                %title(['Force vs. j-index ' fid ', time: ' time], 'FontSize',13);
                title(['Force vs. j-index ' fid], 'FontSize',13);
                
                qstring2 = 'Does the allignment look good?';
                choice2 = questdlg(qstring2,'Question?','Yes','No','No');
            end
            trace{j}.Times_n = Times_n;
            trace{j}.jindex_n_old = jindex_n_old;
            trace{j}.jindex_n = jindex_n;
            trace{j}.ForcepN_n = ForcepN_n;
            
            close(f0);
            
            % Plot figure
            figure(f1);
            set(gca,'FontSize',13,'LineWidth',1.5);
            
            trace{j}.type = 'DNA';
            qstring3 = 'Free DNA or not?';
            choice3 = questdlg(qstring3,'Question?','Yes','No','No');
            
            if strcmp(choice3, 'No')
                trace{j}.type = 'protein';
                

                [peak_lo, peak_force, peak_halfwidth] = find_peaks(movingmean(jindex_n,40), movingmean(ForcepN_n,40));
                
                
                trace{j}.peak_lo = peak_lo;
                trace{j}.peak_force = peak_force;
                trace{j}.peak_halfwidth = peak_halfwidth;
            
            
                
                figure(f1)
                hold all
                plot(peak_lo, peak_force, '+k','MarkerSize',9);
               


                        prompt = {'Comment?'};
                        titlea = 'Add comment';
                        num_line = 1;
                        defaultvals = {'break'};
                        comments = inputdlg(prompt, titlea,num_line,defaultvals);
                        
                        count = count + 1;

                        data3 = [peak_lo;peak_force;peak_halfwidth]; data3 = data3(:)'; 
                        data3 = data3;
                        
                        data4 = [fid comments num2cell(data3)];%

                        disp(data4);
                        if ~isempty(data4)
                            xlswrite('protein.xls', data4, 1, ['A' num2str(count)]);
                        end
                 
            else
                    countdna = countdna + 1;
                    
                        data4 = [fid num2cell( jindex_n(end))];
                    
                    xlswrite('DNA.xls', data4, 1, ['A' num2str(countdna)]);
            end
            
            figure(f1)
            savefig(f1,[fid '_aligned.fig']);
            export_fig(f1,[fid '_aligned.png']);
            
        else
            trace{j}.selected = 'no';
        end
        
        save('data','trace');
    end
    
    %save('data1','trace');
end
end

function [peak_lo, peak_force, peak_halfwidth] = find_peaks(x, y)

index = y>=20 & x <= 4025;
x1 = x(index);
y1 = y(index);

if any(x1) && any(y1)

    djr = 5; jr = -100 : djr : 4025;
    nj = histc(x1, jr);
    jrn = jr + djr/2;

    % finding force peak, location of force peak, disruption force

    n_t = 10;
    index_f = find(nj >= n_t);

    j = 1;
    bin = {};
    bin{1} = index_f(1);
    for i = 1 : length(index_f)-1
        if index_f(i+1) - index_f(i) <= 2.5   % 5 bp to considered as different interactions
            bin{j} = [bin{j} index_f(i+1)];
        else
            j = j+1;
            bin{j} = index_f(i+1);
        end
    end

    L = size(bin,2);

    peak_lo = zeros(1, L);
    peak_force = zeros(1, L);
    peak_halfwidth = zeros(1, L);

    dj_search = 10;
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];

    for i = 1 : L
        pos_temp = mean(jrn(bin{i}));

        index_p_local = abs(x1 - pos_temp) <= dj_search;

        peak_force(i) = max(y1(index_p_local));

        %range = pos_temp - dj_search : 2 : pos_temp + dj_search;
        %n = hist(x1(index_p_local), range);
        %v = range + 1;

        %[xData, yData] = prepareCurveData( v, n );

        %opts.StartPoint = [max(n) pos_temp 3];    
        %[fitresult, gof] = fit( xData, yData, ft, opts );

        %peak_lo(i) = fitresult.b1;
        %peak_halfwidth(i) = fitresult.c1/sqrt(2);
        
        peak_lo(i) = mean(x1(index_p_local));
        peak_halfwidth(i) = std(x1(index_p_local));
    end

    peak_forcerise = peak_lo - 3*peak_halfwidth;

   % figure;
   % hold all
   % plot(x,y);
   % plot(peak_forcerise, peak_force,'+r');
else
    peak_lo = [];  peak_force = []; peak_halfwidth = [];
end
end

function [Times, jindex, ForcepN, fmax, fmean, jend] = Show_each_trace(filename, filename_t, f)
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
    
    %clearvars filename_t delimiter startRow formatSpec fileID_t dataArray ans;
    clearvars  delimiter startRow formatSpec fileID_t dataArray ans;
    
    
    
    %% Filter data based on the AOMSent Magnitude and the length of the sequence
    index1 = find(AOMSentV > 0.1 & ForcepN > 9) ; 
    %% Discard the last np data points
    
    jindex = jindex(index1); ForcepN = ForcepN(index1);Times = Times(index1); 
    jindex = jindex - median(jindex(ForcepN < 11 & ForcepN > 9));
    
    index2 =  ForcepN > 10;
    jindex = jindex(index2); Times = Times(index2); ForcepN = ForcepN(index2);
    
    jindex_s = movingmean(jindex,5); ForcepN_s =   movingmean(ForcepN,5);
    
    %np = 10;
    %Times_n(end-np:end) = []; jindex_n(end-np:end) = []; ForcepN_n(end-np:end) = [];
    
    
    fmax = max(ForcepN_s);
    fmean = mean(ForcepN_s);
    jend= jindex(end-4);
    
%% Ploting
figure(f);
            

plot(jindex_t, forcepN_t,'Color',[0.5, 0.5, 0.5], 'LineWidth',1.5,'LineStyle','-');   

plot(jindex_s,ForcepN_s, 'o-','LineWidth',1,'Color',[0.5, 0.5, 0.5], 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

axis([-100 4200 5 65]);

textp1 = { ['Max force: ' , num2str(fmax) ' pN'] ,['Mean force: ' num2str(fmean) ' pN'] ,['j-index end @: ' num2str(jend) ' bp']};   
text(0, 65, textp1, 'HorizontalAlignment', 'left','VerticalAlignment', 'top','FontSize',11, 'FontWeight','bold', 'BackgroundColor','w', 'EdgeColor','k');

set(gca,'XGrid','on');
title(filename);

%set(gca, 'XTick', -100:200:4600);
xlabel('j-index (bp)', 'FontSize',13); 
ylabel('Force (pN)', 'FontSize',13);
title(['Force vs. j-index ' filename], 'FontSize',13);



end