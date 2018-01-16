%% PURPOSE: Select seasons in which the Madden Julian Oscillation (MJO) is 
% significant conditioned by the local precipitation (Uruguay)
%
%
% This script is part of a thesis to obtain the bachelor degree in
% Atmospheric Sciences - University of the Republic (Uruguay).
%
% Author: Paulina Tedesco
% Tutor: Rafael Terra
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figures of the annual cycle of precipitation in the different
% meteorological station conditioned by the quartiles of the MJO indices
% Index 4, 8 and 9. The significant pentads (5% and 95%) are shaded based
% on a script provided by Alejandra de Vera.
%
% The signs of the difference between the quartiles and the confidence
% intervals are plotted for each meteorological station and each MJO Index.
% 1 means positive and -1 means negative.
%
% The function nanmoving_average.m is used.
% Copyright 2006-2008  Carlos Vargas, nubeobscura@hotmail.com
% $Revision: 3.1 $  $Date: 2008/03/12 18:20:00 $
% 
% The MJO Indices were obtained from the Climate Prediction Center (CPC) 
% of the National Oceanic and Atmospheric Administration (NOAA). See:
% http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_mjo_index/pentad.shtml
%
% Daily precipitation data (1979-2015) was provided by Instituto Uruguayo 
% de Meteorología(INUMET) and Instituto Nacional de Investigación 
% Agropecuaria (INIA).
%
% https://inumet.gub.uy/
% http://www.inia.uy/gras/Clima/Banco-datos-agroclimatico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1) Load the data

close all
clear all

A = csvread('proj_norm_order.csv',1,1); % MJO Indices(idx 4, 8, 9) 
B = csvread('dailyPP.csv',1,3); % Daily precipitation from 5 meteorological stations located in Uruguay.

[ma,na] = size(A);
[mb, nb] = size(B);


%% 2) Transform daily precipitation into pentads (PPpentads)

PPpentads = zeros(mb/5,nb);
j = 1;
for i=1:5:mb
    PPpentads(j,:) = nansum( B(i:i+4, :) );
    j = j+1;
end

%% 3) Create new matrix of pentads of size 73x37 (for precipitation and MJO indices)

PPpentads_resh = reshape(PPpentads,[],(mb/365),nb); % pentads x years x station 
idx_resh = reshape(A,[],(mb/365),na); % pentads x years x station 

%% 4) Annual cycle of precipitation in pentads - (only used in the graphs)

% Moving average of precipitation in pentads
PPpentads_movave = NaN(mb/5,nb);
for j = 1:nb
    [PPpentads_movave(:,j),Nsum] = nanmoving_average(PPpentads(:,j),5);
end

% reshape
PPpentads_movave_resh = reshape(PPpentads_movave,[],(mb/365),nb); % pentads x years x station 

% annual cycle
PPpentads_movave_anual = squeeze(nanmean(PPpentads_movave_resh,2));



%% 5) Copy the first/last 5 rows at the end/beginning 
% This is necessary for creating a matrix of size 73x(37x11).

PPpentads_resh_cat = cat(1, PPpentads_resh(end-4:end,:,:), PPpentads_resh, PPpentads_resh(1:5,:,:));
idx_resh_cat = cat(1, idx_resh(end-4:end,:,:),idx_resh, idx_resh(1:5,:,:));

% Window of 11 pentads
PPpentads_11 = NaN(73,37*11, nb); % 73 pentads x (37x11) years x 5 stations
idx_11 = NaN(73,37*11, na); % 73 pentads x (37x11) years x 3 indices

for i = 6:(73+5)
    for j =1:37
        PPpentads_11(i-5, j*11-10:j*11,:) = permute( PPpentads_resh_cat(i-5:i+5,j,:), [2,1,3] ) ;
        idx_11(i-5, j*11-10:j*11,:) = permute( idx_resh_cat(i-5:i+5,j,:), [2,1,3] ) ;
    end
end

%% 6) sort MJO indices (ascending)

I = NaN(73,407, na); % positions of the sorted indices 
idx_sorted = NaN(73,407, na); % sorted indices

for j=1:na
    [idx_sorted(:,:,j), I(:,:,j)] = sort(idx_11(:,:,j),2);
end

% Define first and last quartiles of the sorted indices
I_Q1 = I(:, 1:102, :); % position of the first 102 sorted indices 
I_Q4 = I(:, end-101:end, :); % position of the last 102 sorted indices


%% 7) q1 and q4 are the annual cycle of the precipitacion pentads in each quartile

q1 = NaN(73,3,nb);
q4 = NaN(73,3,nb);

for i=1:73
    for j=1:na % indices 4, 8, 9
        for k=1:nb % meteorological stations 
            q1(i,j,k) = nanmean( PPpentads_11(i, I_Q1(i, :,j), k), 2 );
            q4(i,j,k) = nanmean( PPpentads_11(i, I_Q4(i, :,j), k), 2 );
        end
    end
end


%% 8) Confidence intervals (Monte Carlo Method)

int = ones(73,nb,1000);

for l = 1:1000
    rand_idx = randi(407,1,102); % randomly select 102 out of 407 positions
    for i = 1:73
        for k = 1:nb
            % average the precipitation in those indices
            int(i,k,l) = nanmean( PPpentads_11(i, rand_idx, k), 2 ); 
        end
    end 
end

% 5% confidence
int5 = quantile(int, 0.05,3); % 0.05 quantil of averaged precipitation

% 95% confidence
int95 = quantile(int, 0.95,3); % 0.95 quantil of averaged precipitation


        
%% 9) Figures

% Meteorological stations
str_stations={'Artigas', 'Paysandú', 'Melo', 'Rocha', 'Las Brujas'};

% MJO Indices
str_idx={'Index4', 'Index8', 'Index9'};


% Plot annual cycles and shade significant regions

for j = 1:na % MJO Index
    for k = 1:nb % Meteorological station
        
        figure
        Q1 = q1(:,j,k); % First quartile
        Q4 = q4(:,j,k); % Last quartile
        clima = PPpentads_movave_anual(:,k); % climatology
        p5 = int5(:,k); % 5% confidence interval
        p95 = int95(:,k); % 5% confidence interval
        
        hold on
        plot(q1(:,j,k),'linewidth', 2)
        hold on
        plot(q4(:,j,k),'r', 'linewidth', 2)
        plot(PPpentads_movave_anual(:,k), 'g', 'linewidth', 2)
        ylim([0,45])
        xlim([0,74])
        plot(int5(:,k), '--', 'Color', [0.4,0.4,0.4], 'linewidth', 2)
        plot(int95(:,k), '--', 'Color', [0.4,0.4,0.4], 'linewidth', 2)

      
        % Shade significant regions
        if mean(Q1)>mean(Q4)

            % first option
            aux = find(Q1>p95);
            trozos=[0;find(diff(aux)>1);length(aux)];
            for i=1:length(trozos)-1
                ini=trozos(i)+1;fin=trozos(i+1);
                fill([aux(ini:fin);aux(fin:-1:ini)],[Q1(aux(ini:fin));p95(aux(fin:-1:ini))],...
                    [.7294 0.8314 0.9569],'LineStyle','none')
            end
            hold on
            aux = find(Q4<p5);
            trozos=[0;find(diff(aux)>1);length(aux)];
            for i=1:length(trozos)-1
                ini=trozos(i)+1;fin=trozos(i+1);
                fill([aux(ini:fin);aux(fin:-1:ini)],[Q4(aux(ini:fin));p5(aux(fin:-1:ini))],...
                    [1 0.8 0.8],'LineStyle','none')
            end

        else
            hold on
            % Second option
            aux = find(Q1<p5);
            trozos=[0;find(diff(aux)>1);length(aux)];
            for i=1:length(trozos)-1
                ini=trozos(i)+1;fin=trozos(i+1);
                fill([aux(ini:fin);aux(fin:-1:ini)],[Q1(aux(ini:fin));p5(aux(fin:-1:ini))],...
                    [.7294 0.8314 0.9569],'LineStyle','none')
            end
            hold on
            aux = find(Q4>p95);
            trozos=[0;find(diff(aux)>1);length(aux)];
            for i=1:length(trozos)-1
                ini=trozos(i)+1;fin=trozos(i+1);
                fill([aux(ini:fin);aux(fin:-1:ini)],[Q4(aux(ini:fin));p95(aux(fin:-1:ini))],...
                    [1 0.8 0.8],'LineStyle','none')
            end

        end
  
        
        hold on
        plot(q1(:,j,k),'linewidth', 2)
        hold on
        plot(q4(:,j,k),'r', 'linewidth', 2)
        plot(PPpentads_movave_anual(:,k),'Color', [0 204 0] ./ 255, 'linewidth', 2)
        ylim([0,45])
        xlim([0,74])
        plot(int5(:,k), '--', 'Color', [0.4,0.4,0.4], 'linewidth', 2)
        plot(int95(:,k), '--', 'Color', [0.4,0.4,0.4], 'linewidth', 2)
        grid
        set(gca, 'XTick', [0 1 2 3 4 5 6 7 8 9 10 11 12]*73/12+1);
        set(gca, 'XTickLabel', {'J','F','M','A','M','J', 'J', 'A', 'S', 'O', 'N', 'D', 'E'}, ...
            'FontSize', 10, 'FontWeight', 'bold') 
        str = [str_stations(k), str_idx(j)];
        title(str, 'FontSize', 12, 'FontWeight', 'bold' )
        legend('Q1', 'Q4', 'Clima', '5%', '95', 1, 'Orientation','horizontal')
        ylabel('Precipitation (mm)', 'FontSize', 12, 'FontWeight', 'bold')
        str = [str_stations(k), str_idx(j)];
        title(str, 'FontSize', 12, 'FontWeight', 'bold' )
    end
end




%%

% Find significant pentads and plot them
t = 1:73; % Pentads
tt = zeros(48,73); % Intersections
sg = zeros(48,73); % Signs
r = 1;

for k = 1:nb
    i = 1;
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for j = 1:na
        subplot(6,2,i)
        a = sign(int5(:,k) - q1(:,j,k))';
        ii = sort([strfind(a,[-1,1]),strfind(a,[1,-1])]);
        plot(a,'linewidth', 2)
        tt(r, 1:length(ii)) = t(ii);
        sg(r,:) = a;
        ylim([-1.5,1.5])
        xlim([0,74])
        hold on
        str = [str_idx(j), 'sg(5 - Q1%)'];
        title(str)
        
        subplot(6,2,i+1)
        a = sign(q1(:,j,k) - int95(:,k))';
        ii = sort([strfind(a,[-1,1]),strfind(a,[1,-1])]);
        plot(a,'linewidth', 2)
        tt(r+1, 1:length(ii)) = t(ii);
        sg(r+1,:) = a;
        ylim([-1.5,1.5])
        xlim([0,74])
        str = [str_idx(j), 'sg(Q1 - 95%)'];
        title(str)
        
        subplot(6,2,i+2)
        a = sign(int5(:,k) - q4(:,j,k))';
        ii = sort([strfind(a,[-1,1]),strfind(a,[1,-1])]);
        plot(a,'linewidth', 2)
        tt(r+2, 1:length(ii)) = t(ii);
        sg(r+2,:) = a;
        ylim([-1.5,1.5])
        xlim([0,74])
        str = [str_idx(j), 'sg(5% - Q4)'];
        title(str)
        
        subplot(6,2,i+3)
        a = sign(q4(:,j,k) - int95(:,k))';
        ii = sort([strfind(a,[-1,1]),strfind(a,[1,-1])]);
        plot(a,'linewidth', 2)
        tt(r+3, 1:length(ii)) = t(ii); 
        sg(r+3,:) = a;
        ylim([-1.5,1.5])
        xlim([0,74])
        str = [str_idx(j), 'sg(Q4 - 95%)'];
        title(str)
        
        i = i + 4;
        r = r + 4;
    end
    
    str = str_stations(k);
    annotation('textbox', [0 0.9 1 0.1], ...
    'String', str, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
end
