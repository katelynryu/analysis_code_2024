
% Suppose you have a table named 'data' with fields (subject, state, fraction, group)
% Example:
data = readtable('umap-combine-wcontrol.csv'); % Read your actual data from a file

% Calculate total fraction spent in each state for each group
data.weighted_fraction = data.Count .* data.Fraction;
groupStateFractions = grpstats(data, {'Group', 'State'},  {'sum', 'mean'}, 'DataVars', {'Count', 'weighted_fraction'});
%groupStateFractions = grpstats(data, {'Group', 'State'}, 'sum', 'DataVars', 'Fraction');

% Pivot the table to have groups as rows and states as columns
pivotTable = unstack(groupStateFractions, {'sum_weighted_fraction', 'mean_Count'}, 'State');
%pivotTable = unstack(groupStateFractions, 'sum_fraction', 'state');

% Extract group names and states
groupNames = pivotTable.Group;
states = pivotTable.Properties.VariableNames(2:end);

% Create stacked bar chart
figure;
bar(groupNames, pivotTable{:, 2:end}, 'stacked');
xlabel('Group');
ylabel('Fraction');
title('Fraction Spent in Each State per Group');
legend(states, 'Location', 'Best');

%% 
% Example data (replace this with your actual data)
% Suppose you have a table named 'data' with fields (subject, state, count, fraction, group, id)
% Example:
data = readtable('umap-combine-wcontrol.csv'); % Read your actual data from a file

% Initialize variables to store aggregated fractions by group and state
uniqueGroups = unique(data.Group);
uniqueStates = unique(data.State);
groupStateFractions = zeros(length(uniqueGroups), length(uniqueStates));
groupStateCounts = zeros(length(uniqueGroups), length(uniqueStates));

% Loop through each row of the table
for i = 1:height(data)
    % Find the index of the group and state for the current row
    groupIndex = find(strcmp(data.Group{i}, uniqueGroups));
    stateIndex = find(data.State(i) == uniqueStates);
    
   groupStateFractions(groupIndex, stateIndex) = groupStateFractions(groupIndex, stateIndex) + data.Fraction(i);
   groupStateCounts(groupIndex, stateIndex) = groupStateCounts(groupIndex, stateIndex) + 1;
end

groupStateAverages = groupStateFractions ./ groupStateCounts;

%% 
groupStateAverages = groupStateAverages([3, 1, 2], :);
%% normalize 
totalFractionPerState = sum(groupStateFractions, 1);

% Normalize the fractions for each group
groupStateFractionsNormalized = groupStateFractions ./ totalFractionPerState;
%%
customColormap = parula(size(groupStateFractions, 1));
figure;
h = bar(groupStateAverages, 'stacked');
%customColors = {[0.2 0.4 0.6], [0.8 0.2 0.1], [0.5 0.7 0.3]};

numColors = 13;

% Create a muted rainbow colormap
rainbowColormap = jet(numColors);
saturationFactor = 0.9; % Adjust saturation factor to mute the colors
rainbowColormap(:, 2) = rainbowColormap(:, 2) * saturationFactor;

% Convert the colormap to cell array of RGB triplets
customColors = mat2cell(rainbowColormap, ones(numColors, 1));

for i = 1:numel(h)
    h(i).FaceColor = customColors{i};
end

xlabel('State');
ylabel('Fraction');
title('Fraction Spent in Each State per Group');
stateLabels = {'unassigned', 'trot', 'slow explore', 'trot', 'rear', 'head grooming', ...
    'walk', 'grooming', 'walk', 'idle', 'stepping', 'fast explore', 'turn'};
legend(stateLabels, 'Location', 'Best');

groupLabels = {'Control (No Bug)', 'Exposure to Bug 1 time/week', 'Exposure to Bug 4 times/week'};
set(gca, 'XTickLabel', groupLabels); % Set x-axis tick labels to state names


%% dwell time of idle 
for i = [1:102, 107:length(dbase)]
    for m=1:2
        idle = dbase(i).umapClustWT(:, 3, m) == 9; 
        idlex = double(idle);

        % Find the indices where the mask changes (from 0 to 1 or from 1 to 0)
        changeIndices = find(diff([0; idlex; 0]));

        % Calculate the lengths (duration) of consecutive segments where the mask is equal to 1
        if(m==1)
            dbase(i).idleLength_m1 = diff(changeIndices);
        else
            dbase(i).idleLength_m2 = diff(changeIndices);
        end
    end
end

%% connected components ver. 
for i = [1:102, 107:length(dbase)]
    for m=1:2
        idle = dbase(i).umapClustWT(:, 3, m) == 9; 
        idlex = double(idle);

        % Find the indices where the mask changes (from 0 to 1 or from 1 to 0)
        CC = bwconncomp(idlex);
        component_lengths = cellfun(@length, CC.PixelIdxList);

        % Calculate the lengths (duration) of consecutive segments where the mask is equal to 1
        if(m==1)
            dbase(i).idleLength_m1CC = component_lengths;
        else
            dbase(i).idleLength_m2CC = component_lengths;
        end
    end

end

%% spatial correlate of idle
Xedges = 20:20:820;
Yedges = 20:20:820;
idle_points = [];
for i =107:length(dbase)%[1:8, 30:35]%[1:102, 107:length(dbase)]
    for m=1:2
        idle = dbase(i).umapClustWT(:, 3, m) == 9; 

        subplot(1, 2, 1);
        idle_points = vertcat(idle_points, dbase(i).mtracks0(idle, 10, :, m)); % Extract x and y positions
        % Create a heatmap for the current cluster
        % 
        % subplot(1, 2, 2);
        % idle_pointsb = dbase(i).smoothedtracks(idle, 1, :);
        % h = histcounts2(idle_points(:, 1, 1), idle_pointsb(:, 1, 2), Xedges, Yedges);
        % imagesc(h);
    end
end

figure;
h = histcounts2(idle_points(:, 1, 1), idle_points(:, 1, 2), Xedges, Yedges, "Normalization", "pdf");
imagesc(h);
colorbar;
xlim([0 31]);
ylim([0 31]);
caxis([0 2*10^-5]);
xlabel("X Position (cm)");
ylabel("Y Position (cm)");


%% idle distance to bug 
Xedges = 0:0.3:30;
Yedges = 0:0.3:30;
idle_points = [];
for i = 1:56%57:102%57:102%102%[1:8, 30:35]%[1:102, 107:length(dbase)]
    for m=1:2
        idle = dbase(i).umapClustWT(:, 3, m) == 9; 
        idle_points = vertcat(idle_points, dbase(i).dist(idle, m) * pixel_size * 0.1); % Extract x and y positions
    end
end

%%
figure;

[counts1, ~] = histcounts(idle_points2, Xedges, "Normalization", "pdf");
plot(Xedges(2:end), counts1, 'LineWidth', 2, 'Color', '#7E2F8E');
hold on;

[counts2, ~] = histcounts(idle_points, Xedges, "Normalization", "pdf");
plot(Xedges(2:end), counts2, 'LineWidth', 2, 'Color', '#77AC30');
hold on;

ylabel("probability");
xlabel("distance from bug (cm)");
legend({'Exposure', 'Additional Exposure'});
%xlim([0 400]);

%% dwell time by group (no bug vs. bug)

idle_bugCC = [];
idle_controlCC = [];

for i=1:102
    idle_bugCC = vertcat(idle_bugCC, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
end

for i=107:length(dbase)
    idle_controlCC = vertcat(idle_controlCC, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
end

%% histogram of idle
figure;
bin_edges = 3:6:600;
h1 = histogram(idle_bugCC(idle_bugCC > 3 & idle_bugCC< 600), bin_edges, "Normalization", "pdf");
hold on;
h2 = histogram(idle_controlCC(idle_controlCC>3 &  idle_controlCC< 600), bin_edges, "Normalization", "pdf");

xlabel('Durations (s)');
ylabel('Frequency');
title('Histogram of Idle Durations');
legend({'with bug', 'control'});
xticklabels([0, 1.25, 2.5, 3.75, 5.5, 6.75, 7.5]);

%%
figure;
bin_edges = 3:6:600;
h1 = histcounts(idle_bugCC(idle_bugCC > 3 & idle_bugCC< 600), bin_edges, "Normalization", "pdf");
h2 = histcounts(idle_controlCC(idle_controlCC>3 &  idle_controlCC< 600), bin_edges, "Normalization", "pdf");

plot(bin_edges(1:end-1), h1, 'LineWidth', 2);
hold on;
plot(bin_edges(1:end-1), h2, 'LineWidth', 2);
hold on;
xlabel('Durations');
ylabel('Frequency');
title('Histogram of Idle Durations (s)');
legend({'with bug', 'control'});
xticklabels([0, 1.25, 2.5, 3.75, 5.5, 6.75, 7.5]);

%% box plot

% Remove outliers from durations1
idle_bug_no_outliers = idle_bugCC(idle_bugCC < quantile(idle_bugCC, 0.75) ...
    + 1.5*iqr(idle_bugCC) & idle_bugCC > quantile(idle_bugCC, 0.25) - 1.5*iqr(idle_bugCC));

% Remove outliers from durations2
idle_control_no_outliers = idle_controlCC(idle_controlCC < quantile(idle_controlCC, 0.75) ...
    + 1.5*iqr(idle_controlCC) & idle_controlCC > quantile(idle_controlCC, 0.25) - 1.5*iqr(idle_controlCC));

idlebno = idle_bugCC(idle_bugCC > 3 & idle_bugCC< 600);
idlecno = idle_controlCC(idle_controlCC>3 &  idle_controlCC< 600);

maxLength = max(length(idlebno), length(idlecno));
idlebno(end+1:maxLength) = NaN;
idlecno(end+1:maxLength) = NaN;
combinedDurations = [idlebno, idlecno];

figure;

% Create side-by-side box plots for each duration vector
boxplot(combinedDurations, 'Labels', {'Bug', 'No Bug'}, 'OutlierSize', 1);

% Add labels and title
xlabel('Groups');
ylabel('Durations (s)');
title('Side-by-Side Box Plots of Idle Durations');

%% ks test 
% Perform the KS test
[h, p, ksstat] = kstest2(idlebno, idlecno, 'Alpha', 0.05);  % Perform KS test with significance level of 0.05

% Display the results
if h == 0
    disp('The null hypothesis cannot be rejected. The two data sets are consistent with each other.');
else
    disp('The null hypothesis is rejected. The two data sets are not consistent with each other.');
end
disp(['p-value: ', num2str(p)]);
disp(['KS statistic: ', num2str(ksstat)]);

%%
% Perform Wilcoxon rank-sum test
[p, h] = ranksum(idle_bugCC, idle_controlCC);

% Display p-value
fprintf('The p-value is: %.9f\n', p);

%% 
for i=[1:102, 107:length(dbase)]
    dbase(i).idleLength_m1CCa = dbase(i).idleLength_m1CC(dbase(i).idleLength_m1CC > 3 & dbase(i).idleLength_m1CC < 600);
    dbase(i).idleLength_m2CCa = dbase(i).idleLength_m2CC(dbase(i).idleLength_m2CC > 3 & dbase(i).idleLength_m2CC < 600);
end

%% evolution of dwell time of idle across weeks 
idle_week3E = [];
idle_week4E = [];
idle_week5E = [];
idle_week6E = [];
idle_week8E = [];
idle_week12E = [];
idle_week16E = [];

for i=57:102
    if(strcmp(dbase(i).week, "3"))
        idle_week3E = vertcat(idle_week3E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    elseif(strcmp(dbase(i).week, "4"))
        idle_week4E = vertcat(idle_week4E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    elseif(strcmp(dbase(i).week, "5"))
        idle_week5E = vertcat(idle_week5E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    elseif(strcmp(dbase(i).week, "6"))
        idle_week6E = vertcat(idle_week6E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    elseif(strcmp(dbase(i).week, "8"))
        idle_week8E = vertcat(idle_week8E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    elseif(strcmp(dbase(i).week, "12"))
        idle_week12E = vertcat(idle_week12E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    else
        idle_week16E = vertcat(idle_week16E, dbase(i).idleLength_m1CC', dbase(i).idleLength_m2CC');
    end
end 

idle_progressionC = [median(idle_week3E), median(idle_week4E), median(idle_week5E), median(idle_week6E), ...
    median(idle_week8E), median(idle_week12E), median(idle_week16E)];
stderrsC = [std(idle_week3E)/sqrt(numel(idle_week3E)), std(idle_week4E)/sqrt(numel(idle_week4E)), ...
    std(idle_week5E)/sqrt(numel(idle_week5E)), std(idle_week6E)/sqrt(numel(idle_week6E)),...
    std(idle_week8E)/sqrt(numel(idle_week8E)), std(idle_week12E)/sqrt(numel(idle_week12E)),...
    std(idle_week16E)/sqrt(numel(idle_week16E))]; 

%%
stderror = stderrsC;
x = 1:numel(idle_progressionC);
y = idle_progressionC;

figure;
errorbar(x, y, stderror, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);

hold on; 

stderror = stderrsC;
y = idle_progressionC; 
errorbar(x, y, stderror, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);


weeklabels = {'3', '4', '5', '6', '8', '12', '16'};
xticklabels(weeklabels);
xlabel("Week");
ylabel("Median Idle Duration (Frames)");
legend({'Exposure', 'No Exposure'});

%% idle based on bug state




%% idle duration (idle/min) over video
for i=[1:102, 107:length(dbase)]
    dbase(i).idlepm(:, 1) = sum(dbase(i).idleLength_m1CCa) / size(dbase(i).time, 1);
    dbase(i).idlepm(:, 2) = sum(dbase(i).idleLength_m2CCa) / size(dbase(i).time, 1);
end

%% evolution of dwell time of idle across weeks 
idle_week3E = [];
idle_week4E = [];
idle_week5E = [];
idle_week6E = [];
idle_week8E = [];
idle_week12E = [];
idle_week16E = [];

for i=1:56
    if(strcmp(dbase(i).week, "3"))
        idle_week3E = vertcat(idle_week3E, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "4"))
        idle_week4E = vertcat(idle_week4E, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "5"))
        idle_week5E = vertcat(idle_week5E, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "6"))
        idle_week6E = vertcat(idle_week6E, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "8"))
        idle_week8E = vertcat(idle_week8E, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "12"))
        idle_week12E = vertcat(idle_week12E, dbase(i).idlepm(:));
    else
        idle_week16E = vertcat(idle_week16E, dbase(i).idlepm(:));
    end
end 

%%
idle_progressionE = [median(idle_week3E), median(idle_week4E), median(idle_week5E), median(idle_week6E), ...
    median(idle_week8E), median(idle_week12E), median(idle_week16E)];
stddevE = [std(idle_week3E), std(idle_week4E), std(idle_week5E), std(idle_week6E), ...
    std(idle_week8E), std(idle_week12E), std(idle_week16E)]; 

%% evolution of dwell time of idle across weeks 
idle_week3C = [];
idle_week4C = [];
idle_week5C = [];
idle_week6C = [];
idle_week8C = [];
idle_week12C = [];
idle_week16C = [];

for i=57:102
    if(strcmp(dbase(i).week, "3"))
        idle_week3C = vertcat(idle_week3C, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "4"))
        idle_week4C = vertcat(idle_week4C, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "5"))
        idle_week5C = vertcat(idle_week5C, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "6"))
        idle_week6C = vertcat(idle_week6C, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "8"))
        idle_week8C = vertcat(idle_week8C, dbase(i).idlepm(:));
    elseif(strcmp(dbase(i).week, "12"))
        idle_week12C = vertcat(idle_week12C, dbase(i).idlepm(:));
    else
        idle_week16C = vertcat(idle_week16C, dbase(i).idlepm(:));
    end
end 

%% 
idle_progressionC = [median(idle_week3C), median(idle_week4C), median(idle_week5C), median(idle_week6C), ...
    median(idle_week8C), median(idle_week12C), median(idle_week16C)];
stddevC = [std(idle_week3C), std(idle_week4C), std(idle_week5C), std(idle_week6C), ...
    std(idle_week8C), std(idle_week12C), std(idle_week16C)];

%% no exposure vs. exposure 
stderror = stddevC;
x = 1:numel(idle_progressionC);
y = idle_progressionC;

figure;
errorbar(x, y, stderror, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);

hold on; 

stderror = stddevE;
y = idle_progressionE; 
errorbar(x, y, stderror, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);


weeklabels = {'3', '4', '5', '6', '8', '12', '16'};
xticklabels(weeklabels);
xlabel("Week");
ylabel("Median Idle Duration (Fraction/min)");
legend({'No Exposure', 'Exposure'});


%% attacker vs non-attacker idle duration


%% ethogram examples 
% 2 mice w/o bug



% 2 mice w bug 

