%% zero out positions for control
for i= 107:length(dbase)
    % mice 
    xmtracks = dbase(i).mtracks(:, :, 1, :); 
    ymtracks = dbase(i).mtracks(:, :, 2, :);
    xmins_m = min(xmtracks);
    xmin_m =  min(xmins_m(xmins_m>=180));
    ymins_m = min(ymtracks);
    ymin_m = min(ymins_m(ymins_m>=300));

    %dbase(i).tracks0 = zeros(size(dbase(i).tracks)); 
    %dbase(i).mtracks0 = zeros(size(dbase(i).mtracks));

    dbase(i).mtracks0(:, :, 1, :) = dbase(i).mtracks(:, :, 1, :) - xmin_m;
    dbase(i).mtracks0(:, :, 2, :) = dbase(i).mtracks(:, :, 2, :) - ymin_m;
end

%% plot trajectories of mice and bug
timeframes = 38400:43200; % 60 s * 80 fps = 4800
xBug = squeeze(dbase(57).smoothedtracks(timeframes, 1, 1));
yBug = squeeze(dbase(57).smoothedtracks(timeframes, 1, 2));

% fill gaps in mouse track
for m = 1:2
    filledtracks = fillgaps(squeeze(dbase(57).mtracks0(timeframes, 10, :, m)));
    dbase(57).smoothedmtracks2(:, :, m) = smoothdata(filledtracks, "movmedian");
end

% xMouse = dbase(indices(i)).smoothedmtracks(:, 1);
% yMouse = dbase(indices(i)).smoothedmtracks(:, 2);

time = 1:size(timeframes, 2);

xMouse1 = dbase(57).smoothedmtracks2(time, 1, 1); % tail base
yMouse1 = dbase(57).smoothedmtracks2(time, 2, 1);

xMouse2 = dbase(57).smoothedmtracks2(time, 1, 2); % tail base
yMouse2 = dbase(57).smoothedmtracks2(time, 2, 2);

% Plot the trajectory of positions over time
figure;
scatter3(xMouse1, yMouse1, time, 10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
%colormap("parula");
hold on;
%p1 = plot(xBug, yBug, 'LineWidth', 2, 'Color', 'r');
scatter3(xMouse2, yMouse2, time, 10, 'filled', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
scatter3(xBug, yBug, time, 10, "filled", 'MarkerFaceColor', [0.9290 0.6940 0.1250]);
%p2 = plot(xMouse, yMouse, 'LineWidth', 2, 'Color', 'b');

zlabel('Time (s)');
xlabel('X Position');
ylabel('Y Position');

newTickLocations = linspace(0, 4800, 7); % 8 ticks from 0 to 4800
newTickLabels = cellstr(num2str(newTickLocations' / 80)); % Divide each tick label by 80
zticks(newTickLocations);
zticklabels(newTickLabels);
daspect([1 1 3]);
xlim([250 850]);
ylim([0 600]);

legend({'mouse 1', 'mouse 2', 'bug'});
%% plot trajectories of two mice no bug
timeframes = 3840:8640; % 60 s * 80 fps = 4800

% fill gaps in mouse track
for m = 1:2
    filledtracks = fillgaps(squeeze(dbase(107).mtracks0(timeframes, 10, :, m)));
    dbase(107).smoothedmtracks2(:, :, m) = smoothdata(filledtracks, "movmedian");
end

% xMouse = dbase(indices(i)).smoothedmtracks(:, 1);
% yMouse = dbase(indices(i)).smoothedmtracks(:, 2);

time = 1:size(timeframes, 2);

xMouse1 = dbase(107).smoothedmtracks2(time, 1, 1); % tail base
yMouse1 = dbase(107).smoothedmtracks2(time, 2, 1);

xMouse2 = dbase(107).smoothedmtracks2(time, 1, 2); % tail base
yMouse2 = dbase(107).smoothedmtracks2(time, 2, 2);

% Plot the trajectory of positions over time
figure;
%colormap("parula");
scatter3(xMouse1, yMouse1, time, 10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
hold on;
%p1 = plot(xBug, yBug, 'LineWidth', 2, 'Color', 'r');
scatter3(xMouse2, yMouse2, time, 10, 'filled', 'MarkerFaceColor', [0.6350 0.0780 0.1840]);
%p2 = plot(xMouse, yMouse, 'LineWidth', 2, 'Color', 'b');

zlabel('Time (s)');
xlabel('X Position');
ylabel('Y Position');

newTickLocations = linspace(0, 4800, 7); % 8 ticks from 0 to 4800
newTickLabels = cellstr(num2str(newTickLocations' / 80)); % Divide each tick label by 80
zticks(newTickLocations);
zticklabels(newTickLabels);
daspect([1 1 3]);
xlim([0 600]);
ylim([0 600]);


legend({'mouse 1', 'mouse 2'});

%% ethograms

figure;
time = 5000:9800;
barh(1:size(time,1), dbase(107).umapClustWT(time, 3, 1), 'stacked');

% for i=numel(customColors)
%     set(get(gca, 'Children'), 'FaceColor', customColors{i});
% end
