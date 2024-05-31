%% zero out positions 
for i= 95:102%101:length(dbase)
    % bug
    xtracks = dbase(i).tracks(:, :, 1); % x bug
    ytracks = dbase(i).tracks(:, :, 2); % y bug
    xmins_b = min(xtracks); % x min bug
    xmin_b = min(xmins_b(xmins_b>=180));
    ymins_b = min(ytracks); % y min bug
    ymin_b = min(ymins_b(ymins_b>=300));

    % mice 
    xmtracks = dbase(i).mtracks(:, :, 1, :); 
    ymtracks = dbase(i).mtracks(:, :, 2, :);
    xmins_m = min(xmtracks);
    xmin_m =  min(xmins_m(xmins_m>=0));
    ymins_m = min(ymtracks);
    ymin_m = min(ymins_m(ymins_m>=0));

    
    xmin = min([xmin_b xmin_m]);
    ymin = min([ymin_b ymin_m]);

    %dbase(i).tracks0 = zeros(size(dbase(i).tracks)); 
    %dbase(i).mtracks0 = zeros(size(dbase(i).mtracks));

    dbase(i).tracks0(:, :, 1) = dbase(i).tracks(:, :, 1) - xmin;
    dbase(i).tracks0(:, :, 2) = dbase(i).tracks(:, :, 2) - ymin;

    dbase(i).mtracks0(:, :, 1, :) = dbase(i).mtracks(:, :, 1, :) - xmin;
    dbase(i).mtracks0(:, :, 2, :) = dbase(i).mtracks(:, :, 2, :) - ymin;
end

%% smooth bug tracks using center of the two nodes

for i=95:102%1:length(dbase)
    nanmeantracks = nanmean(dbase(i).tracks0, 2);
    dbase(i).smoothedtracks = smoothdata(nanmeantracks, "movmedian");
end

%% plot mice position relative to bug 
for i=1:102%107:length(dbase)]
    for m=1:2
        dbase(i).mtracks2b(:, 1, 1, m) = dbase(i).mtracks0(:, 1, 1, m) - dbase(i).smoothedtracks(:, 1, 1); % nose
        dbase(i).mtracks2b(:, 1, 2, m) = dbase(i).mtracks0(:, 1, 2, m) - dbase(i).smoothedtracks(:, 1, 2);
        dbase(i).mtracks2b(:, 2, 1, m) = dbase(i).mtracks0(:, 10, 1, m) - dbase(i).smoothedtracks(:, 1, 1); % tail base
        dbase(i).mtracks2b(:, 2, 2, m) = dbase(i).mtracks0(:, 10, 2, m) - dbase(i).smoothedtracks(:, 1, 2);
    end
end

%% probability histograms of mice position relative to bug

Xedges = -600:40:600;
Yedges = -600:40:600;
for i=1:102%length(dbase)
    for m=1:2
        totalcounts = numel(dbase(i).mtracks2b(:, 1, :, m));
        h = histcounts2(dbase(i).mtracks2b(:, 1, 1, m), dbase(i).mtracks2b(:, 1, 2, m), Xedges, Yedges);
        dbase(i).hprobm2b_n(:, :, m) = h / totalcounts;

        totalcounts = numel(dbase(i).mtracks2b(:, 2, :, m));
        h = histcounts2(dbase(i).mtracks2b(:, 2, 1, m), dbase(i).mtracks2b(:, 2, 2, m), Xedges, Yedges);
        dbase(i).hprobm2b_tb(:, :, m) = h / totalcounts;
    end
end

%%
% attackers
for i=97%9:length(dbase)
    figure('Position', [100, 100, 800, 600]);

    % nose
    hsum = dbase(i).hprobm2b_n(:, :, 1);
    havg = hsum;
    subplot(2, 2, 1);
    imagesc(havg);
    colormap(cm_inferno);
    colorbar;
    yticks(1:2:31);
    xticks(1:2:31);
    yticklabels(Xedges(1:2:end));
    xticklabels(Yedges(1:2:end));
    caxis([0 2.5*10^-3]);
    title("attackers-nose");
    
    hsum = dbase(i).hprobm2b_n(:, :, 2);
    havg = hsum;
    subplot(2, 2, 2);
    imagesc(havg);
    colorbar;
    caxis([0 2.5*10^-3]);
    yticks(1:2:31);
    xticks(1:2:31);
    yticklabels(Xedges(1:2:end));
    xticklabels(Yedges(1:2:end));
    title("non attackers-nose");

    % tail base
    hsum = dbase(i).hprobm2b_tb(:, :, 1); 
    havg = hsum;
    subplot(2, 2, 3);
    imagesc(havg);
    colormap(cm_inferno);
    colorbar;
    yticks(1:2:31);
    xticks(1:2:31);
    yticklabels(Xedges(1:2:end));
    xticklabels(Yedges(1:2:end));
    caxis([0 2.5*10^-3]);
    title("attackers-tail base");
    
    % non-attackers
    hsum = dbase(i).hprobm2b_tb(:, :, 2);
    havg = hsum;
    subplot(2, 2, 4);
    imagesc(havg);
    colorbar;
    caxis([0 2.5*10^-3]);
    yticks(1:2:31);
    xticks(1:2:31);
    yticklabels(Xedges(1:2:end));
    xticklabels(Yedges(1:2:end));
    title("non attackers-tail base");
end 

%%
cm_inferno = inferno(100);
%% average across the four recordings 
Xedges = -600:40:600;
Yedges = -600:40:600;

figure('Position', [100, 100, 800, 600]);

% nose
hsum = dbase(5).hprobm2b_n(:, :, 2) + dbase(6).hprobm2b_n(:, :, 2) + dbase(8).hprobm2b_n(:, :, 1) + dbase(60).hprobm2b_n(:, :, 2);
havg = hsum / 4;
subplot(2, 2, 1);
imagesc(havg);
colormap(cm_inferno);
colorbar;
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end)*0.05);
xticklabels(Yedges(1:2:end)*0.05);
caxis([0 1.5*10^-3]);
title("Attackers - nose");

hsum = dbase(5).hprobm2b_n(:, :, 1) + dbase(6).hprobm2b_n(:, :, 1) + dbase(8).hprobm2b_n(:, :, 2) + dbase(60).hprobm2b_n(:, :, 1);
havg = hsum / 4;
subplot(2, 2, 2);
imagesc(havg);
colorbar;
caxis([0 1.5*10^-3]);
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end)*0.05);
xticklabels(Yedges(1:2:end)*0.05);
title("Non-attackers - nose");

% tail base
hsum = dbase(5).hprobm2b_tb(:, :, 2) + dbase(6).hprobm2b_tb(:, :, 2) + dbase(8).hprobm2b_tb(:, :, 1) + dbase(60).hprobm2b_tb(:, :, 2);
havg = hsum / 4;
subplot(2, 2, 3);
imagesc(havg);
colormap(cm_inferno);
colorbar;
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end)*0.05);
xticklabels(Yedges(1:2:end)*0.05);
caxis([0 1.5*10^-3]);
title("Non-attackers - tail base");

% non-attackers
hsum = dbase(5).hprobm2b_tb(:, :, 1) + dbase(6).hprobm2b_tb(:, :, 1) + dbase(8).hprobm2b_tb(:, :, 2) + dbase(60).hprobm2b_tb(:, :, 1);
havg = hsum / 4;
subplot(2, 2, 4);
imagesc(havg);
colorbar;
caxis([0 1.5*10^-3]);
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end)*0.05);
xticklabels(Yedges(1:2:end)*0.05);
title("non attackers - tail base");

xlabel("Relative X Position from Bug (cm)")
ylabel("Relative Y Position from Bug (cm)")

%% plot trajectory for loop 
timeframes = {[37440:37940], [21300:21800], [52010:52510], [13420:13920]}; 
indices = [5, 6, 8, 60];
attacker = [2, 2, 1, 2];
escaper = [1, 1, 2, 1];

figure('Position', [100, 100, 1200, 900]);
for i=1:4
    % Extract x and y positions from bug vector
    xBug = squeeze(dbase(indices(i)).smoothedtracks(timeframes{i}, 1, 1)); 
    yBug = squeeze(dbase(indices(i)).smoothedtracks(timeframes{i}, 1, 2));

    % fill gaps in mouse track for attacker
    filledtracks = fillgaps(squeeze(dbase(indices(i)).mtracks0(timeframes{i}, 1, :, attacker(i))));
    dbase(indices(i)).smoothedmtracks = smoothdata(filledtracks, "gaussian");
    xMouse1 = dbase(indices(i)).smoothedmtracks(:, 1);
    yMouse1 = dbase(indices(i)).smoothedmtracks(:, 2);

    % fill gaps for escaper
    filledtracks = fillgaps(squeeze(dbase(indices(i)).mtracks0(timeframes{i}, 1, :, escaper(i))));
    dbase(indices(i)).smoothedmtracks = smoothdata(filledtracks, "gaussian");
    xMouse2 = dbase(indices(i)).smoothedmtracks(:, 1);
    yMouse2 = dbase(indices(i)).smoothedmtracks(:, 2);


    time = 1:size(timeframes{i}, 2);

    % Plot the trajectory of positions over time
    subplot(2, 2, i);
    scatter(xBug, yBug, 10, time, "x", 'MarkerFaceAlpha', 0.2);
    colormap("parula");
    hold on; 
    p1 = plot(xBug, yBug, 'LineWidth', 1, 'Color', 'r');
    scatter(xMouse1, yMouse1, 20, time, 'filled', 'MarkerFaceAlpha', 0.2);
    scatter(xMouse2, yMouse2, 20, time, 'filled', 'MarkerFaceAlpha', 0.2);
    
    %p2 = plot(xMouse, yMouse, 'LineWidth', 2, 'Color', 'b');

    % Add rectangle over bug area 
    lastX = xBug(end);
    lastY = yBug(end);

    bugWidth = 20; % pixels (10 mm = 20 pixels) 
    bugLength = 80; 

    if(std(xBug(end-100:end) < std(yBug(end-100:end)))) 
        lastXbug = lastX - bugWidth/2;  
        lastYbug = lastY - bugLength/2;  
        rectangle('Position', [lastXbug, lastYbug, bugWidth, bugLength], ...
        'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 2);
    else
        lastXbug = lastX - bugLength/2;  
        lastYbug = lastY - bugWidth/2;  
        rectangle('Position', [lastXbug, lastYbug, bugLength, bugWidth], ...
        'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 2);
    end
    
    % Add labels and title
    title('Trajectory of Bug and Mouse Positions at Time of Initial Attack');
    xlabel('X Position');
    ylabel('Y Position');
    plots = [p1]%[p1 p2];
    legend(plots, 'bug');%,'mouse');
    xlim([0, 700]);
    ylim([0, 700]);

    cb = colorbar;
    ylabel(cb,'time (s)')
    tickLocations = get(cb, 'Ticks');
    newTickLabels = tickLocations / 80;
    set(cb, 'TickLabels', newTickLabels);
    grid on; 

    old_xticks = get(gca, 'XTick');
    old_yticks = get(gca, 'YTick');
    new_xticks = old_xticks *0.05;
    new_yticks = old_yticks *0.05;
    set(gca, 'XTickLabel', new_xticks);
    set(gca, 'YTickLabel', new_yticks);

end

%saveas(gcf, "attack trajectory 4 examples.png");


%% overlay attacker mouse-bug positions 
timeframes = {[37440:37940], [21300:21800], [52010:52510], [11430:11930]}; 
indices = [5, 6, 8, 60];
attacker = [2, 2, 1, 2];
colors = {'r', 'g', 'b', 'm'};
legends = cell(4, 1);

figure;
for i=1:4
    % fill gaps in mouse track
    filledtracks = fillgaps(squeeze(dbase(indices(i)).mtracks2b(timeframes{i}, 1, :, attacker(i))));
    smoothedm2b = smoothdata(filledtracks, "gaussian");
    xMouse = smoothedm2b(:, 1);
    yMouse = smoothedm2b(:, 2);

    time = 1:size(timeframes{i}, 2);
    
    % Plot the trajectory of positions over time
    scatter(xMouse, yMouse, 50, time, 'filled', 'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.3);
    scatter(xMouse, yMouse, 50, time, 'filled', 'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.3);
    colormap("perula");
    hold on; 
    % p2 = plot(xMouse, yMouse, 'LineWidth', 2, 'Color', colors{1}, ...
    %     'DisplayName', sprintf('Mouse %d', i));

    % Add labels and title
    title('Trajectory of Mouse Positions Relative to Bug Over Time for Attack');
    xlabel('X Position');
    ylabel('Y Position');
    %legends{i} = sprintf('Mouse %d', i);
    %legend(p2, strcat('mouse', i));
    cb = colorbar;
    ylabel(cb,'time (s)')
    tickLocations = get(cb, 'Ticks');
    newTickLabels = tickLocations / 80;
    set(cb, 'TickLabels', newTickLabels);
    grid on;

end
legend('show');
saveas(gcf, "bug-mouse-relative-pos-attack.png")

%% plot trajectory
timeframe5 = 37440:37840;
timeframe6 = 21300:21700;
timeframe8 = 51910:52310;
timeframe60 = 11430:11830;

% Extract x and y positions from vectors A and B
xBug = squeeze(dbase(60).smoothedtracks(timeframe60, 1, 1)); 
yBug = squeeze(dbase(60).smoothedtracks(timeframe60, 1, 2)); 

% fill gaps in mouse track
filledtracks = fillgaps(squeeze(dbase(60).mtracks0(timeframe60, 1, :, 2)));
dbase(60).smoothedmtracks = smoothdata(filledtracks, "movmean");
xMouse = dbase(60).smoothedmtracks(:, 1);
yMouse = dbase(60).smoothedmtracks(:, 2);
% xMouse = squeeze(dbase(5).mtracks0(timeframe5, 1, 1, 2)); % Extract x positions 
% yMouse = squeeze(dbase(5).mtracks0(timeframe5, 1, 2, 2)); % Extract y positions 

time = 1:size(timeframe60, 2);

% Plot the trajectory of positions over time
figure;
scatter(xBug, yBug, 50, time, 'filled');
colormap("parula");
hold on; 
plot(xBug, yBug, 'r');
scatter(xMouse, yMouse, 50, time, 'filled');
plot(xMouse, yMouse, 'b');

% Add labels and title
title('Trajectory of Object Positions Over Time');
xlabel('X Position');
ylabel('Y Position');
legend('Bug', 'Mouse');
colorbar;
grid on; 


%% attack rates over weeks for attacker
attacks = zeros(7, 1);
attacks(1, 1) = dbase(40).hprobm2b_n(15, 15, 2);
attacks(2, 1) = dbase(31).hprobm2b_n(15, 15, 1);
attacks(3, 1) = dbase(24).hprobm2b_n(15, 15, 2);
attacks(4, 1) = dbase(14).hprobm2b_n(15, 15, 2);
attacks(5, 1) = dbase(6).hprobm2b_n(15, 15, 2);
attacks(6, 1) = dbase(45).hprobm2b_n(15, 15, 2);
attacks(7, 1) = dbase(56).hprobm2b_n(15, 15, 1);

figure;
plot(attacks)
xticklabels({'week 3', 'week 4', 'week 5', 'week 6', 'week 8', 'week 12','week 20'})

%% attack average across ELE week 8 males 
attacks = []
for i=21:24
    for m=1:2
        attacks = vertcat(attacks, dbase(i).hprobm2b_n(15, 15, m));
    end
end

disp(mean(attacks));

attacks = []
for i=84:86
    for m=1:2
        attacks = vertcat(attacks, dbase(i).hprobm2b_n(15, 15, m));
    end
end

disp(mean(attacks));

%% plot attack average over weeks
ele = [2.4740e-04, 5.2301e-04, 0.0012,  0.0022, 7.3243e-04];
cont = [5.2808e-04, 4.7093e-04, 0.0011, 4.7093e-04, 2.6259e-04];

figure;
plot(ele);
hold on;
plot(cont);
xticks([1 2 3 4 5]);

%%
array = [];
for i=45:48;
    for m=1:2
        array = vertcat(array, dbase(i).hprobm2b_n(15, 15, m));
    end
end




