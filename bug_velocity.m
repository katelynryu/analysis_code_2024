%%
for i=1:length(dbase) 
    dbase(i).tracksb = getVelocity(dbase(i).SLEAPtracksb(:,1,:), 20, 80, 1/1.97);
    dbase(i).tracksm1 = getVelocity(dbase(i).SLEAPtracksm(:,10,:,1), 20, 80, 1/1.97);
    dbase(i).tracksm2 = getVelocity(dbase(i).SLEAPtracksm(:,10,:,2), 20, 80, 1/1.97);
end 

%%
for i=1:length(dbase)
    histogram(dbase(i).tracksb);
end

%%
v1 = getVelocity(squeeze(dbase(1).tracks(:, 1, :)), 20, 80, 1/1.97);

%%
for i=1:length(dbase) 
   dbase(i).velocity = getVelocity(squeeze(dbase(i).tracks(1:(60*fps),1,:)), 20, 80, 1/1.97);
end

%%
for i=1:length(dbase)
    histogram(dbase(i).velocity, 0:0.1:2, "Normalization", "pdf");
    hold on;
end 

%%
for i=1%:length(dbase)
    plot(timex(1:60*fps), dbase(i).velocity);
    hold on;
    yyaxis right
    plot(timex(1:60*fps), dbase(i).tracks(1:(60*fps), 1, 1));
end

%%
for i=1%:length(dbase)
    x = dbase(i).tracks(1:(60*fps), 1, 1);
    y = dbase(i).tracks(1:(60*fps), 1, 2);
    [counts, edges] = histcounts2(x, y, 50);

    figure
    imagesc(edges(1), edges(2), counts)
    colormap(jet) 
    hold on;
end 

%%
for i=1%:length(dbase)
    x =  x + dbase(i).tracks(1:(60*fps), 1, 1);
    y = y + dbase(i).tracks(1:(60*fps), 1, 2);
end 
binscatter(x, y);

%%
for i=1:length(dbase)
    x = dbase(i).tracks(1:(60*fps), 1, 1);
    y = dbase(i).tracks(1:(60*fps), 1, 2);
    binscatter(x, y);
    hold on;
end 
%%
for i=1:length(dbase)
    x =  x + dbase(i).tracks(:, 1, 1);
    y = y + dbase(i).tracks(:, 1, 2);
    [counts, edges] = histcounts2(x, y, 50);
end 

figure
imagesc(edges(1), edges(2), counts)
colormap(jet)  % choose a colormap
colorbar

%% angular change
ang_vel = getAngularVelocity(dbase(2).tracks(:, :, 1));
ang_v = ang_vel * fps;
ang_sum = cumsum(ang_vel);
(sum(abs(ang_vel) > 0.1))

%% reconstruct bug array
for i=1
    for j=1:length(dbase(i).velocity)
        if (dbase(i).velocity(j,:) > 0.2) % | angular velocity | distance between points? is this already covered)
            dbase(i).SLEAPtracksb(j, :, :) = NaN;
        end
    end
end

%% load movie
vid_name = 'Y:\behdata\2023-06-SocDev\C57Male\week8\OFTsocialgroup-0125-90.h5';
mov = h5read(vid_name, '/pg0', [1,1,1,9445], [Inf, Inf, Inf, 1000]);
for i = 9446:10445
    imshow(mov(:,:,:,i-9445))
    hold on
    plot([dbase(2).tracks(i,1,1), dbase(2).tracks(i,2,1)], [dbase(2).tracks(i,1,2), dbase(2).tracks(i,2,2)]);
    drawnow()
end

% for loop for mouse 1
for i = 9446:10445
    imshow(mov(:,:,:,i-9445))
    hold on
    plot([dbase(2).mtracks(i,5,1,1), dbase(2).mtracks(i,6,1,1)], [dbase(2).mtracks(i,5,2,1), dbase(2).mtracks(i,6,2,1)]);
    drawnow()
end

% mouse 2
for i = 9446:10445
    imshow(mov(:,:,:,i-9445))
    hold on
    plot([dbase(2).mtracks(i,1,1,2), dbase(2).mtracks(i,10,1,2)], [dbase(2).mtracks(i,1,2,2), dbase(2).mtracks(i,10,2,2)]);
    drawnow()
end
%%
for t=1:length(dbase)
    dbase(t).dist = zeros(size(timex, 1), 3);
    for i = 1:2
        for j = 1:size(timex, 1)
        dbase(t).dist(j,i) = returnDist(dbase(t).mtracks(j, 10, :, i),  dbase(1).tracks(j, 1, :)) * pixel_size;
        dbase(t).dist(j,3) = returnDist(dbase(t).mtracks(j, 10, :, 1),  dbase(1).mtracks(j, 10, :, 2)) * pixel_size;
        end
    end
end
%%
bmdistances = calculate_distances_over_time(dbase(1).mtracks(:, 1, :, 1), dbase(1).tracks(:, 1, :));
bmdist = bmdistances / pixel_size;
histogram(bmdist, 100);
hold on

bmdistances2 = calculate_distances_over_time(dbase(5).mtracks(:, 1, :, 1), dbase(5).tracks(:, 1, :));
bmdist2 = bmdistances2 / pixel_size;
%plot(timex, bmdist2);
histogram(bmdist2, 100);

%%
%distnELE = zeros(size(dbase(1).time, 1), 1);
for i=1:8
    bmdistances = calculate_distances_over_time(dbase(i).mtracks(:, 1, :, 1), dbase(i).tracks(:, 1, :));
    dbase(i).bmdistm1 = bmdistances / pixel_size;
end

%%
distnELE2 = horzcat(dbase(1).bmdist2, dbase(2).bmdist2, dbase(3).bmdist2, dbase(4).bmdist2);
distELEm2 = horzcat(dbase(5).bmdistm2, dbase(6).bmdistm2, dbase(7).bmdistm2, dbase(8).bmdistm2);
histogram(distnELE2, 100);
hold on;

histogram(distELE2, 100);

%% attack vs escape
distesc = horzcat(dbase(1).bmdistm2, dbase(2).bmdistm2, dbase(3).bmdistm2, dbase(7).bmdistm2);
distattack = horzcat(dbase(5).bmdistm2, dbase(6).bmdistm2, dbase(4).bmdistm2, dbase(8).bmdistm2);

histogram(distesc, 100);
hold on;
histogram(distattack, 100);

% Compute the histogram
[countsesc, edgesesc] = histcounts(distesc, 100);
[countsatt, edgesatt] = histcounts(distattack, 100);

% Calculate the centers of the bins
bin_centersesc = (edgesesc(1:end-1) + edgesesc(2:end)) / 2;
bin_centersatt = (edgesatt(1:end-1) + edgesatt(2:end)) / 2;

% Plot the line chart
plot(bin_centersesc, countsesc, 'b-', 'LineWidth', 2);
hold on;
plot(bin_centersatt, countsatt, 'r-', 'LineWidth', 2);
hold on;

%%
% Compute the histogram
[counts, edges] = histcounts(distnELE, 100);
[counts2, edges2] = histcounts(distELE, 100);
[countsm2, edgesm2] = histcounts(distELEm2, 100);

% Calculate the centers of the bins
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
bin_centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
bin_centersm2 = (edgesm2(1:end-1) + edgesm2(2:end)) / 2;

% Plot the line chart
plot(bin_centers, counts, 'b-', 'LineWidth', 2);
hold on;
plot(bin_centers2, counts2, 'r-', 'LineWidth', 2);
hold on;
plot(bin_centersm2, countsm2, 'y-', 'LineWidth', 2);

%%
for i=5:8
    bmdistances2 = calculate_distances_over_time(dbase(i).mtracks(:, 1, :, 1), dbase(i).tracks(:, 1, :));
    bmdist2 = bmdistances2 / pixel_size;
    distELE(:, i) = bmdist2;
end

%histogram(distELE, 100);

%%
threshold_distance = 500; % Adjust as needed
is_bout2 = dbase(1).bmdist < threshold_distance;

% Detect bout start and end times
start_times = find(diff([0; is_bout2]) == 1);
end_times = find(diff([is_bout2; 0]) == -1);

% Calculate bout durations
bout_durations2 = end_times - start_times;

%% attack bouts
threshold_distance = 500; % Adjust as needed
is_bout = dbase(5).bmdistm2 < threshold_distance;

% Detect bout start and end times
start_times = find(diff([0; is_bout]) == 1);
end_times = find(diff([is_bout; 0]) == -1);

% Calculate bout durations
bout_durations = end_times - start_times;

disp(bout_durations)

%% attack bouts
threshold_distance = 500; % Adjust as needed
is_bout = dbase(5).bmdistm2 < threshold_distance;

% Detect bout start and end times
start_times = find(diff([0; is_bout]) == 1);
end_times = find(diff([is_bout; 0]) == -1);

% Calculate bout durations
bout_durations = end_times - start_times;

% plot
figure;
plot(1:size(dbase(5).bmdistm2,1), dbase(5).bmdistm2);
hold on;
plot(start_times, dbase(5).bmdistm2(start_times), 'go', 'MarkerSize', 10); % Mark bout start times
plot(end_times, dbase(5).bmdistm2(end_times), 'ro', 'MarkerSize', 10); % Mark bout end times
xlabel('Time');
ylabel('Distance');
title('Distance between points over time with bouts');
legend('Distance', 'Bout Start', 'Bout End');

%% bug velocity vs. distance with mouse
figure;
plot(squeeze(dbase(1).mtracks(:, 10, :, 1)));