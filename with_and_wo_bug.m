%% mouse-mouse distance in bug group
for i=1:102%length(dbase)
    dbase(i).dist = zeros(size(dbase(i).time, 1), 3);
    for m = 1:2
        for j = 1:size(dbase(i).time, 1)
            dbase(i).dist(j,m) = returnDist(dbase(i).mtracks0(j, 10, :, m),  dbase(i).smoothedtracks(j, 1, :)) * pixel_size;
            dbase(i).dist(j,3) = returnDist(dbase(i).mtracks0(j, 10, :, 1),  dbase(i).mtracks0(j, 10, :, 2)) * pixel_size;
        end
    end
end

%% mouse-mouse distance in control group (no bug)
for i=107:length(dbase)
    dbase(i).dist = zeros(size(dbase(i).time, 1), 1);
        for j = 1:size(dbase(i).time, 1)
            dbase(i).dist(j,1) = returnDist(dbase(i).mtracks0(j, 10, :, 1),  dbase(i).mtracks0(j, 10, :, 2)) * pixel_size;
        end
end

%% plot distributions of mouse-mouse distances

mdistbug = [];
mdistcont = [];

for i=1:95
    mdistbug = vertcat(mdistbug, dbase(i).dist(:, 3));
end

for i=107:length(dbase)
    mdistcont = vertcat(mdistcont, dbase(i).dist(:, 1));
end

%% histcounts of mouse mouse 

figure;
bin_edges = 0:4.5:450;
[counts, ~] = histcounts(mdistcont, bin_edges, "Normalization", "probability");
plot(bin_edges(2:end), counts, 'LineWidth', 2);

hold on;
%h1 = histcounts(mdistbug, bin_edges, "Normalization", "probability");
[counts, ~] = histcounts(mdistbug, bin_edges, "Normalization", "probability");
plot(bin_edges(2:end), counts, 'LineWidth', 2);



xlabel('Distances (mm)');
ylabel('Probability');
title('Distribution of Mouse-Mouse Distances');
legend({'Conrol', 'With Bug'});


%% plot histograms
figure;
bin_edges = 0:4.5:450;
h2 = histogram(mdistcont, bin_edges, "Normalization", "probability");
hold on;
h1 = histogram(mdistbug, bin_edges, "Normalization", "probability");


xlabel('Distances (cm)');
ylabel('Probability');
title('Distribution of Mouse-Mouse Distances');
legend({'Control M-M', 'Exposure M-M'});
xlim([0 450]);

old_xticks = get(gca, 'XTick');

% Calculate new tick values (divide by 50)
new_xticks = old_xticks / 10;

% Set new tick labels
set(gca, 'XtickLabel', new_xticks);

%% bug-mouse distance
bm1 = [];
bm2 = [];
m1m2 = [];

for i=1:102
    bm1 = vertcat(bm1, dbase(i).dist(:, 1), dbase(i).dist(:, 2));
    m1m2 = vertcat(m1m2, dbase(i).dist(:, 3));
end

figure;
bin_edges = 0:4:450;
h3 = histogram(m1m2, bin_edges, "Normalization", "pdf", "FaceColor", [0.8500 0.3250 0.0980]);
hold on;
h1 = histogram(bm1, bin_edges, "Normalization", "pdf", "FaceColor", [0.4660 0.6740 0.1880]);

xlabel('Distances (mm)');
ylabel('Frequency');
title('Distribution of Mouse-Mouse Distances');
legend({'Exposure M-M', 'Exposure M-B'});

xlim([0 450]);

old_xticks = get(gca, 'XTick');

% Calculate new tick values (divide by 50)
new_xticks = old_xticks / 10;

% Set new tick labels
set(gca, 'XtickLabel', new_xticks);


%% umaps
idx = [92, 79, 80, 89, 90, 107, 121];
mid = [1, 2, 2, 1, 2, 1, 2];
% 244 m1, 178 m2, 179 m2, 280 m1, 281 m2
umapClustCombined = [];
for i=6:7
    umapClustCombined = vertcat(umapClustCombined, dbase(idx(i)).umapClustWT(:, :, mid(i)));
end

%% 3d colors
tdcolors = cat(3, customColors{:});
tdcolors = squeeze(tdcolors);
%%
figure;
% colors = cellfun(@(x) sprintf('#%02x%02x%02x', ...
%     round(x(1)*255), round(x(2)*255), round(x(3)*255)), customColors, 'UniformOutput', false);

for i = 1:numel(customColors)
    customColors{i} = customColors{i} / max(customColors{i});
end

%%
figure;
gscatter(umapClustCombined(:, 1), umapClustCombined(:, 2), umapClustCombined(:, 3), tdcolors', '.', 10); 