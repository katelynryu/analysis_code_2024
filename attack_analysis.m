%% duration of proximity to bug across all mice

% database number 
for i = 1:102
    for m=1:2 % mouse 
        squeezed = squeeze(dbase(i).mtracks2b(:, 1, :, m)); % noi = 1 for nose
        prox = (abs(squeezed(:, 1)) < 50 & abs(squeezed(:, 2)) < 50);
        prox = double(prox);

        % Find the indices where the mask changes (from 0 to 1 or from 1 to 0)
        CC = bwconncomp(prox);
        component_lengths = cellfun(@length, CC.PixelIdxList);

        % Calculate the lengths (duration) of consecutive segments where the mask is equal to 1
        if(m==1)
            dbase(i).prox_m1CC = component_lengths;
        else
            dbase(i).prox_m2CC = component_lengths;
        end
    end
end

%% create array of contact durations across all mice
allprox = [];
for i = 1:102
    allprox = vertcat(allprox, dbase(i).prox_m1CC', dbase(i).prox_m2CC');
end

%% plot allprox histogram

figure;

Xedges = [0, 3, 11, 41, 50];%2987];
% [counts, edges] = histcounts(allprox, bin_edges, 'Normalization', 'probability');
% %
% % Plot the histogram as a line plot
% plot(edges(1:end-1), counts, 'LineWidth', 2);
% 
% histogram(allprox(allprox < 100), 20);

histogram(allprox, 100);

xlabel('Durations');
ylabel('Frequency');
title('Proximity to Bug Durations');


%% plotting only for four attackers
attackprox = vertcat(dbase(5).prox_m2CC',dbase(6).prox_m2CC', dbase(8).prox_m1CC', dbase(60).prox_m2CC');
figure;
Xedges = 0:5:600;
h1 = histogram(attackprox, Xedges, "Normalization", "pdf");
hold on;

noattackprox = vertcat(dbase(5).prox_m1CC',dbase(6).prox_m1CC', dbase(8).prox_m2CC', dbase(60).prox_m1CC');
h2 = histogram(noattackprox, Xedges, "Normalization", "pdf");
legend({'attackers', 'non-attackers'});

%% clustered tracks for bug
j = 1;
for i=1:102%length(dbase)
    vcond = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    dbase(i).clusters2 = clusterX((cumpoints(j)+1):cumpoints(j+1));
    dbase(i).clustertracksm2b = dbase(i).mtracks2b(vcond, :, :, :);
    j = j+1;
end

%% proximity using cluster tracks

for i = 1:102
    for m=1:2
        squeezed = squeeze(dbase(i).clustertracksm2b(dbase(i).clusters == 1, 2, :, m)); % nose
        prox = (abs(squeezed(:, 1)) < 50 & abs(squeezed(:, 2)) < 50);
        prox = double(prox);

        % Find the indices where the mask changes (from 0 to 1 or from 1 to 0)
        CC = bwconncomp(prox);
        component_lengths = cellfun(@length, CC.PixelIdxList);

        % Calculate the lengths (duration) of consecutive segments where the mask is equal to 1
        if(m==1)
            dbase(i).cluster1prox_m1CC = component_lengths;
        else
            dbase(i).cluster1prox_m2CC = component_lengths;
        end
    end
end
%%
cluster1proxt = [];
cluster2proxt = [];
cluster3proxt = [];

for i=1:102
   cluster1proxt = vertcat(cluster1proxt, dbase(i).cluster1prox_m1CC', dbase(i).cluster1prox_m2CC');
   cluster2proxt = vertcat(cluster2proxt, dbase(i).cluster2prox_m1CC', dbase(i).cluster2prox_m2CC');
   cluster3proxt = vertcat(cluster3proxt, dbase(i).cluster3prox_m1CC', dbase(i).cluster3prox_m2CC');
end

%% plot
figure;
Xedges = 0:10:300;
[counts, ~] = histcounts(cluster1proxt, Xedges, "Normalization", "pdf");
cdf1 = cumsum(counts);

plot(Xedges(2:end), cdf1, 'LineWidth', 2);
hold on;

[counts, ~] = histcounts(cluster2proxt, Xedges, "Normalization", "pdf");
cdf2 = cumsum(counts);

plot(Xedges(2:end), cdf2, 'LineWidth', 2);


[counts, ~] = histcounts(cluster3proxt, Xedges, "Normalization", "pdf");
cdf3 = cumsum(counts);

plot(Xedges(2:end), cdf3, 'LineWidth', 2);

legend({'cluster 1', 'cluster 2', 'cluster 3'});
xticks_new = get(gca, 'XTick');
xticklabels_new = arrayfun(@(x) sprintf('%.2f', x), xticks_new/80, 'UniformOutput', false);

% Set the new tick labels
xticks(xticks_new);
xticklabels(xticklabels_new);
xlabel("duration (s)");
ylabel("probability");

%% significance test
%p = kruskalwallis([cluster1proxn, cluster2proxn, cluster3proxn]);
x = [cluster1proxt; cluster2proxt; cluster3proxt];
grouping = [ones(size(cluster1proxt)); 2*ones(size(cluster2proxt)); 3*ones(size(cluster3proxt))];
[p, ~, stats] = kruskalwallis(x, grouping);
disp(['Kruskal-Wallis p-value: ', num2str(p)]);
c = multcompare(stats, 'display', 'on');

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
disp(tbl);

%% just for attackers
attackprox = vertcat(dbase(5).prox_m2CC',dbase(6).prox_m2CC', dbase(8).prox_m1CC', dbase(60).prox_m2CC');
figure;
Xedges = 0:10:600;
[counts, ~] = histcounts(attackprox, Xedges, "Normalization", "pdf");
cdf1 = cumsum(counts);

plot(Xedges(2:end), cdf1, 'LineWidth', 2, 'Color', "#A2142F");
hold on;


noattackprox = vertcat(dbase(5).prox_m1CC',dbase(6).prox_m1CC', dbase(8).prox_m2CC', dbase(60).prox_m1CC');
[counts, ~] = histcounts(noattackprox, Xedges, "Normalization", "pdf");
cdf2 = cumsum(counts);

plot(Xedges(2:end), cdf2, 'LineWidth', 2, 'Color', "#4DBEEE");


legend({'attackers', 'non-attackers'});
xticks_new = get(gca, 'XTick');
xticklabels_new = arrayfun(@(x) sprintf('%.2f', x), xticks_new/80, 'UniformOutput', false);

% Set the new tick labels
xticks(xticks_new);
xticklabels(xticklabels_new);
xlabel("duration (s)");
ylabel("probability");

%%
attackprox = vertcat(dbase(5).prox_m2CC',dbase(6).prox_m2CC', dbase(8).prox_m1CC', dbase(60).prox_m2CC');
figure;
Xedges = 0:5:300;
[counts, ~] = histcounts(attackprox, Xedges, "Normalization", "pdf");
hold on;
counts_reverse = flip(counts);
log_counts_reverse = log(counts_reverse+1);
cdf1 = cumsum(log_counts_reverse); 

plot(Xedges(2:end), cdf1, 'LineWidth', 2);


noattackprox = vertcat(dbase(5).prox_m1CC',dbase(6).prox_m1CC', dbase(8).prox_m2CC', dbase(60).prox_m1CC');
[counts, ~] = histcounts(noattackprox, Xedges, "Normalization", "pdf");
counts_reverse = flip(counts);
log_counts_reverse = log(counts_reverse+1);
cdf2 = cumsum(log_counts_reverse); 
plot(Xedges(2:end), cdf2, 'LineWidth', 2);

legend({'attackers', 'non-attackers'});
xticks_new = get(gca, 'XTick');
xticklabels_new = arrayfun(@(x) sprintf('%.2f', x), xticks_new/80, 'UniformOutput', false);

% Set the new tick labels
xticks(xticks_new);
xticklabels(xticklabels_new);
xlabel("duration (s)");
ylabel("log frequency");

%% ks test 

[h, p, ksstat] = kstest2(cdf1, cdf2, 'Alpha', 0.05);  % Perform KS test with significance level of 0.05

% Display the results
if h == 0
    disp('The null hypothesis cannot be rejected. The two data sets are consistent with each other.');
else
    disp('The null hypothesis is rejected. The two data sets are not consistent with each other.');
end
disp(['p-value: ', num2str(p)]);
disp(['KS statistic: ', num2str(ksstat)]);

%% map time points of proximity to umap space 
umapClose = [];
umapFar = [];
Xedges = [-4:0.22:7];
Yedges = [-7:0.22:4];
bugMdall = [];
umapcorall = [];

for i = 1:102
    for m=1:2
        bugMd = squeeze(dbase(i).mtracks2b(:, 1, :, m));
        bugMdall = vertcat(bugMdall, bugMd);
        umapcor=dbase(i).umapClustWT(:,1:2, m);
        umapcorall = vertcat(umapcorall, umapcor);
    end
end

%%
figure;
h = histogram2(umapcorall(abs(bugMdall(:,1))<50 & abs(bugMdall(:,2))<50, 1), ...
    umapcorall(abs(bugMdall(:,1))<50 & abs(bugMdall(:,2))<50, 2), Xedges, Yedges, ...
    'DisplayStyle','tile','ShowEmptyBins','on', "Normalization",'pdf');
umapClose = h.Values;
colormap(flipud(gray));
colorbar;
caxis([0 0.1]);

figure;
h = histogram2(umapcorall(abs(bugMdall(:,1))>50 & abs(bugMdall(:,2))>50, 1), ...
    umapcorall(abs(bugMdall(:,1))>50 & abs(bugMdall(:,2))>50, 2), Xedges, Yedges, ...
    'DisplayStyle','tile','ShowEmptyBins','on', "Normalization",'pdf');
umapFar = h.Values;
colormap(flipud(gray));
colorbar;
caxis([0 0.1]);

%%
  umapSub = umapClose - umapFar;
  % umapSub(umapSub == 0) = NaN;
  figure;
  imagesc((umapSub)');
  colormap(redblue);
  %colormap([1 1 1;jet(66);0 0 0])
  caxis([-0.1 0.1]);
  % imagesc(umapFar);
  % figure;
  % imagesc((umapClose-umapFar)');
  % colormap(cool);
  set(gca, 'Ydir', 'normal');
  colorbar;


