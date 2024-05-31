%% calculate distances for bug-m1, bug-m2, m1-m2
for t=1:length(dbase)
    dbase(t).ndist = zeros(size(dbase(t).time,1), 3);
    for i = 1:2
        for j = 1:size(dbase(t).time,1)
        % mouse-bug (tail base)
        dbase(t).ndist(j,i) = returnDist(dbase(t).mtracks(j, 1, :, i),  dbase(1).tracks(j, 1, :)) * pixel_size;
        % mouse-mouse (tail base)
        dbase(t).ndist(j,3) = returnDist(dbase(t).mtracks(j, 1, :, 1),  dbase(1).mtracks(j, 1, :, 2)) * pixel_size;
        end
    end
end
%% plot histograms for each video

for i=4%[1:8, 57:length(dbase)]
    figure
    histogram(dbase(i).dist(:, 1), 100);
    hold on;
    histogram(dbase(i).dist(:, 2), 100);
    hold on;
    histogram(dbase(i).dist(:, 3), 100);
    saveas(gcf, strcat(dbase(i).fileID, "-distance-histogram.png"));
end 

%% combined histogram for each condition/group
dist_ELE_M = [dbase(5).dist; dbase(6).dist; dbase(7).dist; dbase(8).dist];
%zeros(size(dbase(1).time, 1) * 4, 3);
dist_ELE_F = [dbase(1).dist; dbase(2).dist; dbase(3).dist; dbase(4).dist];
dist_Control_M = [dbase(57).dist; dbase(58).dist; dbase(59).dist; dbase(60).dist];
dist_Control_F = [dbase(61).dist; dbase(62).dist; dbase(63).dist; dbase(64).dist];

%% histograms
figure
histogram(dist_ELE_M(:, 1), 100);
hold on; 
histogram(dist_ELE_M(:, 2), 100);
histogram(dist_ELE_M(:, 3), 100);
saveas(gcf, "ELE-M-distance-histogram.png");

figure
histogram(dist_ELE_F, 100);
hold on; 
histogram(dist_ELE_F(:, 2), 100);
histogram(dist_ELE_F(:, 3), 100);
saveas(gcf, "ELE-F-distance-histogram.png");

figure
histogram(dist_Control_M, 100);
hold on; 
histogram(dist_Control_M(:, 2), 100);
histogram(dist_Control_M(:, 3), 100);
saveas(gcf, "Control-M-distance-histogram.png");

figure
histogram(dist_Control_F, 100);
hold on; 
histogram(dist_Control_F(:, 2), 100);
histogram(dist_Control_F(:, 3), 100);
saveas(gcf, "Control-F-distance-histogram.png");

%% plot all at once
figure
histogram(dist_ELE_M(:, 2), 100);
hold on;
histogram(dist_ELE_F(:, 2), 100);
histogram(dist_Control_M(:, 2), 100);
histogram(dist_Control_F(:, 2), 100);

%%
hm1 = histogram(dist_ELE_M(:, 1), 100, "normalization", "pdf");
hm2 = histogram(dist_ELE_M(:, 2), 100, "normalization", "pdf");
hm12 = histogram(dist_ELE_M(:, 3), 100, "normalization", "pdf");

mean_m1 = mean(hm1);
std_dev_m1 = mean(hm1);
ci = 1.96 * std_dev_m1 / sqrt(length(hm1)); % 95% confidence interval
ci_lower = mean_m1 - ci;
ci_upper = mean_m1 + ci;

%%
bins = 0:4:400;
bincounts = 100; 
for i=[1:8, 57:64]
    dbase(i).hval = zeros(bincounts, 3);
    for j = 1:3
        h = histogram(dbase(i).dist(:, j), bins, "Normalization", "pdf");
        dbase(i).hval(:,j) = h.Values(1, :);
    end
end 

%%
bins = 0:4:400;
bincounts = 101; 
for i=[1:8, 57:64]
    dbase(i).ksdensity = zeros(bincounts, 3);
    for j = 1:3
        [f, xi] = ksdensity(dbase(i).dist(:, j), bins);
        dbase(i).ksdensity(:,j) = f;
    end
end 

%% combine pdfs for each group 
bdistELE = [];
mdistELE = [];
bdistControl = [];
mdistControl = [];

%ELE 
for i=1:8
    bdistELE = horzcat(bdistELE, dbase(i).hval(:, 1));
    bdistELE = horzcat(bdistELE, dbase(i).hval(:, 2));
    mdistELE = horzcat(mdistELE, dbase(i).hval(:, 3));
end 

%Control
for i=57:64
    bdistControl = horzcat(bdistControl, dbase(i).hval(:, 1));
    bdistControl = horzcat(bdistControl, dbase(i).hval(:, 2));
    mdistControl = horzcat(mdistControl, dbase(i).hval(:, 3));
end

bdistELE = transpose(bdistELE); 
mdistELE = transpose(mdistELE); 
bdistControl = transpose(bdistControl); 
mdistControl = transpose(mdistControl); 
%% stats 
%ELE
bmeanELE = transpose(mean(bdistELE));
mmeanELE = transpose(mean(mdistELE));
bstddevELE = transpose(std(bdistELE));
mstddevELE = transpose(std(mdistELE));

%Control
bmeanControl = transpose(mean(bdistControl));
mmeanControl = transpose(mean(mdistControl));
bstddevControl = transpose(std(bdistControl));
mstddevControl = transpose(std(mdistControl));

%% confidence intervals

bciELE = 1.96 * bstddevELE / sqrt(length(bdistELE)); % 95% confidence interval
ci_lower = bmeanELE - bciELE;
ci_upper = bmeanELE + bciELE;

%%
x = 1:4:400;                     
xconf = [x, x(end:-1:1)];       
[Data-1.96.*Data_sd, fliplr(Data+1.96.*Data_sd)];
yconf = [bmeanELE-1.96*bciELE, fliplr(bmeanELE+1.96*bciELE)];

figure
p = fill(xconf, yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';           

hold on
plot(x, bmeanELE(:, 1),'r')
hold off

%%
figure
plot(x, bmeanELE(:, 1),'r');
hold on
plot(x, mmeanELE(:, 1), 'b');
plot(x, bmeanControl(:, 1), 'k');
plot(x, mmeanControl(:, 1), 'g');
legend('mouse-bug ELE', 'mouse-mouse ELE', 'mouse-bug Control', 'mouse-mouse Control');

%% males only 
bdistELEM = [];
mdistELEM = [];
bdistControlM = [];
mdistControlM = [];

%ELE 
for i=5:8
    bdistELEM = horzcat(bdistELEM, dbase(i).hval(:, 1));
    bdistELEM = horzcat(bdistELEM, dbase(i).hval(:, 2));
    mdistELEM = horzcat(mdistELEM, dbase(i).hval(:, 3));
end 

%Control
for i=57:60
    bdistControlM = horzcat(bdistControlM, dbase(i).hval(:, 1));
    bdistControlM = horzcat(bdistControlM, dbase(i).hval(:, 2));
    mdistControlM = horzcat(mdistControlM, dbase(i).hval(:, 3));
end

bdistELEM = transpose(bdistELEM); 
mdistELEM = transpose(mdistELEM); 
bdistControlM = transpose(bdistControlM); 
mdistControlM = transpose(mdistControlM); 

%% attackers 
bdistAttack = [];
mdistAttack = [];

for i=[5:6, 8, 60]
    bdistAttack = horzcat(bdistAttack, dbase(i).hval(:, 1));
    bdistAttack = horzcat(bdistAttack, dbase(i).hval(:, 2));
    mdistAttack = horzcat(mdistAttack, dbase(i).hval(:, 3));
end 

bmeanAttack = mean(bdistAttack');

%% tail base attack
bdistAttack4 = [];

for i=[5:6, 60]
    bdistAttack4 = horzcat(bdistAttack4, dbase(i).hval(:, 2));
end 

%bdistAttack4 = horzcat(bdistAttack4, dbase(60).hval(:, 2));
bdistAttack4 = horzcat(bdistAttack4, dbase(8).hval(:, 1));

figure
plot(bdistAttack4);
hold on;
bmeanAttack4 = mean(bdistAttack4');
plot(bmeanAttack4', 'k');

%% nose attack
nbdistAttack4 = [];

for i=[5:6, 60]
    nbdistAttack4 = horzcat(nbdistAttack4, dbase(i).nhval(:, 2));
end 

%bdistAttack4 = horzcat(bdistAttack4, dbase(60).nhval(:, 2));
nbdistAttack4 = horzcat(nbdistAttack4, dbase(8).nhval(:, 1));

figure
plot(nbdistAttack4);
hold on;
nbmeanAttack4 = mean(nbdistAttack4');
plot(nbmeanAttack4', 'k');

%% 
figure
plot(dbase(5).ksdensity(:, 2));
hold on;
plot(dbase(6).ksdensity(:, 2));
plot(dbase(60).ksdensity(:, 2));
plot(dbase(8).ksdensity(:, 1));
%%
bdistNoAttack4 = [];

for i=[5:6, 60]
    bdistNoAttack4 = horzcat(bdistNoAttack4, dbase(i).hval(:, 1));
end 

%bdistAttack4 = horzcat(bdistAttack4, dbase(60).hval(:, 2));
bdistNoAttack4 = horzcat(bdistNoAttack4, dbase(8).hval(:, 2));

figure
plot(bdistNoAttack4);
hold on;
bmeannoAttack4 = mean(nbdistAttack4');
plot(bmeannoAttack4', 'k');

%% 
figure
plot(bdistAttack);
hold on;
plot(bmeanAttack', 'k');
figure
plot(mdistAttack);
%% stats male only 
%ELE
bmeanELEM = transpose(mean(bdistELEM));
mmeanELEM = transpose(mean(mdistELEM));
bstddevELEM = transpose(std(bdistELEM));
mstddevELEM = transpose(std(mdistELEM));

%Control
bmeanControlM = transpose(mean(bdistControlM));
mmeanControlM = transpose(mean(mdistControlM));
bstddevControlM = transpose(std(bdistControlM));
mstddevControlM = transpose(std(mdistControlM));

%% plot male only 

plot(x, bmeanELEM(:, 1),'r');
hold on
plot(x, mmeanELEM(:, 1), 'b');
plot(x, bmeanControlM(:, 1), 'k');
plot(x, mmeanControlM(:, 1), 'g');


%% probability density
bins = 0:4:400;
for i=[1:8, 57:length(dbase)]
    dbase(i).pdf = zeros(3, 4);
    for j = 1:3
        h = histogram(dbase(i).dist(:, j), bins, "Normalization", "pdf");
        bin_counts = h.Values;
        bin_edges = h.BinEdges;
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
        
        mean = sum(bin_centers .* bin_counts); % Mean
        variance = sum((bin_centers - mean).^2 .* bin_counts); % Variance
        std_dev = sqrt(variance); % Standard deviation
        
        % Confidence interval (assuming normal distribution)
        alpha = 0.05; % significance level
        n = sum(bin_counts); % sample size
        z = norminv(1 - alpha/2); % Z-score for two-tailed test

        ci_lower = mean - z * std_dev / sqrt(n);
        ci_upper = mean + z * std_dev / sqrt(n);

        dbase(i).pdf(j, :) = [mean, std_dev, ci_lower, ci_upper];
    end 
end 

%%
for i=[57:length(dbase)]
    dbase(i).pdf(1, :)
end

%%
bin_edges = 0:4:400;
ELE_M_bugm1 = zeros(4, 4);
ELE_M_bugm2 = zeros(4, 4);
ELE_M_m1m2 = zeros(4, 4);

for i=[1:8, 57:length(dbase)]
    dbase(i).pdf = zeros(4, 3);
    for j = 1:3
        h = histogram(dbase(i).dist(:, j), bin_edges, "Normalization", "pdf");
        mean = mean(h);
        std_dev = std(h);
        ci = 1.96 * std_dev / sqrt(length(h)); % 95% confidence interval
        ci_lower = mean - ci;
        ci_upper = mean + ci;
        dbase(i).pdf(:, j) = [mean, std_dev, ci_lower, ci_upper];
    end
end