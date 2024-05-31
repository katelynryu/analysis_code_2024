% constants

timex = dbase(1).time-dbase(1).time(1,1);
timex = size(timex, 1);
%% linear velocity w/o mice in m/s

for i=1:102%length(dbase) 
   if(sum(isnan(dbase(1).tracks(:, 1, :)), "all") < sum(isnan(dbase(1).tracks(:, 2, :)), "all"))
       b = 1;
   else
       b = 2;
   end
   dbase(i).linvel = getVelocity(squeeze(dbase(i).tracks(1:(10*fps),b,:)), 20, 80, 1/1.97);
end

%% plot linear velocities
figure
for i=57:length(dbase)
    plot(dbase(i).linvel);
    hold on;
end 

%% angular velocity w/o mice in rad/s

for i=1:length(dbase)
   if(sum(isnan(dbase(1).tracks(:, 1, :)), "all") < sum(isnan(dbase(1).tracks(:, 2, :)), "all"))
       b = 1;
   else
       b = 2;
   end
   dbase(i).angvel = getAngularVelocity(dbase(i).tracks(1:(10*fps), b, :)) * fps;
end

%% plot angular velocities

figure
for i=1:8%length(dbase)
    plot(dbase(i).angvel);
    hold on;
end 

%% scatterplot linear vs angular velocity 

figure
for i=[1:8, 57:length(dbase)]
    %figure
    scatter(dbase(i).linvel, dbase(i).angvel);
    xlabel("Linear Velocity (m/s)");
    ylabel("Angular Velocity (rad/s)");
    title("linear vs. angular velocity");
    %saveas(gcf, strcat(dbase(i).fileID,'-bugscatter.png'));
end

%% 2D histogram (heatmap) linear vs angular velocity
n = 8; % num vectors 
times = size(dbase(1).linvel, 1);

% combine linear and angular velocities
combined_linvelm = zeros(times, n);
combined_angvelm = zeros(times, n);

for i=1:8
    combined_linvelm(:, i) = dbase(i).linvel;
    combined_angvelm(:, i) = dbase(i).angvel;
end


%% vertical concat version
% combine linear and angular velocities
combined_linvel = [];
combined_angvel = [];

for i=[1:8, 57, 59:length(dbase)]
    combined_linvel = vertcat(combined_linvel, dbase(i).linvel);
    combined_angvel = vertcat(combined_angvel, dbase(i).angvel);
end
%% 2D histogram 

figure
histogram2(combined_linvel, combined_angvel);
xlabel("Linear Velocity (m/s)");
ylabel("Angular Velocity (rad/s)");
title("2D histogram of linear vs. angular velocity");
colorbar;
saveas(gcf, "2D-hist-bug-10s.png");

%% heatmap
figure
combined_linvelm(isnan(combined_linvelm)) = 0;
combined_angvelm(isnan(combined_angvelm)) = 0;
imagesc(combined_linvelm, combined_angvelm);

%% linear velocity w mice (full recording)

for i=1:102%length(dbase)
   if(sum(isnan(dbase(1).tracks(:, 1, :)), "all") < sum(isnan(dbase(1).tracks(:, 2, :)), "all"))
       b = 1;
   else
       b = 2;
   end
   dbase(i).full_linvel = getVelocity(squeeze(dbase(i).tracks(:,b,:)), 20, 80, 1/1.97);
end

%% angular velocity w mice (full recording)

for i=1:102%length(dbase)
   if(sum(isnan(dbase(1).tracks(:, 1, :)), "all") < sum(isnan(dbase(1).tracks(:, 2, :)), "all"))
       b = 1;
   else
       b = 2;
   end
   dbase(i).full_angvel = real(getAngularVelocity(dbase(i).tracks(:, b, :))) * fps;
end


%% 
% combine full linear and angular velocities 
combined_linvel_full = [];
combined_angvel_full = [];

for i=1:102%[1:8, 57, 59:length(dbase)]
    combined_linvel_full = vertcat(combined_linvel_full, dbase(i).full_linvel);
    combined_angvel_full = vertcat(combined_angvel_full, dbase(i).full_angvel);
end

%% %% 
% combine full linear and angular velocities 
combined_linvel_all = [];
combined_angvel_all = [];

for i=1:102%length(dbase)
    combined_linvel_all = vertcat(combined_linvel_all, dbase(i).full_linvel);
    combined_angvel_all = vertcat(combined_angvel_all, dbase(i).full_angvel);
end

%% 2D histogram for full time track
xbins = 0:0.004:0.4;
ybins = -1:0.01:1;

figure
h3 = histogram2(combined_linvel_full, combined_angvel_full, xbins, ybins, "Normalization", "probability")%, ...
    %'DisplayStyle','tile','ShowEmptyBins','on');
xlabel("Linear Velocity (m/s)");
ylabel("Angular Velocity (rad/s)");
title("2D histogram of linear vs. angular velocity");
xlim([0 0.4]);
ylim([-1 1]);
colorbar;
%saveas(gcf, "2D-hist-bug-full.png");
h3_values = h3.Values;


%%
C = imagesc(h3_values);
set(gca, 'YDir', 'normal');
yticks(1:20:101);
xticks(1:20:201);
yticklabels(xbins(1:20:end));
xticklabels(ybins(1:20:end));
ylim([1 102]);
xlim([1 202]);
colormap(parula);
colorbar;
title("Object Linear vs. Angular Velocity");
xlabel("Angular Velocity (rad/s)");
ylabel("Linear Velocity (m/s)");

%% clustering eval
condt = combined_linvel_full < 0.4 & combined_angvel_full > -1 & combined_angvel_full < 1;
meas = horzcat(combined_linvel_full(condt), combined_angvel_full(condt));
X = meas(:, 1:2);

%% 
rng("default");
evaluation = evalclusters(meas, "kmeans", "CalinskiHarabasz", "KList", 1:6); % k=6

%% gmm

rng("default");
k = 3; % Number of GMM components
options = statset('MaxIter',1000);

Sigma = {'diagonal','full'}; % Options for covariance matrix type
nSigma = numel(Sigma);

SharedCovariance = {true,false}; % Indicator for identical or nonidentical covariance matrices
SCtext = {'true','false'};
nSC = numel(SharedCovariance);

%%
d = 100; % Grid length
x1 = linspace(min(X(:,1))-2, max(X(:,1))+2, d);
x2 = linspace(min(X(:,2))-2, max(X(:,2))+2, d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];
%X0 = [0:0.004:0.4 -1:0.01:1];

%% 
figure
threshold = sqrt(chi2inv(0.99,2));
count = 1;
for i = 1:nSigma
    for j = 1:nSC
        gmfit1 = fitgmdist(X,k,'CovarianceType',Sigma{i}, ...
            'SharedCovariance',SharedCovariance{j},'Options',options); % Fitted GMM
        clusterX = cluster(gmfit1,X); % Cluster index 
        mahalDist = mahal(gmfit1,X0); % Distance from each grid point to each GMM component
        % Draw ellipsoids over each GMM component and show clustering result.
        subplot(2,2,count);
        h1 = gscatter(X(:,1),X(:,2),clusterX);
        hold on
            for m = 1:k
                idx = mahalDist(:,m)<=threshold;
                Color = h1(m).Color*0.75 - 0.5*(h1(m).Color - 1);
                h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
                uistack(h2,'bottom');
            end    
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        title(sprintf('Sigma is %s\nSharedCovariance = %s',Sigma{i},SCtext{j}),'FontSize',8)
        legend(h1,{'1','2','3'})
        hold off
        count = count + 1;
    end
end

%%
figure
threshold = sqrt(chi2inv(0.99,2));
count = 1;
for i = 1 % sigma = diagonal
    for j = 2 % shared covariance = false  
        gmfit1 = fitgmdist(X,k,'CovarianceType',Sigma{i}, ...
            'SharedCovariance',SharedCovariance{j},'Options',options); % Fitted GMM
        clusterX = cluster(gmfit1,X); % Cluster index 
        mahalDist = mahal(gmfit1,X0); % Distance from each grid point to each GMM component
        % Draw ellipsoids over each GMM component and show clustering result.
        subplot(2,2,count);
        h1 = gscatter(X(:,1),X(:,2),clusterX);
        hold on
            for m = 1:k
                idx = mahalDist(:,m)<=threshold;
                Color = h1(m).Color*0.75 - 0.5*(h1(m).Color - 1);
                h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
                uistack(h2,'bottom');
            end    
        plot(gmfit1.mu(:,1),gmfit1.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        title(sprintf('Sigma is %s\nSharedCovariance = %s',Sigma{i},SCtext{j}),'FontSize',8)
        legend(h1,{'1','2','3'})
        hold off
        count = count + 1;
    end
end

%% 
points = [0];
for i=[1:8, 57, 59:64]
    condt = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    points = vertcat(points, length(dbase(i).tracks(condt)));
    cumpoints = cumsum(points);
end 

%% clustered tracks for bug
j = 1;
for i=[1:8, 57, 59:64]
    vcond = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    dbase(i).clusters = clusterX((cumpoints(j)+1):cumpoints(j+1));
    dbase(i).clustertracks = dbase(i).tracks(vcond, :, :);
    j = j+1;
end
%%
% Iterate over each cluster
for i=5
    figure('Position', [100, 100, 3000, 400]);
    if(sum(isnan(dbase(i).tracks(:, 1, :)), "all") < sum(isnan(dbase(i).tracks(:, 2, :)), "all"))
       b = 1;
    else
       b = 2;
    end
    for cluster = 1:3
        % Extract points belonging to the current cluster
        cluster_points = dbase(i).clustertracks(dbase(i).clusters == cluster, b, :); % Extract x and y positions
        
        % Create a heatmap for the current cluster
        subplot(1, 3, cluster); % Divide the figure into 3 subplots
        scatter(squeeze(cluster_points(:, 1, 1)), squeeze(cluster_points(:, 1, 2))); % Create heatmap
        
        % Customize the heatmap (optional)
        title(['Cluster ', num2str(cluster)]); % Set title for the subplot
        xlabel('X Position'); % Set label for x-axis
        ylabel('Y Position'); % Set label for y-axis
        colorbar; % Display colorbar
    end
end 

%% plot bug position in arena
% Iterate over each cluster
Xedges = 200:12:800;
Yedges = 350:12:950;
for i=8%[1:8]%, 57, 59:64]
    figure('Position', [100, 100, 3000, 400]);
    if(sum(isnan(dbase(i).tracks(:, 1, :)), "all") < sum(isnan(dbase(i).tracks(:, 2, :)), "all"))
       b = 1;
    else
       b = 2;
    end
    for cluster = 1:3
        % Extract points belonging to the current cluster
        cluster_points = dbase(i).clustertracks(dbase(i).clusters == cluster, b, :); % Extract x and y positions
        subplot(1, 3, cluster);
        % Create a heatmap for the current cluster
        h = histcounts2(cluster_points(:, 1, 1), cluster_points(:, 1, 2), Xedges, Yedges, "Normalization", "pdf");
        imagesc(h);
        % subplot(1, 3, cluster); % Divide the figure into 3 subplots
        % x = cluster_points(:, 1, 1);
        % y = cluster_points(:, 1, 2);
        % [counts, edges] = histcounts2(x, y, 100);% Xedges, Yedges);
        % imagesc(edges(1), edges(2), counts)
        % colormap(jet) 
        
        % Customize the heatmap (optional)
        title(['Cluster ', num2str(cluster)]); % Set title for the subplot
        xlabel('X Bin'); % Set label for x-axis
        ylabel('Y Bin'); % Set label for y-axis
        colorbar; % Display colorbar
        caxis([0 5*10^-5]);
    end
    saveas(gcf, strcat(dbase(i).fileID,"bug position.png"));
end 

%% mouse position in arena by cluster
% clustered tracks for mouse 
j = 1;
for i=[1:8, 57, 59:64]
    vcond = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    dbase(i).clustertracks_m = dbase(i).mtracks(vcond, 10, :, :); % tail base only
    j = j+1;
end

%% plot positions
Xedges = 200:12:800;
Yedges = 350:12:950;
for i=8%[1:8]%, 57, 59:64]
    for m=2%:2 % iterate over mouse
        figure('Position', [100, 100, 3000, 400]);
        for cluster = 1:3
            % Extract points belonging to the current cluster
            cluster_points = dbase(i).clustertracks_m(dbase(i).clusters == cluster, 1, :, m); % Extract x and y positions
            
            % Create a heatmap for the current cluster
            subplot(1, 3, cluster); % Divide the figure into 3 subplots
            h = histcounts2(cluster_points(:, 1, 1), cluster_points(:, 1, 2), Xedges, Yedges, "Normalization", "pdf"); %Xedges, Yedges);
            imagesc(h);
            caxis([0 4*10^-5]);
            % x = cluster_points(:, 1, 1);
            % y = cluster_points(:, 1, 2);
            % [counts, edges] = histcounts2(x, y, 50);
            % imagesc(edges(1), edges(2), counts)
            % colormap(jet) 
            
            % Customize the heatmap (optional)
            title(['Cluster ', num2str(cluster)]); % Set title for the subplot
            xlabel('X Position'); % Set label for x-axis
            ylabel('Y Position'); % Set label for y-axis
            colorbar; % Display colorbar
        end
    end
    %saveas(gcf, strcat(dbase(i).fileID,"-mouse", m, "-position.png"));
end 

%%
for i= 95:length(dbase)
    % bug
    xtracks = dbase(i).tracks(:, :, 1); % x bug
    ytracks = dbase(i).tracks(:, :, 2); % y bug
    xmins_b = min(xtracks); % x min bug
    xmin_b = min(xmins_b(xmins_b>=0));
    ymins_b = min(ytracks); % y min bug
    ymin_b = min(ymins_b(ymins_b>=0));

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

%% probability histogram 

Xedges = 0:40:600;
Yedges = 0:40:600;
for i=65:length(dbase)%[1:8]%, 57, 59:64]
    %figure 
    if(sum(isnan(dbase(i).tracks(:, 1, :)), "all") < sum(isnan(dbase(i).tracks(:, 2, :)), "all"))
       b = 1;
    else
       b = 2;
    end
    % Create a heatmap for the current cluster
    totalcounts = numel(dbase(i).tracks0(:, b, :));
    h = histcounts2(dbase(i).tracks0(:, b, 1), dbase(i).tracks0(:, b, 2), Xedges, Yedges);
    dbase(i).hprob = h / totalcounts;
    % imagesc(dbase(i).hprob);
    % 
    % % Customize the heatmap (optional)
    % xlabel('X Bin'); % Set label for x-axis
    % ylabel('Y Bin'); % Set label for y-axis
    % colorbar; % Display colorbar
    %saveas(gcf, strcat(dbase(i).fileID,"bug position.png"));
end 

%% bug probability histogram across all
hsum = dbase(1).hprob;
count = 1;
for i=2:length(dbase)
    hsum = hsum + dbase(i).hprob; 
    count = count+1; 
end
havg = hsum / count;
figure;
imagesc(havg);
colorbar; 

%% probability histogram mice

Xedges = 0:40:600;
Yedges = 0:40:600;
for i=64:length(dbase)%[1:8]%, 57, 59:64]
    for m=1:2
    %figure 
    % Create a heatmap for the current cluster
    totalcounts = numel(dbase(i).mtracks0(:, 10, :, m));
    h = histcounts2(dbase(i).mtracks0(:, 10, 1, m), dbase(i).mtracks0(:, 10, 2, m), Xedges, Yedges);
    dbase(i).hprobm(:, :, m) = h / totalcounts;

    % imagesc(dbase(i).hprob);
    % 
    % % Customize the heatmap (optional)
    % xlabel('X Bin'); % Set label for x-axis
    % ylabel('Y Bin'); % Set label for y-axis
    % colorbar; % Display colorbar
    %saveas(gcf, strcat(dbase(i).fileID,"bug position.png"));
    end
end 

%% probability histogram across all mice
cm_inferno = inferno(100);
hsum = dbase(1).hprob;
count = 1;
for i=2:length(dbase)
    hsum = hsum + dbase(i).hprob; 
    count = count+1; 
end
havg = hsum / count;
figure('Position', [100, 100, 1500, 500]);
subplot(1, 2, 1);
imagesc(havg);
colormap(cm_inferno);
colorbar; 
caxis([0 6*10^-3]);
title("bug");

hsum = dbase(1).hprobm(:, :, 1) + dbase(1).hprobm(:, :, 2);
count = 1;
for i=2:length(dbase)
    hsum = hsum + dbase(i).hprobm(:, :, 1) + dbase(i).hprobm(:, :, 2); 
    count = count+1; 
end
havg = hsum / (2*count);
subplot(1, 2, 2);
imagesc(havg);
title("all mice");
colorbar; 
caxis([0 6*10^-3]);

%% attackers only - mice
hsum = dbase(5).hprobm(:, :, 1) + dbase(6).hprobm(:, :, 1) + dbase(60).hprobm(:, :, 1) + dbase(8).hprobm(:, :, 2);
havg = hsum / 4;
figure;
imagesc(havg);
colorbar;

%% non attackers
hsum = dbase(5).hprobm(:, :, 2) + dbase(6).hprobm(:, :, 2) + dbase(60).hprobm(:, :, 2) + dbase(8).hprobm(:, :, 1);
havg = hsum / 4;
figure;
imagesc(havg);
colorbar;

%% attackers only - bug
hsum = dbase(5).hprob + dbase(6).hprob + dbase(60).hprob + dbase(8).hprob;
havg = hsum / 4;
figure('Position', [100, 100, 3000, 400]);
subplot(1, 3, 1);
imagesc(havg);
colorbar;
caxis([0 3*10^-3]);
title("bug");

hsum = dbase(5).hprobm(:, :, 1) + dbase(6).hprobm(:, :, 1) + dbase(60).hprobm(:, :, 1) + dbase(8).hprobm(:, :, 2);
havg = hsum / 4;
subplot(1, 3, 2);
imagesc(havg);
colorbar;
caxis([0 3*10^-3]);
title("attackers");

hsum = dbase(5).hprobm(:, :, 2) + dbase(6).hprobm(:, :, 2) + dbase(60).hprobm(:, :, 2) + dbase(8).hprobm(:, :, 1);
havg = hsum / 4;
subplot(1, 3, 3);
imagesc(havg);
colorbar;
caxis([0 3*10^-3]);
title("non attackers");
%% across all controls
hsum = dbase(57).hprob;
count = 1;
for i=58:length(dbase)
    hsum = hsum + dbase(i).hprob; 
    count = count+1; 
end
havg = hsum / count;
figure('Position', [100, 100, 1500, 500]);
subplot(1, 2, 1);
imagesc(havg);
colorbar; 
caxis([0 6*10^-3]);
title("bug");

hsum = dbase(57).hprobm(:, :, 1) + dbase(57).hprobm(:, :, 2);
count = 1;
for i=58:length(dbase)
    hsum = hsum + dbase(i).hprobm(:, :, 1) + dbase(i).hprobm(:, :, 2); 
    count = count+1; 
end
havg = hsum / (2*count);
subplot(1, 2, 2);
imagesc(havg);
title("all controls");
colorbar; 
caxis([0 6*10^-3]);


%% across all ELE
hsum = dbase(1).hprob;
count = 1;
for i=2:56
    hsum = hsum + dbase(i).hprob; 
    count = count+1; 
end
havg = hsum / count;
figure('Position', [100, 100, 1500, 500]);
subplot(1, 2, 1);
imagesc(havg);
colorbar; 
caxis([0 6*10^-3]);
title("bug");

hsum = dbase(1).hprobm(:, :, 1) + dbase(1).hprobm(:, :, 2);
count = 1;
for i=2:56
    hsum = hsum + dbase(i).hprobm(:, :, 1) + dbase(i).hprobm(:, :, 2); 
    count = count+1; 
end
havg = hsum / (2*count);
subplot(1, 2, 2);
imagesc(havg);
title("all ELE");
colorbar; 
caxis([0 6*10^-3]);

%% per group
combined_linvel_ELE = [];
combined_angvel_ELE = [];

for i=1:8
    combined_linvel_ELE = vertcat(combined_linvel_ELE, dbase(i).full_linvel);
    combined_angvel_ELE = vertcat(combined_angvel_ELE, dbase(i).full_angvel);
end

figure
histogram2(combined_linvel_ELE, combined_angvel_ELE);
xlabel("Linear Velocity (m/s)");
ylabel("Angular Velocity (rad/s)");
title("2D histogram of linear vs. angular velocity");
colorbar;
xlim([0 0.5]);
ylim([-1 1]);
saveas(gcf, "2D-hist-bug-ELE.png");

%% find negative values
for i=1:length(dbase)
    if(find(dbase(i).mtracks0 < 0) > 0)
        disp(num2str(i));
    end
end

%% smooth bug tracks
for i=95:length(dbase)
    nanmeantracks = nanmean(dbase(i).tracks0, 2);
    dbase(i).smoothedtracks = smoothdata(nanmeantracks, "movmedian");
end
%% plot mice position relative to bug 
for i=95:length(dbase)
    for m=1:2
        dbase(i).mtracks2b(:, 1, 1, m) = dbase(i).mtracks0(:, 1, 1, m) - dbase(i).smoothedtracks(:, 1, 1); 
        dbase(i).mtracks2b(:, 1, 2, m) = dbase(i).mtracks0(:, 1, 2, m) - dbase(i).smoothedtracks(:, 1, 2);
        dbase(i).mtracks2b(:, 2, 1, m) = dbase(i).mtracks0(:, 10, 1, m) - dbase(i).smoothedtracks(:, 1, 1); % tail base
        dbase(i).mtracks2b(:, 2, 2, m) = dbase(i).mtracks0(:, 10, 2, m) - dbase(i).smoothedtracks(:, 1, 2);
    end
end
%% 
Xedges = -600:40:600;
Yedges = -600:40:600;
for i=95:length(dbase)
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
for i=95:length(dbase)
    figure('Position', [100, 100, 800, 600]);

    % nose
    hsum = dbase(i).hprobm2b_n(:, :, 2); %+ dbase(6).hprobm2b(:, :, 1) + dbase(60).hprobm2b(:, :, 1) + dbase(8).hprobm2b(:, :, 2);
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
    
    hsum = dbase(i).hprobm2b_n(:, :, 1); %+ dbase(6).hprobm2b(:, :, 2) + dbase(60).hprobm2b(:, :, 2) + dbase(8).hprobm2b(:, :, 1);
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
    hsum = dbase(i).hprobm2b_tb(:, :, 2); %+ dbase(6).hprobm2b(:, :, 1) + dbase(60).hprobm2b(:, :, 1) + dbase(8).hprobm2b(:, :, 2);
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
    hsum = dbase(i).hprobm2b_tb(:, :, 1); %+ dbase(6).hprobm2b(:, :, 2) + dbase(60).hprobm2b(:, :, 2) + dbase(8).hprobm2b(:, :, 1);
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

%% average across
figure('Position', [100, 100, 800, 600]);

% nose
hsum = dbase(95).hprobm2b_n(:, :, 2) + dbase(96).hprobm2b_n(:, :, 2) + dbase(97).hprobm2b_n(:, :, 2) + dbase(98).hprobm2b_n(:, :, 2);
havg = hsum / 4;
subplot(2, 2, 1);
imagesc(havg);
colormap(cm_inferno);
colorbar;
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end));
xticklabels(Yedges(1:2:end));
caxis([0 1.5*10^-3]);
title("attackers-nose");

hsum = dbase(95).hprobm2b_n(:, :, 1) + dbase(96).hprobm2b_n(:, :, 1) + dbase(97).hprobm2b_n(:, :, 1) + dbase(98).hprobm2b_n(:, :, 1);
havg = hsum / 4;
subplot(2, 2, 2);
imagesc(havg);
colorbar;
caxis([0 1.5*10^-3]);
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end));
xticklabels(Yedges(1:2:end));
title("non attackers-nose");

% tail base
hsum = dbase(95).hprobm2b_tb(:, :, 2) + dbase(96).hprobm2b_tb(:, :, 2) + dbase(97).hprobm2b_tb(:, :, 2) + dbase(98).hprobm2b_tb(:, :, 2);
havg = hsum / 4;
subplot(2, 2, 3);
imagesc(havg);
colormap(cm_inferno);
colorbar;
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end));
xticklabels(Yedges(1:2:end));
caxis([0 1.5*10^-3]);
title("attackers-tail base");

% non-attackers
hsum = dbase(95).hprobm2b_tb(:, :, 1) + dbase(96).hprobm2b_tb(:, :, 1) + dbase(97).hprobm2b_tb(:, :, 1) + dbase(98).hprobm2b_tb(:, :, 1);
havg = hsum / 4;
subplot(2, 2, 4);
imagesc(havg);
colorbar;
caxis([0 1.5*10^-3]);
yticks(1:2:31);
xticks(1:2:31);
yticklabels(Xedges(1:2:end));
xticklabels(Yedges(1:2:end));
title("non attackers-tail base");
