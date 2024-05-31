
%% 
points = [0];
for i=1:102%length(dbase)%[1:8, 57, 59:64]
    condt = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    points = vertcat(points, length(dbase(i).tracks(condt)));
    cumpoints = cumsum(points);
end 

%% clustered tracks for bug
j = 1;
for i=1:102%length(dbase)
    vcond = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    dbase(i).clusters = clusterX((cumpoints(j)+1):cumpoints(j+1));
    dbase(i).clustertracks = dbase(i).tracks0(vcond, :, :);
    j = j+1;
end

%% plot bug position in arena
% Iterate over each cluster
Xedges = 0:16:800;
Yedges = 0:16:800;
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
        h = histcounts2(cluster_points(:, 1, 1), cluster_points(:, 1, 2), Xedges, Yedges);
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
    end
    %saveas(gcf, strcat(dbase(i).fileID,"bug position.png"));
end 
%% mouse position in arena by cluster
% clustered tracks for mouse 
j = 1;
for i=1:102%length(dbase)%[1:8, 57, 59:64]
    vcond = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    dbase(i).clustertracks_m2 = dbase(i).mtracks0(dbase(dbase(i).tracks0(vcond, :, :), 1, :, :);
    %dbase(i).clustertracks_m = dbase(i).mtracks0(vcond, 10, :, :); % tail base only
    j = j+1;
end

%% plot positions
Xedges = 0:16:800;
Yedges = 0:16:800;

for i=8%[1:8]%, 57, 59:64]
    for m=1:2 % iterate over mouse
        figure('Position', [100, 100, 3000, 400]);
        for cluster = 1:3
            % Extract points belonging to the current cluster
            cluster_points = dbase(i).clustertracks_m(dbase(i).clusters == cluster, 1, :, m); % Extract x and y positions
            
            % Create a heatmap for the current cluster
            subplot(1, 3, cluster); % Divide the figure into 3 subplots
            h = histcounts2(cluster_points(:, 1, 1), cluster_points(:, 1, 2), Xedges, Yedges);
            imagesc(h);
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


%% plot positions of mice without bug 
Xedges = 30:12:630;
Yedges = 30:12:630;

for i=109%[1:8]%, 57, 59:64]
    for m=1%:2 % iterate over mouse
        figure; %('Position', [100, 100, 3000, 400]);
        %for cluster = 1:3
            % Extract points belonging to the current cluster
            cluster_points = dbase(i).mtracks0(:, 1, :, m); % Extract x and y positions
            
            % Create a heatmap for the current cluster
            %subplot(1, 3, cluster); % Divide the figure into 3 subplots
            h = histcounts2(cluster_points(:, 1, 1), cluster_points(:, 1, 2), Xedges, Yedges, "Normalization", "pdf");
            imagesc(h);
            
            % Customize the heatmap (optional)
            title(['Cluster ', num2str(cluster)]); % Set title for the subplot
            xlabel('X Position'); % Set label for x-axis
            ylabel('Y Position'); % Set label for y-axis
            colorbar; % Display colorbar
            caxis([0 4*10^-5]);
  
    end
    %saveas(gcf, strcat(dbase(i).fileID,"-mouse", m, "-position.png"));
end 

%% plot mice position in arena when idle 



%% 
umapcor=dbase(5).umapClustWT(:,1:2,2);
bugMd = squeeze(dbase(5).mtracks2b(:,1,:,1));

h=histogram2(umapcor(bugMd(:,1)>50 & bugMd(:,2)>50,1),umapcor(bugMd(:,1)>50 & bugMd(:,2)>50,2),50,...
    'DisplayStyle','tile','ShowEmptyBins','on', "Normalization",'pdf');
umapFar=h.Values;

h = histogram2(umapcor(bugMd(:,1)<50 & bugMd(:,2)<50,1),umapcor(bugMd(:,1)<50 & bugMd(:,2)<50,2),50,...
    'DisplayStyle','tile','ShowEmptyBins','on', "Normalization","pdf");
umapClose=h.Values;
scatter(umapcor(bugMd(:,1)<50 & bugMd(:,2)<50,1),umapcor(bugMd(:,1)<50 & bugMd(:,2)<50,2), 10, [],"filled")