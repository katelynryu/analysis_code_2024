%% generate csv files for umap
for i=1:length(dbase)
    for m=1:2
        umap = dbase(i).umapClustWT(:, :, m);
        csvname = strcat(dbase(i).fileID, "-m", num2str(m), "-umap.csv");
        writematrix(umap, csvname);
    end
end


%% mouse position in arena by cluster
% clustered tracks for mouse 
j = 1;
for i=1:length(dbase)%[1:8, 57, 59:64]
    vcond = dbase(i).full_linvel < 0.4 & dbase(i).full_angvel > -1 & dbase(i).full_angvel < 1;
    dbase(i).clustertracks_m = dbase(i).mtracks0(vcond, 1, :, :); % nose
    j = j+1;
end

%% plot positions
for i=3%[1:8]%, 57, 59:64]
    for m=1:2 % iterate over mouse
        figure('Position', [100, 100, 3000, 400]);
        for cluster = 1:3
            % Extract points belonging to the current cluster
            cluster_points = dbase(i).clustertracks_m(dbase(i).clusters == cluster, 1, :, m); % Extract x and y positions
            
            % Create a heatmap for the current cluster
            subplot(1, 3, cluster); % Divide the figure into 3 subplots
            h = histcounts2(cluster_points(:, 1, 1), cluster_points(:, 1, 2), 100); %Xedges, Yedges);
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