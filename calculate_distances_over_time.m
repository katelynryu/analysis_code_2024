function distances = calculate_distances_over_time(points1, points2)
    % points1 and points2 should be 3-dimensional arrays with dimensions (time, noi, 2)

    % Get the number of time steps
    num_time_steps = size(points1, 1);

    % Initialize array to store distances over time
    distances = zeros(num_time_steps, 1);

    % Calculate distances at each time step
    for t = 1:num_time_steps
        % Extract x and y coordinates from each point at time t
        x1 = points1(t, 1, 1);
        y1 = points1(t, 1, 2);

        x2 = points2(t, 1, 1);
        y2 = points2(t, 1, 2);

        % Calculate the Euclidean distance between the points at time t
        distances(t) = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    end
end
