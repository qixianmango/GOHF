function [points, num_points, radii, total_points, group_ids] = arrange_concentric_circles(center, num_circles, min_radius, spacing)

    points = [];
    total_points = zeros(1, num_circles);
    radii = zeros(1, num_circles);
    radii(1) = min_radius;

    for i = 2:num_circles
        radii(i) = radii(i-1) + max(spacing(i-1), spacing(i));
    end

    prev_points = 0;
    group_ids = [];

    for i = 1:num_circles
        radius = radii(i);
        circumference = 2 * pi * radius;
        num_points = floor(circumference / spacing(i));
        delta_angle = 2*pi / num_points;

        if i > 1
            offset = pi / (prev_points + num_points);
        else
            offset = 0;
        end

        for j = 1:num_points
            angle = (j-1) * delta_angle + offset;
            x = center(1) + radius * cos(angle);
            y = center(2) + radius * sin(angle);
            points = [points; x y];
            group_ids = [group_ids, i];
        end

        prev_points = num_points;
    end

    distances = sqrt(sum((points - [0, 0]).^2, 2));
    points(distances > 350, :) = [];
    total_points = size(points, 1);
    group_ids = group_ids(1:total_points);
end