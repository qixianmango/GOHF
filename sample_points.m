function points_inside_rectangle = sample_points(rect_corners, num_points)

    v1 = rect_corners(1, :);
    v2 = rect_corners(2, :);
    v3 = rect_corners(3, :);
    v4 = rect_corners(4, :);

    edge1 = v2 - v1;
    edge2 = v4 - v1;
    area = norm(cross(edge1, edge2));
    increment = sqrt(area / num_points);
    length1 = norm(edge1);
    length2 = norm(edge2);
    num_points_edge1 = round(length1 / increment);
    num_points_edge2 = round(length2 / increment);
    points_inside_rectangle = zeros(num_points_edge1 * num_points_edge2, 3);
    k = 1;
    for i = 1:num_points_edge1
        for j = 1:num_points_edge2
            u = (i - 0.5) * increment / length1;
            v = (j - 0.5) * increment / length2;
            point = v1 + u * edge1 + v * edge2;
            points_inside_rectangle(k, :) = point;
            k = k + 1;
        end
    end
end