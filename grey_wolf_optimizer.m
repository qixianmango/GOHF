function grey_wolf_optimizer()
    collector_tower = [0, 0, 80];
    init_circle_num = 75;
    min_r = 100;
    width = 6;
    dist = width + 5;

    num_wolves = 5;
    max_iter = 50;
    dim = 2;
    lb = [3, 3];
    ub = [8, 8];
    positions = lb + (ub - lb) .* rand(num_wolves, dim);

    function unit_energy = calculate_objective(install_height, width, collector_tower, init_circle_num, min_r, dist)
        height = install_height;
        [~, data_len] = generate_concentric_circles(collector_tower, init_circle_num, min_r, dist);
        [~, unit_energy, ~] = evaluate_performance(repmat(install_height, data_len, 1), repmat(height, data_len, 1), ...
            repmat(width, data_len, 1), init_circle_num, min_r, dist, collector_tower);
    end

    total_calculations = 0;
    for iter = 1:max_iter
       
        fitness = zeros(num_wolves, 1);
        for i = 1:num_wolves
            install_height = positions(i, 1);
            width = positions(i, 2);
            dist = width + 5;
            fitness(i) = calculate_objective(install_height, width, collector_tower, init_circle_num, min_r, dist);
            total_calculations = total_calculations + 1;
        end

        [sorted_fitness, idx] = sort(fitness, 'descend');
        alpha_pos = positions(idx(1), :);
        beta_pos = positions(idx(2), :);
        delta_pos = positions(idx(3), :);

        a = 2 - iter * (2 / max_iter);
        for i = 1:num_wolves
            if i ~= idx(1) && i ~= idx(2) && i ~= idx(3)
                A1 = 2 * a * rand(1, dim) - a;
                C1 = 2 * rand(1, dim);
                A2 = 2 * a * rand(1, dim) - a;
                C2 = 2 * rand(1, dim);
                A3 = 2 * a * rand(1, dim) - a;
                C3 = 2 * rand(1, dim);

                D_alpha = abs(C1 .* alpha_pos - positions(i, :));
                D_beta = abs(C2 .* beta_pos - positions(i, :));
                D_delta = abs(C3 .* delta_pos - positions(i, :));

                X1 = alpha_pos - A1 .* D_alpha;
                X2 = beta_pos - A2 .* D_beta;
                X3 = delta_pos - A3 .* D_delta;
                new_position = (X1 + X2 + X3) / 3;

                new_position = max(new_position, lb);
                new_position = min(new_position, ub);
                positions(i, :) = new_position;
            end
        end

    end

    install_height = alpha_pos(1);
    width = alpha_pos(2);
    height = install_height;
    dist = width + 5;
    [~, data_len] = generate_concentric_circles(collector_tower, init_circle_num, min_r, dist);
    [final_energy, final_unit_energy, final_eta] = evaluate_performance(repmat(install_height, data_len, 1), repmat(height, data_len, 1), ...
        repmat(width, data_len, 1), init_circle_num, min_r, dist, collector_tower);

end

function [points, nums, circle_r, point_num, group] = generate_concentric_circles(center, circle_num, min_r, dist)
    points = [];
    point_num = zeros(1, circle_num);
    circle_r = zeros(1, circle_num);
    circle_r(1) = min_r;
    for i = 2:circle_num
        circle_r(i) = circle_r(i-1) + max(dist(i-1), dist(i));
    end
    prev_num_points = 0;
    group = [];
    
    for i = 1:circle_num
        R = circle_r(i);
        circumference = 2 * pi * R;
        num_points = floor(circumference / dist(i));
        delta_theta = 2*pi / num_points;
        
        if i > 1
            offset = pi / (prev_num_points + num_points);
        else
            offset = 0;
        end

        for j = 1:num_points
            theta = (j-1) * delta_theta + offset;
            x = center(1) + R * cos(theta);
            y = center(2) + R * sin(theta);
            points = [points; x y];
            group = [group, i];
        end
    
        prev_num_points = num_points;
    end
    
    distances = sqrt(sum((points - [0, 0]).^2, 2));
    points(distances > 350, :) = [];
    nums = size(points, 1);
    group = group(1:nums);
end

function [sum_energy, unit_energy, avg_eta, results_table] = evaluate_performance(install_heights, heights, widths, circle_num, min_r, dist, collector_tower)
    positions = generate_concentric_circles(collector_tower, circle_num, min_r, dist);
    positions = [positions, install_heights];

    dates = {'2025-1-21', '2025-2-21', '2025-3-21', '2025-4-21', '2025-5-21', '2025-6-21', '2025-7-21', '2025-8-21','2025-9-21', '2025-10-21','2025-11-21', '2025-12-21'};
    times = {'9:00:00', '10:30:00', '12:00:00', '13:30:00', '15:00:00'};
    datetimeMatrix = repmat(datetime('now'), 12, 5);
    for i = 1:length(dates)
        for j = 1:length(times)
            datetimeString = [dates{i} ' ' times{j}];
            datetimeMatrix(i, j) = datetime(datetimeString, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        end
    end
    
    num_mirrors = size(positions, 1);
    num_entries = num_mirrors * length(dates) * length(times);
    results_table = table('Size', [num_entries, 5], 'VariableTypes', {'datetime', 'double', 'double', 'double', 'double'}, 'VariableNames', {'DateTime', 'MirrorID', 'OpticalEfficiency', 'TruncationEfficiency', 'OutputPower'});
    idx = 1;

    total_energy = 0;
    total_optical_efficiency = 0;
    for j = 1:length(dates)
        for i = 1:length(times)
            results = compute_optical_efficiencies(datetimeMatrix(j, i), positions, heights, widths, collector_tower);
            for m = 1:num_mirrors
                results_table.DateTime(idx) = datetimeMatrix(j, i);
                results_table.MirrorID(idx) = m;
                results_table.OpticalEfficiency(idx) = results(m).optical_efficiency;
                results_table.TruncationEfficiency(idx) = results(m).truncation_efficiency;
                results_table.OutputPower(idx) = results(m).output_power;
                total_energy = total_energy + results(m).output_power;
                total_optical_efficiency = total_optical_efficiency + results(m).optical_efficiency;
                idx = idx + 1;
            end
        end
    end

    total_energy = total_energy / (length(dates) * length(times));
    avg_eta = total_optical_efficiency / (num_mirrors * length(dates) * length(times));
    total_area = sum(widths .* heights);
    unit_energy = total_energy / total_area;
end

function [results] = compute_optical_efficiencies(date, mirror_positions, mirror_heights, mirror_widths, tower_position)
    num_mirrors = size(mirror_positions, 1);
    [sun_elevation, sun_azimuth, sun_direction] = calculate_sun_position(98.5, 39.4, date);
    sun_direction = real(sun_direction);
    mirror_normals = zeros(num_mirrors, 3);
    mirror_corners = zeros(4, 3, num_mirrors);

    for i = 1:num_mirrors
        mirror_normals(i, :) = derive_mirror_normal(mirror_positions(i, :), sun_direction, tower_position);
        [mirror_corners(1, :, i), mirror_corners(2, :, i), mirror_corners(3, :, i), mirror_corners(4, :, i)] = ...
            define_mirror_corners(mirror_positions(i, :), mirror_heights(i), mirror_widths(i), mirror_normals(i, :));
    end

    shadow_efficiency = zeros(num_mirrors, 1);
    truncation_efficiency = zeros(num_mirrors, 1);
    cosine_efficiency = zeros(num_mirrors, 1);
    atmospheric_transmittance = zeros(num_mirrors, 1);

    for i = 1:num_mirrors
        [shadow_efficiency(i), truncation_efficiency(i), cosine_efficiency(i), atmospheric_transmittance(i)] = ...
            evaluate_efficiency(mirror_positions, i, sun_direction, mirror_normals(i, :), tower_position, mirror_corners);
    end

    mirror_reflectivity = 0.92;
    solar_constant = 1.366;
    altitude = 3;
    a = 0.4237 - 0.00821 * (6 - altitude)^2;
    b = 0.5055 + 0.00595 * (6.5 - altitude)^2;
    c = 0.2711 + 0.01858 * (2.5 - altitude)^2;
    direct_normal_irradiance = solar_constant * (a + b * exp(-c / sind(sun_elevation)));

    results = struct('optical_efficiency', [], 'truncation_efficiency', [], 'output_power', []);

    for i = 1:num_mirrors
        optical_efficiency = shadow_efficiency(i) * truncation_efficiency(i) * cosine_efficiency(i) * atmospheric_transmittance(i) * mirror_reflectivity;
        output_power = direct_normal_irradiance * mirror_heights(i) * mirror_widths(i) * optical_efficiency;
        results(i).optical_efficiency = optical_efficiency;
        results(i).truncation_efficiency = truncation_efficiency(i);
        results(i).output_power = output_power;
    end
end

function [sun_elevation, sun_azimuth, sun_direction] = calculate_sun_position(lon_deg, lat_deg, date_time)
    lon_rad = deg2rad(lon_deg);
    lat_rad = deg2rad(lat_deg);
    spring_equinox = datetime(year(date_time), 3, 21, 0, 0, 0);
    D = days(date_time - spring_equinox);

    sin_delta = sin(2*pi*D / 365) * sin(2*pi*23.45 / 360);
    delta_rad = asin(sin_delta);

    local_time_hours = hour(date_time) + minute(date_time) / 60;

    w = (local_time_hours - 12) * (pi / 12);
    sin_as = sin(lat_rad) * sin(delta_rad) + cos(lat_rad) * cos(delta_rad) * cos(w);
    sun_elevation_rad = asin(sin_as);

    cos_sun_azimuth = (sin_delta - sin_as * sin(lat_rad)) / (cos(sun_elevation_rad) * cos(lat_rad));
    sun_azimuth_rad = acos(cos_sun_azimuth);

    x = cos(sun_elevation_rad) * sin(sun_azimuth_rad);
    y = cos(sun_elevation_rad) * cos(sun_azimuth_rad);
    z = sin(sun_elevation_rad);
    sun_direction = [x, y, z];

    sun_elevation = rad2deg(sun_elevation_rad);
    sun_azimuth = rad2deg(sun_azimuth_rad);
end

function mirror_normal = derive_mirror_normal(mirror_position, sun_direction, tower_position)
    reflection_direction = tower_position - mirror_position;
    reflection_direction = reflection_direction / norm(reflection_direction);
    mirror_normal = sun_direction + reflection_direction;
    mirror_normal = mirror_normal / norm(mirror_normal);
end

function mirror_corners = define_mirror_corners(mirror_position, mirror_height, mirror_width, mirror_normal)
    up_vector = [0, 0, 1];
    side_vector1 = cross(mirror_normal, up_vector);
    side_vector1 = side_vector1 / norm(side_vector1);
    side_vector2 = cross(mirror_normal, side_vector1);
    side_vector2 = side_vector2 / norm(side_vector2);

    mirror_corners = zeros(4, 3);
    mirror_corners(1, :) = mirror_position - (mirror_width / 2) * side_vector1 - (mirror_height / 2) * side_vector2;
    mirror_corners(2, :) = mirror_position + (mirror_width / 2) * side_vector1 - (mirror_height / 2) * side_vector2;
    mirror_corners(3, :) = mirror_position + (mirror_width / 2) * side_vector1 + (mirror_height / 2) * side_vector2;
    mirror_corners(4, :) = mirror_position - (mirror_width / 2) * side_vector1 + (mirror_height / 2) * side_vector2;
end

function [shadow_efficiency, truncation_efficiency, cosine_efficiency, atmospheric_transmittance] = ...
    evaluate_efficiency(mirror_positions, mirror_index, sun_direction, mirror_normal, tower_position, mirror_corners)
    num_points = 20;
    nearby_distance = 20;
    num_epochs = 1;

    target_mirror = mirror_positions(mirror_index, 1:3);
    distances = sqrt(sum((target_mirror - mirror_positions).^2, 2));
    nearby_indices = distances <= nearby_distance & distances > 0;
    nearby_mirror_corners = mirror_corners(:, :, nearby_indices);
    points = sample_points(mirror_corners(:, :, mirror_index), num_points);
    counter_not_blocked = 0;
    counter_reach_collector = 0;

    for e = 1:num_epochs
        for i = 1:size(points, 1)
            blocked = false;
            new_direction = perturb_direction(sun_direction, 4.65e-3);
            reflection_direction = reflect_ray(-sun_direction, mirror_normal);
            blocked = detect_cylinder_overlap(points(i, :), new_direction, [tower_position(1), tower_position(2), 0], 84, 3.5);
            for j = 1:size(nearby_mirror_corners, 3)
                if blocked
                    break;
                end
                blocked = blocked | detect_rectangle_overlap(points(i, :), new_direction, nearby_mirror_corners(:, :, j));
                blocked = blocked | detect_rectangle_overlap(points(i, :), reflection_direction, nearby_mirror_corners(:, :, j));
            end
            if ~blocked
                counter_not_blocked = counter_not_blocked + 1;
                if detect_cylinder_overlap(points(i, :), reflection_direction, [tower_position(1), tower_position(2), 0], 84, 3.5)
                    counter_reach_collector = counter_reach_collector + 1;
                end
            end
        end
    end

    shadow_efficiency = counter_not_blocked / (num_epochs * size(points, 1));

    if counter_not_blocked == 0
        truncation_efficiency = 0;
    else
        truncation_efficiency = counter_reach_collector / counter_not_blocked;
    end

    cosine_efficiency = cos(0.5 * acos(dot(real(sun_direction), real(reflection_direction))));
    assert(isreal(cosine_efficiency));
    d_HR = norm(tower_position - target_mirror);
    atmospheric_transmittance = 0.99321 - 0.0001176 * d_HR + 1.97e-8 * d_HR^2;
end

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

function new_direction = perturb_direction(original_direction, max_offset_angle)
    sigma = 1/3 * atan(max_offset_angle);
    original_direction = original_direction / norm(original_direction);
    if original_direction(1) ~= 0 || original_direction(2) ~= 0
        temp_vector = [0; 0; 1];
    else
        temp_vector = [1; 0; 0];
    end

    basis1 = cross(original_direction, temp_vector);
    basis1 = basis1 / norm(basis1);

    basis2 = cross(original_direction, basis1);
    basis2 = basis2 / norm(basis2);

    rand_point = sigma * randn(2, 1);

    point = (basis1 * rand_point(1, :) + basis2 * rand_point(2, :))';
    new_direction = original_direction + point';
    new_direction = new_direction / norm(new_direction);
end

function reflection = reflect_ray(input, mirror_normal)
    reflection = input - 2 * dot(input, mirror_normal) * mirror_normal;
    reflection = reflection / norm(reflection);
end

function is_intersect = detect_cylinder_overlap(ray_origin, ray_direction, cylinder_center, cylinder_height, cylinder_radius)
    t1 = (cylinder_center(3) - ray_origin(3)) / ray_direction(3);
    t2 = (cylinder_center(3) + cylinder_height - ray_origin(3)) / ray_direction(3);

    p1 = ray_origin + t1 .* ray_direction;
    p2 = ray_origin + t2 .* ray_direction;

    if t1 > 0 && (norm(p1(1:2) - cylinder_center(1:2)) <= cylinder_radius)
        is_intersect = false;
        return;
    elseif t2 > 0 && (norm(p2(1:2) - cylinder_center(1:2)) <= cylinder_radius)
        is_intersect = true;
        return;
    end

    d = ray_origin(1:2) - cylinder_center(1:2);
    A = ray_direction(1)^2 + ray_direction(2)^2;
    B = 2 * (ray_direction(1)*d(1) + ray_direction(2)*d(2));
    C = d(1)^2 + d(2)^2 - cylinder_radius^2;

    discriminant = B^2 - 4*A*C;

    if discriminant >= 0
        t3 = (-B + sqrt(discriminant)) / (2*A);
        t4 = (-B - sqrt(discriminant)) / (2*A);

        p3 = ray_origin + t3 .* ray_direction;
        p4 = ray_origin + t4 .* ray_direction;

        if t3 > 0 && (cylinder_center(3) <= p3(3) && p3(3) <= cylinder_center(3) + cylinder_height)
            is_intersect = true;
            return;
        elseif t4 > 0 && (cylinder_center(3) <= p4(3) && p4(3) <= cylinder_center(3) + cylinder_height)
            is_intersect = true;
            return;
        end
    end

    is_intersect = false;
end

function is_intersect = detect_rectangle_overlap(ray_origin, ray_direction, rect_corners)
    edge1 = rect_corners(2, :) - rect_corners(1, :);
    edge2 = rect_corners(3, :) - rect_corners(1, :);
    normal = cross(edge1, edge2);
    denom = dot(normal, ray_direction);
    t = dot(rect_corners(1,:) - ray_origin, normal) / denom;

    if t < 0
        is_intersect = false;
        return;
    end

    intersection_point = ray_origin + t * ray_direction;

    AB = rect_corners(2,:) - rect_corners(1,:);
    AM = intersection_point - rect_corners(1,:);
    BC = rect_corners(3,:) - rect_corners(2,:);
    BM = intersection_point - rect_corners(2,:);

    if (dot(AB, AM) >= 0) && (dot(AB, AM) <= dot(AB, AB)) && ...
       (dot(BC, BM) >= 0) && (dot(BC, BM) <= dot(BC, BC))
        is_intersect = true;
    else
        is_intersect = false;
    end
end