function demo()
    tower_position = [0, 0, 80];
    min_radius = 100;
    max_radius = 350;
    num_circles = 75;
    spacing = 5.75 + 5;
    a0 = 5.75;
    a1 = 2.875;
    population_size = 50;
    max_generations = 100;
    mutation_rate = 0.01;
    crossover_rate = 0.8;
    cone_angle = deg2rad(4.65);
    num_rays = 1000;

    [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths, radii] = ...
        generate_heliostat_layout(tower_position, num_circles, min_radius, spacing, a0, a1);

    population = initialize_population(population_size, mirror_positions, mirror_install_heights, mirror_heights, mirror_widths);

    best_fitness = -Inf;
    best_individual = [];
    for generation = 1:max_generations
        fitness = evaluate_population(population, tower_position, mirror_positions, radii, cone_angle, num_rays);

        selected_population = selection(population, fitness, crossover_rate);

        crossed_population = crossover(selected_population, crossover_rate);

        mutated_population = mutation(crossed_population, mutation_rate);

        population = mutated_population;

        [current_best_fitness, best_index] = max(fitness);
        if current_best_fitness > best_fitness
            best_fitness = current_best_fitness;
            best_individual = population(best_index, :);
        end
    end
end

function [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths, radii] = ...
    generate_heliostat_layout(tower_position, num_circles, min_radius, spacing, a0, a1)
    [points, ~, radii, ~, group] = generate_concentric_circles(tower_position, num_circles, min_radius, spacing);
    mirror_positions = points;
    mirror_install_heights = a1 * (radii(group) - min_radius) + a0;
    mirror_heights = a0 * (radii(group) - min_radius) + a0;
    mirror_widths = a0 * (radii(group) - min_radius) + a0;
end

function population = initialize_population(population_size, mirror_positions, mirror_install_heights, mirror_heights, mirror_widths)
    num_mirrors = size(mirror_positions, 1);
    population = zeros(population_size, 3 * num_mirrors);
    for i = 1:population_size
        for j = 1:num_mirrors
            population(i, 3*j-2) = mirror_install_heights(j) + 0.1 * randn;
            population(i, 3*j-1) = mirror_heights(j) + 0.1 * randn;
            population(i, 3*j) = mirror_widths(j) + 0.1 * randn;
        end
    end
end

function fitness = evaluate_population(population, tower_position, mirror_positions, radii, cone_angle, num_rays)
    num_individuals = size(population, 1);
    fitness = zeros(num_individuals, 1);
    for i = 1:num_individuals
        [mirror_install_heights, mirror_heights, mirror_widths] = extract_params(population(i, :), size(population, 2) / 3);
        [total_energy, ~, avg_efficiency] = compute_optical_efficiency(mirror_positions, ...
            mirror_install_heights, mirror_heights, mirror_widths, tower_position, cone_angle, num_rays);
        fitness(i) = avg_efficiency;
    end
end

function selected_population = selection(population, fitness, crossover_rate)
    num_individuals = size(population, 1);
    selected_population = zeros(size(population));
    for i = 1:num_individuals
        if rand < crossover_rate
            [idx1, idx2] = datasample(1:num_individuals, 2, 'Replace', false);
            selected_population(i, :) = population(idx1, :);
            selected_population(i, :) = crossover_individual(selected_population(i, :), population(idx2, :));
        else
            selected_population(i, :) = population(i, :);
        end
    end
end

function crossed_population = crossover(selected_population, crossover_rate)
    num_individuals = size(selected_population, 1);
    crossed_population = zeros(size(selected_population));
    for i = 1:num_individuals
        if rand < crossover_rate
            [idx1, idx2] = datasample(1:num_individuals, 2, 'Replace', false);
            crossed_population(i, :) = crossover_individual(selected_population(idx1, :), selected_population(idx2, :));
        else
            crossed_population(i, :) = selected_population(i, :);
        end
    end
end

function mutated_population = mutation(population, mutation_rate)
    num_individuals = size(population, 1);
    num_params = size(population, 2);
    mutated_population = population;
    for i = 1:num_individuals
        for j = 1:num_params
            if rand < mutation_rate
                mutated_population(i, j) = mutated_population(i, j) + 0.1 * randn;
            end
        end
    end
end

function [install_heights, heights, widths] = extract_params(individual, num_mirrors)
    install_heights = individual(1:3:end);
    heights = individual(2:3:end);
    widths = individual(3:3:end);
end

function crossed_individual = crossover_individual(individual1, individual2)
    num_params = length(individual1);
    crossover_point = randi(num_params);
    crossed_individual = [individual1(1:crossover_point), individual2(crossover_point+1:end)];
end

function [total_energy, unit_energy, avg_efficiency] = compute_optical_efficiency(mirror_positions, ...
    mirror_install_heights, mirror_heights, mirror_widths, tower_position, cone_angle, num_rays)
    num_mirrors = size(mirror_positions, 1);
    shadow_efficiency = zeros(num_mirrors, 1);
    truncation_efficiency = zeros(num_mirrors, 1);
    cosine_efficiency = zeros(num_mirrors, 1);
    atmospheric_transmittance = zeros(num_mirrors, 1);

    date_time = datetime('now');
    [sun_elevation, sun_azimuth, sun_direction] = compute_sun_position(98.5, 39.4, date_time);

    for i = 1:num_mirrors
        mirror_position = mirror_positions(i, :);
        mirror_install_height = mirror_install_heights(i);
        mirror_height = mirror_heights(i);
        mirror_width = mirror_widths(i);

        mirror_normal = derive_mirror_normal(mirror_position, sun_direction, tower_position);

        [discrete_rays, ray_energies] = discretize_light_cone(sun_direction, cone_angle, num_rays);

        shadow_blocked = zeros(num_rays, 1);
        for j = 1:num_rays
            ray_origin = mirror_position + mirror_install_height * [0, 0, 1];
            ray_direction = discrete_rays(j, :);
            for k = 1:num_mirrors
                if k ~= i
                    mirror_corners = define_mirror_corners(mirror_positions(k, :), mirror_heights(k), mirror_widths(k), mirror_normal);
                    if detect_rectangle_overlap(ray_origin, ray_direction, mirror_corners)
                        shadow_blocked(j) = 1;
                        break;
                    end
                end
            end
        end

        shadow_efficiency(i) = 1 - sum(shadow_blocked) / num_rays;

        truncation_blocked = zeros(num_rays, 1);
        for j = 1:num_rays
            if ~shadow_blocked(j)
                ray_origin = mirror_position + mirror_install_height * [0, 0, 1];
                ray_direction = discrete_rays(j, :);
                if detect_cylinder_overlap(ray_origin, ray_direction, tower_position, 8, 3.5)
                    truncation_blocked(j) = 1;
                end
            end
        end

        truncation_efficiency(i) = sum(truncation_blocked) / num_rays;

        cosine_efficiency(i) = dot(mirror_normal, sun_direction);

        d_HR = norm(tower_position - mirror_position);
        atmospheric_transmittance(i) = 0.99321 - 0.0001176 * d_HR + 1.97e-8 * d_HR^2;
    end

    total_energy = 0;
    for i = 1:num_mirrors
        optical_efficiency = shadow_efficiency(i) * truncation_efficiency(i) * cosine_efficiency(i) * ...
            atmospheric_transmittance(i) * 0.92;
        total_energy = total_energy + mirror_heights(i) * mirror_widths(i) * optical_efficiency;
    end

    total_area = sum(mirror_heights .* mirror_widths);
    unit_energy = total_energy / total_area;

    avg_efficiency = mean(shadow_efficiency .* truncation_efficiency .* cosine_efficiency .* atmospheric_transmittance);
end

function [sun_elevation, sun_azimuth, sun_direction] = compute_sun_position(lon_deg, lat_deg, date_time)
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

function is_intersect = detect_rectangle_overlap(ray_origin, ray_direction, mirror_corners)
    edge1 = mirror_corners(2, :) - mirror_corners(1, :);
    edge2 = mirror_corners(3, :) - mirror_corners(1, :);
    normal = cross(edge1, edge2);
    denom = dot(normal, ray_direction);
    t = dot(mirror_corners(1,:) - ray_origin, normal) / denom;

    if t < 0
        is_intersect = false;
        return;
    end

    intersection_point = ray_origin + t * ray_direction;

    AB = mirror_corners(2,:) - mirror_corners(1,:);
    AM = intersection_point - mirror_corners(1,:);
    BC = mirror_corners(3,:) - mirror_corners(2,:);
    BM = intersection_point - mirror_corners(2,:);

    if (dot(AB, AM) >= 0) && (dot(AB, AM) <= dot(AB, AB)) && ...
       (dot(BC, BM) >= 0) && (dot(BC, BM) <= dot(BC, BC))
        is_intersect = true;
    else
        is_intersect = false;
    end
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