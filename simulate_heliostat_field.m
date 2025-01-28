function simulate_heliostat_field()
    tower_position = [0, 0, 80];
    num_circles = 75;
    min_radius = 100;
    spacing = 5.75 + 5;
    a0 = 5.75;
    a1 = 2.875;
    cone_angle = deg2rad(4.65);
    num_rays = 1000;
    date_time = datetime('now');

    [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths, radii] = ...
        generate_heliostat_layout(tower_position, num_circles, min_radius, spacing, a0, a1);

    [sun_elevation, sun_azimuth, sun_direction] = compute_sun_position(98.5, 39.4, date_time);

    [total_energy, unit_energy, avg_efficiency] = compute_optical_efficiency(mirror_positions, ...
        mirror_install_heights, mirror_heights, mirror_widths, tower_position, sun_direction, cone_angle, num_rays);

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

function [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths, radii] = ...
    generate_heliostat_layout(tower_position, num_circles, min_radius, spacing, a0, a1)
    [points, ~, radii, ~, group] = generate_concentric_circles(tower_position, num_circles, min_radius, spacing);
    mirror_positions = points;
    mirror_install_heights = a1 * (radii(group) - min_radius) + a0;
    mirror_heights = a0 * (radii(group) - min_radius) + a0;
    mirror_widths = a0 * (radii(group) - min_radius) + a0;
end

function [discrete_rays, ray_energies] = discretize_light_cone(center_ray, cone_angle, num_rays)
    center_ray = center_ray / norm(center_ray);
    temp_vector = [0, 0, 1];
    if isequal(center_ray, temp_vector) || isequal(center_ray, -temp_vector)
        temp_vector = [1, 0, 0];
    end

    basis1 = cross(center_ray, temp_vector);
    basis1 = basis1 / norm(basis1);
    basis2 = cross(center_ray, basis1);
    basis2 = basis2 / norm(basis2);

    theta = randn(num_rays, 1) * cone_angle / 3;
    phi = randn(num_rays, 1) * cone_angle / 3;

    discrete_rays = zeros(num_rays, 3);
    for i = 1:num_rays
        discrete_rays(i, :) = center_ray + theta(i) * basis1 + phi(i) * basis2;
        discrete_rays(i, :) = discrete_rays(i, :) / norm(discrete_rays(i, :));
    end

    sigma = cone_angle / 3;
    ray_energies = exp(-0.5 * (theta.^2 + phi.^2) / sigma^2);
    ray_energies = ray_energies / sum(ray_energies);
end

function [total_energy, unit_energy, avg_efficiency] = compute_optical_efficiency(mirror_positions, ...
    mirror_install_heights, mirror_heights, mirror_widths, tower_position, sun_direction, cone_angle, num_rays)
    num_mirrors = size(mirror_positions, 1);
    shadow_efficiency = zeros(num_mirrors, 1);
    truncation_efficiency = zeros(num_mirrors, 1);
    cosine_efficiency = zeros(num_mirrors, 1);
    atmospheric_transmittance = zeros(num_mirrors, 1);

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