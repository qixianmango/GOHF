function [results] = compute_mirror_performance(date, mirror_pos, mirror_heights, mirror_widths, tower_loc)

    num_mirrors = size(mirror_pos, 1);
    [sun_angle, sun_dir, sunlight_vec] = determine_sun_position(98.5, 39.4, date);
    sunlight_vec = real(sunlight_vec);
    mirror_normals = zeros(num_mirrors, 3);
    mirror_corners = zeros(4, 3, num_mirrors);

    for i = 1:num_mirrors
        mirror_normals(i, :) = derive_mirror_normal(mirror_pos(i, :), sunlight_vec, tower_loc);
        [mirror_corners(1, :, i), mirror_corners(2, :, i), mirror_corners(3, :, i), mirror_corners(4, :, i)] = ...
            define_mirror_corners(mirror_pos(i, :), mirror_heights(i), mirror_widths(i), mirror_normals(i, :));
    end

    shadow_eff = zeros(num_mirrors, 1);
    trunc_eff = zeros(num_mirrors, 1);
    cosine_eff = zeros(num_mirrors, 1);
    atm_trans = zeros(num_mirrors, 1);

    for i = 1:num_mirrors
        [shadow_eff(i), trunc_eff(i), cosine_eff(i), atm_trans(i)] = ...
            evaluate_efficiency(mirror_pos, i, sunlight_vec, mirror_normals(i, :), tower_loc, mirror_corners);
    end

    mirror_reflectivity = 0.92;
    solar_constant = 1.366;
    altitude = 3;
    a = 0.4237 - 0.00821 * (6 - altitude)^2;
    b = 0.5055 + 0.00595 * (6.5 - altitude)^2;
    c = 0.2711 + 0.01858 * (2.5 - altitude)^2;
    dni = solar_constant * (a + b * exp(-c / sind(sun_angle)));

    results = struct('optical_efficiency', [], 'truncation_efficiency', [], 'output_power', []);

    for i = 1:num_mirrors
        optical_efficiency = shadow_eff(i) * trunc_eff(i) * cosine_eff(i) * atm_trans(i) * mirror_reflectivity;
        output_power = dni * mirror_heights(i) * mirror_widths(i) * optical_efficiency;
        results(i).optical_efficiency = optical_efficiency;
        results(i).truncation_efficiency = trunc_eff(i);
        results(i).output_power = output_power;
    end
end