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
        discrete_rays(i, :) = discrete_rays(i, :) / norm(discrete_rays(i, :)); % 归一化
    end

    sigma = cone_angle / 3;
    ray_energies = exp(-0.5 * (theta.^2 + phi.^2) / sigma^2);
    ray_energies = ray_energies / sum(ray_energies); 
end