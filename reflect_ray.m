function reflection = reflect_ray(input, mirror_normal)

    reflection = input - 2 * dot(input, mirror_normal) * mirror_normal;
    reflection = reflection / norm(reflection);
end