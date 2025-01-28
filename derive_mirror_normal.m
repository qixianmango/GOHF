function [mirror_normal] = derive_mirror_normal(target_mirror, sunlight_dir, tower_loc)

    reflection_dir = tower_loc - target_mirror;
    assert(reflection_dir(3) * sunlight_dir(3) >= 0);
    sunlight_dir = sunlight_dir / norm(sunlight_dir);
    reflection_dir = reflection_dir / norm(reflection_dir);
    mirror_normal = sunlight_dir + reflection_dir;
    mirror_normal = mirror_normal / norm(mirror_normal);
end