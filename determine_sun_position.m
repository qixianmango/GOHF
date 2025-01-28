function [sun_elevation, sun_azimuth, sunlight_dir] = determine_sun_position(lon_deg, lat_deg, date_time)

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
    sunlight_dir = [x, y, z];

    sun_elevation = rad2deg(sun_elevation_rad);
    sun_azimuth = rad2deg(sun_azimuth_rad);
end