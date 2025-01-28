function optimize_heliostat_field()
    tower_position = [0, 0, 80];
    min_radius = 100;
    max_radius = 350;
    num_circles = 75;
    spacing = 5.75 + 5;
    num_mirrors = 1745;
    a0 = 5.75;
    a1 = 2.875;
    population_size = 50;
    max_generations = 100;
    mutation_rate = 0.01;
    crossover_rate = 0.8;

    [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths, radii] = ...
        initialize_heliostat_field(tower_position, num_circles, min_radius, spacing, a0, a1);

    population = initialize_population(population_size, mirror_positions, mirror_install_heights, mirror_heights, mirror_widths);

    for generation = 1:max_generations
        fitness = evaluate_population(population, tower_position, radii);

        selected_population = selection(population, fitness, crossover_rate);

        crossed_population = crossover(selected_population, crossover_rate);

        mutated_population = mutation(crossed_population, mutation_rate);

        population = mutated_population;

        [best_fitness, best_index] = max(fitness);
        best_individual = population(best_index, :);
        
        tower_position = update_tower_position(best_individual, tower_position);
        [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths] = ...
            update_heliostat_params(best_individual, mirror_positions, mirror_install_heights, mirror_heights, mirror_widths);
    end

end

function [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths, radii] = ...
    initialize_heliostat_field(tower_position, num_circles, min_radius, spacing, a0, a1)
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

function fitness = evaluate_population(population, tower_position, radii)
    num_individuals = size(population, 1);
    fitness = zeros(num_individuals, 1);
    for i = 1:num_individuals
        [mirror_install_heights, mirror_heights, mirror_widths] = ...
            extract_params(population(i, :), size(population, 2) / 3);
        [valid, ~, ~, ~] = validate_configuration(mirror_install_heights, mirror_heights, mirror_widths, radii);
        if valid
            [total_energy, ~, avg_efficiency] = evaluate_performance(mirror_install_heights, mirror_heights, mirror_widths, ...
                size(population, 2) / 3, min(radii), radii, tower_position);
            fitness(i) = avg_efficiency;
        else
            fitness(i) = 0;
        end
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

function tower_position = update_tower_position(best_individual, tower_position)
    tower_position(1) = tower_position(1) + 0.1 * randn;
    tower_position(2) = tower_position(2) + 0.1 * randn;
end

function [mirror_positions, mirror_install_heights, mirror_heights, mirror_widths] = ...
    update_heliostat_params(best_individual, mirror_positions, mirror_install_heights, mirror_heights, mirror_widths)
    num_mirrors = size(mirror_positions, 1);
    for i = 1:num_mirrors
        mirror_install_heights(i) = best_individual(3*i-2);
        mirror_heights(i) = best_individual(3*i-1);
        mirror_widths(i) = best_individual(3*i);
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