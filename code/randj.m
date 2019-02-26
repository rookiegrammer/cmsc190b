function position = randj(center, space_min, space_max, nearness_pressure)
	rand_arr = rand(length(center),1).*2-1
    rand_min = abs(min(rand_arr, 0)).^nearness_pressure;
    rand_max = max(rand_arr, 0).^nearness_pressure;

    position = center - (center - space_min) .* rand_min + (-center + space_max) .* rand_max;
end
