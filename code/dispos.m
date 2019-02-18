function mp = dispos(position, spec)
    len = length(position);
    if len == 1
        mp = plot(position(1), 0, spec);
    elseif len == 2
        mp = plot(position(1), position(2), spec);
    elseif len == 3
        mp = plot3(position(1), position(2), position(3), spec);
    else
        matcl = num2cell(position);
        fprintf('- Position (');
        fprintf('%g, ', matcl{1:end});
        fprintf(')\n');
        mp = 0;
    end    
end

