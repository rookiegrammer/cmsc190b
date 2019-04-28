function [params, ind, mindata] = paramfind(plength, filename)
    % default 11
    data = csvread(filename);
    dlength = size(data, 2);
    fitdata = data(:, (plength+1):dlength);
    pardata = data(:, 1:plength);
    
    avedata = mean(fitdata, 2);
    
    mindata = min(avedata);
    ind = find(avedata==mindata);
    params = pardata(ind,:);
end

