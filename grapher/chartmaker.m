function chartmaker(filename, filename2)

    data = csvread(filename);
    [iters, itmin, itmax] = getiters(data);
    
    diff = (itmax - itmin)./2;
    ave = (itmax + itmin)./2;
    
    data2 = csvread(filename2);
    [iters2, itmin2, itmax2] = getiters(data2);
    
    diff2 = (itmax2 - itmin2)./2;
    ave2 = (itmax2 + itmin2)./2;
    
    errorbar(iters', ave', diff', '.');
    hold on;
    errorbar(iters2', ave2', diff2', '.');
    hold off;
    figure();
    plot(iters', itmin');
    hold on;
    plot(iters2', itmin2');
    hold off;
end

function [iters, itmin, itmax] = getiters(data)
    col1 = data(:,1);
    iters = unique(col1);
    idata = [];
    
    for i=1:length(iters)
        itern = iters(i);
        indic = find(col1 == itern);
        for j=1:length(indic)
            idata(i,j) = data(indic(j),2);
        end
    end
    
    itmin = min(idata,[],2);
    itmax = max(idata,[],2);

end