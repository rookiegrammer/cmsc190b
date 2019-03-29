function chartmaker2(filename, filename2)
    
    data1 = csvread(filename);
    data2 = csvread(filename2);
    
    iters = 0:(size(data1,1)-1);
    
    itmin1 = min(data1,[],2);
    itmax1 = max(data1,[],2);
    
    itmin2 = min(data2,[],2);
    itmax2 = max(data2,[],2);
    
    itdmin = max(itmin1, itmin2);
    itdmax = min(itmax1, itmax2);
    
    diff1 = (itmax1 - itmin1)./2;
    ave1 = (itmax1 + itmin1)./2;
    
    diff2 = (itmax2 - itmin2)./2;
    ave2 = (itmax2 + itmin2)./2;
    
    diffd = (itdmax - itdmin)./2;
    aved = (itdmax + itdmin)./2;
    
    figure();
    
    errorbar(iters', ave1', diff1', '.', 'CapSize',5,'MarkerSize', 0.001,'LineWidth', 1.3);
    hold on;
    errorbar(iters', ave2', diff2', '.', 'CapSize',5,'MarkerSize', 0.001,'LineWidth', 1.3);
    errorbar(iters', aved', diffd', 'k.', 'CapSize',0,'MarkerSize', 0.001,'LineWidth', 1.3);
    legend('LOA', 'iLOA', 'Overlap');
    
    print(['out/' filename '-bars.png'], '-dpng');
    hold off;
    
    figure();
    plot(iters', itmin1','LineWidth', 1.3);
    hold on;
    plot(iters', itmin2','LineWidth', 1.3);
    legend('LOA', 'iLOA');
    
    print(['out/' filename '-fit.png'], '-dpng');
    hold off;
    
end