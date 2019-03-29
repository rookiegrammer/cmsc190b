dinfo = dir('loa-*.csv');
for K = 1 : length(dinfo)
  thisfilename = dinfo(K).name;  %just the name
  chartmaker2(thisfilename, ['i' thisfilename]);
end
close all