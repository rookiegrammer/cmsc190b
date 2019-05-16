clc
format shortE

dinfo = dir('conv-test-*.csv');
for K = 1 : length(dinfo)
  thisfilename = dinfo(K).name  %just the name
  matrx = csvread(thisfilename);
  succ = find(matrx(:,5)==1);
  rati = length(succ)*100/length(matrx)
  mm = mean(matrx(succ,:))
  
end
close all