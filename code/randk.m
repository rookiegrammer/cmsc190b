function random = randk(limit, pop, k)
  random = zeros(1,pop);
  for i=1:pop
    newnum = ceil(rand()^k*limit);
    while ismember(newnum, random)
        newnum = ceil(rand()^k*limit);
    end
    random(1,i)= newnum;
  end
end
