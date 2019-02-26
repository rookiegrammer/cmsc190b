function random = randk(limit, pop, k) {
  random = ceil(rand(pop,1).^k.*limit);
}
