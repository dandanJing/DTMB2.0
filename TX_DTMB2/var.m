function result = var(x)
temp = mean(x);
energy = (x - temp).*conj(x-temp);
result = sum(energy)/length(x);
