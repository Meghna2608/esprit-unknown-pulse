function [seq,err] = tau_k_gen(min_sep,K)
% K pulse locations with consecutive separations between min_sep and 1.1*min_sep.
rng(ceil(min_sep*K));

err = 1;
max_iter = 10;
iter = 0;

while err == 1 && iter < max_iter

    iter = iter + 1;

    seperation = min_sep + [0 rand(1,K-1)*(min_sep/10)];
    seperation = seperation(randperm(K-1));
    seq = [0 cumsum(seperation)];

    tmpmat = abs(repmat(seq',[1 K]) - repmat(seq,[K 1]));
    tmpmat = min(tmpmat,1-tmpmat) + diag(ones(K,1));

    if min(tmpmat(:)) <= min_sep
        err = 0;
    end
end

if iter == max_iter
    error('tau_k_gen failed')
end

end
