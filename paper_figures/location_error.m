function error = location_error(true_loc,est_loc,T)
% Matching distance: min over permutations of max_k |tau_k - tau_hat_pi(k)| on the torus [0,T).

tmpmat = perms(est_loc);
tmperr = abs(repmat(true_loc,size(tmpmat,1),1) - tmpmat);
tmperr = min(tmperr, T-tmperr);
maxerr = max(tmperr,[],2);
error = min(maxerr);

end
