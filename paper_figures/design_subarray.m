function [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,mode)
% Sub-array design.  Use as  Y_hat = Y_tilde(sel,:), then
% Y_hat(idx1,:) is sub-array 1 and Y_hat(idx2,:) is sub-array 2.
%
% mode 1: maximum overlapping sub-arrays. sel = 1:2M (consecutive grid points),
%         sub-array 2 = sub-array 1 shifted by one grid step.
% mode 2: random doublets on the grid (0:2*M_tilde-1)/(2T): M points chosen at
%         random from the odd indices (l/T), each paired with the next index
%         (shift 1/(2T)). mc_idx sets the random seed.

switch mode

    case 1 % maximum overlapping
        sel1 = 1:(2*M-1);
        sel2 = sel1+1;

    case 2 % random doublets
        rng(M*mc_idx);
        tmpv = randperm(M_tilde);
        sel1 = 2*tmpv(1:M)-1;
        sel2 = sel1 + 1;

end

sel = union(sel1,sel2);
[~,idx1,~] = intersect(sel,sel1);
[~,idx2,~] = intersect(sel,sel2);

end
