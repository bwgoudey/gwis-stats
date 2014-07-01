function prevPerm = getPrevPerm(ctControl, ctCase)
% generates uint32 array of values that define prevalence permutation order
% for a probe pair
% ie if prevPerm = 786321459, then genotype 9 in the contingency table for
% the pair has the highest prevalence followed by genotype 5 then 4, 1, 2 and
% so on.
    if (~isequal(size(ctControl), size(ctCase)))
        err = MException('getPrevPerm:wrongSize', ...
                 'control and case contingency tables must be the same size');
        throw(err);
    end

    minGenoSum = 0.00001; % avoid division by 0

    % find prevalence
    %prev = ctCase ./ max(minGenoSum, ctControl + ctCase);
    %bwgoudey - Use relative risk the same as Adam implemented in gaimBM    
    eps=1e-8;
    t0=sum(ctControl(1, :), 2);
    t1=sum(ctCase(1, :), 2);
    p0=max(eps, ctControl/t0);
    p1=max(eps, ctCase/t1);
    prev  = p1./p0;
    
    [a, prev] = sort(prev, 2, 'descend'); 

    % convert into a single number
    [nP, nV] = size(ctControl);
    if (nV > 9)
        err = MException('getPrevPerm:badEncoding', ...
                'Encoding fails for more than 9 genotypes');
        throw(err);
    end
    prevPerm = repmat((10*ones(1,nV)).^(0:nV - 1), nP, 1);
    prevPerm = uint32(sum(prevPerm .* prev, 2));
end
