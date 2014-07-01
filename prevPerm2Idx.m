function prevIdx = prevPerm2Idx(prevPerm, nV)
% recover indices of cells in prevalence permutation given a uint32 value
% representing the ordering 
% see getPrevPerm
    if(~exist('nV', 'var'))
        nV=ceil(log10(double(prevPerm(1))));
    end
    if (nV > 9)
        err = MException('prevPerm2Idx:badEncoding', ...
                         'Encoding fails for more than 9 genotypes');
        throw(err);
    end

    prevIdx = (10*ones(1,nV)) .^ (0:nV-1);
    prevIdx = mod(cell2mat(arrayfun(@(z) idivide(uint32(prevPerm), z), prevIdx, ...
                            'uniformoutput', false)), 10);
end
