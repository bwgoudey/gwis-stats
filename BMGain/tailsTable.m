function table = tailsTable(method, LowUp, p, xIn, t, dp, dx)
% generates tables of binomial cdf tails
% CALL: table = tailsTable(method,lowUp,p,x, t, dp,dx)
%   method = 
%     ='binomDirect' - generates exhuastively the whole binomial distribution
%                      for all values for 0:x for logLow 
%                      and for all values x:t for logUp 
%   LowUp - mode for tail generation, 
%         = 'Low' lower tail (sum from 0 up to limitUp > 0)
%         = 'Up' upper tail (sum from t down to limitLo < t)
%         = 'LowUp' both
%   p   - vector of slected proportions, if =[], then 0:dp:1 will be used  
%   xIn - vector of selected counts, if = [], then 0:dx:t will be used
%   t   - integer, the total number of samples
%   dp  - used if p = [], can be empty otherwise
%   dx  - used if x = [], can be empty otherwise
% binomBeta:
%           biLow(p) = t {t-1 \choose x}\int_0^{1-p} s^{t-x-1}(1-s)^x ds
%         = t {t-1 \choose x} F^{-t}\int_0^{1-p} r^{t-x-1}(F-r)^x dr
% table = 
%     logLow: npVals * nxVals 
%      logUp: []
%      xVals: 100
%      pVals: [11x1 double]
 

    if (nargin ~= 7)
        err = MException('tailsTable:wrongNumInputs', ...
                         'Wrong number of input arguments.');
        throw(err);
    end

    if (isempty(p))
        if (isempty(dp))
            err = MException('tailsTable:emptyProb', ...
                             'One of p or dp must be non-empty');
            throw(err);
        elseif(~isscalar(dp))
            err = MException('tailsTable:nonScalarProb', ...
                             'dp should be a positive scalar');
            throw(err);
        elseif (dp <= 0)
            err = MException('tailsTable:negativeProb', ...
                             'dp should be positive');
            throw(err);
        end
    else
        if (~(size(p, 2)==1 && size(p,1)>0))%~iscolumn(p)
            p = p';
        end
    end

    if (isempty(xIn))
        if (isempty(dx))
            err = MException('tailsTable:emptyCount', ...
                             'One of xIn or dx must be non-empty');
            throw(err);
        elseif(~isscalar(dx))
            err = MException('tailsTable:nonScalarCount', ...
                             'dx should be a positive scalar');
            throw(err);
        elseif (dx <= 0)
            err = MException('tailsTable:negativeCount', ...
                             'dx should be positive');
            throw(err);
        end
    else
        if (~isrow(xIn))
            xIn = xIn';
        end
    end

    table.logLow = [];
    table.logUp = [];
    negINF = -Inf;
    switch method
        case 'binomDirect'
            if isempty(p)
                % set basic vector of values of r, excluding extremes,
                % i.e. 0 and F
                p = (0: dp: 1)';
            end
            np = length(p);
            logp = log(p);
            logq = log(1-p);
            
            logp(logp == -Inf) = negINF;
            logq(logq == -Inf) = negINF;
            
            xVals = 0:t; % maximal range of values of x to be used
            switch LowUp      
                case 'Low'
                    if (~isempty(xIn))
                        xVals = 0:max(xIn);
                    end
                case 'Up'
                    if (~isempty(xIn))
                        xVals = min(xIn):t;
                    end
                case 'LowUp'
                    % do nothing here - already done
                otherwise
                    err = MException('tailsTable:unknownMode', ...
                         ['Mode ' LowUp ' unrecognised.']);
                    throw(err);
            end

            nx = length(xVals);

            if ismember(LowUp, {'Low' 'LowUp'})
                table.logLow  = NaN * zeros(np, nx);
                table.logLow(:, 1) = t * logq;
                for  ix = 2 : nx
                    x = xVals(ix);
                    table.logLow(:, ix) = logSumLog(table.logLow (:, ix-1), ...
                                    log_nCk(t, x) + x * logp + (t-x) * logq);
                end
                % recalculate last column in case logq = -inf
                % because -inf * 0 == NaN
                if (xVals(nx) == t)
                    x = t;
                    table.logLow(:, nx) = logSumLog(table.logLow (:, nx-1), ...
                                    log_nCk(t, x) + x * logp);
                end
            end

            if ismember(LowUp, {'Up' 'LowUp'})
                table.logUp = NaN * zeros(np,nx);
                table.logUp(:, nx) = t * logp;
                for  ix = nx-1 : -1 :  1
                    x = xVals(ix);
                    table.logUp(:, ix) = logSumLog(table.logUp(:, ix+1), ...
                                    log_nCk(t, x) + x * logp + (t-x) * logq);
                end
                % recalculate first column in case logp = -inf
                % because -inf * 0 == NaN
                if (xVals(1) == 0)
                    x = 0;
                    table.logUp(:, ix) = logSumLog(table.logUp(:, ix+1), ...
                                    log_nCk(t, x) + (t-x) * logq);
                end
            end
            
            [a, pos0] = ismember(0, p);
            pos0 = pos0(pos0>0);
            [a, pos1] = ismember(1, p);
            pos1 = pos1(pos1>0);
           
            if ~isempty(xIn)
                [a, loc] = ismember(xIn, xVals);
                xVals = xIn;
            else
                [a, loc] = ismember(0:dx:t, xVals);
                xVals = 0:dx:t;
            end

            if ismember(LowUp, {'Low' 'LowUp'})
                table.logLow  = table.logLow(:, loc);
                table.logLow(pos1, xVals < t) = negINF;
                table.logLow(pos1, xVals == t) = 0;
            end

            if ismember(LowUp, {'Up' 'LowUp'})
                table.logUp  = table.logUp(:, loc);
                table.logUp(pos0, xVals > 0) = negINF;
                table.logUp(pos0, xVals == 0) = 0;
            end
           
            table.xVals = xVals;
            table.pVals = p;
        case 'normalApprox'
            if isempty(p)
                % set basic vector of values of r, excluding extremes,
                % i.e. 0 and F
                p = (0: dp: 1)';
            end
            np = length(p);
            logp = log(p);
            logq = log(1-p);
            logt = log(t);
            calcNorm = @(x) -0.5*(logt + logp + logq + log(2*pi)) - ((x-t*p).^2)./(2*t*p.*(1-p));

            logp(logp == -Inf) = negINF;
            logq(logq == -Inf) = negINF;
            
            xVals = 0:t; % maximal range of values of x to be used
            switch LowUp      
                case 'Low'
                    if (~isempty(xIn))
                        xVals = 0:max(xIn);
                    end
                case 'Up'
                    if (~isempty(xIn))
                        xVals = min(xIn):t;
                    end
                case 'LowUp'
                    % do nothing here - already done
                otherwise
                    err = MException('tailsTable:unknownMode', ...
                         ['Mode ' LowUp ' unrecognised.']);
                    throw(err);
            end

            nx = length(xVals);
            %%
            if ismember(LowUp, {'Low' 'LowUp'})
                table.logLow  = NaN * zeros(np, nx);
                table.logLow(:, 1) = calcNorm(xVals(1));
                for  ix = 2 : nx
                    x = xVals(ix);
                    table.logLow(:, ix) = logSumLog(table.logLow (:, ix-1), ...
                                    calcNorm(x));
                end
                % recalculate last column in case logq = -inf
                % because -inf * 0 == NaN
%                 if (xVals(nx) == t)
%                     x = t;
%                     table.logLow(:, nx) = logSumLog(table.logLow (:, nx-1), ...
%                                     log_nCk(t, x) + x * logp);
%                 end
            end
            %%
            if ismember(LowUp, {'Up' 'LowUp'})
                table.logUp = NaN * zeros(np,nx);
                table.logUp(:, nx) = calcNorm(xVals(nx));
                for  ix = nx-1 : -1 :  1
                    x = xVals(ix);
                    table.logUp(:, ix) = logSumLog(table.logUp(:, ix+1), ...
                                    calcNorm(x));
                end
                % recalculate first column in case logp = -inf
                % because -inf * 0 == NaN
%                 if (xVals(1) == 0)
%                     x = 0;
%                     table.logUp(:, ix) = logSumLog(table.logUp(:, ix+1), ...
%                                     log_nCk(t, x) + (t-x) * logq);
%                 end
            end
            
            [a, pos0] = ismember(0, p);
            pos0 = pos0(pos0>0);
            [a, pos1] = ismember(1, p);
            pos1 = pos1(pos1>0);
           
            if ~isempty(xIn)
                [a, loc] = ismember(xIn, xVals);
                xVals = xIn;
            else
                [a, loc] = ismember(0:dx:t, xVals);
                xVals = 0:dx:t;
            end

            if ismember(LowUp, {'Low' 'LowUp'})
                table.logLow  = table.logLow(:, loc);
                table.logLow(pos1, xVals < t) = negINF;
                table.logLow(pos1, xVals == t) = 0;
            end

            if ismember(LowUp, {'Up' 'LowUp'})
                table.logUp  = table.logUp(:, loc);
                table.logUp(pos0, xVals > 0) = negINF;
                table.logUp(pos0, xVals == 0) = 0;
            end
           
            table.xVals = xVals;
            table.pVals = p;            
        otherwise
            err = MException('tailsTable:badMethod', ['Unkown method ' method]);
            throw(err);
    end
end
