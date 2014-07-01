function p2 = NWboundary(p1,param)
% call: p2 = NWboundary(p1,param)
% computes NW boundary curve
% param contains different settings of the paramters
%   param.H0Type = 'ratio';
%   param.rho0 = 1;
%   
%   param.HType = 'diff';       % for 'difference'
%   param.delta0 = 1;
%   
%   param.HType = 'pieceLin'    % for 'piecewiseLinear';
%   param.p1 = [0 : 0.2:1].^2;
%   param.p2 = [0 : 0.2:1] ;
% written: Adam Kowalczyk, ,14 June, 2012


switch param.H0Type
    case 'ratio'
        p2 = param.rho0*p1;
    case 'diff'
        p2 = p1 + param.delta0;
    case 'pieceLin'
        p2 = NaN*p1;
        n = length(p1);
        m = length(param.p1);
        for i = 1 : n
            j = nnz(param.p1 <= p1(i));  % this means param.p1(j) <= p1(i)  < param.p1(j+1) <= 
            if  j ==m
                p2(i) = param.p1(j);
            else
                p2(i) = param.p2(j+1)*( p1(i) - param.p1(j)) + param.p2(j) * ( param.p1(j+1) - p1(i)) ;
                p2(i) = p2(i)/ (param.p1(j+1) - param.p1(j));
            end
        end
    otherwise
        error
        
end
if p1 < 0 | p1>1
   error( 'values of p1 out f range')
end
if p2 < 0 | p2>1
   error( 'values of p2 out f range')
end