function v = log_nCk(n,k)
% approximation of log(nchoosek) based on Stirling approximation
% log n! = n log n  - n + 0.5 * log(2*pi*n) + 1/12/n

if (n == k || k == 0)
    v = 0;
elseif (k == n-1 || k == 1)
    v = log(n);
else
    v = k * log(n/k) + (n-k) * log(n/(n-k)) - 0.5 * log(2*pi*k*(n-k)/n) ;
    v = v +  (1/n - 1/(n-k)-1/k)/12;
end

