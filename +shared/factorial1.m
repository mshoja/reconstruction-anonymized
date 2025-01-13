function res = factorial1(n)
    %speed up small factorials
    f = [1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800 87178291200 1307674368000 20922789888000 355687428096000];
    if all(n < length(f))
        res = f(n);
    else
        res = factorial(n);
    end
end
