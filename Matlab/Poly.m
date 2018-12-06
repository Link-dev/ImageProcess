function result=Poly(cos, n)
    
    syms x;
    y = (x^2-1)^n;
    dy = diff(y, x, n);
    
    x = cos;
    
    result = 1.0/2^n * factorial(n) * subs(dy);
    
end