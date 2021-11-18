using Printf

function newton(f, df, p0, n_max, rel_tol; verbose = true)
    
    converged = false;
    p = p0;
    p_old = p0;

    for i in 1:n_max

        p = p_old - f(p_old)/df(p_old);
        
        if verbose
            @printf(" %d: p = %.15g, f(p) = %g\n", i, p, f(p));
        end

        
        if (i>1)
            if abs(p-p_old)/abs(p)< rel_tol
                converged = true;
                break
            end
        end

        p_old = p;

    end
    
    if !converged
        @printf("ERROR: Did not converge after %d iterations\n", n_max);
    end

    return p
    
end

function bisection(f, a, b, n_max, rel_tol; verbose = true)
    
    converged = false;
    p_old = 0;
    p = 0;
    for i in 1:n_max

        p = 0.5 * (a+b)
        
        if verbose
            @printf(" %d: a = %g, b = %g, p = %g, f(p) = %g\n", i, a, b, p, f(p));
        end

        if ( f(a) * f(p)<=0)
            b = p;
        else
            a = p
        end
        
        if (i>1)
            if abs(p-p_old)/abs(p)< rel_tol
                converged = true;
                break
            end
        end

        if(abs(f(p))==0)
            converged = true;
            break
        end
        p_old = p;

    end
    
    if !converged
        @printf("ERROR: Did not converge after %d iterations\n", n_max);
    end

    return p
    
end

function secant(f, p0, p1, n_max, rel_tol; verbose = true)
    
    converged = false;
    
    p = p0;
    for i in 1:n_max

        p = p1 - f(p1) * (p1-p0)/(f(p1)-f(p0));
        
        if verbose
            @printf(" %d: p = %.12g, f(p) = %g\n", i, p, f(p));
        end

        
        if (i>1)
            if abs(p-p1)/abs(p1)< rel_tol
                converged = true;
                break
            end
        end
        p0 = p1;
        p1 = p;

    end
    
    if !converged
        @printf("ERROR: Did not converge after %d iterations\n", n_max);
    end

    return p
    
end