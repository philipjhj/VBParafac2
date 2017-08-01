function [a,b] = hp_update_gamma(a,b,alpha,ln_alpha)
    b=1./b;
    aold = a-1;
    while any((a-aold)/abs(a) > 1e-9)
        aold = a;
        a = a.*exp(-(psi(a)-log(a)+log(alpha)-ln_alpha)./(a.*psi(1,a)-1));
    end
    
    b = a./alpha;
    
    b=1./b;
    assert(all(a>0));
    assert(all(b>0));
end