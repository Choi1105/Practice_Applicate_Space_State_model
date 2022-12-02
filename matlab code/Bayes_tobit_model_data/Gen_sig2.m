function sig2 = Gen_sig2(Y, X, v_0, d_0, beta)

   % Parameters
    v_1 = v_0 + rows(Y);
    d_1 = d_0 + (Y - X*beta)'*(Y - X*beta);
    
    %Sampling and Save
    sig2 = randig(v_1/2,d_1/2,1,1);
    
end