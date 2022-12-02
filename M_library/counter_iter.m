function counter_iter(iter, freq)

if nargin < 2
    freq = 100;
end

    [~, resid] = minresid(iter,freq);
        if resid == 0
             clc
            disp(['Current iteration is ',num2str(iter)]);
        end
        
end