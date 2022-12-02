function [residual, g1, g2] = static_resid_g1_g2(T, y, x, params, T_flag)
% function [residual, g1, g2] = static_resid_g1_g2(T, y, x, params, T_flag)
%
% Wrapper function automatically created by Dynare
%

    if T_flag
        T = FinalExam.static_g2_tt(T, y, x, params);
    end
    [residual, g1] = FinalExam.static_resid_g1(T, y, x, params, false);
    g2       = FinalExam.static_g2(T, y, x, params, false);

end
