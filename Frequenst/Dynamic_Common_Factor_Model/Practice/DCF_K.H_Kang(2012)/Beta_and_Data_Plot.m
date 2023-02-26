% Figure
% Index 
i = 1:rows(Beta_ttm); 

% Filtered value
tiledlayout(3,3)
for z = 1:3
nexttile
plot(i, Beta_ttm(:,z) ,'k', i, Beta_LB(:,z), 'b:', i, Beta_UB(:,z),'r:','LineWidth',1.5)
legend('Common Factor', 'Low Band', 'High Band');
title('Filtered Common Factor and Confidence Interval');
end

% Smoothed value
for z = 1:3
nexttile
plot(i, Beta_tTm(:,z) ,'k', i, Beta_LB_SM(:,z), 'b:', i, Beta_UB_SM(:,z),'r:','LineWidth',1.5)
legend('Common Factor', 'Low Band', 'High Band');
title('Smoothed Common Factor and Confidence Interval');
end

% Data Value
nexttile
plot(data(:,10), 'k', LineWidth=1.5);
legend('M120');
title('M120 Data Plot');
nexttile
plot(data(:,10)-data(:,1), 'k', LineWidth=1.5);
legend('M120-M3');
title('(M120-M3) Data Plot');
nexttile
plot(((data(:,6)*2)-(data(:,1) + data(:,10))), 'k', LineWidth=1.5);
legend('(M24*2)-(M3+M120)');
title('(M24*2)-(M3+M120) Data Plot');