% CHEME 5440, Prelim 1
% Problem 2

AMP_concen=[0.000; 0.055; 0.093; 0.181; 0.405; 0.990]; % mM 
Rate=[3.003; 6.302; 29.761; 52.002; 60.306; 68.653]; % microM/hr
err=[0.59; 1.20; 5.7; 10.2; 11.8; 13.3]; % 95% confidence for exp data 

Rate_mM_hr=zeros(6,1); 
for i=1:6
Rate_mM_hr(i,1)=Rate(i,1)*0.001; %convert to micromolar/hr to millimolar/hr
err_mM_hr(i,1)=err(i,1)*0.001;   % convert error to mM/hr as well...
end
% y-data is Rate_mM_hr; x-data is AMP_concen 
% concen unit: mM, time unit: hrs

% Given Parameter Values 
kcat=1440; % hr-1 (0.4 s-1) 
E1=0.00012; % mM (0.12 microM), concen of PFK in solution
F6P=0.1; % mM
k_F6P=0.11; % mM
ATP=2.3; % mM
k_ATP=0.42; % mM

% Researched Parameter Values (from Nissler et al., 1983)
K2=0.52; % nM
n2=1.2; % dimensionless

% Compute kinetic limit of PFK 
r = kcat*E1*(F6P/(k_F6P+F6P))*(ATP/(k_ATP+ATP)); % kinetic limit in mM/hr 

% W1 and W2 determined via data-fit: 
a=0.05595; % W1
b=6.705;   % W2

% remember to include error-bars in graph!
Predicted_Rate=zeros(6,1); 
for j=1:6
Predicted_Rate(j,1)=0.0696*(a+(b*(((AMP_concen(j,1)/0.52)^1.2)/(1+((AMP_concen(j,1)/0.52)^1.2))))/(1+a+(b*(((AMP_concen(j,1)/0.52)^1.2)/(1+((AMP_concen(j,1)/0.52)^1.2))))));
end

% Plot experimental data (with error bars) and model data 
figure(1)
set(gca,'DefaultTextFontSize',26)
hold on 
errorbar(AMP_concen,Rate_mM_hr,err_mM_hr,'ko') % plot with error bars
hold on 
plot(AMP_concen,Predicted_Rate,'g*')
legend("Measured Rate", "Predicted Rate") %move legend to bottom right corner
xlabel("3'-5'-AMP Concentration (mM)")
ylabel("Rate (mM/hr)")
hold off 
