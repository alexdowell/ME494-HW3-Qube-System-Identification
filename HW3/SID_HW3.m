load sweep_5v_300s
time = data(1,:)';
volt = data(2,:)';
ang = data(3,:)';
vel = data(4,:)';
acc = deriv(vel,.01);
load qube_data_multistep
time2 = data(1,:)';
volt2 = data(2,:)';
ang2 = data(3,:)';
vel2 = data(4,:)';
acc2 = deriv(vel2,.01);
x = [time, volt, ang, vel];
T_hat =(x'*x)\x'*acc;
Y_hat = x*T_hat;
x2 = [time2, volt2, ang2, vel2];
T_hat2 =(x2'*x2)\x2'*acc2;
Y_hat2 = x2*T_hat2;
disp('and the estimated parameters are:')
disp('frequency sweep')
T_hat
disp('step input')
T_hat2
v = (acc - Y_hat);
s_sq = sum(v.^2)/(length(v) -length(T_hat));
T_hat_confid = 2 * sqrt(s_sq) * diag((x'*x)^-1);
disp('and the confidence intervals are (frequency sweep):')
T_hat_confid
v2 = (acc2 - Y_hat2);
s_sq2 = sum(v2.^2)/(length(v2) -length(T_hat2));
T_hat_confid2 = 2 * sqrt(s_sq2) * diag((x2'*x2)^-1);
disp('and the confidence intervals are (step-input):')
T_hat_confid2
percent_T_hat_confid = (T_hat_confid./T_hat)*100*100;
disp('and the confidence intervals percent are (frequency sweep):')
percent_T_hat_confid
percent_T_hat_confid2 = (T_hat_confid2./T_hat2)*100*100;
disp('and the confidence intervals percent are (step-input):')
percent_T_hat_confid2
param = {'A', 'B', 'C', 'D'};
param = param';
param2 = {'A', 'B', 'C', 'D'};
param2 = param2';
disp('(frequency sweep) and (step-input)')
T = table(param,T_hat,T_hat_confid, percent_T_hat_confid)
T2 = table(param2,T_hat2,T_hat_confid2, percent_T_hat_confid2)
R = acc - Y_hat; % #1
R2 = acc2 - Y_hat2; % #1
disp('model structure for the frequency sweep is:') % #1
disp('Y_hat = .0255*Time + 173.271*Voltage -.0243*Angle - 5.1593*Velocity') % #1
disp('model structure for the step input is:') % #1
disp('Y_hat = -297.625*Time + 406.9984*Voltage - 3.0379*Angle - 10.1334*Velocity') % #1
% plotting %
figure(1) % #1
plot(time, Y_hat, time, acc, time, Y_hat + 2 * sqrt(s_sq),'--r', time, Y_hat - 2 * sqrt(s_sq), '--r') % #1
title('measured and modeled angular acceleration (frequency sweep) vs. time')
xlabel('time(s)')
ylabel('angular acceleration(rad/s^2)')
legend('measured','modeled')
figure(2) % #1
plot(time2, Y_hat2, time2, acc2, time2, Y_hat2 + 2 * sqrt(s_sq2),'--r', time2, Y_hat2 - 2 * sqrt(s_sq2),'--r')
% #1
title('measured and modeled angular acceleration (step-input) vs. time')
xlabel('time(s)')
ylabel('angular acceleration(rad/s^2)')
legend('measured','modeled')
figure(3) % #2
plot(time,R,'.') % #1
title('(frequency sweep)residual vs. time')
xlabel('time(s)')
ylabel('residual')
figure(4) % #2
plot(Y_hat, R, '.') % #2
title('(frequency sweep)residual vs. Yhat')
xlabel('Yhat(rad/s^2)')
ylabel('residual')
figure(5) % #2
plot(time2, R2,'.') % #2
title('(step-input)residual vs. time')
xlabel('time(s)')
ylabel('residual')
figure(6) % #2
plot(Y_hat2, R2, '.') % #2
title('(step-input)residual vs. Yhat')
xlabel('Yhat(rad/s^2)')
ylabel('residual')
disp('The frequency sweep R squared value is:')
v = acc - Y_hat;
R_sq = (T_hat'*x'*acc - length(acc)*mean(acc)^2) / (acc'*acc - length(acc) * mean(acc)^2) % #1
disp('The step input R squared value is:')
v2 = acc2 - Y_hat2;
R_sq2 = (T_hat2'*x2'*acc2 - length(acc2)*mean(acc2)^2) / (acc2'*acc2 - length(acc2) * mean(acc2)^2) % #1
sjj_time = 0;
sjj_volt = 0;
sjj_ang = 0;
sjj_vel = 0;
sjj_acc = 0;
qj_time = [];
qj_volt = [];
qj_ang = [];
qj_vel = [];
qj_acc = [];
for i = 1:length(time)
 sjj_time = sjj_time + (time(i) - mean(time))^2;
 sjj_volt = sjj_volt + (volt(i) - mean(volt))^2;
 sjj_ang = sjj_ang + (ang(i) - mean(ang))^2;
 sjj_vel = sjj_vel + (vel(i) - mean(vel))^2;
 sjj_acc = sjj_acc + (acc(i) - mean(acc))^2;
end
qj_time = (time - mean(time))/sqrt(sjj_time);
qj_volt = (volt - mean(volt))/sqrt(sjj_volt);
qj_ang = (ang - mean(ang))/sqrt(sjj_ang);
qj_vel = (vel - mean(vel))/sqrt(sjj_vel);
qj_acc = (acc - mean(acc))/sqrt(sjj_acc);
xr = [qj_time, qj_volt, qj_ang, qj_vel];
T_hatnorm =(xr'*xr)\xr'*qj_acc % #3
sjj_time2 = 0;
sjj_volt2 = 0;
sjj_ang2 = 0;
sjj_vel2 = 0;
sjj_acc2 = 0;
qj_time2 = [];
qj_volt2 = [];
qj_ang2 = [];
qj_vel2 = [];
qj_acc2 = [];
for i = 1:length(time2)
 sjj_time2 = sjj_time2 + (time2(i) - mean(time2))^2;
 sjj_volt2 = sjj_volt2 + (volt2(i) - mean(volt2))^2;
 sjj_ang2 = sjj_ang2 + (ang2(i) - mean(ang2))^2;
 sjj_vel2 = sjj_vel2 + (vel2(i) - mean(vel2))^2;
 sjj_acc2 = sjj_acc2 + (acc2(i) - mean(acc2))^2;
end
qj_time2 = (time2 - mean(time2))/sqrt(sjj_time2);
qj_volt2 = (volt2 - mean(volt2))/sqrt(sjj_volt2);
qj_ang2 = (ang2 - mean(ang2))/sqrt(sjj_ang2);
qj_vel2 = (vel2 - mean(vel2))/sqrt(sjj_vel2);
qj_acc2 = (acc2 - mean(acc2))/sqrt(sjj_acc2);
xr2 = [qj_time2, qj_volt2, qj_ang2, qj_vel2];
T_hatnorm2 = (xr2'*xr2)\xr2'*qj_acc2 % #3
corr = xr'*xr % #4
VIF = diag(corr^-1);
r2 = 1-1./VIF % #4
d2 = svd(corr).^2;
cond_ind = (max(d2)./d2) % #4
corr2 = xr2'*xr2 % #4
VIF2 = diag(corr2^-1);
r22 = 1-1./VIF2 % #4
d22 = svd(corr2).^2;
cond_ind2 = (max(d22)./d22) % #4
Y_hatval = x2*T_hat; % validating frequency sweep model using step input data %
Y_hat2val = x*T_hat2; % validating step input model using frequency sweep data %
figure(7) % #5
plot(time2, acc2, time2, Y_hatval) % #5
title('measured and modeled angular acceleration (frequency sweep) vs. time')
xlabel('time(s)')
ylabel('angular acceleration(rad/s^2)')
legend('measured','modeled')
figure(8) % #5
plot(time, acc, time, Y_hat2val) % #5
title('measured and modeled angular acceleration (step-input) vs. time')
xlabel('time(s)')
ylabel('angular acceleration(rad/s^2)')
legend('measured','modeled')
disp('The frequency sweep model validation using step input data R squared value is:')
v = acc - Y_hat;
R_sq_val = (T_hat'*x2'*acc2 - length(acc2)*mean(acc2)^2) / (acc2'*acc2 - length(acc2) * mean(acc2)^2)
% #5
disp('The step input model validation using frequency sweep data R squared value is:')
v2 = acc2 - Y_hat2;
R_sq2_val = (T_hat2'*x'*acc - length(acc)*mean(acc)^2) / (acc'*acc - length(acc) * mean(acc)^2) % #5
%6 m1cmd, m2cmd, roll_rate, and the bias
load quadroll
m1cmd = quad.m1cmd;
m2cmd = quad.m2cmd;
roll_rate = quad.roll_rate;
bias = ones(6999,1);
roll_acc = deriv(quad.roll_rate,.01);
sjj_m1cmd = 0;
sjj_m2cmd = 0;
sjj_roll_rate = 0;
sjj_bias = 0;
sjj_roll_acc = 0;
qj_m1cmd = [];
qj_m2cmd = [];
qj_roll_rate = [];
qj_bias = [];
qj_roll_acc = [];
for i = 1:length(m2cmd)
 sjj_m1cmd = sjj_m1cmd + (m1cmd(i) - mean(m1cmd))^2;
 sjj_m2cmd = sjj_m2cmd + (m2cmd(i) - mean(m2cmd))^2;
 sjj_roll_rate = sjj_roll_rate + (roll_rate(i) - mean(roll_rate))^2;
 sjj_bias = sjj_bias + (bias(i) - mean(bias))^2;
 sjj_roll_acc = sjj_roll_acc + (roll_acc(i) - mean(roll_acc))^2;
end
qj_m1cmd = (m1cmd - mean(m1cmd))/sqrt(sjj_m1cmd);
qj_m2cmd = (m2cmd - mean(m2cmd))/sqrt(sjj_m2cmd);
qj_roll_rate = (roll_rate - mean(roll_rate))/sqrt(sjj_roll_rate);
qj_bias = (bias - mean(bias))/sqrt(sjj_bias);
qj_roll_acc = (roll_acc - mean(roll_acc))/sqrt(sjj_roll_acc);
xr3 = [qj_m1cmd, qj_m2cmd, qj_roll_rate]; %, qj_bias];
T_hatnorm3 = (xr3'*xr3)\xr3'*qj_roll_acc;
corr3 = xr3'*xr3 % #6
%corr3(isnan(corr3))=0;
VIF3 = diag(corr3^-1);
r23 = 1-1./VIF3 % #6
d23 = svd(corr3).^2;
cond_ind3 = (max(d23)./d23) % #6
% Cross Correlation
[Xa,Ya,D] = alignsignals(volt,acc);
dt = 0.01;
adj_time = time + D*dt;
figure(9);
hold on;
plot(time,volt,'DisplayName','Output');
ylabel('Voltage (V)');
yyaxis right
plot(adj_time,volt,'DisplayName','Input');
ylabel('Voltage (V)');
xlabel('Time (sec)');
xlim([0,time(length(time))]);
title('Frequency Sweep Input-Output Angular Velocity Lag vs Time');
legend('location','best');
[Xa1,Ya1,D1] = alignsignals(vel,acc);
D1
dt = 0.01;
adj_time1 = time + D1*dt;
figure(10);
hold on;
plot(time,vel,'DisplayName','Output');
ylabel('angular velocity Command Input (rad/s)');
yyaxis right
plot(adj_time1,vel,'DisplayName','Input');
ylabel('angular velocity Command Input (rad/s)');
xlabel('Time (s)');
xlim([0,time(length(time))]);
title('Frequency Sweep Input-Output angular velocity Lag vs Time');
legend('location','best');
[Xa,Ya,D] = alignsignals(volt2,acc2);
dt = 0.01;
adj_time = time2 + D*dt;
figure(11);
hold on;
plot(time2,volt2,'DisplayName','Output');
ylabel('Voltage (V)');
yyaxis right
plot(adj_time,volt2,'DisplayName','Input');
ylabel('Voltage (V)');
xlabel('Time (sec)');
xlim([0,time2(length(time2))]);
title('Step input Input-Output Angular Velocity Lag vs Time');
legend('location','best');
[Xa1,Ya1,D1] = alignsignals(vel2,acc2);
dt = 0.01;
adj_time1 = time2 + D1*dt;
D1
figure(12);
hold on;
plot(time2,vel2,'DisplayName','Output');
ylabel('angular velocity Command Input (rad/s)');
yyaxis right
plot(adj_time1,vel2,'DisplayName','Input');
ylabel('angular velocity Command Input (rad/s)');
xlabel('Time (s)');
xlim([0,time2(length(time2))]);
title('step input Input-Output angular velocity Lag vs Time');
legend('location','best');