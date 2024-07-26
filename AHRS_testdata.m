clc
clear
close all

load('phone_imu_data_1ms_with_acc.mat')
Ts = 0.001;
figure
subplot(3,1,1)
plot(time,Accel')
title('Phone sensor Data')
legend('x_a_c_c_e_l','y_a_c_c_e_l','z_a_c_c_e_l')
ylabel('m/s^2')
subplot(3,1,2)
plot(time,Gyro')
legend('x_g_y_r_o','y_g_y_r_o','z_g_y_r_o')
ylabel('rad/sec')
subplot(3,1,3)
plot(time,Mag')
legend('x_m_a_g','y_m_a_g','z_m_a_g')
ylabel('value')

A = [1, 0, 0, -Ts,   0,  0;
     0, 1, 0,   0, -Ts,  0;
     0, 0, 1,   0,   0,-Ts;
     0, 0, 0,   1,   0,  0;
     0, 0, 0,   0,   1,  0;
     0, 0, 0,   0,   0,  1];
B = [Ts, 0, 0;
     0, Ts, 0;
     0,  0,Ts;
     0,  0, 0;
     0,  0, 0;
     0,  0, 0];
C = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0];

Q = diag([0.5,0.5,0.5,0.05,0.05,0.05]);
R = 100*eye(3);
P = 1e-3*eye(1);
xplus = zeros(6,1);

for k=1:numel(time)                 

    P = A*P*A' + Q;
    
    S = C*P*C' + R; % Estimate error 
    K = (P*C')/S; % Kalman gain
    
    %% meas
    phi_ = -atan(Accel(k,2)/sqrt(Accel(k,1)^2+Accel(k,3)^2));
    theta_ = atan(Accel(k,1)/sqrt(Accel(k,2)^2+Accel(k,3)^2));
    Zt = [phi_;
          theta_;
          atan2(Mag(k,3) * sin(phi_) - Mag(k,2) * cos(phi_), Mag(k,1) * cos(theta_) + Mag(k,2) * sin(theta_) * sin(phi_) + Mag(k,3) * sin(theta_) * cos(phi_))  ];% atan2(-Mag(2),Mag(1))
    meas(k,:,:) = Zt;
    if k==2; xplus = [Zt;0;0;0];
    else
    omega = double(Gyro(k,:)');
    x_ = A*xplus+B*omega;
    innov = Zt - C*x_ ;
    xplus = x_ + K*innov;
    end
    states(k,:,:) = xplus;
    
    % Update the error covariance 
    
    P = (eye(6) - K*C)*P;
   
end
figure
subplot(3,1,1)
plot(time,rad2deg(states(:,1)));hold on
plot(time,rad2deg(meas(:,1)))
title("Roll")
legend('estimation','measurement')
ylabel("Deg")
subplot(3,1,2)
plot(time,rad2deg(states(:,2)));hold on
plot(time,rad2deg(meas(:,2)))
title("Pitch")
legend('estimation','measurement')
ylabel("Deg")
subplot(3,1,3)
plot(time,rad2deg(states(:,3)));hold on
plot(time,rad2deg(meas(:,3)))
title("Yaw")
legend('estimation','measurement')
ylabel("Deg")
