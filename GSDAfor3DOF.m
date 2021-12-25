function GSDAfor3DOF()
% Mahmoud Muhammad Yahaya's 3DOF  model implementation solve via NLS algo.
% mahmoudpd@gmail.com
clear;
t0 = 0;
tf = 20; % [0, 20s]
g = 0.05; % sampling period
Kmax = (tf - t0)/g; % 

theta0 = [0, pi/3, pi/2];
thetamah = theta0;
for k = 1:Kmax
    tk = k*g;
    rd = [1.5 + 0.4*sin(pi*tk/5), sqrt(3)/2 + 0.4*sin(pi*tk/5 +   pi/3)]'; % ring
    thetak = GSDA(rd, thetamah(end,:)');
    thetamah = [thetamah; thetak']; % a matrix of joint angular vectors
end

x0 = zeros(Kmax+1,3);
x1 = [cos(thetamah(:,1)), sin(thetamah(:,1))];
x2 = [cos(thetamah(:,1)) + cos(thetamah(:,1) + thetamah(:,2)), sin(thetamah(:,1)) + sin(thetamah(:,1) + thetamah(:,2))];
x3 = [cos(thetamah(:,1)) + cos(thetamah(:,1) + thetamah(:,2)) + cos(thetamah(:,1) + thetamah(:,2) + thetamah(:,3)), sin(thetamah(:,1)) + sin(thetamah(:,1) + thetamah(:,2)) + sin(thetamah(:,1) + thetamah(:,2)+ thetamah(:,3))];

figure
for k = 1:Kmax
    plot([x0(k,1); x1(k,1)],[x0(k,2); x1(k,2)],'r');
    grid on
    hold on
    plot([x1(k,1); x2(k,1)],[x1(k,2); x2(k,2)],'g');
    grid on
    hold on
    plot([x2(k,1); x3(k,1)],[x2(k,2); x3(k,2)],'b');
    grid on
end
% plot that draw the Lassajous curve at the end- effector
rdz = [];
for k = 0:Kmax
    tk = k*g;
 rd = [1.5 + 0.4*sin(pi*tk/5), sqrt(3)/2 + 0.4*sin(pi*tk/5 + pi/3)]; % ring
 rdz = [rdz; rd];
end
plot(rdz(:,1),rdz(:,2))
grid on
    
figure
plot(rdz(:,1),rdz(:,2),'ob', x3(:,1),x3(:,2),'-.r');
legend('Desired path','Actual trajectory')
grid on

figure
%subplot(122)
semilogy(abs(x3(:,2) - rdz(:,2)))
axis([0,200,0,1])
legend('\rm e_2')
grid on

figure
%subplot(121)
semilogy(abs(x3(:,1) - rdz(:,1)))
axis([0,200,0,1])
legend('\rm e_1')
grid on

