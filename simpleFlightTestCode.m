%% simple flight test code

clear
clc

u_exit = 1620;

thrust = 3310;

m_dot = thrust/u_exit;

throat_r = 0.01;

exit_r = sqrt(8)*throat_r;

body_r = 0.1;

cross_area = pi*body_r*body_r;

Cd = 1.0;

mass = 180;

empty_mass = 80;

dt = 0.1;

y(1) = 1;

y_dot(1) = 0;

tburnout = 0;
i = 1;
while y(i) >= 0
    rho = exp(-y(i)/8500);
    q = 0.625*rho*y_dot(i)^2;
    f = - sign(y_dot(i))*Cd*cross_area*q;
    if (tburnout < 1)
        f = f+thrust;
        mass = mass - m_dot*dt;
        if (mass < empty_mass)
            tburnout = i*dt;
        end
    end
    a = f/mass-9.8;
    y(i+1) = y(i) + y_dot(i)*dt + a*dt*dt*0.5;
    y_dot(i+1) = y_dot(i) + a*dt;
    i = i+1;
end

figure

plot(y./1600);