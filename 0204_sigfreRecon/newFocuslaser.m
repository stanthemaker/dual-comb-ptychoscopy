current = 44:2:90;
power= [-2 -1.23 -0.55 0.03 0.47 1.02 1.41 1.76 2.13 2.4 2.71 3 3.2 3.5 3.79 4 4.22 4.40 4.59 4.76 4.93 5.09 5.25 5.4];

figure;
plot(current, power);


xlabel("current (mA)");
ylabel("power (dBm)");

set(gca,"Fontsize",18);
grid on;
