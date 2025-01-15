

tspan = linspace(0,100);

x0 = ones(12,1);

[t, y] = ode45(@(t,y)ODE(t,y), tspan, x0);


subplot(3,1,1)
plot(t,y(:,1:3:10),"-")
title("Transcriptional states")
legend(["TopoI" "Gyrase" "Fis" "cspA"])

subplot(3,1,2)
plot(t, y(:, 2:3:11), "-")
title("Protein states")
legend(["TopoI" "Gyrase" "Fis" "cspA"])

subplot(3,1,3)
plot(t, y(:, 3:3:12), "-")
title("Supercoiling state")

legend(["TopoI" "Gyrase" "Fis" "cspA"])