

tspan = linspace(0,10);

x0 = ones(15,1);

[t, y] = ode45(@(t,y)ODE(t,y), tspan, x0);


subplot(3,1,1)
plot(t,y(:,1:3:13),"-")
title("Transcriptional states")
legend(["topA" "gyrA" "fis" "cspA" "hns"])

subplot(3,1,2)
plot(t, y(:, 2:3:14), "-")
title("Protein states")
legend(["TopoI" "Gyrase" "Fis" "CspA" "H-NS"])

subplot(3,1,3)
plot(t, y(:, 3:3:15)-1, "-")
title("Supercoiling state")

xlabel("Time")
legend(["topA" "gyrA" "fis" "cspA" "hns"])