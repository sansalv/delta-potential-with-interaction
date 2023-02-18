# Numerical simulation of a quantum physical system

## System and progress:

Two particles in a box with Dirac's delta potential and harmonic quadratic interaction. The Hamiltonian (in natural units) is
$$H(x_1, x_2) = \frac12 p^2_1+\frac12 p^2_2+v_1\delta(x_1)+v_2\delta(x_2)+V_\text{int}(x_1-x_2)$$
where the interaction potential is harmonic quadratic potential (attractive or repulsive) between the two particles
$$V_\text{int}(x_1-x_2) = \frac12 k (x_1-x_2)^2.$$
Few years ago, I calculated the analytical result for a particle in a box with Dirac's delta potential (which is the pdf in Finnish). Without the interaction, two particles separate to same solutions, obliously. But with quadratic interactions, analytical result seems non-trivial to calculate.

This is why me and @Jesse Huhtala (PhD student from theoretical physics department at the University of Turku) worked on this numerical simulation of system. Energy modes are quite simple to prosess. These energy eigenstate solutions are plotted. We are also interested in the time dependend Gaussian wave packet resonances but they are less trivial to simulate. The code simulating these states are not here.
