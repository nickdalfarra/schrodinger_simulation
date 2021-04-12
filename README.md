# schrodinger_simulation
Simulate and compare numerical methods on the 1D Schrödinger equation. Group project for a course on Numerical Methods for PDEs.
Contributors: Luke Sianchuk, María Preciado, Nick Dal Farra.

Each folder contains Leapfrog, Crank-Nicolson, and classic RK4 solvers for their respective problems. Most of the solvers can easily be applied to other problems by substituting the quantum potential function and changing the parameters x, N, dt, and t_end.

Stability conditions are not calculated automatically and therefore dt must be sufficiently small for choices of dx in the Leapfrog and RK4 methods. The RK4 method requires particularly small choices of dt. Crank-Nicolson is stable for all choices of dt and dx, which is a common theme for implicit numerical schemes.

