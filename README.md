# PS-Methods_1D-Vibrating-Beams
Solutions to Inverse Problems in the field of Mechanical Vibration: Utilizing Pseudo-Spectral Method for the Dynamic Analysis of One-Dimensional Non-Ideal Cantilevers &amp; Fixed Beams.

## Abstract
In this paper, we address the free vibration of non-ideal cantilevers and fixed beams in
one dimension. According to the Euler-Bernoulli theory, a 4th order ODE models the free
vibration of beams. The ODE is refined to an eigenvalue problem, where computing the natural
frequencies become evident. The analytic approach for getting the natural frequencies is proven
to be lengthy and problem-specific. Alternatively, we compute the natural frequencies of the
beams using forward solvers that utilize the Pseudo-spectral method based on the Chebyshev
polynomials of the first kind. The solvers are designed to output the natural frequencies of
the beams subject to non-ideal boundary conditions (i.e supplying decay parameter(s) kL/kR
as input(s)). Inverse solvers for the equivalent inverse problems are also designed, where the
damage parameter(s) to be computed after supplying the natural frequencies of the beams as
inputs. All solvers are written in MATLAB [3] language, and the computational approach is
conclusively accurate. The codes can be modified to suit other beam configurations as well.
