# Simulation Files Of Aircraft Control Allocation
the book "Aircraft Control Allocation" by Wayne Durham, Kenneth A. Bordignon and Roger Beck, which provide this files on [the companion site ](https://www.wiley.com//legacy/wileychi/durham2/). this repos. use this files to developed some new method about control allocation. The master branch of this repos. is default files in the simulation, but the other branch is change by [mengchaoheng](https://github.com/mengchaoheng).

## DESCRIPTION

Aircraft Control Allocation addresses the problem of allocating redundant flight controls. It provides introductory material on flight dynamics and control to provide the context, and then describes in detail the geometry of the problem. The book includes a large section on solution methods, including ‘Banks’ method’, a previously unpublished procedure. Generalized inverses are also discussed at length. There is an introductory section on linear programming solutions, as well as an extensive and comprehensive appendix on linear programming formulations and solutions. Discrete-time or ‘frame-wise’ allocation is described, including rate-limiting, nonlinear data, and preferred solutions.

## Key features:

* Written by pioneers in the field of control allocation

* Comprehensive explanation and discussion of the major control-allocation solution methods

* Extensive treatment of linear programming solutions to control allocation

* A companion web site contains the code of a MATLAB/Simulink light simulation with modules that incorporate all of the major solution methods

* Includes examples based on actual aircraft

* The book is a vital reference for researchers and practitioners working in aircraft control, as well as graduate students in aerospace engineering.


## USING THE LINEAR SIMULATION
the detail of the usage of the simulation can be find on `SimulationQuickStartGuide.pdf`, and appendix B of this book git some review too.

## Use as lib
For test LP control allocation method, run `test_LPwrap`, other test file can use like this.

## Note 
For aircraft simulation, the reader maybe need to learn about the ADMIRE, but for study control allocation method, just use the algorithm Implements files `xx_wrap` as function, then create and run `test_xx` is enough.
