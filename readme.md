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

## Test report
For test control allocation method of the book, run `test_ACA.m` of [control_allocation](https://github.com/mengchaoheng/control_allocation.git) , other test file can use like this.

we Use the moments in all directions on the unit sphere as output and test the output, or we can use the moments input from a fly data (load fly.mat).

* LPwrap
    1. Testing using hover flight data:
        - Effectors jitter appears when LPmethod=0, 1. 
        - Effectors saturation when LPmethod=2
        - warning when LPmethod=0, 1, 4.
    2. Tests in any unit direction：
        - All tests passed but have warning: "Warning:
        - Matrix is singular within working precision" when LPmethod=0, 1, 4.
* CGIwrap:
    - All tests passed

* DAwrap:
    - All tests passed
## Change log
1. add readme.md
2. fix bug of `LPwrap.m`,  the diff is:
```sh
diff --git a/LPwrap.m b/LPwrap.m
index a7d0281..1c223f0 100644
--- a/LPwrap.m
+++ b/LPwrap.m
@@ -88,9 +88,9 @@
 wu=ones(m_act,1);
 emax=ones(k,1)*2e1;
 itlim=5e2;
-lam=0.01;
+lam=1;
 eMax=emax;
-w=0.01*wu;
+w=wu;
 
 switch LPmethod
     case 0
@@ -307,7 +307,7 @@
 %
 indn = 1:n;
 numzer = length(find(eyd==0));
-inBi = [indn(eyd>0) n+indn(eyd<0) (2*n)+(1:numzer)];
+inBi = [indn(eyd>0) n+indn(eyd<0) (2*n):( (2*n)-1+numzer )];
 
 e = true(2*n+m,1);
 
@@ -798,9 +798,9 @@
 
 if errout ~=0  % Construct an incorrect solution to accompany error flags
     xout = zeros(m+1,1);
-    indv = inB1<=(m+1);
+    indv = indB1<=(m+1);
     xout(inB1(indv)) = y1(indv);
-    xout(~e1(1:m+1)) = -xout(~e1(1:m+1))+h(~e1(1:m+1));
+    xout(~e1(1:m+1)) = -xout(~e(1:m+1))+h(~e1(1:m+1));
     
 else  % No Error continue to solve problem
     
@@ -949,7 +949,7 @@
 if errout ~=0  % Construct an incorrect solution to accompany error flags
     xout = zeros(m,1);
     xout(inB1(1:m)) = y1(1:m);
-    xout(~e1(1:m)) = -xout(~e1(1:m))+h(~e1(1:m));
+    xout(~e1(1:m)) = -xout(~e(1:m))+h(~e1(1:m));
     
 else  % No Error continue to solve problem
     
@@ -1087,7 +1087,7 @@
 %
 indn = 1:n;
 numzer = length(find(eyp==0));
-inBi = [indn(eyp>0) n+indn(eyp<0) (2*n)+(1:numzer)];
+inBi = [indn(eyp>0) n+indn(eyp<0) (2*n):( (2*n)-1+numzer )];
 e = true(2*(n+m),1);
 
 %Solve using Bounded Revised Simplex

```
## Note 
For aircraft simulation, the reader maybe need to learn about the ADMIRE, but for study control allocation method, just use the algorithm Implements files `xx_wrap` as function, then create and run `test_xx` is enough.
