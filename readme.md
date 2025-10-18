# Simulation Files Of Aircraft Control Allocation
the book "Aircraft Control Allocation" by Wayne Durham, Kenneth A. Bordignon and Roger Beck, which provide this files on [the companion site ](https://www.wiley.com//legacy/wileychi/durham2/). this repos. use this files to developed some new method about control allocation. The master branch of this repos. is default files in the simulation.



## Change log

By comparing master with fd4b5.
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

3. update `yout[n,1]  = Optimal output variable` to `yout[m,1]  = Optimal output variable`

4. add LPwrap_par.m for  add parallel computing support.
## Note 
For aircraft simulation, the reader maybe need to learn about the ADMIRE, but for study control allocation method, just use the algorithm Implements files `xx_wrap` as function, then create and run `test_xx` is enough.


 
