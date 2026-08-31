restart

--this parametrizes the neuromanifold, as a pair of quadratic forms with a common real linear factor, as in Section 4.2:
S = QQ[a_(1,1)..a_(2,2),v_1,v_2];
S[x_1,x_2];
alpha_1 = a_(1,1)*x_1+a_(1,2)*x_2;
alpha_2 = a_(2,1)*x_1+a_(2,2)*x_2;
nu = v_1*x_1+v_2*x_2;
C = (last coefficients (alpha_1*nu))||(last coefficients (alpha_2*nu));

--this computes the defining equation of the neurovariety:
R = QQ[y_1..y_6];
F = map(S,R,transpose sub(C,S));
eq = (gens kernel F)_(0,0);

--now we compute the parametrization of the neuromanifold in its original coordinates, using v and w and q1 and q2 as at the beginning of Section 4.2, to compute the data-dependent distance from it::
T = QQ[y_1..y_6,x_1,x_2,v_1,v_2,w_1,w_2];
mons = sub(first coefficients (alpha_1*nu), T);
q_1 = (mons * transpose (matrix {{y_1,y_2,y_3}}))_(0,0);
q_2 = (mons * transpose (matrix {{y_4,y_5,y_6}}))_(0,0);
subV = {x_1 => v_1, x_2 => v_2};
subW = {x_1 => w_1, x_2 => w_2};
eq1 = v_1*sub(q_1,subV) + v_2*sub(q_2,subV) + v_1*sub(q_1,subW) + v_2*sub(q_2,subW)
eq2 = w_1*sub(q_1,subW) + w_2*sub(q_2,subW) + w_1*sub(q_1,subV) + w_2*sub(q_2,subV)


---------------
--now we compute the regular critical points of the MSE loss for generic data of cardinality s=1
s = 1;
r = 2*s;
k = 6-r;
inpVars = {v_1,v_2,w_1,w_2};
SdataInp = apply(s, i -> apply(4, j -> inpVars#j => random(QQ)));
SdataOut = apply(s, i -> apply(2, j -> random(QQ)));
dist = sum apply(s, i -> (sub(eq1, SdataInp#i)-SdataOut#i#0)^2+(sub(eq2, SdataInp#i)-SdataOut#i#1)^2);

--this computes the regular critical points:
J = jacobian eq;
AJ = (jacobian sub(dist,R)) | J;
I = ideal(eq) + minors(2,AJ);
Sing = radical ideal singularLocus ideal eq;
Crit = saturate(I, Sing);
dec = decompose Crit;
dec/dim, dec/degree
--this is a threefold of degree 4

--it consists of zero loss solutions:
sub(dist, R/(first dec))
---------------

---------------
--now we compute the regular critical points of the MSE loss for generic data of cardinality s=2
s = 2;
r = 2*s;
k = 6-r;
inpVars = {v_1,v_2,w_1,w_2};
SdataInp = apply(s, i -> apply(4, j -> inpVars#j => random(QQ)));
SdataOut = apply(s, i -> apply(2, j -> random(QQ)));
dist = sum apply(s, i -> (sub(eq1, SdataInp#i)-SdataOut#i#0)^2+(sub(eq2, SdataInp#i)-SdataOut#i#1)^2);

--this computes the regular critical points:
J = jacobian eq;
AJ = (jacobian sub(dist,R)) | J;
I = ideal(eq) + minors(2,AJ);
Sing = radical ideal singularLocus ideal eq;
Crit = saturate(I, Sing);
dec = decompose Crit;
dec/dim, dec/degree
--this is 2 curves of degree 2

--the first consists of zero loss solutions:
sub(dist, R/(first dec))

--now we compute the ramification locus and check that it contains the other curve:
Q = jacobian transpose jacobian sub(dist, R);
K = gens kernel Q;
RamCond = (transpose K) * J;
Ram1 = ideal(eq) + ideal {RamCond_0_0, RamCond_0_1};
Ram = saturate(Ram1, Sing);
dim Ram
isSubset(Ram,last dec)
radical ideal singularLocus last dec

--we can see that there are only 3 independent affine linear equations in the ideal defining this curve:
transpose mingens last dec
--this means that the affine-linear span of the curve is a 3-dimensional affine space
--an irreducible quadratic curve only spans a 2-dimensional space, hence the curve must factor into 2 skew lines
---------------

---------------
--now we compute the regular critical points of the MSE loss for generic data of cardinality s=3
s = 3;
r = 2*s;
k = 6-r;
inpVars = {v_1,v_2,w_1,w_2};
SdataInp = apply(s, i -> apply(4, j -> inpVars#j => random(QQ)));
SdataOut = apply(s, i -> apply(2, j -> random(QQ)));
dist = sum apply(s, i -> (sub(eq1, SdataInp#i)-SdataOut#i#0)^2+(sub(eq2, SdataInp#i)-SdataOut#i#1)^2);

--this computes the regular critical points:
J = jacobian eq;
AJ = (jacobian sub(dist,R)) | J;
I = ideal(eq) + minors(2,AJ);
Sing = radical ideal singularLocus ideal eq;
Crit = saturate(I, Sing);
dec = decompose Crit;
dec/dim, dec/degree
--this is 14 points
---------------

---------------
--now we compute the regular critical points *on the singular locus of the neurovariety* for the MSE loss for generic data of cardinality s=3
s = 3;
r = 2*s;
k = 6-r;
inpVars = {v_1,v_2,w_1,w_2};
SdataInp = apply(s, i -> apply(4, j -> inpVars#j => random(QQ)));
SdataOut = apply(s, i -> apply(2, j -> random(QQ)));
dist = sum apply(s, i -> (sub(eq1, SdataInp#i)-SdataOut#i#0)^2+(sub(eq2, SdataInp#i)-SdataOut#i#1)^2);

--this computes the regular critical points:
Sing = radical ideal singularLocus ideal eq;
Sing2 = radical ideal singularLocus Sing;
J = jacobian Sing;
AJ = (jacobian sub(dist,R)) | J;
I = Sing + minors(3,AJ);
Crit = saturate(I, Sing2);
dec = decompose Crit;
dec/dim, dec/degree
--this is 10 points
---------------

---------------
--now we compute the regular critical points *on the singular locus of the neurovariety* for the MSE loss for generic data of cardinality s=2
s = 2;
r = 2*s;
k = 6-r;
inpVars = {v_1,v_2,w_1,w_2};
SdataInp = apply(s, i -> apply(4, j -> inpVars#j => random(QQ)));
SdataOut = apply(s, i -> apply(2, j -> random(QQ)));
dist = sum apply(s, i -> (sub(eq1, SdataInp#i)-SdataOut#i#0)^2+(sub(eq2, SdataInp#i)-SdataOut#i#1)^2);

--this computes the regular critical points:
Sing = radical ideal singularLocus ideal eq;
Sing2 = radical ideal singularLocus Sing;
J = jacobian Sing;
AJ = (jacobian sub(dist,R)) | J;
I = Sing + minors(3,AJ);
Crit = saturate(I, Sing2);
dec = decompose Crit;
dec/dim, dec/degree
--this is a curve of degree 6

--the loss is not 0 along the curve:
sub(dist, R/(first dec))

Q = jacobian transpose jacobian sub(dist, R);
K = gens kernel Q;
RamCond = (transpose K) * J;
rank sub(RamCond, R/Sing)
Ram1 = Sing + minors(2,RamCond);
Ram = saturate(Ram1, Sing2);
dim Ram, degree Ram
isSubset(Ram,first dec)
radical (Crit+Ram)
--the curve does not intersect the ramification locus
---------------

---------------
--now we compute the regular critical points *on the singular locus of the neurovariety* for the MSE loss for generic data of cardinality s=1
s = 1;
r = 2*s;
k = 6-r;
inpVars = {v_1,v_2,w_1,w_2};
SdataInp = apply(s, i -> apply(4, j -> inpVars#j => random(QQ)));
SdataOut = apply(s, i -> apply(2, j -> random(QQ)));
dist = sum apply(s, i -> (sub(eq1, SdataInp#i)-SdataOut#i#0)^2+(sub(eq2, SdataInp#i)-SdataOut#i#1)^2);

--this computes the regular critical points:
Sing = radical ideal singularLocus ideal eq;
Sing2 = radical ideal singularLocus Sing;
J = jacobian Sing;
AJ = (jacobian sub(dist,R)) | J;
I = Sing + minors(3,AJ);
Crit = saturate(I, Sing2);
dec = decompose Crit;
dec/dim, dec/degree
--this is 2 planes

--the second consists of zero loss solutions:
sub(dist, R/(last dec))

--the loss is constant along the first plane:
sub(dist, R/(first dec))

--the first plane is contained in the ramification locus:
Q = jacobian transpose jacobian sub(dist, R);
K = gens kernel Q;
RamCond = (transpose K) * J;
Ram1 = Sing + minors(2,RamCond);
Ram = saturate(Ram1, Sing2);
dim Ram
isSubset(Ram,first dec)
---------------


