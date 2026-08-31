restart

--this parametrizes the neuromanifold, as pairs of quadratic forms with a common real linear factor, as in Section 4.2:
S = QQ[a_(1,1)..a_(2,2),v_1,v_2]
S[x_1,x_2]
alpha_1 = a_(1,1)*x_1+a_(1,2)*x_2
alpha_2 = a_(2,1)*x_1+a_(2,2)*x_2
nu = v_1*x_1+v_2*x_2
C = (last coefficients (alpha_1*nu))||(last coefficients (alpha_2*nu))

--this computes the vanishing equation of the neurovariety:
R = QQ[y_1..y_6]
F = map(S,R,transpose sub(C,S))
eq = (gens kernel F)_(0,0)

--this computes the singular locus of the neurovariety:
Sing = radical ideal singularLocus ideal eq

--the singular locus of the singular locus is just the origin:
Sing2 = radical ideal singularLocus Sing

--now we compute the regular critical points for generic quadrics Q with k-dimensional kernel
----- k=0 -----
k = 0
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian eq
AJ = (jacobian dist) | J

I = ideal(eq) + minors(2,AJ);
Crit = saturate(I,Sing);
dec = decompose Crit;
dec/dim, dec/degree
-- this is 14 points
---------------

---------------
----- k=1 -----
k = 1
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian eq
AJ = (jacobian dist) | J

I = ideal(eq) + minors(2,AJ);
Crit = saturate(I,Sing);
dec = decompose Crit;
dec/dim, dec/degree
-- this is 4+10 points

--the first 4 are zero loss solutions:
sub(dist, R/(first dec))

--the other 10 lie on the ramification locus:
M = J | A
I = ideal(eq) + minors(r+1,M);
RamL = saturate(I,Sing);
isSubset(RamL, last dec)
---------------

---------------
----- k=2 -----
k = 2
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian eq
AJ = (jacobian dist) | J

I = ideal(eq) + minors(2,AJ);
Crit = saturate(I,Sing);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a curve of degree 4 plus 4 points

--the curve consists of zero loss solutions:
sub(dist, R/(first dec))

--the other 4 points lie on the ramification locus:
M = J | A
I = ideal(eq) + minors(r+1,M);
RamL = saturate(I,Sing);
isSubset(RamL, last dec)
---------------

---------------
----- k=3 -----
k = 3
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian eq
AJ = (jacobian dist) | J

I = ideal(eq) + minors(2,AJ);
Crit = saturate(I,Sing);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a surface of degree 4

--it consists of zero loss solutions:
sub(dist, R/(first dec))

--the ramification locus is empty:
M = J | A
I = ideal(eq) + minors(r+1,M);
RamL = saturate(I,Sing)
---------------

---------------
----- k=4 -----
k = 4
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian eq
AJ = (jacobian dist) | J

I = ideal(eq) + minors(2,AJ);
Crit = saturate(I,Sing);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a threefold of degree 4

--it consists of zero loss solutions:
sub(dist, R/(first dec))

--the ramification locus is empty:
M = J | A
I = ideal(eq) + minors(r+1,M);
RamL = saturate(I,Sing)
---------------

---------------
----- k=5 -----
k = 5
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian eq
AJ = (jacobian dist) | J

I = ideal(eq) + minors(2,AJ);
Crit = saturate(I,Sing);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a fourfold of degree 4

--it consists of zero loss solutions:
sub(dist, R/(first dec))

--the ramification locus is empty:
M = J | A
I = ideal(eq) + minors(r+1,M);
RamL = saturate(I,Sing)
---------------

--now we compute the regular critical points *on the singular locus of the neurovariety* for generic quadrics Q with k-dimensional kernel
----- k=0 -----
k = 0
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian Sing
AJ = (jacobian dist) | J

I = Sing + minors(3,AJ);
Crit = saturate(I,Sing2);
dec = decompose Crit;
dec/dim, dec/degree
-- this is 10 points
---------------

---------------
----- k=1 -----
k = 1
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian Sing
AJ = (jacobian dist) | J

I = Sing + minors(3,AJ);
Crit = saturate(I,Sing2);
dec = decompose Crit;
dec/dim, dec/degree
-- this is 10 points
---------------

---------------
----- k=2 -----
k = 2
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian Sing
AJ = (jacobian dist) | J

I = Sing + minors(3,AJ);
Crit = saturate(I,Sing2);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a 3+7 points

--the first 3 points are zero loss solutions:
sub(dist, R/(first dec))

--the other 7 points lie on the ramification locus:
M = J | A
I = Sing + minors(r+2,M);
RamL = saturate(I,Sing2);
isSubset(RamL, last dec)
---------------

---------------
----- k=3 -----
k = 3
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian Sing
AJ = (jacobian dist) | J

I = Sing + minors(3,AJ);
Crit = saturate(I,Sing2);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a curve of degree 3 plus 3 points

--the curve consists of zero loss solutions:
sub(dist, R/(first dec))

--the 3 points lie on the ramification locus:
M = J | A
I = Sing + minors(r+2,M);
RamL = saturate(I,Sing2);
isSubset(RamL, last dec)
---------------

---------------
----- k=4 -----
k = 4
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian Sing
AJ = (jacobian dist) | J

I = Sing + minors(3,AJ);
Crit = saturate(I,Sing2);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a surface of degree 3

--it consists of zero loss solutions:
sub(dist, R/(first dec))

--the ramification locus is empty:
M = J | A
I = Sing + minors(r+2,M);
RamL = saturate(I,Sing2)
---------------

---------------
----- k=5 -----
k = 5
r = 6-k
A = matrix apply(6, i -> apply(r, j -> random(QQ)))
rank A == r
Q = A*transpose(A)

u = apply(6, i -> random(QQ))
c = gens R
v = matrix{c-u}

dist = (v*Q*transpose(v))_(0,0)

J = jacobian Sing
AJ = (jacobian dist) | J

I = Sing + minors(3,AJ);
Crit = saturate(I,Sing2);
dec = decompose Crit;
dec/dim, dec/degree
-- this is a threefold of degree 3

--it consists of zero loss solutions:
sub(dist, R/(first dec))

--the ramification locus is empty:
M = J | A
I = Sing + minors(r+2,M);
RamL = saturate(I,Sing2)
---------------






