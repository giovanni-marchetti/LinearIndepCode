m=2;
n=3;
d=m*n;
r=1;
S=1;
FF=QQ;
R=FF[x_{1,1}..x_{m,n}];
-- Construct the determinantal variety X of mxn matrices of rank r.
M=matrix{{x_{1,1}..x_{1,n}}};
for i from 2 to m do M=M||matrix{{x_{i,1}..x_{i,n}}};
I=minors(r+1,M);
-- Take a random point u
u = random(FF^1,FF^d);
q={};
for i from 1 to m do ( for i from 1 to n do (if i<S+1 then q=q|{1} else q=q|{0}))
Q=diagonalMatrix(q);
-- Construct the ideal D corresponding to C(u)
j = numgens I;
W=FF[x_{1,1}..x_{m,n},y_1..y_j];
X=matrix{{x_{1,1}..x_{m,n}}};
L=((X-u)*Q*transpose(X-u))_(0,0);
f=map(W,R);
finv=map(R,W);
for i from 1 to j do L=L-y_i*f(I_(i-1));
J=minors(1,diff(vars W, L));
D=finv(eliminate(toList(y_1..y_j), J));
-- Verify
D==I+D
-- Obtain the components
primes = minimalPrimes(D);
primes/dim
primes/degree


m=2;
n=3;
d=m*n;
r=1;
S=2;
FF=QQ;
R=FF[x_{1,1}..x_{m,n}];
-- Construct the determinantal variety X of mxn matrices of rank r.
M=matrix{{x_{1,1}..x_{1,n}}};
for i from 2 to m do M=M||matrix{{x_{i,1}..x_{i,n}}};
I=minors(r+1,M);
-- Take a random point u
u = random(FF^1,FF^d);
q={};
for i from 1 to m do ( for i from 1 to n do (if i<S+1 then q=q|{1} else q=q|{0}))
Q=diagonalMatrix(q);
-- Construct the ideal D corresponding to C(u)
j = numgens I;
W=FF[x_{1,1}..x_{m,n},y_1..y_j];
X=matrix{{x_{1,1}..x_{m,n}}};
L=((X-u)*Q*transpose(X-u))_(0,0);
f=map(W,R);
finv=map(R,W);
for i from 1 to j do L=L-y_i*f(I_(i-1));
J=minors(1,diff(vars W, L));
D=finv(eliminate(toList(y_1..y_j), J));
-- Verify
D==I+D
-- Obtain the components
primes = minimalPrimes(D);
primes/dim
primes/degree


m=2;
n=3;
d=m*n;
r=1;
S=3;
FF=QQ;
R=FF[x_{1,1}..x_{m,n}];
-- Construct the determinantal variety X of mxn matrices of rank r.
M=matrix{{x_{1,1}..x_{1,n}}};
for i from 2 to m do M=M||matrix{{x_{i,1}..x_{i,n}}};
I=minors(r+1,M);
-- Take a random point u
u = random(FF^1,FF^d);
q={};
for i from 1 to m do ( for i from 1 to n do (if i<S+1 then q=q|{1} else q=q|{0}))
Q=diagonalMatrix(q);
-- Construct the ideal D corresponding to C(u)
j = numgens I;
W=FF[x_{1,1}..x_{m,n},y_1..y_j];
X=matrix{{x_{1,1}..x_{m,n}}};
L=((X-u)*Q*transpose(X-u))_(0,0);
f=map(W,R);
finv=map(R,W);
for i from 1 to j do L=L-y_i*f(I_(i-1));
J=minors(1,diff(vars W, L));
D=finv(eliminate(toList(y_1..y_j), J));
-- Verify
D==I+D
-- Obtain the components
primes = minimalPrimes(D);
primes/dim
primes/degree



