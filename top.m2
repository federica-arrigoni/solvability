--- Macaulay2 Functions
--- Tomas Pajdla
--- pajdla@gmail.com
--- 2018-11-27
---
--- gb with F4: https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.12/share/doc/Macaulay2/Macaulay2Doc/html/_groebner__Basis.html
---
--- Settings
needsPackage("MapleInterface");
printWidth = 150; -- long print lines
printingPrecision = 3; -- readable printing precision
TOP_PRECRR = 128; -- implicit precision toRRp precision
-- Functions and constants
E3 = matrix{{1,0,0},{0,1,0},{0,0,1}} 
toQQ = (x,n) -> ((round(x*10^n))_QQ)/(10^n) -- convert RR number to QQ number with precision n digits after the decimal point
toRRp = (x) -> toRR(TOP_PRECRR,x) -- preset precision toRR
mat = matrix -- matrix shortcut
inv = inverse -- matrix shortcut
rdim = A -> numRows(A) -- number of rows
cdim = A -> numColumns(A) -- number of columns
dims = A -> (rdim(A),cdim(A)) -- matrix size
rnum = A -> rdim(A)
cnum = A -> cdim(A)
ones = (m,n) -> trn(mat({toList(m:1)}))*mat({toList(n:1)}) -- matrix of ones
zeros = (m,n) ->0*ones(m,n) -- matrix of zeros
vo = A -> trn(matrix({flatten(entries(trn(A)))})) -- matrix vectorization
tra = A -> trace(A) -- trace
trn = A -> transpose(A) -- transposition shortcut
dia = x -> diagonalMatrix(x) -- diagonal matrix shortcut
eye = n -> dia(ones(n,1))
diag = A -> (n:=min(dims(A)); apply(toList(0..n-1),i->A_(i,i))) -- extract the diagonal from amatrix
wide = A -> (if rdim(A)>cdim(A) then trn(A) else A) -- make matrices wide
tall = A -> (if cdim(A)>rdim(A) then trn(A) else A) -- make matrices tall
fliplr = A -> A_(toList(reverse(0..cdim(A)-1))) -- flip matrix A horizontally
flipud = A -> A^(toList(reverse(0..rdim(A)-1))) -- flip matrix A verically
xx = x -> (if (class(x)===Matrix) then mat{{0,-x_(2,0),x_(1,0)},{x_(2,0),0,-x_(0,0)},{-x_(1,0),x_(0,0),0}}
                                  else mat{{0,-x_(2),x_(1)},{x_(2),0,-x_(0)},{-x_(1),x_(0),0}}) -- 3x3 antisymmetric marix
ll = x -> mat{{-x_(1,2)},{x_(0,2)},{-x_(0,1)}} -- extract line coordinates from the skew symmetric matrix
elm  = (A,B) -> (local a; local b; m:=matrix({apply(flatten(entries(A)),flatten(entries(B)),(a,b)->a*b)}); trn(reshape((ring(m))^(cdim(A)),(ring(m))^(rdim(A)),m)) ) -- elemenetwise matrix division
matf = (A,f) -> ( m:=matrix({apply(flatten(entries(A)),f)}); trn(reshape((ring(m))^(cdim(A)),(ring(m))^(rdim(A)),m)) ) -- apply function to a matrix
mmf = (A,B,f) -> (local a; local b; m:=matrix({apply(flatten(entries(A)),flatten(entries(B)),(a,b)->f(a,b))}); trn(reshape((ring(m))^(cdim(A)),(ring(m))^(rdim(A)),m)) ) -- apply elementwise function to a pair of matrices
coln = A -> ones(1,rdim(A))*elm(A,A) -- matrix column norms
nul1 = (A) -> (local i; a:=wide(A); m:=trn(mat({apply(toList(0..cnum(a)-1),i->(-1)^i*det(submatrix'(a,,{i})))})); if rdim(A)>cdim(A) then trn(m) else m )
adj = A -> (local i; local j; matrix table(cnum(A),cnum(A),(i,j)->((-1)^(i+j))*det submatrix'(A,{j},{i}))) -- adjugate matrix
M2L = A -> (local i; apply(cnum(A),i->A_{i})) -- matrix to list of column vectors
L2M = L -> matrix({L}) -- list of column vectors to matrix
ent = A -> flatten(entries(trn(A))) -- list of elements of a matrix
Minors = (n,A) -> ent(gens(minors(n,A))) -- List of n x n minors of matrix A 
maxMinors = A -> ent(gens(minors(min(dims(A)),A))) -- list of maxial minors of matrix A 
maxMinorVec = A -> (X:=tall(A); trn(mat({apply(subsets(rdim(X),cdim(X)),i->det(X^i))}))) -- vectors of maximal minors
eqnmaxabs = e -> apply(e,e->1/max(ent(gens(content(e)))/abs)*e) -- normalize coeffients of a list of polynomials to be in [-1,1]
spy = A -> matf(A,x->(if x==0 then 0 else 1))
--- List operations
LIntersect = (A,B) -> (local x; select(A,x->member(x,B)))
LMinus = (A,B) -> (local x; select(A,x->not(member(x,B))))
LZip = (A,B) -> apply(A,B,(a,b)->{a,b})
---
primpart = f -> lcm((ent(gens(content(f))))/denominator)*f -- primitive part of a polynomial in ZZ         

--- GB computation via Maple
-- gbF4Maple = e -> (p1:=toString(e);
--                   p2:=replace("\\}",")",replace("\\{","tdeg(",toString(support(ideal(e))))); 
--                   callMaple(p1,p2,"with(Groebner);returnvalue:=Basis(placeholder1,placeholder2);")
--                  )
--
gbF4Maple = f -> (
    redStr:={"[\\(,\\)\\{\\}]",""};
    resStr:={{"\\(","\\("},{"\\)","\\)"},{"\\{","\\{"},{"\\}","\\}"}};
    vi:=support(ideal(f)); -- M2 unknowns 
    vis:=apply(vi,v->toString(v)); -- M2 unknowns strings
    vos:=apply(vis,v->replace(redStr#0,redStr#1,v)); -- Maple unknowns strings
    fvis:=apply(vis,v->fold((q,r)->replace(r#0,r#1,q),v,resStr)); -- generator unknowns string regular expressions
    fs:=apply(f,f->toString(f)); -- generator strings    	   
    fos:=apply(fs,v->fold((q,r)->replace(r#0,r#1,q),v,LZip(fvis,vos))); -- generators in Maple variables	  
    p1:=toString(fos); -- list of generators in Maple format				  
    p2:=replace("\\}",")",replace("\\{","tdeg(",toString(vos))); -- monomial ordering in Maple format
    gos:=callMapleRawOut(p1,p2,"with(Groebner);returnvalue:=Basis(placeholder1,placeholder2);"); -- GB in Maple format
    gs:=fold((q,r)->replace(r#0,r#1,q),gos,LZip(vos,vis)); -- GB string in M2 format
    value(gs) -- GB in current ring
)	     

--- Transform I to the same I generate by GB computated via Maple, F3 in Macaulay2 or by gb to achieve maximal efficiency
IGeneratedByGBofI = I->(i := sequence(I); if #i==0 then ( -- unit tests
	Rn := QQ[x1,x2,MonomialOrder=>GRevLex];
	Fg := {x1^2+x2^2-2,((3/5)*x1)^2+((4/5)*x2)^2-1};
	Id := ideal(Fg);
	Jd := IGeneratedByGBofI(Id);
	Id == Jd
	) 
    -- function
    else (
  	mpf:=false;
	J:={};
	F:=flatten(entries(gens((I)))); -- list of generators
	R:=ring(F#0);
	try (J=gbF4Maple(F)) then (
	       if instance(J,Symbol) then (
		  mpf=true;
		  print("IGeneratedByGBofI: gbF4Maple empty output, check variable names - going via native Macaulay2 groebnerBasis/gb");
	       ) else (mpf=false;
		       J = matrix {J};
	              );
	    ) else (mpf=true;	    	    
	            print("IGeneratedByGBofI: gbF4Maple call failed going via native Macaulay2 groebnerBasis/gb");
		   ); 
	if (mpf) then (
	    if isQuotientOf(ZZ,coefficientRing(ring(I))) or coefficientRing(ring(I)) === ZZ then (
	        print("IGeneratedByGBofI: groebnerBasis: Strategy=>F4");		
		J=groebnerBasis(I, Strategy=>"F4");
		) else (
		        print("IGeneratedByGBofI: gb");
		        J=generators(gb(I)));
	    );
	J = sub(J,R); -- fix case J=|1| to get J = |1_R|
	forceGB(J); 
	ideal(J)
	)
    )
-- Solve 0-dim system of algebraic equations with QQ coefficients over CC
-- Solving polynomial equations via multiplicative matrix of a random pomynomial in R/I of a radical ideal I for a monomial basis b
-- returns {sols, vars, equation residuals, evaluated monomials, basis b}
SolvePolyEqByEigV = {DebugLevel => 0} >> opt -> X ->
(
    X = sequence(X); -- variable input must work even for single input
    F := X#0; -- equations
    R := ring(F#0); -- ring
    -- handle variable number of input arhumants
    r := {}; if #X>3 then r = X#3;
    b := {}; if #X>2 then b = X#2;
    f := 0_R; if #X>1 then f = X#1;
    -- 
    if opt.DebugLevel>0 then print("Equations F = " | toString(F));
    if f==0 then (f = ((mat(randomMutableMatrix(1,cdim(vars(R)),0.0,100))*trn(vars(R)))_(0,0))_R + 1_R;) else (f = promote(f,R););-- random linear multiplication polynomial;
    if opt.DebugLevel>0 then print("Multiplier f = " | toString(f));
    -- f must not be constant
    if (support(f)=={}) then error("f must not be constant") else if opt.DebugLevel>0 then print("OK - f not constant");
    -- the idel
    I := ideal(F);
    -- Compute the Groebner basis and install it to the list of known GBs
    J := IGeneratedByGBofI(I);
    -- Find the dimension and the degree of V(I) and check that the dimension=0
    dm := dim(J);
    dg := degree(J);
    if (dm=!=0) then error("dim I(F) must be 0") else if opt.DebugLevel>0 then print("OK - dim(I) = 0, deg(I) = " | toString(dg));
    mJ := rsort(unique(flatten(apply(ent(gens(J)),f->ent(monomials(f)))))); -- monmials of J
    s := {{},{},{}}; vr := {}; vs := {}; re := {}; 
    if all(flatten(mJ/degree),d->d<2) then ( -- linear system
	(B,A) := coefficients(gens(J),Monomials => mat({mJ})); A = trn(A); -- A x = 0 system
	a := -A_{cdim(A)-1}; -- right hand side
	A = A_{0..cdim(A)-2}; -- matrix of the system
	vs = solve(A,a); -- solve th elinear syste A x = a
	vr = mJ_{0..#mJ-2}; -- unknowns
	re = apply({ent(vs)},v->apply(F,f->sub(f,apply(vr,v,(b,v)->b=>v)))); -- equation residuals
	s = {ent(vs),vr,re,a2h(vs),trn(B)} -- costruct solutions	
    ) else (
        -- Factor ring
	A = R/J;
	use(R);
	-- The standard monomial basis of A
	B = sub(sort(basis(A),MonomialOrder=>Descending),R);
	if b=={} then (b = B;) else (b = mat({b}););
	if opt.DebugLevel>0 then print("Basis b = " | toString(b));
	-- Find thransformation from B to b (b may have more elements than B)
	C := (coefficients(b%J,Monomials => B))_1;
	if numcols B <= numcols b then (iC := trn(C)*inv(C*trn(C));) else (iC = inv(trn(C)*C)*trn(C););
	if eye(cdim(C))=!=sub(iC*C,ZZ) then (error("Bad change matrix from b to B");) else if opt.DebugLevel>0 then (print("OK - the change matrix from b to B is good"););
	-- Reduce f multiple of B by J
	if r=={} then (r = f*B;) else (r = mat({r}););
	if opt.DebugLevel>0 then print("Reduced pols r = " | toString(r));
	-- Multiplicaton matrix of y in R/I wrt B
	MfB := ((coefficients(r%I,Monomials => B))_1); -- wrt the standard monomial basis B
	Mf := iC*MfB*C; -- wrt the the given basis b
	if opt.DebugLevel>1 then print("Mf = " | toString(Mf));
	-- Check numerical solutions if computing in QQ
	if (coefficientRing(R)===QQ) then (
	    MR := matf(sub(Mf,QQ),x->toRRp(x)); -- to inexact real number
	    if opt.DebugLevel>1 then print("Mf_RR = " | toString(MR));
	    (e,v) := eigenvectors(trn(MR)); -- eigenvalues can be computer only in inexact complex numbers
	    if cdim(b)>1 and last(ent(b))==1 then ( -- normalize monomial vectors
		v = mmf(v,ones(cdim(v),1)*v^{rdim(v)-1},(x,y)->x/y); -- normalize to get evaluation of B in the eigenvectors 
		ix := positions(ent(b),m->(degree(m))#0==1); -- linear monomials
		if set((ent(b))_ix) === set(ent(vars(R))) then ( -- all variabes present, then recover solutions
		    vix := v^ix; -- values of variables
		    vr = ent(trn(b_ix)); -- variables 
		    vs = trn(entries(vix)); -- solution list 
		    re = apply(vs,v->apply(F,f->sub(f,apply(vr,v,(b,v)->b=>v)))); -- equation residuals
		    s = {vs,vr,re}; -- solution structure
		);
	    );
	    s = join(s,{v,trn(b)}) -- return solutions and evauations on the basis
	) 
        else (error("Ring =!= QQ");)
    )
)
--- Perspective projection
h2a = x -> (local a; (trn(mat(apply(entries(trn(x)),a->a/a_(-1)))))^{0..(rdim(x)-2)}) -- dividing the columns of matrices by their last coorinate
a2h = x -> x || ones(1,cdim(x)) -- add row of ones to a matrix
q2R = v -> matrix{{v_(1,0)^2+v_(0,0)^2-v_(2,0)^2-v_(3,0)^2,2*v_(1,0)*v_(2,0)-2*v_(3,0)*v_(0,0),2*v_(1,0)*v_(3,0)+2*v_(2,0)*v_(0,0)},{2*v_(1,0)*v_(2,0)+2*v_(3,0)*v_(0,0),v_(2,0)^2+v_(0,0)^2-v_(1,0)^2-v_(3,0)^2,2*v_(2,0)*v_(3,0)-2*v_(1,0)*v_(0,0)},{2*v_(1,0)*v_(3,0)-2*v_(2,0)*v_(0,0),2*v_(2,0)*v_(3,0)+2*v_(1,0)*v_(0,0),v_(3,0)^2+v_(0,0)^2-v_(1,0)^2-v_(2,0)^2}} -- unit quaternion to rotation matrix
C2R = c -> matrix{{c_(0,0)^2-c_(1,0)^2-c_(2,0)^2+1, 2*c_(0,0)*c_(1,0)+2*c_(2,0), 2*c_(0,0)*c_(2,0)-2*c_(1,0)}, {2*c_(0,0)*c_(1,0)-2*c_(2,0), -c_(0,0)^2+c_(1,0)^2-c_(2,0)^2+1, 2*c_(1,0)*c_(2,0)+2*c_(0,0)}, {2*c_(0,0)*c_(2,0)+2*c_(1,0), 2*c_(1,0)*c_(2,0)-2*c_(0,0), -c_(0,0)^2-c_(1,0)^2+c_(2,0)^2+1}} -- Non-normalized Cayley parameterization of rotations
c2R = c -> (1/(1+c_(0,0)^2+c_(1,0)^2+c_(2,0)^2))*C2R(c) -- Caley parameterization of rotations
R2c = R -> ll(inv(id_((ring(R))^3)+R)*(id_((ring(R))^3)-R)) -- From rotation R to the Cayley parameterization
--- Radial distortion
ru2u = (u,k) -> u || (ones(1,cdim(u))+k*coln(u)) -- radial undistortion - 1p divivision model
u2ru = (v,k) -> ( -- radial distortion - 1p division model
    local a; local b;
    vv := substitute(v,QQ);
    ru2:= matf(coln(vv),toRRp);
    rd := 1/2*(-2*k*ru2+ones(rdim(ru2),cdim(ru2))-matf(-4*k*ru2+ones(rdim(ru2),cdim(ru2)),sqrt))*(1/k_RR)^2;
    rd  = mmf(rd,ru2,(a,b)->sqrt(a/b));
    s  := mmf(rd,ru2,(a,b)->a/sqrt(b));
    elm(vv,ones(rdim(vv),1)*s)
    )
-- Gauss-Jordan elimination
gjelim = A -> (
    x := local x;
    T := ring(A)[x_1 .. x_(cdim(A))]; 
    X := trn(vars(T));
    G := ideal(sub(A,T)*X); 
    G  = mingens(gb(G)); 
    G  = sort(G,MonomialOrder=>Descending);
    (y,B) := coefficients(G,Monomials=>trn(X));
    sub(trn(B),ring(A)))
-- Remove fractions by row mutiplications
RowDeFract = A -> (
    r := M2L(trn(A)); -- list of row vectprs
    r  = apply(r,r->lcm((ent(r))/denominator)*r);
    trn(L2M(r))
    )
--- Plucker coordinates
-- Plucker Matrix from Plucker coordinates or generators
PLM = X -> (x:=tall(X); 
            if dims(x)==(3,1) then xx(x) 
            else if dims(x)==(4,2) then x_{0}*trn(x_{1})-x_{1}*trn(x_{0})
	                           else mat({{       0,-x_(0,0),-x_(1,0),-x_(2,0)},
		    	                     { x_(0,0),       0,-x_(3,0),-x_(4,0)},
			                     { x_(1,0), x_(3,0),       0,-x_(5,0)},
			                     { x_(2,0), x_(4,0), x_(5,0),       0}}))
-- Plucker coordinates from Plucker matrix or 2x4 or 4x2 matrix
PLC = X -> (if (dims(X)==(3,1) or dims(X)==(1,3)) then tall(X)
            else if dims(X)==(3,3) then gens(kernel(X))
            else if dims(X)==(4,4) then trn(mat({{X_(1,0),X_(2,0),X_(3,0),X_(2,1),X_(3,1),X_(3,2)}})) 
                                   else (-maxMinorVec(X))^{0,1,3,2,4,5})
-- Dual to Plucker matrix
--  sub((1/(gens(radical(minors(4,L))))_(0,0))*promote(adj(L),frac(R)),R)
-- -trn(adj(L^{2,3}_{2,3})) |  trn(adj(L^{0,1}_{2,3})) || trn(adj(L^{2,3}_{0,1})) | -trn(adj(L^{0,1}_{0,1}))
PLD = L -> (if (dims(L)==(3,3) or dims(L)==(3,1)) then L
            else if dims(L)==(4,4) then mat({{ L_(2,2), L_(3,2), L_(1,3),-L_(1,2)},
                                             { L_(2,3), L_(3,3),-L_(0,3), L_(0,2)},
		                             { L_(3,1),-L_(3,0), L_(0,0), L_(1,0)},
		                             {-L_(2,1), L_(2,0), L_(0,1), L_(1,1)}})
                              else PLC(PLD(PLM(L))))
-- Plucker Klein Quadric from Plucker matrx or Plucker coordinates
PKQ = L -> (if dims(L)==(4,4) then (PLD(L)*L)_(0,0)
                              else (PLD(PLM(L))*PLM(L))_(0,0))
-- Plucker coordinate of two incident lines must satisfy the following incidence condition
PLCI = (K,L) -> (trn(L)*fliplr(diagonalMatrix(ring(K),{1,-1,1,1,-1,1}))*K)_(0,0)
-- Robotics
DhVars = n -> c_1..c_n | s_1..s_n | p_1..p_n | q_1..q_n | a_1..a_n | d_1..d_n | (r_1,r_2,r_3,t_1,t_2,t_3) -- Dennavit-Hartenbergh variables for n axes
DhIds = n -> flatten(transpose(toList(apply(1..n,i->{c_i^2+s_i^2-1,p_i^2+q_i^2-1})))) -- The i-th DH motion matrix
DhM = i -> matrix({{c_i,-s_i*p_i, s_i*q_i,a_i*c_i},
 	           {s_i, c_i*p_i,-c_i*q_i,a_i*s_i},
 		   {  0,     q_i,     p_i,    d_i},
		   {  0,       0,       0,      1}})
DhMi = i-> (M:=DhM(i); trn(M^{0..2}_{0..2})  | -(trn(M^{0..2}_{0..2})*M^{0..2}_{3..3}) || matrix({{0_Rng,0,0,1}})) -- The inverse of the i-th DH motion matrix
DhFKT = n -> fold(apply(toList(1..n),i->DhM(i)),(x,y)->x*y) -- DH Forward kinematic task motion matrix for n motion axes 
t2c = t->(1-t^2)/(1+t^2) -- a rational parameterization of the unit circle
t2s = t->(2*t)/(1+t^2)
DhRandPars = (n,mx) -> (r:=apply(toList(1..n),i->apply(toList(1..4),i->random(mx))); -- random DH parameters for n motion axes
                        flatten(transpose(apply(toList(1..n),i->{c_i=>t2c(r#(i-1)#0),s_i=>t2s(r#(i-1)#0),p_i=>t2c(r#(i-1)#1),q_i=>t2s(r#(i-1)#1),a_i=>r#(i-1)#2,d_i=>r#(i-1)#3})))) 
DhEqns = (n,Mh) -> ent((fold(apply(toList(1..floor(n/2)),i->DhMi(i)),(x,y)->y*x)*Mh-fold(apply(toList((floor(n/2)+1)..n),i->DhM(i)),(x,y)->x*y))^{0..2}) | (DhIds(n))_{0..(n-1)} -- DH equations
end
-- Plucker unit tests
restart
clearAll
load "top.m2"
R = QQ[x_1..x_4,y_1..y_4,L01,L02,L03,L12,L13,L23,p_(1,1)..p_(3,4)]
X = genericMatrix(R,x_1,4,1) | genericMatrix(R,y_1,4,1)
PM = PLM(X)
PLM(trn(X))-PM
PC = PLC(X)
PLC(trn(X))-PC
PLM(PC)-PM
PLC(PM),PC
PLD(PM)*X
PLM(PLD(PC))*X
PLD(PM)*PM
PM = PLM(mat({{L01,L02,L03,L12,L13,L23}}))
det(PM)-((mingens(minors(1,PLD(PM)*PM)))_(0,0))^2
PKQ(PLC(PM))
PKQ(PLM(X))
PKQ(PLC(X))
trn(PLD(PLC(X)))*PLC(X)
PKQ(PM)-PKQ(PLC(PM))
-- in 2D
X = genericMatrix(R,x_1,3,1)
PM = PLM(X)
PLM(trn(X))-PM
PC = PLC(X)
PLC(trn(X))-PC
PLM(PC)-PM
PLC(PM)-PC
PLD(PM)*X
PLM(PLD(PC))*X
PLD(PM)*PM


