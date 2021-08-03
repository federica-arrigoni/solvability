
------------------------------------------------------
--------------   GROBNER BASIS solver  ---------------
------------------------------------------------------
solveI = (I,R,inQQ) -> (

J:={};
-- compute Grobner basis
if inQQ then 
   time J = IGeneratedByGBofI I 
else (
time J = groebnerBasis(I, Strategy=>"F4"); -- assumes working in a prime field
);
J = sub(J,R); -- fix case J=|1| to get J = |1_R|
forceGB J; 
I = ideal J;

dm:={}; dg:={}; -- initialize variables
time dm = dim I;
time dg = degree I;

-- solve equations
time IG:=gens gb I;

-- prepare output
O:=new MutableList from {};
O#0=(dm,dg);
O#1=IG;
O

)



------------------------------------------------------
--------------   JACOBIAN check  ---------------------
------------------------------------------------------
jacobianCheck = (I,R,nc,mL,verbose) -> (

-- Jacobian check

Jc:={};
time Jc = jacobian(I); -- 4 s
sol:=ones(1,nc) | zeros(1,4*mL) | -ones(1,mL);
sv:= apply(gens R,ent sol,(x,v)->x=>v); -- fabricate a solution
g:= ent gens I;
v:= apply(g,g->sub(g,sv)); 

-- must be all zero for a fabricated solution
if verbose then print(all(v,v->v==0)); 

Jv:= sub(Jc,sv); -- evaluate J at the fabricated solution

r:= rank(Jv); -- the rank of the jacobian

r==numgens R -- MUST be true for a finite # of solutions

)



------------------------------------------------------
--------------   NECESSARY condition  ----------------
------------------------------------------------------
necessarySolvable = (G) -> (

-- create variables about the input graph
V:= vertexSet G;
E:= edges G;
n:= length V; -- number of vertices
m:= length E; -- number of edges 

-- check necessary conditions: the graph is NOT solvable if a condition is not satisfied  
NC1:= 7*m >= 11*n-15; -- number of edges
NC2:= isConnected G; -- connected
NC3:= not (vertexConnectivity G == 1); -- biconnected 

NC1 and NC2 and NC3

)



------------------------------------------------------
--------   CYCLE-CONSISTENCY equations  --------------
------------------------------------------------------
cycleConsistency = (R,C,BL,VL,u,c,b,verbose) -> (

nc:=#C;

I:= ideal(); -- empty ideal is in Z by default
promote(I,R);
I4:= id_(R^4); -- the identity matrix

-- equations for all cycles
p:=0;
while p <=nc-1 do 
(
cp:=C#p; -- current cycle

if verbose then (
print("Current cycle = ");
print(cp);
);

J:=I4;
SC:=1; -- initialization
for l from 0 to length(cp)-2 do
(
alpha:=cp#(l); -- current edge in GL
beta:=cp#(l+1);
k:=getEdge(BL,alpha,beta); -- corresponding edge in BL: used to index variables u 
i:=getCommonNodeLine(VL,alpha,beta); -- node in G: used to index camera centres 
sc:=(1+( trn(c_{i})*trn(u^{k}) )_(0,0) )*(1-edgeSign(alpha,beta))/2 + (1+edgeSign(alpha,beta))/2;
J=J*(sc*I4+edgeSign(alpha,beta)*c_{i}*u^{k});
SC=SC*sc;
); 
I=I+ideal ent( SC*b_p*I4 - J ); 
p=p+1;
);

I
)



------------------------------------------------------
--------   PATH-INVARIANCE equations  ----------------
------------------------------------------------------
pathInvariance = (R,C,BL,VL,u,c,b,verbose) -> (

nc:=#C;

I:= ideal(); -- empty ideal is in Z by default
promote(I,R);
I4:= id_(R^4); -- the identity matrix

-- equations for all cycles
p:=0;
while p <=nc-1 do 
(
cp:=C#p; -- current cycle
ind:=ceiling((length cp)/2);
Ca:=cp_{0..ind-1};
Cb:=reverse(cp_{ind-1..(length cp)-1});

if verbose then (
print("Current cycle = ");
print(cp); print(Ca); print(Cb);
);

Ja:=I4; -- initialization
Jb:=I4; -- initialization
SCa:=1; -- initialization
SCb:=1; -- initialization

for l from 0 to length(Ca)-2 do
(
alpha:=Ca#(l); -- (alpha,beta) is the current edge in the cycle in LG
beta:=Ca#(l+1);
k:=getEdge(BL,alpha,beta); -- corresponding edge in BL: used to index variables u 
i:=getCommonNodeLine(VL,alpha,beta); -- node in G: used to index camera centres 
sc:=(1+( trn(c_{i})*trn(u^{k}) )_(0,0) )*(1-edgeSign(alpha,beta))/2 + (1+edgeSign(alpha,beta))/2;
Ja=Ja*(sc*I4+edgeSign(alpha,beta)*c_{i}*u^{k});
SCa=SCa*sc;
); 

for l from 0 to length(Cb)-2 do
(
alpha:=Cb#(l); -- (alpha,beta) is the current edge in the cycle in LG
beta:=Cb#(l+1);
k:=getEdge(BL,alpha,beta); -- corresponding edge in BL: used to index variables u 
i:=getCommonNodeLine(VL,alpha,beta); -- node in G: used to index camera centres 
sc:=(1+( trn(c_{i})*trn(u^{k}) )_(0,0) )*(1-edgeSign(alpha,beta))/2 + (1+edgeSign(alpha,beta))/2;
Jb=Jb*(sc*I4+edgeSign(alpha,beta)*c_{i}*u^{k});
SCb=SCb*sc;
); 

I=I+ideal ent( SCa*b_p*Jb - SCb*Ja ); 
p=p+1;
);

I
)



------------------------------------------------------
--------  Invertible Matrices Constraint -------------
------------------------------------------------------
constraintInvertible = (R,mL,BL,VL,c,u,z) -> (

I:= ideal(); -- empty ideal is in Z by default
promote(I,R);
I4:= id_(R^4); -- the identity matrix

-- add auxiliary equations that remove non-invertible matrices
for k from 0 to mL-1 do
(
ind:=getNodes(BL,k); -- nodes corresponding to edge k
alpha:=ind#0; -- left endpoint
beta:=ind#1; -- right endpoint
i:=getCommonNodeLine(VL,alpha,beta); -- node in G: used to index camera centres 
I=I+ideal(z_k*det(I4+c_{i}*u^{k})+1); -- enforce matrices to be invertible
);

I
)




------------------------------------------------------
--------  SOLVABILITY Equations ----------------------
------------------------------------------------------
solvabilityEq = (G,useSparse,usePath,allCycles,useLex,verbose,inQQ) -> (

if verbose then gbTrace=3 else gbTrace=0;

-- create variables about the input graph
V := vertexSet G;
E := edges G;
--A := adjacencyMatrix G;
--B := incidenceMatrix G;
n := length V; -- number of vertices
m := length E; -- number of edges 

-- compute the line graph
L := lineGraph G;
VL := vertexSet L;
EL := edges L; 
BL := incidenceMatrix L;
AL := adjacencyMatrix L;
nL := length VL; -- number of vertices
mL := length EL; -- number of edges 

-- compute a set of cycles
C:={};
if allCycles then
(
C=getCycles(graph(AL)); -- ALL cycles
)
else
(
C =fundCycleBasis(L); -- compute a fundamental cycle basis
);
nc:=#C; -- number of cycles

if verbose then (
ll:={}; -- length of cycles
for k from 0 to nc-1 do (
ll=append(ll,length(C#k));
);
print("Length of cycles = "); print(ll); 
);

-- the ambient field
FF:={};
if inQQ then FF=QQ else FF = ZZ/30011; 

-- known variables (random camera centers):
c:={};
if useSparse then (
if n<4 then c=random(FF^4,FF^n) else c=(id_(FF^4) | random(FF^4,FF^(n-4)));)
else c = random(FF^4,FF^n);

if verbose then (
print("Camera centers = "); print(c);
print("Rank of camera centers = "); print(rank c);  -- check generic camera centres
);

-- count unknowns:  
-- one scale for each cycle in LG
-- one 4-vector for each edge in LG 
-- one auxiliary variable for each edge in LG (remove invertible matrices)

R:={};
if useLex then 
(
R= FF[b_0..b_(nc-1),x_(0)..x_(4*mL-1),z_0..z_(mL-1),MonomialOrder=>Lex];
)
else 
(
R= FF[b_0..b_(nc-1),x_(0)..x_(4*mL-1),z_0..z_(mL-1)]; 
);
-- R:= FF[b_0..b_(nc-1),x_(0)..x_(4*mL-1),z_0..z_(mL-1)];
u:=genericMatrix(R,(vars R)_(0,nc),mL,4); 
I4:= id_(R^4); -- the identity matrix

-- construct equations by path invariance
I:={};
if usePath then 
(
I = pathInvariance(R,C,BL,VL,u,c,b,verbose);
)
else
(
I = cycleConsistency(R,C,BL,VL,u,c,b,verbose);
);

if verbose then (
print("Number of Equations = "); print(16*nc); -- number of equations
print("Number of Unknowns = "); print(4*mL+nc); -- number of unknowns
);

-- add auxiliary equations to impose invertible matrices
I = I + constraintInvertible(R,mL,BL,VL,c,u,z);

if verbose then (
print("Number of auxiliary variables/equations = "); print(mL); 
);

-- prepare output
O:=new MutableList from {};
O#0=I;
O#1=R;
O#2=mL;
O#3=nc;
O

)



------------------------------------------------------
------------------  SOLVABILITY ----------------------
------------------------------------------------------
solvability = (G,useSparse,usePath,allCycles,useLex,verbose,inQQ) ->(

O:=solvabilityEq(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);
I:=O#0;
R:=O#1;
mL:=O#2;
nc:=O#3;

-- compute Grobner basis
oo:=solveI(I,R,inQQ);
print("Is the graph solvable?");
print(oo#0#0==0 and oo#0#1==1);

join(O,oo)

)



------------------------------------------------------
------------------  FINITE SOLVABILITY ---------------
------------------------------------------------------
finiteSolvability = (G,useSparse,usePath,allCycles,useLex,verbose,inQQ) ->(

O:=solvabilityEq(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);
I:=O#0;
R:=O#1;
mL:=O#2;
nc:=O#3;

-- check finite number of solutions
oj:=jacobianCheck(I,R,nc,mL,verbose);
print("Is there a finite number of solutions?");
print(oj);

append(O,oj)

)



------------------------------------------------------
--------  SOLVABILITY + SATURATION Equations ---------
------------------------------------------------------
solvabilityEqSaturation = (G,useSparse,usePath,allCycles,useLex,verbose,inQQ) -> (

if verbose then gbTrace=3 else gbTrace=0;

-- create variables about the input graph
V := vertexSet G;
E := edges G;
--A := adjacencyMatrix G;
--B := incidenceMatrix G;
n := length V; -- number of vertices
m := length E; -- number of edges 

-- compute the line graph
L := lineGraph G;
VL := vertexSet L;
EL := edges L; 
BL := incidenceMatrix L;
AL := adjacencyMatrix L;
nL := length VL; -- number of vertices
mL := length EL; -- number of edges 

-- compute a set of cycles
C:={};
if allCycles then
(
C=getCycles(graph(AL)); -- ALL cycles
)
else
(
C =fundCycleBasis(L); -- compute a fundamental cycle basis
);
nc:=#C; -- number of cycles

if verbose then (
ll:={}; -- length of cycles
for k from 0 to nc-1 do (
ll=append(ll,length(C#k));
);
print("Length of cycles = "); print(ll); 
);

-- the ambient field
FF:={};
if inQQ then FF=QQ else FF = ZZ/30011; 

-- known variables (random camera centers):
c:={};
if useSparse then (
if n<4 then c=random(FF^4,FF^n) else c=(id_(FF^4) | random(FF^4,FF^(n-4)));)
else c = random(FF^4,FF^n);

if verbose then (
print("Camera centers = "); print(c);
print("Rank of camera centers = "); print(rank c);  -- check generic camera centres
);

-- count unknowns:  
-- one scale for each cycle in LG
-- one 4-vector for each edge in LG 
-- one auxiliary variable for each edge in LG (remove invertible matrices)

R:={};
if useLex then 
(
R= FF[b_0..b_(nc-1),x_(0)..x_(4*mL-1),MonomialOrder=>Lex];
)
else 
(
R= FF[b_0..b_(nc-1),x_(0)..x_(4*mL-1)];
);
u:=genericMatrix(R,(vars R)_(0,nc),mL,4); 
I4:= id_(R^4); -- the identity matrix

-- construct equations by path invariance
I:={};
if usePath then 
(
I = pathInvariance(R,C,BL,VL,u,c,b,verbose);
)
else
(
I = cycleConsistency(R,C,BL,VL,u,c,b,verbose);
);

if verbose then (
print("Number of Equations = "); print(16*nc); -- number of equations
print("Number of Unknowns = "); print(4*mL+nc); -- number of unknowns
);

-- solve equations BEFORE saturation
oo:=solveI(I,R,inQQ);
if verbose then (
dm:=oo#0#0;
dg:=oo#0#1;
print("Degree = "); print(dg); -- degree
print("Dimension = "); print(dm); -- dimension
);

I=ideal(oo#1);

if verbose then print("Performing saturation...");

-- enforce matrices to be invertible
for k from 0 to mL-1 do
(
ind:=getNodes(BL,k); -- nodes corresponding to edge k
alpha:=ind#0; -- left endpoint
beta:=ind#1; -- right endpoint
i:=getCommonNodeLine(VL,alpha,beta); -- node in G: used to index camera centres 
I=saturate(I,det(I4+c_{i}*u^{k}));
);

-- prepare output
O:=new MutableList from {};
O#0=I;
O#1=R;
O#2=mL;
O#3=nc;
O

)



------------------------------------------------------
------------------  SOLVABILITY + SATURATION ---------
------------------------------------------------------
solvabilitySaturation = (G,useSparse,usePath,allCycles,useLex,verbose,inQQ) ->(

O:=solvabilityEqSaturation(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);
I:=O#0;
R:=O#1;
mL:=O#2;
nc:=O#3;

-- compute Grobner basis
oo:=solveI(I,R,inQQ);

if verbose then (
dm:=oo#0#0;
dg:=oo#0#1;
print("After saturation:");
print("Degree = "); print(dg); -- degree
print("Dimension = "); print(dm); -- dimension
);

print("Is the graph solvable?");
print(oo#0#0==0 and oo#0#1==1);

join(O,oo)

)
