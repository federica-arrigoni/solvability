
-- Viewing Graph Solvability demo

clearAll
load "top.m2"
needsPackage "Matroids" 

load "graphUtils.m2"
load "solvabilityUtils.m2"

-- load EXAMPLES

load "examples_minimal.m2" -- minimal graphs up to 25 nodes
--load "examples.m2" -- other graphs

-- choose a graph
G=G20;

-- true if necessary conditions are satisfied
necessarySolvable(G) 

-- parameters (default values are recommended)
useSparse=true; -- TRUE: the first 4 centres form the identity matrix 
usePath=true; -- TRUE: use path-invariance (instead of cycle-consistency)
allCycles=false; -- FALSE: consider a fundamental cycle basis (instead of all cycles)
verbose=false;
inQQ=false; -- FALSE: work in a prime field
useLex=false; -- FALSE: use grevlex order

-- check finite solvability (necessary condition for solvability)
-- Of=finiteSolvability(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);


-- check solvability
time O=solvability(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);

-- I=O#0; -- ideal
-- R=O#1; -- ambient ring
-- mL=O#2; -- number of edges in the line graph
-- nc=O#3; -- number of cycles in the line graph
-- Ig=O#5 -- generators of the Grobner Basis
-- (dm,dg)=O#4 -- (dimension,degree)


-- check solvability with saturation (slow)
-- Os=solvabilitySaturation(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);


