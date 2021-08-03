              
-- Viewing Graph Solvability demo

clearAll
load "top.m2"
needsPackage "Matroids" 
load "graphUtils.m2"
load "solvabilityUtils.m2"

-- load DATA
--load "./ICCV21Data/DataPumpkin.txt"
--load "./ICCV21Data/DataAlcatrazCourtyard.txt"
--load "./ICCV21Data/DataBuddahTooth.txt"
--load "./ICCV21Data/DataSkansenKronan.txt"
--load "./ICCV21Data/DataTsarNikolaiI.txt"
--load "./ICCV21Data/DataAlamo.txt"
--load "./ICCV21Data/DataEllis_Island.txt"
--load "./ICCV21Data/DataGendarmenmarkt.txt"
--load "./ICCV21Data/DataMadrid_Metropolis.txt"
--load "./ICCV21Data/DataMontreal_Notre_Dame.txt"
--load "./ICCV21Data/DataNotre_Dame.txt"
--load "./ICCV21Data/DataNYC_Library.txt"
--load "./ICCV21Data/DataPiazza_del_Popolo.txt"
--load "./ICCV21Data/DataPiccadilly.txt"
--load "./ICCV21Data/DataRoman_Forum.txt"
--load "./ICCV21Data/DataTower_of_London.txt"
--load "./ICCV21Data/DataTrafalgar.txt"
--load "./ICCV21Data/DataUnion_Square.txt"
--load "./ICCV21Data/DataVienna_Cathedral.txt"
--load "./ICCV21Data/DataYorkminster.txt"
load "./ICCV21Data/DataQuad.txt"

-- parameters
useSparse=true; -- TRUE: the first 4 centres form the identity matrix 
usePath=true; -- TRUE: use path-invariance (instead of cycle-consistency)
allCycles=false; -- FALSE: consider a fundamental cycle basis (instead of all cycles)
verbose=false;
inQQ=false; -- FALSE: work in a prime field
useLex=false; -- FALSE: use grevlex order

ll={}; gunsol={}; gsol={}; gchord={}; gfin={};
nnCount=0; -- graphs satisfying necessary conditions
ccCount=0; -- chordal graphs
fnCount=0; -- finite solvable graphs
ssCount=0; -- solvable graphs
--length(GG)-1
for k from 0 to length(GG)-1 do (
nn={}; solvable={}; cc={}; dm={}; dg={};
print(k);
G=GG#k;
nn=necessarySolvable(G); -- true if necessary conditions are satisfied
if nn then 
(
cc=isChordal(G); -- check if the graph is chordal (triangulated): sufficient condition
if cc then (solvable=true; dm=0; dg=1;)
else
(
-- check solvability
O=solvability(G,useSparse,usePath,allCycles,useLex,verbose,inQQ);
(dm,dg)=O#4; -- (dimension,degree)
print("Degree = "); print(dg); -- degree
print("Dimension = "); print(dm); -- dimension
solvable=(dm==0 and dg==1);
--fn=dg; -- number of solutions (if finite solvable)
);
) else (
solvable=false;
print("error: no necessary condition satisfied");
print(k);
);
if nn then nnCount=nnCount+1;
if solvable then (ssCount=ssCount+1; gsol=append(gsol,k);) else gunsol=append(gunsol,k);
if cc then (ccCount=ccCount+1; gchord=append(gchord,k););
if dm==0 then fnCount=fnCount+1;
ll=append(ll,{nn,cc,(dm,dg),solvable});
if (not(solvable) and (dm==0)) then gfin=append(gfin, k);
)

print("Total number of graphs = "); print(length(GG)); 
print("Graphs satisfying necessary conditions = "); print(nnCount); 
print("Chordal graphs = "); print(ccCount); 
print("Finite solvable graphs = "); print(fnCount); 
print("Solvable graphs = "); print(ssCount); 
print("Unsolvable graphs = "); print(nnCount-ssCount); 

print("Chordal graphs = "); gchord -- chordal graphs
print("Solvable graphs = "); gsol 
print("Unsolvable graphs = "); gunsol 
print("Finite solvable but not solvable = "); gfin -- finite solvable and not solvable


