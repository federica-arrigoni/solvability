
-- input: incidence matrix B, nodes i and j
-- output: index k corresponding to the edge in the incidence matrix
getEdge = (B,i,j) -> (
rowList:=entries (B^{i}+B^{j});
ind:=positions(rowList#0,l->l==2);
ind#0
)

-- input: incidence matrix B, index k corresponding to one edge in the incidence matrix
-- output: nodes i and j that are endpoints of the edge (as list)
getNodes = (B,k) -> (
-- select one column (=edge) of the incidence matrix
col:=B_{k};
colList:=entries transpose col;
ind:=positions(colList#0,l->l>0);
i:=ind#0; -- 0-based
j:=ind#1; -- 0-based 
{i,j} -- output
)


edgeSign = (alpha,beta) ->(
if alpha>beta then -1 else 1
)


-- input: vertex set VL of the line graph, edge (i,j) in the line graph
-- output: node k that is endpoint of two edges in the original graph 
-- (corresponding to the two nodes in the line graph) 
getCommonNodeLine = (VL,i,j) -> (
ni:=toList(VL#i); -- node i in GL corresponds to an edge in G
nj:=toList(VL#j); -- node j in GL corresponds to an edge in G
nall:=append(append(ni,nj#0),nj#1); -- the desired vertex appears twice
k:=commonest(nall);
k#0
)



fundCycleBasis = (L) -> (

VL:= vertexSet L;
EL:= edges L; 
AL:= adjacencyMatrix L;
nL:= length VL; -- number of vertices
mL:= length EL; -- number of edges 

-- compute spanning forest for the line graph (=spanning tree, the graph is connected)
F:= spanningForest L;
AF:= adjacencyMatrix F; 
BF:= incidenceMatrix graph(AL-AF); -- edges not belonging to the tree

-- compute cycles from spanning tree (Fundamental Cycle Basis)
C:=new MutableList from {};
for k from 0 to mL-nL+1-1 do ( -- edges not belonging to the tree
ind:=getNodes(BF,k); -- nodes corresponding to edge k
i:=ind#0; -- left endpoint
j:=ind#1; -- right endpoint
Cij:=getCycles addEdge(graph(AF),set {i,j}); -- this is a list: we known that it has a unique element (unique cycle)
C#k=Cij#0; -- current cycle
); 

C

)


 
graphOneVertex = (G) -> (
V:= vertexSet G;
E:= edges G;
n:= length V; -- number of vertices
m:= length E; -- number of edges 

Eo:={};
-- add edges in the first graph
for k from 0 to m-1 do (
Eo=append(Eo,toList(E#k));
);
Eo=append(Eo,{V#0,n});
Eo=append(Eo,{V#1,n});

graph(Eo)
)


mergeGraphs = (G1,G2) -> (

V1 := vertexSet G1;
E1 := edges G1;
n1 := length V1; -- number of vertices
m1 := length E1; -- number of edges 

V2 := vertexSet G2;
E2 := edges G2;
n2 := length V2; -- number of vertices
m2 := length E2; -- number of edges 

-- indentify two edges
e1:=toList(E1#0);
v1l:=e1#0;
v1r:=e1#1;
e2:=toList(E2#0);
v2l:=e2#0;
v2r:=e2#1;

E:={};

-- add edges in the first graph
for k from 0 to m1-1 do (
E=append(E,toList(E1#k));
);

-- add edges in the second graph
for k from 1 to m2-1 do (
e:=toList(E2#k); -- current edge
vl:=e#0; -- left endpoint
vr:=e#1; -- right endpoint
ind1:={};
ind2:={};
if vl==v2l then ind1=v1l else (
if vl==v2r then ind1=v1r else ind1=vl+n1;
);
if vr==v2l then ind2=v1l else (
if vr==v2r then ind2=v1r else ind2=vr+n1;
);
E=append(E,{ind1,ind2});
);

G:=graph(E);
A:= adjacencyMatrix G;
graph(A)

)








