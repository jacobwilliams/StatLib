FUNCTION EdgeOperation(VAR Graph : GraphType; m, n : NodeType;
                       NodeSet : CliqueType; Add : boolean) : integer;

{Algorithm AS 247.1  Appl. Statist. (1989) Vol.38, No.2}

VAR
    exist      : boolean;
    error      : integer;
    i          : integer;

BEGIN
exist := false;
WITH Graph DO
    FOR i := 1 TO NoOfCliques DO
        IF ([m,n] <= cliques[i]) THEN
            exist := true;
error := 0;

{ Check for invalid inputs.}
IF NOT (m IN NodeSet) OR NOT (n IN NodeSet) THEN
    m := n;
IF m = n THEN
    error := 1
ELSE
    IF exist AND Add THEN
        error := 2
    ELSE
        IF NOT exist AND NOT Add THEN
            error := 3
        ELSE
            {If input ok then perform operation.}
            IF Add THEN
                AddEdge(Graph,m,n)
            ELSE
                DeleteEdge(Graph,m,n);
EdgeOperation := error;
END; { of EdgeOperation }



 PROCEDURE AddEdge (VAR Graph : GraphType;
         m, n : NodeType);

 {Algorithm AS 247.2  Appl. Statist. (1989) Vol.38, No.2}

    { This procedure replaces the cliques of the current graph held in the 
      structure cliques of type GraphType by the cliques of the graph      
      obtained by adding edge m-n.}

  VAR
   mSubgraph, nSubgraph, mnSubgraph : GraphType;
   ThisClique, mClique, nClique, mnClique : CliqueType;
   i, j, k, NoOfCliques : CliqueIndex;
   mCliqueOK, nCliqueOK, mnCliqueOK : ARRAY[1..Cliquemax] OF Boolean;

 BEGIN
  NoOfCliques := Graph.NoOfCliques;
  Graph.NoOfCliques := 0;
  mSubgraph.NoOfCliques := 0;
  nSubgraph.NoOfCliques := 0;
  mnSubgraph.NoOfCliques := 0;

    { 1. Scan through the cliques of the current graph in Graph. if a
      clique contains neither node of the edge, copy it as a clique of the
      new graph back to Graph. If it contains one of the nodes copy it as
      one of the temporary cliques to mSubgraph or nSubgraph setting the   
      associated flag true.}

  FOR i := 1 TO NoOfCliques DO
   BEGIN
    ThisClique := Graph.Cliques[i];
    IF m IN ThisClique THEN
     WITH mSubgraph DO
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := ThisClique;
       mCliqueOK[NoOfCliques] := true
      END
    ELSE IF n IN ThisClique THEN
     WITH nSubgraph DO
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := ThisClique;
       nCliqueOK[NoOfCliques] := true
      END
    ELSE
     WITH Graph DO
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := ThisClique
      END
   END;

    { 2.  Compare each clique from mSubgraph with each clique from         
      nSubgraph in turn, generate the clique created from the intersection 
      of the pair combined with the nodes m and n (It will, most often, be 
      simply the edge m-n) and put it into mnSubgraph, setting the         
      associated flag true. Compare both of the creating cliques with the  
      created one and if either is a subset of it, set the associated flag 
      false.  Compare the created one with all previous created ones;  if  
      it is a subset of any  of them remove it; if any of them is a subset 
      of it set their associated flags false. }

  FOR i := 1 TO mSubgraph.NoOfCliques DO
   BEGIN
    mClique := mSubgraph.Cliques[i];
    FOR j := 1 TO nSubgraph.NoOfCliques DO
     BEGIN
      nClique := nSubgraph.Cliques[j];
      WITH mnSubgraph DO
       BEGIN
        mnClique := mClique * nClique + [m, n];
        IF mClique <= mnClique THEN
         mCliqueOK[i] := false;
        IF nClique <= mnClique THEN
         nCliqueOK[j] := false;
        NoOfCliques := NoOfCliques + 1;
        Cliques[NoOfCliques] := mnClique;
        mnCliqueOK[NoOfCliques] := true;
        k := 1;
        WHILE NOT (mnClique <= Cliques[k]) DO
         BEGIN
          IF mnClique >= Cliques[k] THEN
           mnCliqueOK[k] := false;
          k := k + 1
         END;
        IF k <> NoOfCliques THEN
         NoOfCliques := NoOfCliques - 1
       END
     END
   END;

    { 3.  Scan mnSubgraph, mSubgraph and nSubgraph in turn, copying to    
      Graph each clique whose associated flag is still set true.}

  WITH Graph DO
   BEGIN
    FOR i := 1 TO mnSubgraph.NoOfCliques DO
     IF mnCliqueOK[i] THEN
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := mnSubgraph.Cliques[i]
      END;
    FOR i := 1 TO mSubgraph.NoOfCliques DO
     IF mCliqueOK[i] THEN
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := mSubgraph.Cliques[i]
      END;
    FOR i := 1 TO nSubgraph.NoOfCliques DO
     IF nCliqueOK[i] THEN
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := nSubgraph.Cliques[i]
      END
   END
 END;{of AddEdge }


 PROCEDURE DeleteEdge (VAR Graph : GraphType;
         m, n : NodeType);

 {Algorithm AS 247.3  Appl. Statist. (1989) Vol.38, No.2}

    { This procedure replaces the cliques of the current graph held in the 
      structure cliques of type GraphType by the cliques of the graph      
      obtained by deleting edge m-n.                                      }

  VAR
   Subgraph : GraphType;
   ThisClique : CliqueType;
   i, j, NoOfCliques, buffer : CliqueIndex;

 BEGIN
  NoOfCliques := Graph.NoOfCliques;
  Graph.NoOfCliques := 0;
  Subgraph.NoOfCliques := 0;

    { 1. Scan through the cliques of the current graph held in Graph. For
      each clique containing the edge m-n, store in Subgraph the two
      (possibly redundant) cliques obtained by deleting each node of the   
      edge. Copy each clique not containing the edge back to Graph as a    
      clique  of the new graph.                                           }

  FOR i := 1 TO NoOfCliques DO
   BEGIN
    ThisClique := Graph.Cliques[i];
    IF [m, n] <= ThisClique THEN
     WITH Subgraph DO
      BEGIN
       NoOfCliques := NoOfCliques + 2;
       Cliques[NoOfCliques - 1] := ThisClique - [m];
       Cliques[NoOfCliques] := ThisClique - [n];
      END
    ELSE
     WITH Graph DO
      BEGIN
       NoOfCliques := NoOfCliques + 1;
       Cliques[NoOfCliques] := ThisClique
      END
   END;

    { 2. Compare each clique in Subgraph with those directly placed in     
      Graph.  If it is a subset of any clique in Graph it is ignored,      
      otherwise it is added to Graph. The searching process uses a         
      variation of the technique of adding the item as a sentinel at the
      end of the sequence being searched to simplify the loop termination. 
      Here, since Graph is being extended but the searching takes place    
      only over its initial part, a buffer is left at the interface to be
      filled in just before exit.  }

  WITH Graph DO
   BEGIN
    NoOfCliques := NoOfCliques + 1;
    buffer := NoOfCliques;
    FOR i := 1 TO Subgraph.NoOfCliques DO
     BEGIN
      ThisClique := Subgraph.Cliques[i];
      Cliques[buffer] := ThisClique;
      j := 1;
      WHILE NOT (ThisClique <= Cliques[j]) DO
       j := j + 1;
      IF j = buffer THEN
       BEGIN
        NoOfCliques := NoOfCliques + 1;
        Cliques[NoOfCliques] := ThisClique
       END
     END;
    Cliques[buffer] := Cliques[NoOfCliques];
    NoOfCliques := NoOfCliques - 1
   END
 END;{of DeleteEdge }
