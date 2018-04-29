softPower = 18
adjacencyA1 = adjacency(t(datExprA1g),power=softPower,type="signed");
diag(adjacencyA1)=0
TOMA1 = TOMsimilarity(adjacencyA1, TOMType="signed")
dissTOMA1 = 1 - TOMA1
geneTreeA1 = flashClust(as.dist(dissTOMA1), method="ward")

adjacencyA2 = adjacency(t(datExprA2g),power=softPower,type="signed");
diag(adjacencyA2)=0
TOMA2 = TOMsimilarity(adjacencyA2, TOMType="signed")
dissTOMA2 = 1 - TOMA2
geneTreeA2 = flashClust(as.dist(dissTOMA2), method="ward") 
