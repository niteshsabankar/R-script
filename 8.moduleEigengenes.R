PCs1A = moduleEigengenes(t(datExprA1g), colors=moduleColors)
ME_1A = PCs1A$eigengenes

PCs2A = moduleEigengenes(t(datExprA2g), colors=moduleColors)
ME_2A = PCs2A$eigengenes

# Keep module Eigengenes of only 10 highly preserved modules
ME_1A_preserved <- ME_1A[, c(2, 10, 26, 9, 24, 3, 33, 25, 1, 12)]
ME_2A_preserved <- ME_2A[, c(2, 10, 26, 9, 24, 3, 33, 25, 1, 12)]
