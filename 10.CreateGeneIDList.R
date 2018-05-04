# Variable Clist contains 10 highly preserved modules
# Requires A1 and A2 directories already presernt in 'present working directory'

modulesA1 <- datExprA1g
modulesA1$module <- moduleColors

modulesA2 <- datExprA2g
modulesA2$module <- moduleColors

clist <- c("blue", "darkred", "saddlebrown", "darkorange", "red", "brown", "violet", "royalblue", "black", "green")
for (i in clist) {
	
   i <- noquote(i)
   A1path <- paste("A1/", paste(i, "A1", sep=""), sep="")
   A2path <- paste("A2/", paste(i, "A2", sep=""), sep="")
   write.table(rownames(modulesA1[grep(paste("^", i, "$", sep = ""),modulesA1$module),]), file=A1path, col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(rownames(modulesA2[grep(paste("^", i, "$", sep = ""),modulesA2$module),]), file=A2path, col.names=FALSE, row.names=FALSE, quote=FALSE)
}
