
# conda activate R

infile <- 'resources/GTEx_CPM/CPM_full.RData'
outfold <- 'resources/GTEx_CPM/'

load(infile)

for (n in names(CPM.all)) { 
    print(n)
    outfile <- paste(outfold, gsub(' ','_',n), ".txt" , sep = "")
    write.table(data.frame(CPM.all[[n]]),  file = outfile, quote = FALSE, sep ="\t" ,row.names = TRUE, col.names = TRUE)
}
