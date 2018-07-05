setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("~/Dropbox/Classes/BCB_II/BCB_II_HW_5/")
# setwd("c:/Users/Schuyler/Dropbox/Classes/BCB_II/BCB_II_HW_5/")
sapply(list.files("functions",pattern="*.R", full.name =TRUE),source,.GlobalEnv)

# only works easily on unix machines, if using windows 
# https://github.com/cjgb/rPython-win
# or run the Python function outside of R.
python.load("functions/fastq_to_fasta_function.py")
python.load("functions/fastq_adapter_barcode_trimmer.py")

fastq_file <- "files/SRR2241783_2.noN.fastq"
fasta_file <- python.call("convert_fastq", fastq_file)

# fasta_file <- "files/SRR2241783_2.fasta"

trim_lengths <- EM_trim_calc(fasta_file)
saveRDS(trim_lengths, "files/trim_lengths.Rds")
trim_lengths <- readRDS("files/trim_lengths.Rds")
# this doesn't work as-is yet because I cheated by removing some reads in the EM algorithm.. I haven't figured out yet how to deal with them...
python.call("trim_fastq", fastq_file, trim_lengths$MLE_trim[[2]])

sink("MLE_file.txt")
print(trim_lengths$delta)
print(trim_lengths$q[[2]])
print(trim_lengths$q[[1]])
sink()

sink("prediction_file.txt")
print(trim_lengths$MLE_trim[,c(1,3,4)])
sink()

write.table(trim_lengths$MLE_trim[,c(1,3,4)], file="prediction_file.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)
