# Combine EpiNano results into a single table of methylation probabilities and
# a single table of coverage.  Only sites that had parent-specific transcripts
# in both crosses.

# Targets list and sanity checks ####
targets <- read.table("src/targets.txt", sep = "\t", header = TRUE)
targets <- targets[targets$Use,]
targets <- targets[order(targets$Cross, targets$Aggression),]
targets$Name[c(1:3, 6)] <- c("T3-Ag1:A2C", "T3-Ag1:A7C", "T3-Ag1:A9C", "T3-Ag2:A1T")
str(targets)

myfiles <- list.files("results/epinano/2022-07-06/", pattern = ".*genomic\\.csv$")
myfiles
# [1] "1Tb.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"      
# [2] "2C.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"       
# [3] "2T.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"       
# [4] "3T.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"       
# [5] "3Tb.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"      
# [6] "5C.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"       
# [7] "6C.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"       
# [8] "T2Ag1B8T.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv" 
# [9] "T3Ag1A2C.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv" 
# [10] "T3Ag1A7Ct.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"
# [11] "T3Ag1A9Ct.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"
# [12] "T3Ag2A1Tx.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv"

targets$File <- paste0(targets$Sample,
                       ".modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.genomic.csv")

setequal(myfiles, targets$File) # TRUE

testfile <- read.csv(paste0("results/epinano/2022-07-06/", myfiles[1]))
str(testfile)
table(testfile$Parent)
#  Both SRR14654188 SRR14654189 
# 17261       22590       28894

# Read all 12 files and filter ####
site_tables <-
  sapply(targets$File,
         function(x){
           y <- read.csv(paste0("results/epinano/2022-07-06/", x))
           y <- y[y$Parent != "Both",]
         },
         simplify = FALSE)

sapply(site_tables, function(x) table(x$Parent), simplify = FALSE)
head(site_tables[[12]])
table(site_tables[[12]]$Parent)

# Get all sites for each cross
cross_A_sites <- unique(do.call(rbind, lapply(site_tables[1:6],
                                              function(x) x[,c("Transcript", "Gene", "Chromosome",
                                                               "Genomic_pos", "Genomic_strand")])))

head(cross_A_sites)

cross_B_sites <- unique(do.call(rbind, lapply(site_tables[7:12],
                                              function(x) x[,c("Transcript", "Genomic_pos")])))
# 44,635 for cross A and 101,751 for cross B

sites_both_crosses <-
  cross_A_sites[paste(cross_A_sites$Transcript, cross_A_sites$Genomic_pos) %in%
                  paste(cross_B_sites$Transcript, cross_B_sites$Genomic_pos),]
# 35,782 sites

sites_paste <- paste(sites_both_crosses$Transcript, sites_both_crosses$Genomic_pos)

# Build matrices ####
mycolnames <- paste0(rep(targets$Name, each = 2), c("_father", "_mother"))

probM_mat <- matrix(NA_real_, nrow = length(sites_paste), ncol = length(mycolnames),
                    dimnames = list(sites_paste, mycolnames))

cov_mat <- matrix(0L, nrow = length(sites_paste), ncol = length(mycolnames),
                    dimnames = list(sites_paste, mycolnames))

for(i in seq_along(site_tables)){
  tab <- site_tables[[i]]
    if(targets$Cross[i] == "A"){
    mother <- "SRR14654188"
    father <- "SRR14654189"
  } else {
    mother <- "SRR14654186"
    father <- "SRR14654187"
  }
  tab_mother <- tab[tab$Parent == mother,]
  tab_father <- tab[tab$Parent == father,]
  tab_mother$sites <- paste(tab_mother$Transcript, tab_mother$Genomic_pos)
  tab_father$sites <- paste(tab_father$Transcript, tab_father$Genomic_pos)
  sites_mother_common <- intersect(sites_paste, tab_mother$sites)
  sites_father_common <- intersect(sites_paste, tab_father$sites)
  rows1m <- match(sites_mother_common, sites_paste)
  rows2m <- match(sites_mother_common, tab_mother$sites)
  rows1f <- match(sites_father_common, sites_paste)
  rows2f <- match(sites_father_common, tab_father$sites)
  probM_mat[rows1m, paste0(targets$Name[i], "_mother")] <-
    tab_mother$ProbM[rows2m]
  probM_mat[rows1f, paste0(targets$Name[i], "_father")] <-
    tab_father$ProbM[rows2f]
  cov_mat[rows1m, paste0(targets$Name[i], "_mother")] <-
    tab_mother$cov[rows2m]
  cov_mat[rows1f, paste0(targets$Name[i], "_father")] <-
    tab_father$cov[rows2f]
}

# Export ####
write.table(cbind(sites_both_crosses, probM_mat), sep = "\t", row.names = FALSE,
            file = "results/epinano/2022-07-06/probM_matrix_2022-08-03.txt")
write.table(cbind(sites_both_crosses, cov_mat), sep = "\t", row.names = FALSE,
            file = "results/epinano/2022-07-06/coverage_matrix_2022-08-03.txt")
