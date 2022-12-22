# Convert transcript coordinates to genomic coordinates in EpiNano output
# Example use:
# Rscript update_coordinates_batch.R results/epinano/2022-07-06/2C.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.csv 2022-06-14 B

library(GenomicRanges)
library(GenomicFeatures)
library(magrittr)
library(stringi)
library(Biostrings)

# Parameter setup
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
flairdate <- args[2]
cross <- args[3]

outfile <- sub("\\.csv$", ".genomic.csv", infile)

rangesloc <- paste0("results/parent_transcriptomes/", flairdate, "/")

stopifnot(cross %in% c("A", "B"))

if(cross == "A"){
  mother <- "SRR14654188"
  father <- "SRR14654189"
}
if(cross == "B"){
  mother <- "SRR14654186"
  father <- "SRR14654187"
}

gtffile <- paste0("results/flair/", flairdate, "/All_samples_collapse.isoforms.gtf")

fafile <- paste0("results/flair/", flairdate, "/All_samples_collapse.isoforms.fa")

# function to find indels
find_indels <- function(irl, irl_alt){
  txpts <- names(irl)
  out <- vector(mode = "list", length = length(txpts))
  names(out) <- txpts
  for(i in seq_along(txpts)){
    wd1 <- width(irl[[i]])
    wd2 <- width(irl_alt[[i]])
    wd_change <- wd2 - wd1
    keep <- wd_change != 0
    res <- IRanges(start = start(irl[[i]][keep]),
                   end = end(irl[[i]][keep]),
                   insert = wd_change[keep],
                   alt = as.character(irl_alt[[i]][keep]))
    names(res) <- names(irl[[i]])[keep]
    out[[i]] <- res
  }
  return(out)
}

# Import transcriptome sequence
txptome_seq <- readDNAStringSet(fafile)

# Import information on variant location within transcripts
load(paste0(rangesloc, mother, "_ranges_debug.RData"))

mother_indels <- find_indels(irl, irl_alt)

load(paste0(rangesloc, father, "_ranges_debug.RData"))

father_indels <- find_indels(irl, irl_alt)

# Import transcript annotation
gtf <- rtracklayer::import(gtffile)
gtf$transcript_id_fixed <- sub("_N[CW]$", "", gtf$transcript_id)

# Import EpiNano M6A calls
en_tab <- read.csv(infile)
en_tab$Parent <- sub("-.*$", "", en_tab$X.Ref)
temp <- sub("[^\\-]*-", "", en_tab$X.Ref)

# Regular expressions for splitting up transcript labels
readidpat <- "^[[:xdigit:]]{8}-([[:xdigit:]]{4}-){3}[[:xdigit:]]{12}(-[[:digit:]])?"
txptpat <- "^[NX][MR]_[[:digit:]]+\\.[[:digit:]]+(-[[:digit:]])?"
pseudogene_pat <- "^gene[[:digit:]]+(-[[:digit:]])?"

unknowngene_pat <- "N[CW]_[[:digit:]]+\\.[[:digit:]]+:[[:digit:]]+$"
knowngene_pat <- "GeneID:[[:digit:]]+$"
viralgene_pat <- "(DW|ABP|IAP|SB|CBP|VD)V(:[[:digit::]+)?$"
rRNA_pat <- "l-rRNA$"

en_tab$Transcript <- sub(unknowngene_pat, "", temp) %>%
  sub(knowngene_pat, "", .) %>%
  sub(viralgene_pat, "", .) %>%
  sub(rRNA_pat, "", .) %>%
  sub("_$", "", .)

en_tab$Gene <- sub(readidpat, "", temp) %>%
  sub(txptpat, "", .) %>%
  sub(pseudogene_pat, "", .) %>%
  sub("^_", "", .)

# Translate position in modified txpt to position in original txpt
en_tab$txpt_pos <- NA_real_

txpt_lookup_mother <- sub("_N[CW]$", "", names(mother_indels))
txpt_lookup_father <- sub("_N[CW]$", "", names(father_indels))

for(i in seq_len(nrow(en_tab))){
  if(en_tab$Parent[i] %in% c(mother, "Both")){
    theseranges <- mother_indels[[match(en_tab$Transcript[i], txpt_lookup_mother)]]
  } else {
    theseranges <- father_indels[[match(en_tab$Transcript[i], txpt_lookup_father)]]
  }
  if(length(theseranges) == 0){
    en_tab$txpt_pos[i] <- en_tab$pos[i]
    next
  }
  theseranges <- theseranges[order(start(theseranges))]
  mcols(theseranges)$end_mod <-
    end(theseranges) + cumsum(mcols(theseranges)$insert)
  mcols(theseranges)$start_mod <-
    mcols(theseranges)$end_mod - width(theseranges) - mcols(theseranges)$insert + 1L
  theseranges <- theseranges[mcols(theseranges)$start_mod < en_tab$pos[i]]
  # Modify position based on indels before it
  en_tab$txpt_pos[i] <- en_tab$pos[i] - sum(mcols(theseranges)$insert)
  # Is the position inside an indel? Get position by matching 5-mer.
  theseranges <- theseranges[mcols(theseranges)$end_mod > en_tab$pos[i]]
  if(length(theseranges) > 0){
    stopifnot(length(theseranges) == 1)
    altseq <- mcols(theseranges)$alt
    refseq <- as.character(txptome_seq[[paste(en_tab$Transcript[i], en_tab$Gene[i], sep = "_")]][theseranges])
    st5 <- en_tab$pos[i] - mcols(theseranges)$start_mod - 1L # start pos of 5-mer
    fivemer <- 
      substring(altseq, st5, st5 + 4)
    mtchR <- stri_locate_all_fixed(refseq, fivemer, opts_fixed=stri_opts_fixed(overlap = TRUE))[[1]][,1]
    # Continue with adjusting position only if there is a match
    if(!is.na(mtchR[1])){
      # If there are multiple matches, try to file the right one
      if(length(mtchR) > 1){
        mtchA <- stri_locate_all_fixed(altseq, fivemer, opts_fixed=stri_opts_fixed(overlap = TRUE))[[1]][,1]
        if(length(mtchR) == length(mtchA) && st5 %in% mtchA){
          mtchR <- mtchR[mtchA == st5]
        }
      }
      # If match is ambiguous, pick the closest one
      mtchR <- mtchR[which.min(abs(mtchR - st5))]
      # Remove the indel width and add back in width based on 5-mer position
      en_tab$txpt_pos[i] <- en_tab$txpt_pos[i] + mcols(theseranges)$insert +
        mtchR - st5
    }
  }
}

# Add kmers to table

txptome <- readDNAStringSet(paste0(rangesloc,
                                   "Cross", cross ,"_transcriptome.fa"))

en_tab$kmer <- sapply(seq_len(nrow(en_tab)),
                      function(i){
                        txptome[[en_tab$X.Ref[i]]][(en_tab$pos[i] - 2):(en_tab$pos[i] + 2)] %>%
                          as.character()
                      })

# Transcript coordinates to map back to genome
txpt_coords <- GRanges(en_tab$Transcript,
                       IRanges(en_tab$txpt_pos, width = 1),
                       strand = "*")

# Set up GRangesList from GTF
gtf1 <- gtf[gtf$type == "exon"]

gtf_list <- split(gtf1, gtf1$transcript_id_fixed)

# Map to genome
genome_coords <- mapFromTranscripts(txpt_coords, gtf_list)

stopifnot(identical(genome_coords$xHits, seq_len(nrow(en_tab)))) # TRUE

en_tab$Chromosome <- seqnames(genome_coords) %>% as.character()
en_tab$Genomic_pos <- start(genome_coords)
en_tab$Genomic_strand <- strand(genome_coords) %>% as.character()

# export
write.csv(en_tab, file = outfile, row.names = FALSE)
