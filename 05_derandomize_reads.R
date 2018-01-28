library(tidyverse)
library(stringdist)

# Get the keystreams from PDF and prepend 0
keystreams <- c("0002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221",
                "0202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022",
                "0221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100",
                "0100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020")

# Function to run provided perl script
runPerlscript <- function(DNA, keystream){
  command <- str_c("./keystream_randomization.pl ",
                   DNA,
                   " ",
                   keystream,
                   sep = "")
  
  result <- system(command, intern = T)
  
  return(result)
}

seqs <- read_tsv("DNA/IDcheck_out.tsv") %>%
  mutate(mod4 = index %% 4) %>%
  mutate(keystream = keystreams[mod4 + 1]) %>% 
  mutate(decodedDNA = map2_chr(code, keystream, ~ runPerlscript(.x, .y))) %>%
  mutate(decodedDNA = if_else(orientation == "REVERSE", # Reverse complement of reverse orientation
         stringi::stri_reverse(chartr("ACTG", "TGAC", decodedDNA)),
         decodedDNA))

seqs %>%
  select(seqname, decodedDNA) %>%
  write_tsv("DNA/all_decoded_seqs.tsv", col_names = F)

# Plot
seqs %>%
  ggplot(aes(x = ID)) +
  geom_bar()

seqs %>%
  select(ID, index, count, decodedDNA ) %>%
  write_tsv("seqs.tsv")


