# mamba create -n gr bioconductor-genomicranges "r-base>=4.0" \
# r-devtools "r-essentials>=4.0" r-dplyr r-furrr r-purrr r-foreach \
# r-doParallel r-data.table

# this version parses the peak information more carefully, looking at peak and summit together
print('Loading libraries')
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(furrr)
  library(foreach)
  library(doParallel)
  library(rtracklayer)
  library(purrr)
  library(GenomicRanges)
})


source_bed = './ATAC_peak_files/schMedS3h1_full_peak.tab'
liftover_dir = './schMedS3h1_based_full/'

sps = c(
  'schPol2',
  'schLug1',
  'schNov1'
)
atac_files = list(
  './ATAC_peak_files/schPol2_full_peak.tab',
  './ATAC_peak_files/schLug1_full_peak.tab',
  './ATAC_peak_files/schNov1_full_peak.tab'
)
names(atac_files) <- sps

################################################################################
# print what files will be used

print(paste0('working on these species:', sps))
print(paste0('source_bed = ', source_bed))
print(paste0('liftover_dir = ', liftover_dir))
print(paste0('ATAC-Seq peaks in the target species will be taken from these files: ', atac_files))

################################################################################

# this can be used as the main result file with the accounting of what happend to the peak
df_smed_peaks <- read.table(source_bed)
names(df_smed_peaks) <- c('chrom', 'start', 'end', 'name')

hed = c('chrom', 'start', 'end', 'name')

################################################################################

## settings

################################################################################

# this can be used as the main result file with the accounting of what happend to the peak
paste0('Reading summit files')
summit_liftover <- lapply(sps, function(sp){
  # load the summit lifovers into a list
  paste0('Reading: ', liftover_dir, sp, "_halLiftover_schMedS3h1_summits_50bp.bed")
  s1 <- fread(paste0(liftover_dir, sp, "_halLiftover_schMedS3h1_summits_50bp.bed"))
  names(s1) <- hed
  s1$stage <- 'summit'
  s1$species <- sp
  return(s1)
})
names(summit_liftover) = sps

hed = c('chrom', 'start', 'end', 'name')

span_threshold = 1000
# switch to run on subset for testing
nelements <- nrow(df_smed_peaks)
threads = 120 


options(future.globals.maxSize= 209715200000)

################################################################################

## big loop

################################################################################

print(paste0('Starting main iteration'))

smed_based_result <- list()

for (sp in sps) {
  print(paste0('Working on species: ', sp))
  ################################################################################
  ## generate result df
  res_df <- df_smed_peaks
  res_df$species <- sp
  
  ################################################################################
  ## get liftover bedfile
  print(paste0('Getting liftover bedfile: ', liftover_dir, sp, '_halLiftover_schMedS3h1_full_peak.bed'))
  bedfile <- paste0(liftover_dir, sp, '_halLiftover_schMedS3h1_full_peak.bed')
  # make empty genome for loading
  genome <- Seqinfo(genome = NA_character_)
  # read the bedfile as ranges object
  gr1 <- import(bedfile , genome = genome)
  ################################################################################
  ## annotate CRE with no hit in raw liftover
  res_df$rawhit_type <- ifelse(res_df$name %in% gr1$name, 'rawHit', 'no rawHit')
  # Split the GRanges object by the 'name' column
  gr_split <- split(gr1, gr1$name)
  ################################################################################
  ## use map and granges to merge each gr range by name
  # Merge ranges that are 100bp apart using `reduce`
  # filter short liftovers that remain
  plan(multicore, workers = threads)
  #plan("future::multisession")
  print(paste0('number of query elements with a liftover: ', length(gr_split)))
  print(paste0('Number of threads for parallel computation: ', threads))

  print(paste0('Split liftover into ranges for each query peak.'))
    # can be a future map
    lst <- future_map(1:length(gr_split), .progress = TRUE,
                      .options = furrr_options(stdout = FALSE, conditions = character()), function(i) {

      # split into range for each liftover peak
      gr <- gr_split[[i]] 
      # TODO consider the implications of the filter merge order
      # filter out short liftovers AFTER merging below
      #  <- gr[width(gr) > 21,]
      # intersect ranges
      ans <- GenomicRanges::reduce(gr, ignore.strand = TRUE, with.revmap = TRUE, 
                    min.gapwidth = 101L)
      # simplify the name of the merged segments to avoid duplications
      mcols(ans) <- aggregate(gr, mcols(ans)$revmap,
                              name.distinct = unstrsplit(unique(name), ","), drop = FALSE)
      # return the output
      ans <- ans[width(ans) > 21,]
      return(ans)
    })

  print(paste0('Split list created.'))
  # turn off multicore
  #plan(sequential)
  # keep the names for the next step. It gets lost of genes without a hit after filtering
  names(lst) <- names(gr_split)[1:length(gr_split)]
  
  
  ################################################################################
  ## bookkeeping after filtering
  print(paste0('Categorize liftover distribution for each query.'))
  res <- future_map(1:length(lst), .progress = TRUE,
                    .options = furrr_options(stdout = FALSE, conditions = character())
                    , function(j){
    current_gr <- lst[[j]]
    peak = names(lst)[j]
    current_gr_df <- as.data.frame(current_gr)
    # number of liftover hits
    nhits = nrow(current_gr_df)
    
    # if there are not hits
    if (nhits == 0) {
      hit_type = 'no hit after filtering'
    } else {
      # if there is a hit then
      # categorize the hit type
      sum <- current_gr_df |> 
        group_by(seqnames) |> 
        summarise(
          span = max(c(start, end)) - min(c(start, end)),
          nhits = n()
        )
      
      # Unclear why there are several answers
      hit_type <- unique(case_when(
        nrow(current_gr_df) == 0 ~ 'no_hit',
        nrow(current_gr_df) == 1 ~ 'one_hit',
        length(unique(sum$seqnames)) > 1 ~ 'fragmented several chromosomes',
        sum$span < span_threshold ~ paste0('fragmented span < ', span_threshold ),
        sum$span > span_threshold ~ paste0('fragmented span > ', span_threshold )
      ))
    }
    
    return(data.frame(name = peak, nhits = nhits, hit_type = hit_type))
    
  })


  print(paste0('Categorization done.'))
  res <- do.call(rbind, res)
  res_df <- left_join(res_df, res, by = 'name')


    ################################################################################
    ## intersect with the summit liftover and ATACSeq data and book keeping
    print(paste0('Intersect with the summit liftover and ATACSeq data'))

    gr_summit <- summit_liftover[sp]

    bed <- atac_files[ names(atac_files) == sp ][[1]]
    atac_gr <- import(bed, format = 'bed')
    
  res <- future_map(1:length(lst), 
                    .options = furrr_options(stdout = FALSE, conditions = character()),
                    .progress = TRUE,  function(j){
      current_gr <- lst[[j]]
      peak = names(lst)[j]
      # check if there is a summit liftover
      if (! peak %in% gr_summit[[1]]$name) {
        summit_type = 'no_liftover'
        # if there is at least one intersect 
        ans <- current_gr
        mcols(ans)$summit_OL_count <- rep(0, length(ans))
      } else {
        df <- gr_summit[[1]]
        df <- df[df$name == peak,]
        # bedtools_intersect("-a a.bed -b b.bed")
        current_summit_gr <- makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
        
        # merge close summits
        summit <- GenomicRanges::reduce(current_summit_gr, ignore.strand = TRUE, with.revmap = TRUE, 
                      min.gapwidth = 11L)
        # simplify the name of the merged segments to avoid duplications
        mcols(summit) <- aggregate(current_summit_gr, mcols(summit)$revmap,
                                name.distinct = unstrsplit(unique(name), ","), drop = FALSE)
        
        ans <- current_gr
        mcols(ans)$summit_OL_count <- countOverlaps(current_gr, summit, ignore.strand = TRUE)
        cnt <- sum(ans$summit_OL_count)
        # get the type of overlap
        summit_type = case_when(
          cnt == 0 ~ 'no_OL',
          cnt == 1 ~ 'one_OL',
          cnt > 1 ~ 'many_OL',
        )
      }
        # annotate the overlaps 
        # return peaks with OL
        ans2 <- ans
        mcols(ans2)$atac_OL_count <- countOverlaps(ans, atac_gr, ignore.strand = TRUE)
        cnt <- sum(ans$atac_OL_count)
        atac_type = case_when(
          cnt == 0 ~ 'no_OL',
          cnt == 1 ~ 'one_OL',
          cnt > 1 ~ 'many_OL',
        )

      
      return(list(
        df = data.frame(name = peak, summit_type = summit_type, atac_type = atac_type),
        bed = ans2))
  })

  print(paste0('Preparing output for ', sp))
  # todo do a nice conversion of the bed files
    dfs  <- lapply(res, '[[', 1)
    beds <- lapply(res, '[[', 2)


  # write the bed file
  print('Merging coordinate files')
  bed_df <- do.call(rbind, lapply(beds, data.frame))
  print('Done merging coordinate files.')
  df <- do.call(rbind, dfs)
    
    res_df <- left_join(res_df, df, by = 'name')

    ################################################################################
    ## output the result to list
    smed_based_result[[sp]] <- list(type = res_df, bed = bed_df)
}


################################################################################

## save output

################################################################################
print('Saving results.')
saveRDS(smed_based_result, 'schMedS3h1_based_full_result.rds')

