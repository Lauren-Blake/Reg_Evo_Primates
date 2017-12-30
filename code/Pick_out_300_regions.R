y <- c(100,300,500,530,560,570,590,1000,1200,1300,2000)
y2 <- c(y[1],y[-length(y)])
distance <- y-y2

# annotate regions
region_index <- which(distance > 300)
region_array <- data.frame(starting=c(1,region_index),
                           ending=c(region_index-1,length(y)))

# take out region with less than 3 CpGs
region_array_subset <- region_array[with(region_array, which((ending-starting+1)>2)), ] 

region_list <- vector("list",nrow(region_array_subset))

# Pull the average methylation for each region
for (index in 1:length(region_list)) {
  region_tmp <- unlist(region_array_subset[index,])
  region_tmp2 <- c(region_tmp[1]:region_tmp[2])
  region_list[[index]] <- mean(y[region_tmp2])
}



# for (index3 in 1:length(region_list)) {
#   region_index <- distance[]
#   for (index in seq_along(region_index)) {
#     region_dist <- distance[1:(region_index[index]-1)]
#     
#     for (index2 in seq_along(region_dist)) {
#       dd <- distance[index2]  
#       if (dd < 300) {
#         region_list[[index3]] <- c(region,index2)
#       } 
#     }
#   }
# }



# Read the orth CpG file

orth_cpg <- read.table("../data/dfCovnoX_hg19.bed")

y <- c(orth_cpg$V2)
y2 <- c(y[1],y[-length(y)])
distance <- y-y2
