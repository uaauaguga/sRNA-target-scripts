#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-t", "--tree", required=TRUE, help="input tree file")
parser$add_argument("-i", "--input", required=TRUE, help="input pvalues of tips")
parser$add_argument("-o", "--output", required=TRUE, help="output reconstruction")
parser$add_argument("-m", "--method", default="GLS", help="reconstruction method to use")
args <- parser$parse_args()


suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("phytools"))
epsilon <- 0.000001 # prevent numeric error in probit/logit

message(paste('load tree from ',args$tree ,' ...'))
tree <- ape::read.tree(args$tree)
# drop edges with tip length of zero
zero.length.edge <- tree$edge[tree$edge.length==0,]
tip.ids.to.drop <- zero.length.edge[duplicated(zero.length.edge[,1]),2]
tree <- drop.tip(tree,tree$tip.label[tip.ids.to.drop])

# load scores of pairwise sciring 
message(paste('load scores from ',args$input ,' ...'))
pvalues <- list()
fin  <- file(args$input, open = "r")

header <- readLines(fin, n = 1, warn = FALSE)
while (length(oneLine <- readLines(fin, n = 1, warn = FALSE)) > 0) {
  fields <- unlist(strsplit(oneLine, "\t"))
  seq.id.1 <- fields[1]
  seq.id.1 <- unlist(strsplit(seq.id.1,":"))
  sRNA.id <- paste(seq.id.1[1:length(seq.id.1)-1],collapse=":")
  genome.id <- seq.id.1[length(seq.id.1)]
  seq.id.2 <- fields[4] 
  seq.id.2 <- unlist(strsplit(seq.id.2,":"))
  target.id <- paste(seq.id.2[1:length(seq.id.2)-1],collapse=":")
  #RF00001-5S_rRNA::NC_016810.1:3588291-3588407(-):GCF_000006945.2 3       16      WP_000002386.1:GCF_000006945.2  219     232     -7.57165        0.4210257338506649      7       8
  pvalue <- as.numeric(fields[8])
  if( is.null(pvalues[[sRNA.id]])){
    pvalues[[sRNA.id]] <- list()
  }
  if( is.null(pvalues[[sRNA.id]][[target.id]])){
    pvalues[[sRNA.id]][[target.id]] <- c()
  }
  pvalues[[sRNA.id]][[target.id]][genome.id] <- pvalue
} 

close(fin)

# performance ancestral score reconstruction
# sRNA -> target -> genome,score
cat("genome.id","sRNA.id", "target.id","score","branch.length","\n",sep="\t",file=args$output,append=FALSE)
for (sRNA.id in names(pvalues)){
  #print(sRNA.id)
  for(target.id in names(pvalues[[sRNA.id]])){
    #print(target.id)
    scores <-  pvalues[[sRNA.id]][[target.id]]
    scores[scores<epsilon] <- epsilon
    scores[scores>(1-epsilon)] <- 1 - epsilon
    #tip.scores <- log((1-scores)/scores) # logit transformation
    tip.scores <- qnorm(1-scores) # probit transformation
    tip.scores <- tip.scores[!is.na(tip.scores)]
    if (length(scores)>=3){
      message(paste("processing",sRNA.id,target.id,"..."))
      genome.ids <- names(tip.scores)
      genome.ids <- genome.ids[genome.ids %in% tree$tip.label]
      subtree <- drop.tip(tree, tree$tip.label[-match(genome.ids, tree$tip.label)])
      subtree <- midpoint.root(subtree)   
      tip.scores <- tip.scores[subtree$tip.label]
      tryCatch({
          if (args$method == "GLS"){
           res <- ace(tip.scores, subtree, method="GLS",corStruct=corBrownian(1, subtree))
          }else{
          res <- ace(tip.scores, subtree, method="pic")
          }
         },error=function(c){
         message( paste("Wrong with",sRNA.id,target.id,collapse=" "))
         print(c)
     })
    subtree <- reorder(subtree,order="pruningwise")
    cumulative.branch.lengths <- list()
    ancestral.scores <- list()
    for(i in seq(nrow(subtree$edge))){
      parent <- subtree$edge[i,1]
      child <- subtree$edge[i,2]
      parent.s <- as.character(parent)
      child.s <- as.character(child)
      edge.length <- subtree$edge.length[i]   
      if(child<=length(subtree$tip.label)){
        # tip        
        cumulative.branch.lengths[[child.s]] <- list(c(0))
        tmp <- tip.scores[subtree$tip.label[child]]
        #names(tmp) <- NULL
        ancestral.scores[[child.s]] <- list(c(tmp))
      }
      # should consider two child
      if(! parent.s %in% names(cumulative.branch.lengths)){
        cumulative.branch.lengths[[parent.s]] <- list()
        ancestral.scores[[parent.s]] <- list()
      }
    for(i in seq(length(cumulative.branch.lengths[[child.s]]))){      
        last.branch.length <- cumulative.branch.lengths[[child.s]][[i]]
        last.branch.length <- last.branch.length[length(last.branch.length)]
        cumulative.branch.lengths[[parent.s]][[length(cumulative.branch.lengths[[parent.s]])+1]] <- c(cumulative.branch.lengths[[child.s]][[i]], last.branch.length + edge.length)
        
    }  
    cumulative.branch.lengths[[child.s]] <- NULL  # delete child info as they will be never used 
    for(i in seq(length(ancestral.scores[[child.s]]))){
        tmp <- res$ace[parent.s]
        #names(tmp) <- NULL
        ancestral.scores[[parent.s]][[length(ancestral.scores[[parent.s]])+1]] <- c(ancestral.scores[[child.s]][[i]],tmp)
    }
    ancestral.scores[[child.s]] <- NULL # delete child info as they will be never used 
    }
    # save results of current pair
    root.node <- setdiff(subtree$edge[,1],subtree$edge[,2])
    root.node <- as.character(root.node)
    for(i in seq(length(ancestral.scores[[root.node]]))){
      genome.id <- names(ancestral.scores[[root.node]][[i]])[1]
      current.ancestral.scores <- paste(as.character(round(ancestral.scores[[root.node]][[i]],5)),collapse=",")
      current.cumulative.branch.lengths <- paste(as.character(round(cumulative.branch.lengths[[root.node]][[i]],5)),collapse=",")
      line <- sprintf("%s\t%s\t%s\t%s\t%s\n",genome.id,sRNA.id,target.id,current.ancestral.scores,current.cumulative.branch.lengths)    
      cat(line,file=args$output,append=TRUE)
     }
   }
  }
}







