#'@title doEnrich
#'@description given a list of genes, this function combines topological properties of the disease ontology structure for enrichment analysis.
#'@author Haixiu Yang
#'@param interestGenes a vector of gene IDs.
#'@param test One of "fisherTest","hypergeomTest","binomTest","chisqTest" and "logoddTest" statistical model. Default is hypergeomTest.
#'@param method One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr" and "none",for P value correction.
#'@param m Set the maximum number of ancestor layers for ontology enrichment. Default is layer 1.
#'@param maxGsize indicates that doterms with more annotation genes than maxGsize are ignored, and the P value of these doterms is set to 1.
#'@param minGsize indicates that doterms with less annotation genes than minGsize are ignored, and the P value of these doterms is set to 1.
#'@param traditional a logical variable, TRUE for traditional enrichment analysis, FALSE for enrichment analysis with weights. Default is FALSE.
#'@param delta Set the threshold of nodes, if the p value of doterm is greater than delta, the nodes are not significant, and these nodes are not weighted.
#'@param result_do Receives the file output by the wrireResult function, which is used to visually display the enrichment results (without running the enrichment operation again).
#'@param penalize Logical value, whether to add a penalty to the node.Adding a penalty will look for nodes with more branches.
#'@return A \code{enrich} instance.
#'@importFrom magrittr `%>%`
#'@importFrom purrr map pwalk map_dbl walk2 map2 map_int
#'@importFrom dplyr mutate filter select arrange
#'@importFrom stringr str_c
#'@importFrom BiocGenerics intersect
#'@export
#'@examples
#'#The enrichment results were obtained by using demo.data
#'demo.data=c(1636,351,102,2932,3077,348,4137,54209,5663,5328,23621,3416,3553)
#'demo_result<-doEnrich(interestGenes=demo.data)
#'
#'#setting the penalty to FALSE, the algorithm can mitigate the extent of reduced nodes' weights
#'penalF_demo<-doEnrich(interestGenes=demo.data,penalize = FALSE)
#'
#'#Statistical models and P-value correction can be set
#'demo_result2<-doEnrich(demo.data,test="hypergeomTest",method="holm")
#'
#'#Using the traditional enrichment analysis method.
#'Tradition_demo<-doEnrich(demo.data,traditional=TRUE)
#'
#'#Between the number of genes in minGsize and maxGsize doterm enrichment analysis
#'demo_result3<-doEnrich(demo.data,minGsize=5,maxGsize=500)
#'
#'#Draw from wrireResult output files
#'#Firstly, read the wrireResult output file,using the following two lines
#'data<-read.delim(file.path(system.file("examples", package = "EnrichDO"),"result.txt"))
#'doEnrich(result_do = data)
#'#then, Use the drawing function you need
#'drawGraphViz(enrich)    #Tree diagram
#'drawPointGraph(enrich)  #Bubble diagram
#'drawBarGraph(enrich)    #Bar plot



#主要函数
doEnrich <- function(interestGenes,test="hypergeomTest",method="BH",m=1,maxGsize=5000,minGsize=5,traditional=FALSE,delta=0.01,result_do=NULL,penalize=TRUE) {
  if(!is.null(result_do)){
	  TermStruct(result_do=result_do)
  return("Now you can use the drawing function")
  }
  #初始化
  init(traditional)
  interestGenes<-intersect(interestGenes,dotermgenes)
  enrich <<- enrich %>%
    mutate(cg.arr=map(gene.arr, intersect, interestGenes)) %>%
    mutate(cg.len=map_int(cg.arr, length)) %>%
    mutate(ig.len=length(interestGenes))
  #2024-1-18提取出和兴趣集因有交集,且满足minGsize和maxGsize的DOID
  current_enrich<-filter(enrich,cg.len!=0,gene.len>=minGsize,gene.len<=maxGsize)
  DOID_C<<-current_enrich$DOID

  #2024-1-18
  if(traditional==TRUE){
    pwalk(current_enrich,function(DOID,p,gene.arr,gene.w,...){
      p<-Test(test,interestGenes ,gene.arr,gene.w)
      enrich[enrich$DOID==DOID,]$p<<-p
    })
  }else{
    #2024-1-18
    enrich$child.arr<<-sapply(enrich$child.arr,function(x){intersect(x,DOID_C)})
    enrich$child.len<<-sapply(enrich$child.arr,function(x){length(x)})

    #从叶子节点向父节点逐级计算
    for(i in max(enrich$level):m) {
      #获取当前层
      #2024-1-18,选择满足DOID_C，minGsize,maxGsize的DOID
      current.level.terms <- enrich %>% filter(level==i,DOID%in%DOID_C)
      #2024-1-18记录每层的节点数和注释基因数
      level_DOID<-length(current.level.terms$DOID)
      level_gene<-length(unique(unlist(current.level.terms$gene.arr)))
      cat(str_c("LEVEL: ", i,"\t",level_DOID," nodes\t",level_gene," genes to be scored\n"))
      #迭代每一行
      pwalk(current.level.terms, function(DOID, gene.arr, gene.w, parent.arr, child.arr, ...) {
        #计算term相似度
        computeTermSig(interestGenes, i, DOID, gene.arr, gene.w, parent.arr, child.arr,test ,traditional,delta,penalize)
      })

    }
  }

  enrich <<- enrich %>% arrange(p)
  enrich <<- mutate(enrich, p.adjust=p.adjust(p, method=method))
  print(method)
  print(test)
  return(enrich)
}
