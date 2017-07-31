#' Create A Manhattan Plot
#' 
#' This function intakes GWAS summary statsitcs and produces a Manhattan Plot.
#' @param sumstat_data_path Summary statistics data path
#' @param clump_data_path Plink clumped data path
#' @param output_filepath Where you will output the Manhattan Plot
#' @param rsID The title of your SNP ID data column in GWAS data (include quotation marks)
#' @param Chr The title of your Chromosome data column in GWAS data (include quotation marks)
#' @param P.value The title of your p-value data column in GWAS data (include quotation marks)
#' @param position The title of your base pair position data column in GWAS data (include quotation marks)
#' @param clump_rsID The title of your SNP ID data column in clumped data (include quotation marks)
#' @param clump_Chr The title of your Chromosome data column in clumped data (include quotation marks)
#' @param clump_P.value The title of your p-value data column in clumped data (include quotation marks)
#' @param clump_position The title of your base pair position data column in clumped data (include quotation marks)
#' @param delimiter The way your data is separated (often with "/t" or ",")
#' @param title What you would like the title to be
#' @export
#' @examples
#' ssgac_manhattan(sumstat_data_path = "path/sumstata.txt" , 
#' output_filepath = "outpath", 
#' rsID="snpid", Chr="chr", 
#' P.value = "pval", position="bppos",
#' clump_rsID="snpid", clump_Chr="chr", 
#' clump_P.value = "P", clump_position="pos",
#' delimiter = "/t", 
#' title="ssgac_manhattanplot"))
#' 
#' 
## Manhattan Plot Function, modified from Richard's code by Tuan, Tushar and Rosie
## Created on July 24th, 2017 and updated on July 27th, 2017

##Change title.main, input file & ... Remember that it is case-sensitive, check that headers match code! E.g. CHR vs Chr
##Check column headers and contents with read.table(nrow = 2, <file>)

ssgac_manhattan <-function(sumstat_data_path, clump_data_path, output_filepath, rsID, Chr, P.value, position, clump_rsID, clump_Chr, clump_P.value, clump_position, delimiter, title) {
  
##Input validation
  if (missing(sumstat_data_path))
    stop("Need to specify summary statistics filepath and filename to be plotted, e.g. \"../INPUT/sumstats1.txt\" .")

  if (missing(clump_data_path))
    stop("Need to specify clump data filepath and filename to be plotted, e.g. \"../INPUT/sumstats1.txt\" .")
  
  if (missing(rsID))
    stop("Need to specify header of rsID column within double quotes, e.g.\"rsID\"")
  if (is.numeric(rsID))
    stop("Need to specify the header of the rsID column as a character string.")
  
  if (missing(Chr))
    stop("Need to specify header of Chromosome column within double quotes, e.g.\"chr\"")
  if (is.numeric(Chr))
    stop("Need to specify the header of the Chromosome column as a character string.")
  
  if (missing(P.value))
    stop("Need to specify header of p-value column within double quotes, e.g. \"Pval\"")
  if (is.numeric(P.value))
    stop("Need to specify the header of the p-value column as a character string.")
  
  if (missing(position))
    stop("Need to specify header of position column within double quotes, e.g. \"bppos\"")
  if (is.numeric(position))
    stop("Need to specify the header of the position column as a character string.")
  
  if (missing(clump_rsID))
    stop("Need to specify header of the clumped rsID column within double quotes, e.g.\"rsID\"")
  if (is.numeric(clump_rsID))
    stop("Need to specify the header of the clumped rsID column as a character string.")
  
  if (missing(clump_Chr))
    stop("Need to specify header of the clumped Chromosome column within double quotes, e.g.\"chr\"")
  if (is.numeric(clump_Chr))
    stop("Need to specify the header of the clumped Chromosome column as a character string.")

  if (missing(clump_P.value))
    stop("Need to specify header of the clumped p-value column within double quotes, e.g. \"Pval\"")
  if (is.numeric(clump_P.value))
    stop("Need to specify the header of the p-value column as a character string.")
  
  if (missing(clump_position))
    stop("Need to specify header of the clumped position column within double quotes, e.g. \"bppos\"")
  if (is.numeric(clump_position))
    stop("Need to specify the header of the clumped position column as a character string.")
  
  if (missing(delimiter))
    stop("Need to specify the field delimiter within double quotes, e.g. \",\" or \" \".")
  
  if (missing(output_filepath))
    stop("Need to specify the output filepath and prefix withing double quotes, e.g. \"../OUTPUT/qq_plot_GWAS.jpeg\".")
  
  if (missing(title))
    stop("Need title for the Manhattan plot, e.g. \"My Manhattan Plot\"")
  
## Open necessary ggplot2 library 
require(ggplot2)
library(ggplot2)

## Read GWAS sumstat file and read in only needed columns (decrease run time by ~35% in test)
header_row_temp_var = read.table(sumstat_data_path, sep=delimiter, header = TRUE, nrow=1)
temp_col_num <- rep(0,4)
temp_varname <- c(rsID,position,P.value,Chr)
temp_col_num <-which(colnames(header_row_temp_var) %in% temp_varname)

gwas_colclass_variable <- rep(0, length(header_row_temp_var))
for(i in 1:length(header_row_temp_var)){
  if(i == temp_col_num[1]) {
    gwas_colclass_variable[i] <- c("character")
  } 
  else if (i %in% temp_col_num) {
    gwas_colclass_variable[i] <- c("numeric")
  } else {gwas_colclass_variable[i] <- c("NULL")}
}
GWASDATA = read.table(sumstat_data_path, sep=delimiter, colClasses = gwas_colclass_variable, header = TRUE)

#Read in plink.clumped file and read in only needed columns
header_row_temp_var = read.table(clump_data_path, sep=delimiter, header = TRUE, nrow=1)
temp_col_num <- rep(0,4)
temp_varname <- c(clump_rsID,clump_position,clump_P.value,clump_Chr)
temp_col_num <-which(colnames(header_row_temp_var) %in% temp_varname)

clump_colclass_variable <- rep(0, length(header_row_temp_var))
for(i in 1:length(header_row_temp_var)){
  if(i == temp_col_num[1]) {
    clump_colclass_variable[i] <- c("character")
  } 
  else if (i %in% temp_col_num) {
    clump_colclass_variable[i] <- c("numeric")
  }
  else {clump_colclass_variable[i] <- c("NULL")}
}
clump_data = read.table(clump_data_path, sep=delimiter, colClasses = clump_colclass_variable, header = TRUE)


##Loop for calculation of BP plotting positions conditional on CHR used in -> ##Plot GWAS data
nchr <- length(unique(GWASDATA[[Chr]])) #Count number of chr in data
nchr2 <- length(unique(clump_data[[clump_Chr]])) #Count number of chr in clump_data
plottemp <- rep(0, nchr) #Create temp variable

##Calculation of position for chromosome labels
bptemp <- rep(0, nchr) #create empty temp vector
bppos <- rep(0, nchr)  #creat empty vector for calculated positions

#loop
for(i in 1:nchr) {
  #finds the highest BP for each chromosome
  plottemp[i] <- max(GWASDATA[[position]][GWASDATA[[Chr]] == i]) 
  #if chr==1 then the positions are the original bp positions
  if(i==1) {GWASDATA$plotpos <- GWASDATA[[position]] #for gwasdata (SNP) plotting position
  bppos[i] <- plottemp[i]/2} #For Chromosome label position
  #ifelse then we add the sum of the highest BP of all previous chr and add it to the position for actual chr
  else {GWASDATA$plotpos[GWASDATA[[Chr]]==i] <- GWASDATA[[position]][GWASDATA[[Chr]]==i] + sum(plottemp[1:i-1]) #for gwasdata (SNP) plotting position
  bppos[i] <- (plottemp[i]/2 + sum(plottemp[1:i-1]))} #For Chromosome label position
}

#loop for labels in clump_data see -> ##Add rsID labels for significant AND ingwasdataendent observations from clump_data
######## Below Commented out by Rosie for now I don't have clumped data
for(i in 1:nchr2){
  current_chr <- unique(clump_data[[clump_Chr]])[i]
  ifelse(current_chr==1, clump_data$plotpos[clump_data[[clump_Chr]] == 1] <- clump_data[[clump_position]][clump_data[[clump_Chr]] == 1], #if chr==1 then the positions are the original bp positions
         clump_data$plotpos[clump_data[[clump_Chr]] == current_chr] <- clump_data[[clump_position]][clump_data[[clump_Chr]] == current_chr] + sum(plottemp[1:current_chr-1]) #ifelse then we add the sum of the highest BP of all previous chr and add it to the position for actual chr
  )
}

##Loop for calculation of lead snp getting different pch (point type)
#create cph = 19 // col = "royalblue1" vectors
GWASDATA$pointtype <- rep(19, length(GWASDATA[[P.value]]))
GWASDATA$colortype <- rep("royalblue1", length(GWASDATA[[P.value]]))
GWASDATA$cex_size <- rep(3, length(GWASDATA[[P.value]]))
#fill color vector with alternating color
GWASDATA$colortype <- c("royalblue1", "skyblue1", "royalblue1", "skyblue1", "royalblue1", 
                        "skyblue1", "royalblue1", "skyblue1", "royalblue1", "skyblue1", 
                        "royalblue1", "skyblue1", "royalblue1", "skyblue1", "royalblue1", 
                        "skyblue1", "royalblue1", "skyblue1", "royalblue1", "skyblue1", 
                        "royalblue1", "skyblue1") [GWASDATA[[Chr]]]

#loop for changing pch and color for lead snp and saving to vector pointtype / colortype
for(i in 1:length(unique(clump_data[[clump_rsID]]))){
  temp_clumppos <- clump_data[[clump_position]][i]
  GWASDATA$pointtype[GWASDATA[[position]] == temp_clumppos] <- 4
  GWASDATA$colortype[GWASDATA[[position]] == temp_clumppos] <- "red"
  GWASDATA$cex_size[GWASDATA[[position]] == temp_clumppos] <- 5
}

## Create variables and values that will be used in ggplot2
GWASDATA$y<- -log10(GWASDATA[[P.value]])
clump_data$yclump <- -log10(clump_data[[clump_P.value]])
clump_data$xclump <- -log10(clump_data[[clump_position]])
clump_data$cross  <- 4
# Reference line
hline1 <- data.frame(yi=5)
hline2 <- data.frame(zi=-log10(5*10^-8))
# Title string
plottitle <- unlist(strsplit(title_var, split='.', fixed=TRUE))[1]
# y axis value range
limityaxis <- max(-log10(GWASDATA[[P.value]])) + 5

# calculate mean chisq
z_temp <- qnorm(GWASDATA$P.value/2)  #before plotting

# add legends for SNPs and pvalues
leadsnp <- paste(length(unique(clump_data[[clump_rsID]])))
snplegend <- paste("Lead SNPs"," (", leadsnp, " total)", sep="")
pval5 <- paste("p-value = ",paste(5),"x",paste(10^-8))
pval8 <- paste("p-value = ",paste(1),"x",paste(10^-5))

# add annotate text for chi-squared
GWASDATA$z_temp <- qnorm(GWASDATA[[P.value]]/2)  # calculate mean chisq
summary(GWASDATA$z_temp)
chix <- max(GWASDATA$plotpos[GWASDATA[[Chr]] == 2])
chiy <- max(-log10(GWASDATA[[P.value]])) + 3
chisymbol <- paste(expression(chi))
chiresult <- mean(GWASDATA$z_temp^2)
chilegend <- sprintf("\"Mean\"~chi^2 == %0.2f", chiresult)

## ggplot 2
man <- ggplot(GWASDATA, aes(plotpos, y)) 
man <- man + geom_point(data = GWASDATA, shape = GWASDATA$pointtype, 
                        color = GWASDATA$colortype, size = GWASDATA$cex_size, 
                        stroke=0.5,show.legend = FALSE) 
man <- man + geom_point(data = clump_data, aes(x=plotpos,y=yclump), color="red", shape=4,
                        stroke=2, size=5, show.legend = FALSE) 
# for showing point (red cross) legends
man <- man + geom_point(data=clump_data,aes(color=snplegend), size=-100, stroke=2, show.legend=TRUE)
# for showing horizontal reference lines
man <- man + geom_hline(data=hline1, aes(yintercept = yi), linetype="dotted", size=1.5, show.legend = FALSE) 
man <- man + geom_hline(data=hline2, aes(yintercept = zi), linetype="longdash", size=1.5, show.legend = FALSE) 
# for showing line legends
man <- man + geom_hline(data=hline2, aes(yintercept = zi, color=pval8), linetype="longdash", size=0, show.legend = TRUE)
man <- man + geom_hline(data=hline1, aes(yintercept = yi, color=pval5), linetype="dashed", size=0, show.legend = TRUE)
# Specify legend name (null), color, shape, size and line type
man <- man + scale_color_manual("",values=c("red","black","black"),guide=guide_legend(override.aes = list(shape=c(4,NA,NA),linetype=c(0,3,5), size=1)))
# Specify axis labels, titles  
man <- man + scale_y_continuous(expand = c(0,0), limits = c(0, limityaxis)) 
man <- man + scale_x_continuous(expand = c(0,0), 
                          breaks = bppos[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)],
                          label = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,'\n16',17,'\n18',19,'\n20',21,'\n22')) 
man <- man + labs(x = "Chromosome", y = expression(paste(-log[10],"(",italic(p),"-value)")))
# Add chi-squared annotation
man <- man + annotate("text", x = chix, y = chiy, label = as.character(chilegend), size=9, parse=T)
# Title and format options
man <- man + ggtitle(plottitle) 
man <- man + theme(plot.title = element_text(hjust=0.5, vjust=-5, lineheight=3, face="bold", color="black", size=30),
             axis.ticks.x=element_line(c(rep(0,15),1,0,1,0,1,0,1)), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.line = element_line(colour = "black"),
             legend.justification=c(1,1), legend.position=c(1,1),
             axis.text.y = element_text(size=20, face="bold", color="black"), axis.title.y = element_text(size=30, face="bold", hjust=1), 
             axis.text.x = element_text(size=20, face="bold", color="black"), axis.title.x=element_text(size=30, face="bold", margin = margin(t=-5)),
             legend.text=element_text(size=25),
             legend.key.size=unit(0.8,"cm"),
             legend.key = element_rect(colour = 'white', fill = 'white'))

ggsave(title, plot=man, device="jpeg", path=output_filepath, unit="cm", dpi=1000, width = 84.5, height = 20, limitsize=FALSE) # 26.67cm=converted from 2000 pixels in Tuan's code

}