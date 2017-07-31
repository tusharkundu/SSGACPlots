#' Create A QQ Plot
#' 
#' This function intakes GWAS summary statsitcs and produces a QQ Plot.
#' @param sumstat_data_path Summary statistics data path
#' @param output_filepath Where you will output the QQ Plot
#' @param P.value The title of your p-value data column (include quotation marks)
#' @param delimiter The way your data is separated (often with "/t" or ",")
#' @param title What you would like the title to be
#' @export
#' @examples 
#' ssgac_qq(sumstat_data_path = "path/sumstata.txt" , 
#' output_filepath = "outpath", 
#' P.value = "pval", delimiter = "/t", 
#' title="ssgac_qqplot"))
#' 
###############################
##QQ-plot function definition##
###############################

ssgac_qq <- function(sumstat_data_path, output_filepath, P.value, delimiter, title){
  
  ## Install ggplot2
  require(ggplot2)
  library(ggplot2)
  
  ##Read in argument to delimiter variable
  delimiter_variable <- delimiter_var
  ##Read in argument to p-value column header variable
  P.value_variable <-  P.value_var
  
  ##Input validation
  if (missing(sumstat_data_path))
    stop("Need to specify sumstats filepath and filename to be plotted, e.g. \"../INPUT/sumstats1.txt\" .")
  
  if (missing(P.value_variable))
    stop("Need to specify header of p-value column within double quotes, e.g. \"P\"")
  if (is.numeric(P.value_variable))
    stop("Need to specify the header of the p-value column as a character string.")
  
  if (missing(delimiter))
    stop("Need to specify the field delimiter within double quotes, e.g. \",\" or \" \".")
  
  if (missing(output_filepath))
    stop("Need to specify the output filepath and prefix withing double quotes, e.g. \"../OUTPUT/qq_plot_GWAS.jpeg\".")
  
  ##Print executing...
  print("Executing QQ-plot function...")
  
  ##Read in first row of sumstats and then check which field (column) corresponds to the p-value header
  header_row_temp_var = read.table(file = sumstat_data_path, 
                                   sep="\t",
                                   header = FALSE,
                                   nrow=1)
  pval_column_number <- which(header_row_temp_var == P.value_variable)
  
  colclass_variable <- rep(0, length(header_row_temp_var))
  for(i in 1:length(header_row_temp_var)){
    if(i == pval_column_number) {
      colclass_variable[i] <- c("numeric")
    } else {colclass_variable[i] <- c("NULL")}
  }
  
  ## Create results for axis maximum values
  pval_data_temp <- unlist(read.table(file = sumstat_data_path, skip = 1, header = FALSE, colClasses = colclass_variable, sep = delimiter_variable, na.strings = "NA", dec = ".", stringsAsFactors = FALSE, blank.lines.skip = TRUE, check.names = FALSE))
  
  # obs <- readfile; p-values only
  obs <- -log10(pval_data_temp) ## read in your p-values,
  N <- length(pval_data_temp) ## count p-values
  
  # Create generic rounding function
  round.choose <- function(x, round.val, dir = 1) {
    if(dir == 1) {  ##ROUND UP
      x + (round.val - x %% round.val)
    } else {
      if(dir == 0) {  ##ROUND DOWN
        x - (x %% round.val)
      }
    }
  }
  
  # Find max -log10(p-value)
  MAX <- max(obs, na.rm = TRUE)
  # Make sure to round to the same incremental steps as the axis labels further down.
  if(MAX <= 20) {MAX_round <- (round.choose(MAX, 2, 0) + 2)} else if (MAX > 20 && MAX <= 50)
  {MAX_round <- (round.choose(MAX, 4, 0) + 4)} else if (MAX > 50 && MAX <= 100) 
  {MAX_round <- (round.choose(MAX, 8, 0) + 8)} else if(MAX > 100 && MAX <= 200) 
  {MAX_round <- (round.choose(MAX, 10, 0) + 10)} else if(MAX > 200) 
  {MAX_round <- (round.choose(MAX, 25, 0) + 25)}
  
  # Specify the max value for annotation and axis
  ylim <- MAX_round + 3
  ytext <- MAX_round + 1
  
  
  # create empty vectors for the confidence intervals
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  ## Reference from: http://www.gettinggeneticsdone.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
  
  for(i in 1:N){
    c95[i] <- qbeta(0.975,i,N-i+1) 
    c05[i] <- qbeta(0.025,i,N-i+1)
  }
  c95 <- rev(-log(c95,10))
  c05 <- rev(-log(c05,10))
  
  print(paste("Starting plot of ", N, " variants..."), sep = "")
  
  ## Read data for plotting
  dataname <- read.table(file = sumstat_data_path, 
                         sep=delimiter,  
                         header = T, colClasses = colclass_variable, na.strings = "NA", 
                         dec = ".", stringsAsFactors = FALSE, blank.lines.skip = TRUE, 
                         check.names = FALSE)
  
  
  # Create variable on x axis
  N_data <- length(dataname[,1])
  null <- -log(1:N_data/N_data,10)
  null <- sort(null)
  # Create variable on y axis
  column <- which(colnames(dataname) %in% P.value)
  pvalname <- -log10(dataname[,column])
  pvalname <- sort(pvalname)
  # Create data to be used in ggplot2, which works best with direct variables in dataset
  newdata <- cbind(dataname,null,pvalname,c95,c05)
  
  # Add annotate text for lambda
  newdata$z_temp <- qnorm(newdata[[P.value]]/2)  # calculate mean chisq
  lamx <- max(newdata$null)
  lamy <- ytext
  result <- median(newdata$z_temp^2)/0.4549364
  legend <- sprintf("\"GWAS\"~lambda == %0.2f", result)
  
  # Specify filename and titles
  filename <- paste("QQ_", current_file_var, "_lambda.jpeg", sep="")
  ytitle <- paste("Observed ","-log10","(P-value)")
  xtitle <- paste("Expected ","-log10","(P-value)")
  
  # Good to go
  qq <- ggplot(newdata, aes(null, pvalname), color="royalblue") 
  qq <- qq + geom_point(shape=19, color="royalblue") 
  qq <- qq + geom_ribbon(aes(ymin=c05,ymax=c95), fill="gray")  
  qq <- qq + geom_abline(intercept=0, linetype="dashed") 
  qq <- qq + labs(x = xtitle, y = ytitle)
  qq <- qq + ggtitle(title) 
  qq <- qq + annotate("text", x = lamx, y = lamy, label = legend, size=9, parse=T)
  qq <- qq + scale_y_continuous(limits=c(0,ylim)) 
  qq <- qq + scale_x_continuous(limits=c(0,ylim)) 
  qq <- qq + theme(plot.title = element_text(hjust=0.5, vjust=-5, lineheight=3, face="bold", color="black", size=30),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text.y = element_text(size=20, face="bold", color="black"), axis.title.y = element_text(size=25, face="bold"), 
                   axis.text.x = element_text(size=20, face="bold", color="black"), axis.title.x=element_text(size=25, face="bold", margin = margin(t=5)))
  ggsave(filename, plot=qq, device="jpeg", path=output_filepath, units = "cm", width = 26.67, height = 26.67) # converted from 2000 pixels in Tuan's code
}
