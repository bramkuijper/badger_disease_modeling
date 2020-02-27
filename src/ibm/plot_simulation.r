#!/usr/bin/env Rscript

library("ggplot2")
library("gridExtra")

args <- commandArgs(trailingOnly=T)


if (length(args) < 1)
{
    print("provide a simulation file name")
    stop()
}


# find out where the parameter listing starts
# so that we can read in the data part of the file 
# without having it messed up by the subsequent parameter listing
find_out_param_line <- function(filename) {

    f <- readLines(filename)

    # make a reverse sequence
    seqq <- seq(length(f),1,-1)

    # go through each line in the data file and find first line
    # where data is printed (i.e., a line which starts with a digit)
    for (line_i in seqq)
    {
        print(f[[line_i]])
        print(line_i)
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }

    return(NA)
}

parameter_row <- find_out_param_line(args[1])

if (is.na(parameter_row))
{
    print("cannot find data...")
    stop()
}

# read in data frame of corresponding simulation
the.data <- read.table(args[1], header=T, nrow=parameter_row - 1, sep=";")

str(the.data)

p1.a <- ggplot(data=the.data
        ,aes(x=time)) +
            geom_line(aes(y = N_I, colour="Infected")) +
            geom_line(aes(y = N_S, colour="Susceptible")) +
            geom_line(aes(y = N_L, colour="Latent")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("N")

p1.b <- ggplot(data=the.data
        ,aes(x=time)) +
            geom_line(aes(y = meanv, colour="virulence")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("v")

big_plot <- arrangeGrob(p1.a, p1.b,nrow=2,ncol=1)


the.base.name <- basename(args[1])

output_file_name <- paste(
        "graph_"
        ,the.base.name
        ,".pdf"
        ,sep="")

ggsave(output_file_name, big_plot, height = 25)

