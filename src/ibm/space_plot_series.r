#!/usr/bin/env Rscript

library("lattice")
library("stringr")
library("colorRamps")

args <- commandArgs(trailingOnly=T)

time.col <- "time"

ffmpeg.exe <- "ffmpeg"

if (length(args) < 1)
{
    print("no files provided, meh...")
}


for (file in args) 
{
    # check whether this is indeed the desired file name
    # sim_023142334_space
    if (length(grep("sim_.*_space",file)) < 1)
    {
        next
    }

    dat <- read.table(file,sep=";",header=T)

    scale <- seq(0, max(dat$N), max(dat$N) / 50)

    gen.u <- sort(unique(dat[,time.col]))

    pad.width <- str_length(as.character(max(gen.u)))

    file.list <- c()

    for (gen.i in gen.u)
    {
        filename <- paste(basename(file)
                ,"_"
                ,str_pad(as.character(gen.i)
                        ,width=pad.width
                        ,side="left"
                        ,pad="0")
                ,".png"
                ,sep="")

        # append file to list of all files
        file.list <- c(file.list,filename)

        sub.dat <- dat[
                        dat[,time.col] == gen.i
                        & dat[,"sex"] == "F"
                        & dat[,"inf_state"] == 2,]

        png(file = filename)
        print(
                levelplot(
                        N ~ row * column
                        ,data=sub.dat
                        ,col.regions=matlab.like
                        ,xlab="X coordinate"
                        ,ylab="Y coordinate"
                        ,at=scale
                        )
                )
        dev.off()
    }

    video.file <- paste(file,"_video.avi",sep="") 
    image.file.pattern <- paste(file,"_*.png",sep="")

    ffmpeg.command <- paste( 
            ffmpeg.exe # executable
            ,"-f image2" # image2 format, which means that you have specify images using globs
            ,"-pattern_type glob"
            ,"-i",paste("\"",image.file.pattern,"\"",sep="")
            ,video.file)

    system(ffmpeg.command)

    print(image.file.pattern)
    unlink(image.file.pattern)
}
