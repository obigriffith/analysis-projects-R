#!/usr/bin/env /csb/home/obig/bin/R64/R-2.11.1/bin/Rscript

args=(commandArgs(TRUE))
args.all=commandArgs(trailingOnly = FALSE)
CELfile=args[1]

file.arg.name="--file="
R_script=sub(file.arg.name, "", args.all[grep(file.arg.name, args.all)])
R_script_name=strsplit(R_script,"/")[[1]][length(strsplit(R_script,"/")[[1]])]
script_dir=getwd()
script_path=paste(script_dir,R_script_name,sep="/")

sink("test.tex")

cat ("\\documentclass{article}\n")
cat ("\\usepackage{amsmath}\n")
cat ("\\usepackage{amscd}\n")
cat ("\\usepackage[tableposition=top]{caption}\n")
cat ("\\usepackage{ifthen}\n")
cat ("\\usepackage[utf8]{inputenc}\n")
cat ("\\usepackage[hmargin=3cm,vmargin=3cm]{geometry}\n")
cat ("\\usepackage{hyperref}\n")
cat ("\\hypersetup{\n")
cat ("   colorlinks = true,\n")
cat ("   linkcolor = blue}\n")

cat ("\\begin{document}\n")

cat ("\\title{Treatment Response Prediction Report}\n")
cat ("\\author{Anneleen Daemon and Obi Griffith (PI: Gray/Spellman)}\n")
cat ("\\maketitle\n")

cat ("\\subsubsection*{Sample Information:}\n")
cat ("CEL File: ",CELfile,"\n")
cat ("R script: ",script_path,"\n")


cat ("\\subsubsection*{Treatment Response Predictions:}\n")

cat ("\\begin{tabular}{| l || l | l | c | c | c |}\n")
cat ("\\hline\n")
cat ("1 & 2 & 3 & 1 & 2 & 3 \\\\\n")
cat ("4 & 5 & 6 & 1 & 2 & 3 \\\\\n")
cat ("7 & 8 & 9 & 1 & 2 & 3 \\\\\n")
cat ("\\hline\n")
cat ("\\end{tabular}\n")

cat ("\\subsubsection*{Recommended Treatment:}\n")

cat ("\\end{document}\n")

sink()

system(paste('pdflatex',"test.tex")) 


