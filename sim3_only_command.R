# Command for running the script that allows us to go offline
system("nohup Rscript sim3_only_offline_script.R > sim3_log.txt 2>&1 & echo $!")

#For printing some of the lines of the log
cat(tail(readLines("sim3_log.txt"), 20), sep = "\n")

#substitute the PID for the number 480954 here - prints the usage of memory, CPU etc.
system("ps -p 557995 -o pid,etime,%cpu,%mem,cmd")
