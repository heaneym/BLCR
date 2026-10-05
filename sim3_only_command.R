# Command for running the script that allows us to go offline
system("nohup Rscript sim3_only_offline_script.R > sim3_log.txt 2>&1 & echo $!")


# how many replicates done for each sample size
sapply(c("N500", "N300", "N150"), function(nm) {
  f <- paste0("sim3_progress_", nm, ".log")
  if (file.exists(f)) length(readLines(f)) else 0L
})

# last few lines of a progress file
tail(readLines("./sim3_progress_N500.log"), 5)
tail(readLines("./sim3_progress_N300.log"), 5)
tail(readLines("./sim3_progress_N150.log"), 5)

#substitute the PID for the number 480954 here - prints the usage of memory, CPU etc.
system("ps -p 22744 -o pid,etime,%cpu,%mem,cmd")
