# Command for running the script that allows us to go offline
system("nohup Rscript sim3_only_offline_script.R > sim3_log.txt 2>&1 & echo $!")
