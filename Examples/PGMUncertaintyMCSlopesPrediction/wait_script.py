import psutil 
import time 
import subprocess
import sys
sys.path.append("..")
import master_script
pid = 27590 
while(1): 
    if not psutil.pid_exists(pid): 
        master_script.start_sim()
        break
    else: 
        print("PID 27590 is still running")
        time.sleep(5) 