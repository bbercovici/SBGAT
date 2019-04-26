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
<<<<<<< HEAD
    	print("PID 27590 is still running")
=======
        print("PID 27590 is still running")
>>>>>>> f4f95d47ac214ccae5e0ba8bbc53e8c900d5005a
        time.sleep(5) 