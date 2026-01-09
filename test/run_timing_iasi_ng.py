import subprocess
import os
import re
import sys

def run_tests():
    tests = [
        "forward",
        "k_matrix",
        "adjoint",
        "tangent_linear"
    ]
    
    thread_counts = [1, 2, 4, 8, 16]
    sensor_id = "iasi-ng_metop-sg-a1"
    
    print(f"{'Test':<20} {'Threads':<10} {'Time (s)':<15}")
    print("-" * 45)
    
    for test in tests:
        # Assuming running from build/test
        executable = f"../bin/test_{test}_test_OMPoverChannels"
        
        for threads in thread_counts:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(threads)
            
            cmd = [executable, sensor_id]
            
            try:
                result = subprocess.run(
                    cmd,
                    env=env,
                    capture_output=True,
                    text=True,
                    cwd="." 
                )
                
                output = result.stdout + result.stderr
                
                if result.returncode != 0:
                    print(f"{test:<20} {threads:<10} FAILED (RC={result.returncode})")
                    continue
                
                match = re.search(r"CRTM_\w+ wall time \(s\):\s+([\d\.]+)", output)
                if match:
                    time_s = match.group(1)
                    print(f"{test:<20} {threads:<10} {time_s:<15}")
                else:
                    print(f"{test:<20} {threads:<10} Time not found")
                    
            except Exception as e:
                print(f"{test:<20} {threads:<10} ERROR: {e}")

if __name__ == "__main__":
    run_tests()
