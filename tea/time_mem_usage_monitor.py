import sys
import os
import tracemalloc
import time

rounds = 25

def execute_grass(json_path):
    time_list = []
    mem_list = []
    for i in range(rounds):
        tracemalloc.start()
        start_time = time.time()

        os.system(f'python ../regridder/regridder.py {json_path}')

        exec_time = (time.time() - start_time)
        exec_memory = tracemalloc.get_traced_memory()[1]
        time_list.append(exec_time)
        mem_list.append(exec_memory)
        tracemalloc.stop()
    #get avg time and memory spent
    avg_time = sum(time_list) / len(time_list)
    avg_mem = sum(mem_list) / len(mem_list)

    print('Avg time taken for regridder execution', round( avg_time, 3 ), ' seconds')
    print('Avg memory taken for regridder execution: ', round( avg_mem, 3 ), ' KiB')

def execute_lears(json_path):
    
    time_list = []
    mem_list = []
    for i in range(rounds):
        tracemalloc.start()
        start_time = time.time()

        os.system(f'python ../lears-mini/lears-mini.py {json_path}')

        exec_time = (time.time() - start_time)
        exec_memory = tracemalloc.get_traced_memory()[1]
        time_list.append(exec_time)
        mem_list.append(exec_memory)
        tracemalloc.stop()
    #get avg time and memory spent
    avg_time = sum(time_list) / len(time_list)
    avg_mem = sum(mem_list) / len(mem_list)

    print('Avg time taken for Lears execution', round( avg_time, 3 ), ' seconds')
    print('Avg memory taken for Lears execution: ', round( avg_mem, 3 ), ' KiB')

if __name__ == "__main__":
    execute_lears(sys.argv[1])
    execute_grass(sys.argv[1])
    