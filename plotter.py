import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd

par_results = []
CSR2_results = []
CSR3_results = []

with open('to_plot/1138_bus.mtx_results.txt', 'r') as f:
    text = f.read()
    text_list = eval(text)
    par_results.append(text_list[0])
    CSR2_results.append(text_list[1])
    CSR3_results.append(text_list[2])

with open('to_plot/bcspwr08.mtx_results.txt', 'r') as f:
    text = f.read()
    text_list = eval(text)
    par_results.append(text_list[0])
    CSR2_results.append(text_list[1])
    CSR3_results.append(text_list[2])

with open('to_plot/bcspwr10.mtx_results.txt', 'r') as f:
    text = f.read()
    text_list = eval(text)
    par_results.append(text_list[0])
    CSR2_results.append(text_list[1])
    CSR3_results.append(text_list[2])

with open('to_plot/bcsstk25.mtx_results.txt', 'r') as f:
    text = f.read()
    text_list = eval(text)
    par_results.append(text_list[0])
    CSR2_results.append(text_list[1])
    CSR3_results.append(text_list[2])

with open('to_plot/illc1850.mtx_results.txt', 'r') as f:
    text = f.read()
    text_list = eval(text)
    par_results.append(text_list[0])
    CSR2_results.append(text_list[1])
    CSR3_results.append(text_list[2])


barWidth = 0.25
fig = plt.subplots(figsize =(12, 8)) 

br1 = np.arange(len(par_results)) 
br2 = [x + barWidth for x in br1] 
br3 = [x + barWidth for x in br2] 

plt.bar(br1, par_results, color ='b', width = barWidth, 
        edgecolor ='grey', label ='BASE PARALLEL SPEEDUP')
plt.bar(br2, CSR2_results, color ='y', width = barWidth, 
        edgecolor ='grey', label ='CSR-2 PARALLEL SPEEDUP')
plt.bar(br3, CSR3_results, color ='r', width = barWidth, 
        edgecolor ='grey', label ='CSR-3 PARALLEL SPEEDUP')

plt.xlabel('Matrixes', fontweight ='bold', fontsize = 15) 
plt.ylabel('Speedup %', fontweight ='bold', fontsize = 15) 
plt.xticks([r + barWidth for r in range(len(par_results))], 
        ['1138_bus', 'bcspwr08', 'bcspwr10', 'bcsstk25', 'illc1850'])

plt.legend()
plt.savefig('speedup_plot.png', bbox_inches='tight')