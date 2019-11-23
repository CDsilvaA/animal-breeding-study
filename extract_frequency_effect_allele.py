'''
Developers: [
Thiago Roberto do Prado <trprado@outlook.com https://github.com/trprado>,
 Cherlynn Silva <cdnsprado@gmail.com https://github.com/CDsilvaA>
 ]

Example of data
1:-0.263323 2: 0.112048    
1: 0.064858 2:-0.086547    
1: 0.596754 2:-0.658579    
1: 0.251089 2:-0.063362    
1: 0.201582 2:-0.025194 3:-0.478914   
1:-0.304997 2: 0.114748    
1:-0.789275 2: 0.542554    
1: 0.003843 2:-0.003385    
1: 0.163325 2:-0.585013    
1:-0.015535 2: 0.005458    
1:-0.050356 2: 0.051015    
1: 0.101408 2:-0.059016    
1:-0.010441 2: 0.008478    
1: 0.178348 2:-0.053310    
1: 0.124572 2:-0.073397    
1:-0.113552 2:-0.781778 3: 0.894587 4: 0.306456 

Objective: Extract information from the frequencies or effects of the alleles in the QMSim output. 
The extracted information will be saved in a matrix where each column corresponds to the allele, missing data will be completed with zero.

The user must change the name of the input and output files manually.
'''


text = []
with open('ef30', 'r') as file:
    for line in file.readlines():
        text.append(line.replace(': ', ':').strip())

data = []
for line in text:
    aux1, aux2, aux3, aux4 = 0, 0, 0, 0
    
    for val in line.split():
        if val.startswith('1:'):
            aux1 = val[2:]
        elif val.startswith('2:'):
            aux2 = val[2:]
        elif val.startswith('3:'):
            aux3 = val[2:]
        elif val.startswith('4:'):
            aux4 = val[2:]
    
    data.append({
        '1': aux1,
        '2': aux2,
        '3': aux3,
        '4': aux4
    })

import pandas as pd

df = pd.DataFrame(data)
df.to_csv('efeitos_rep30.csv', index=False)
