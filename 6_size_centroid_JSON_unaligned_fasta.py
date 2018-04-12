#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 09:49:52 2018

@author: jordanzehr
"""

import json
import itertools as it
#import random

## change the number (3 spots) at the front of the 3 commands based on the data input

with open('6_clone.json','r') as input_file:
    data = json.load(input_file)

with open('6_size_20_centroid_unaligned_fasta.txt', 'w') as output_file:
    for i, item in enumerate(it.chain.from_iterable(data)):
#        value1 = random.randint(1, 1000)
#        value2 = random.randint(1001, 2000)
#        value3 = random.randint(2001, 3000)
        if int(item["size"]) > 20:
            #print(j)
            output_file.write(''.join(
                    ['>seq' + str(i) + str('_6_') + str(item["size"]) + str(item["centroid"]).replace(':','.')]
            ))
#str(value1) + str(value2) + str(value3)
# str('11')
#str(item["size"])            
#str('.size:') +
output_file.close()