'''
@name: mocap_blanks_removal_v1.py
@description: Remove special blanks in all csv files from the motion capture
@author: Vu Phan
@date: 10 Feb 2020
'''

import csv
import glob
import time
import os
from greeting import *

# Define macros
MIN_row_length = 5 # a row contains id - time - x - y - z (older version)
PATH = '*.csv'
FILE_NAMES = glob.glob(PATH) # get all csv files in the folder
row_length = 0; # initial row length

# Greetings
greetings()

start_time = time.time()

file_count = 0;
for file_name in FILE_NAMES:
    try:
        temp_file = [] # temporary list to store file's content

        ''' Open csv file and copy it into a list (i.e., temp_file) '''
        with open(file_name) as in_file:
            for row in csv.reader(in_file):
                temp_file.append(row) # add every row into the list
        
        ''' List processing '''
        # Remove empty entries
        for i in range(len(temp_file)):
            # print(i) # DEBUG
            ''' Get row_length '''
            row_length = 0;
            for k in range(len(temp_file[i])):
                if temp_file[i][k] != '':
                    row_length = row_length + 1
                else:
                    pass # do nothing
            # Adjust row_length to pad values of 1000
            if row_length < MIN_row_length:
                row_length = MIN_row_length
            else:
                pass # do nothing
            ''' Remove empty space '''
            while len(temp_file[i]) > row_length:
                temp_file[i].remove('')

        # Replace remaining empty with 1000
        for i in range(len(temp_file)):
            for j in range(len(temp_file[i])):
                if temp_file[i][j] == '':
                    temp_file[i][j] = '1000'
                else:
                    pass # do nothing

        ''' Overwrite the csv file '''
        with open(file_name, 'w') as out_file:
            out_file.truncate() # delete all the content
            for row in temp_file:
                for i in range(len(row)):
                    if i < (len(row) - 1):
                        out_file.write(row[i] + ',')
                    else:
                        out_file.write(row[i])

                out_file.write('\n') # EOL

        # Print file count
        file_count = file_count + 1
        print('File ' + str(file_count) + ' done!')
    except Exception as inst:
        print("Cannot process file " + str(file_count) + " because " + str(inst))
        print("Wrong data in position " + str(i + 1)  + "\n")

print('Processing time = ' + str(time.time() - start_time))
print('Done!')

os.system("pause")
