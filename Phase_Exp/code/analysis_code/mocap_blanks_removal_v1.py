'''
@name: mocap_blanks_removal_v1.py
@description: Remove special blanks in all csv files from the motion capture
@author: Vu Phan
@date: 10 Feb 2020
'''

import csv
import glob
import time

# Define macros
ROW_LENGTH = 5 # a row contains id - time - x - y - z
PATH = '*.csv'
FILE_NAMES = glob.glob(PATH) # get all csv files in the folder

print('Processing ...')
start_time = time.time()

for file_name in FILE_NAMES:
    temp_file = [] # temporary list to store file's content

    ''' Open csv file and copy it into a list (i.e., temp_file) '''
    with open(file_name) as in_file:
        for row in csv.reader(in_file):
            temp_file.append(row) # add every row into the list

    ''' List processing '''
    # Remove empty entries
    for i in range(len(temp_file)):
        while len(temp_file[i]) > ROW_LENGTH:
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

print('Processing time = ' + str(time.time() - start_time))
print('Done!')
