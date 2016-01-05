import csv
import numpy as np

def write_data(file_name, file_ext, data):
    file_full_name = file_name + '.' + file_ext
    ## if the data of data is 3D, otherwise m, n without p
    m, n, p = data.shape
    data = data.reshape(m, n*p)
    print data
    with open(file_full_name, 'a+') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        for z_index in range(m):
            writer.writerow(data[z_index])

def read_data(file_name, file_ext):
    file_full_name = file_name + '.' + file_ext
    with open(file_full_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', lineterminator='\n')
        ## get the dimension of the data of data
        z_index = 0
        for row in reader:
            z_index += 1
    data = np.zeros((z_index, len(row)))
    
    with open(file_full_name, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', lineterminator='\n')
        z_index = 0
        for row in reader:
            xy_plane = np.zeros(len(row))
            for i in range(len(row)):                
                xy_plane[i] = eval(row[i])
            data[z_index] = xy_plane
            z_index += 1
    return data

## test codes here
temperature = np.arange(4*5*6).reshape(6, 5, 4)
print temperature
temperature = temperature.transpose((2, 0, 1))
print temperature
write_data('temperature', 'csv', temperature)
temperature = read_data('temperature' , 'csv')
print temperature
