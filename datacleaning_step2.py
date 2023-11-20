import numpy as np
from scipy import loadtxt
from os.path import join
import os
import glob
import re
import shutil
import pandas as pd

if __name__ == "__main__":
    cwd = os.getcwd()

    list_file = []

    vil_used = np.loadtxt('./vil_used.txt', dtype='int', delimiter=',')

    file_temp = open('./network_list.txt', 'r')
    network_txt = file_temp.read()
    network_list = network_txt.split(',')
    file_temp.close()

    # os.chdir('./2010-0760_Data/Data/Raw_csv')

    print('Scanning files...')

    dataset = pd.read_csv('./processed/hh_cov_cleaned.csv', delimiter=',')

    for file in glob.glob('*.csv'):
        print('Loading files...' + str(file), end='\r')
        list_file.append(file)

    print('Scan successful. Total: ' + str(len(list_file)))
    print('Analyzer starts')

    os.chdir(cwd)
    if os.path.exists('./output'):
        shutil.rmtree('./output')

    os.mkdir('./output')

    for counter, vil_num in enumerate(vil_used):
        print('Working on village #' + str(vil_num))

        os.mkdir('./output/' + str(vil_num))

        # Start with separating the village covariate from the original file
        village_covariate = dataset.loc[dataset['village'] == vil_num]

        # Save to a csv file
        village_covariate.to_csv('./output/vil_covariate' + str(vil_num) + '.csv', sep=',', index_label=False)

        key = np.loadtxt('./datav4.0/Data/1. Network Data/Adjacency Matrix Keys/key_HH_vilno_' + str(vil_num) + '.csv', dtype='int', delimiter=',')
        population = key[-1]  # an original population of each village can be deduced by the last number in the key

        entry = village_covariate.loc[:, 'HHnum_in_village']
        entry = entry.array
        population_surveyed = entry.shape[0]

        microfinance = np.loadtxt('./datav4.0/Matlab Replication/India Networks/MF' + str(vil_num) + '.csv', dtype='int', delimiter=',')

        for i in range(population):
            if i + 1 != key[i]:
                key = np.insert(key, i, i + 1)
                microfinance = np.insert(microfinance, i, 0)

        if microfinance.shape[0] != population:
            print('The take-up dimension is incorrect')
            break

        mf_status = 0  # A switch to indicate whether the MF take-up data is processed

        # Make adjacency matrices only with the surveyed individuals
        # The matrix key provided by the original authors is not relevant anymore. Use only the household ID.
        for network_ind in range(12):
            temp_matrix = np.zeros([population, population])
            network_name = network_list[network_ind]
            print(str(network_name), end='\r')
            survey_raw = np.loadtxt('./2010-0760_Data/Data/Raw_csv/' + str(network_name) + str(vil_num) + '.csv', dtype='int',
                                    delimiter=',')
            length = survey_raw.shape[0]
            width = survey_raw.shape[1]

            # Extract the household number from each entry
            for i in range(length):
                respondent = survey_raw[i, 0]
                respondent = str(respondent)
                respondent_hh = int(respondent[-5:-2])

                for j in range(1, width):
                    response = survey_raw[i, j]
                    if response != 0 and response != 7777777 and response != 5555555 and response != 9999999: # Ignore invalid responses
                        response = str(response)
                        hh_num_response = int(response[-5:-2])

                        if hh_num_response <= population:
                            temp_matrix[respondent_hh - 1, hh_num_response - 1] = 1
                        else:
                            print("Flag: response not in range (respondent:", respondent_hh, ", response:", hh_num_response, "network:", network_name, ")")


            # Make a new adjacency matrix for only those surveyed

            for i in range(population)[::-1]:
                if i + 1 not in entry:
                    temp_matrix = np.delete(np.delete(temp_matrix, i, 0), i, 1)
                    if mf_status == 0:
                        microfinance = np.delete(microfinance, i, 0)

            mf_status = 1  # indicate once it's done
            if temp_matrix.shape != (population_surveyed, population_surveyed):
                print("The output adjacency matrix dimension is incorrect")
                break

            if microfinance.shape[0] != population_surveyed:
                print("The MF output matrix dimension is incorrect")
                break

            np.savetxt('./output/MF_outcome' + str(vil_num) + '.csv', microfinance, fmt='%1.0f', delimiter=',')

            temp_matrix = temp_matrix - np.diagflat(np.diag(temp_matrix))
            filename = str('./output/' + str(vil_num) + '/' + str(network_name) + str(vil_num) + '.csv')
            np.savetxt(filename, temp_matrix, fmt='%1.0f', delimiter=',')

print('Complete!')
