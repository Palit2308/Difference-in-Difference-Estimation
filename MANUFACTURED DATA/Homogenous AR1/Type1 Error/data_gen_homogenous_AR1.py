import numpy as np
import pandas as pd

def generate_staggered_law_ar1_data(N, T, rho, num_individuals, mean=0, std_dev=1):
    
    white_noise = np.random.normal(mean, std_dev, size=(N, num_individuals, T)) # generatingk random white noise for each individual

    data = np.zeros((N, num_individuals, T))  # Initializing the array to store the data

    # Including state and time effects 

    alphas = np.random.normal(0,1, size=N)
    betas = np.random.normal(0,1, size = T)

    # Generating the AR(1) process data for each individual

    for i in range(N):
        alpha = alphas[i]
        for j in range(num_individuals):
            for t in range(T):
                beta = betas[t]
                if t == 0:
                    data[i, j, t] = alpha + beta + white_noise[i, j, t]
                else:
                    data[i, j, t] = alpha + beta + rho * data[i, j, t - 1] + white_noise[i, j, t]

    
    reshaped_data = data.reshape((N * num_individuals, T))  # Reshaping the data array for easier DataFrame creation

 
    df = pd.DataFrame(reshaped_data, columns=[f'{t}' for t in range(T)])    # Creating a DataFrame with column names as time periods

    
    df['state'] = np.repeat(np.arange(1, N + 1), num_individuals)  # Add a new 'state' column with repeated state values

    
    df['individual'] = np.tile(np.arange(1, num_individuals + 1), N)  # Add a new 'individual' column with repeated individual values


    melted_df = pd.melt(df, id_vars=['state', 'individual'], var_name='time', value_name= 'value') # to take the data in the long format

   
    melted_df['time'] = melted_df['time'].astype(int)   # Converting the 'time' column to int


    data = melted_df.copy()

   
    return data