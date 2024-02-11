import numpy as np
import pandas as pd

def generate_staggered_law_ar1_hetero_data(N, T, num_individuals, mean=0, std_dev=1):
    
    white_noise = np.random.normal(mean, std_dev, size=(N, num_individuals, T))  # Generating random white noise for each individual

    data = np.zeros((N, num_individuals, T))  # Initializing the array to store the data

    # state, time effects and rho for each state

    rhos = np.random.uniform(0.2,0.9, size = N)
    alphas = np.random.normal(0,1, size=N)
    betas = np.random.normal(0,1, size = T)

    # Generate the AR(1) process data for each individual
    
    for i in range(N):
        alpha = alphas[i]
        rho = rhos[i]
        for j in range(num_individuals):
            for t in range(T):
                beta = betas[t]
                if t == 0:
                    data[i, j, t] = alpha + beta + white_noise[i, j, t]
                else:
                    data[i, j, t] = alpha + beta + rho * data[i, j, t - 1] + white_noise[i, j, t]

 
    reshaped_data = data.reshape((N * num_individuals, T))     # Reshaping the data array for easier DataFrame creation

    
    df = pd.DataFrame(reshaped_data, columns=[f'{t}' for t in range(T)])  # Create a DataFrame with column names as time periods

    df['state'] = np.repeat(np.arange(1, N + 1), num_individuals)  # Add a new 'state' column with repeated state values

    df['individual'] = np.tile(np.arange(1, num_individuals + 1), N)  # Add a new 'individual' column with repeated individual values

    melted_df = pd.melt(df, id_vars=['state', 'individual'], var_name='time', value_name= 'value')

    melted_df['time'] = melted_df['time'].astype(int)  # Convert the 'time' column to int


    data = melted_df.copy()

    data['time'] = data['time'].astype(int)
    
    state_dummies = pd.get_dummies(data['state'], prefix='state', drop_first = True)  # Create state dummy variables

    state_dummies = state_dummies.astype(int)

    time_dummies = pd.get_dummies(data['time'].astype(int), prefix='time', drop_first = True) # Create time dummy variables

    time_dummies = time_dummies.astype(int)

    data = pd.concat([data, state_dummies, time_dummies], axis=1)
   
    return data