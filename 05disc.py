import pandas as pd

df = pd.read_csv('03data_bbknn_b/all.txt', sep='\t')
# Discretize the values in df
df_discretized = df.copy()
df_discretized = df_discretized.drop(columns=["tissue"])

for column in df_discretized.columns[2:]: # Exclude the '@name' column
    quantiles = df_discretized[column].quantile([0.001, 0.75])
    q0 = quantiles[0.001]+1.0e-10
    q75 = quantiles[0.75]+1.0e-10
    # Handle cases where q0 and q75 are the same
    if q0 == q75:
        df_discretized[column] = pd.cut(df_discretized[column],
                                        bins=[-float('inf'), q0, float('inf')],
                                        labels=[0, 1],  # Use two labels if only two bins are possible
                                        right=False,
                                        include_lowest=True)
    else:
        df_discretized[column] = pd.cut(df_discretized[column],
                                        bins=[-float('inf'), q0, q75, float('inf')],
                                        labels=[0, 1, 2],
                                        right=False,
                                        include_lowest=True)


# Identify columns with no variation (all values are the same)
cols_to_drop = [col for col in df_discretized.columns if df_discretized[col].nunique() == 1]

# Print the list of columns to be dropped
print("Columns with no variation:")
print(cols_to_drop)

# Drop the identified columns from the DataFrame
df_discretized_filtered = df_discretized.drop(columns=cols_to_drop)


# Create a copy to avoid modifying the original DataFrame
df_transformed = df_discretized_filtered.copy()

# Create lists to store new columns
new_cols_0_or_greater = {}
new_cols_is_2 = {}

for col in df_transformed.columns[1:]: # Exclude the '@name' column
    # Column for 0 or >0
    new_cols_0_or_greater[f'{col}'] = (df_transformed[col] > 0).astype(int)

    # Column for value 2 (only if value 2 exists in the column)
    #if 2 in df_transformed[col].unique():
    #    new_cols_is_2[f'{col}_high'] = (df_transformed[col] == 2).astype(int)

# Create DataFrames from the new columns
df_0_or_greater = pd.DataFrame(new_cols_0_or_greater, index=df_transformed.index)
df_is_2 = pd.DataFrame(new_cols_is_2, index=df_transformed.index)


# Concatenate the original DataFrame (excluding the original numeric columns) with the new columns
df_transformed = pd.concat([df_transformed[['@name']], df_0_or_greater, df_is_2], axis=1)


# Identify columns with no variation in the transformed DataFrame
cols_to_drop_transformed = [col for col in df_transformed.columns if df_transformed[col].nunique() == 1]

# Print the list of columns to be dropped
print("Columns with no variation after transformation:")
print(cols_to_drop_transformed)

# Drop the identified columns
df_transformed_filtered = df_transformed.drop(columns=cols_to_drop_transformed)
df_transformed_filtered["tissue"]=df["tissue"]

df_transformed_filtered['@name'] = df_transformed_filtered['@name'].apply(lambda x: '|'.join([x.split('|')[0], f'{int(x.split("|")[1][:-1]):02d}{x.split("|")[1][-1]}', x.split('|')[2]]))

# sort
df_transformed_filtered = df_transformed_filtered.sort_values(by="@name")


df_transformed_filtered.to_csv('03data_bbknn_b/all_disc.txt',sep="\t",index = False)




