import pandas as pd


y = [2, 4, 1]
df = pd.DataFrame()
print(df)
df[1] = [3, 2, 1]
df[2] = [7, 3, 9]
print(df)
for i, row in df.iterrows():
    print(tuple(row))
