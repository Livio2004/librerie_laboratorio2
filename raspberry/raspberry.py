import pandas as pd 
import matplotlib.pyplot as plt

data = pd.read_csv('data.txt')
y_data = data.iloc[:]
x_data = list(range(1, len(data) + 1))
plt.figure(figsize=(10, 6))
plt.plot(x_data, y_data, label='LDR Data', color='blue')
plt.xlabel('time')
plt.ylabel('LDR Value')
plt.title('LDR Data over Time')
plt.legend()
plt.grid()
plt.show()