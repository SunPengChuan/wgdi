import numpy as np
import pandas as pd



class block_matirx():
    def dist_matirx(self,data)
        dm = np.zeros(shape=(len(data),len(data)))
        dm.fill(np.nan)
        for index1, row1 in data.iterrows():
            for index2, row2 in data.iterrows():
                if index2>=index1:
                        continue
                    a, b= row1[['y1','x1',]].values,row2[['y1','x1']].values
                    dist[index1][index2]=min(abs(a[1]-b[0]),abs(a[0]-b[1]))
        for chr2, group in data.groupby('chr2'):
            for index1, row1 in group.iterrows():
                for index2, row2 in group.iterrows():
                    if index2<=index1:
                        continue
                    a, b = row1[['y2', 'x2', ]].values, row2[['y2', 'x2']].values
                    dist[index1][index2] = min(abs(a[1] - b[0]), abs(a[0] - b[1]))
        df2=pd.DataFrame(dist)
        df2=(df2.max().max()-df2)/df2.max().max()