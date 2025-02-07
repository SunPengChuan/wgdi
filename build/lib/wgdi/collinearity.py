import numpy as np
import pandas as pd


class collinearity:
    def __init__(self, options, points):
        # Default values
        self.gap_penalty = -1
        self.over_length = 0
        self.mg1 = 40
        self.mg2 = 40
        self.pvalue = 1
        self.over_gap = 3
        self.points = points
        self.p_value = 0
        self.coverage_ratio = 0.8
        
        # Set user-defined options
        for k, v in options:
            setattr(self, str(k), v)

        # Initialize grading and mg values
        self.grading = [50, 40, 25] if not hasattr(self, 'grading') else [int(k) for k in self.grading.split(',')]
        self.mg1, self.mg2 = [40, 40] if not hasattr(self, 'mg') else [int(k) for k in self.mg.split(',')]

        # Convert string values to floats
        self.pvalue = float(self.pvalue)
        self.coverage_ratio = float(self.coverage_ratio)

    def get_matrix(self):
        """Initialize the matrix for the collinearity points."""
        self.points['usedtimes1'] = 0
        self.points['usedtimes2'] = 0
        self.points['times'] = 1
        self.points['score1'] = self.points['grading']
        self.points['score2'] = self.points['grading']
        self.points['path1'] = self.points.index.to_numpy().reshape(len(self.points), 1).tolist()
        self.points['path2'] = self.points['path1']
        self.points_init = self.points.copy()
        self.mat_points = self.points

    def run(self):
        """Run the main collinearity processing."""
        self.get_matrix()
        self.score_matrix()
        data = []

        # Process points for maxPath in the positive direction
        points1 = self.points[['loc1', 'loc2', 'score1', 'path1', 'usedtimes1']].sort_values(by=['score1'], ascending=False)
        points1.drop(index=points1[points1['usedtimes1'] < 1].index, inplace=True)
        points1.columns = ['loc1', 'loc2', 'score', 'path', 'usedtimes']
        
        while (self.over_length >= self.over_gap or len(points1) >= self.over_gap):
            if self.max_path(points1):
                if self.p_value > self.pvalue:
                    continue
                data.append([self.path, self.p_value, self.score])

        # Process points for maxPath in the negative direction
        points2 = self.points[['loc1', 'loc2', 'score2', 'path2', 'usedtimes2']].sort_values(by=['score2'], ascending=False)
        points2.drop(index=points2[points2['usedtimes2'] < 1].index, inplace=True)
        points2.columns = ['loc1', 'loc2', 'score', 'path', 'usedtimes']

        while (self.over_length >= self.over_gap) or (len(points2) >= self.over_gap):
            if self.max_path(points2):
                if self.p_value > self.pvalue:
                    continue
                data.append([self.path, self.p_value, self.score])

        return data

    def score_matrix(self):
        """Calculate the scoring matrix for the points."""
        for index, row, col in self.points[['loc1', 'loc2']].itertuples():
            # Get points within a certain range
            points = self.points[(self.points['loc1'] > row) & 
                                 (self.points['loc2'] > col) & 
                                 (self.points['loc1'] < row + self.mg1) & 
                                 (self.points['loc2'] < col + self.mg2)]
            
            row_i_old, gap = row, self.mg2
            for index_ij, row_i, col_j, grading in points[['loc1', 'loc2', 'grading']].itertuples():
                if col_j - col > gap and row_i > row_i_old:
                    break
                score = grading + (row_i - row + col_j - col) * self.gap_penalty
                score1 = score + self.points.at[index, 'score1']
                if score > 0 and self.points.at[index_ij, 'score1'] < score1:
                    self.points.at[index_ij, 'score1'] = score1
                    self.points.at[index, 'usedtimes1'] += 1
                    self.points.at[index_ij, 'usedtimes1'] += 1
                    self.points.at[index_ij, 'path1'] = self.points.at[index, 'path1'] + [index_ij]
                    gap = min(col_j - col, gap)
                    row_i_old = row_i

        # Reverse processing to handle negative direction
        points_reverse = self.points.sort_values(by=['loc1', 'loc2'], ascending=[False, True])
        for index, row, col in points_reverse[['loc1', 'loc2']].itertuples():
            points = points_reverse[(points_reverse['loc1'] < row) & 
                                    (points_reverse['loc2'] > col) & 
                                    (points_reverse['loc1'] > row - self.mg1) & 
                                    (points_reverse['loc2'] < col + self.mg2)]
            
            row_i_old, gap = row, self.mg2
            for index_ij, row_i, col_j, grading in points[['loc1', 'loc2', 'grading']].itertuples():
                if col_j - col > gap and row_i < row_i_old:
                    break
                score = grading + (row - row_i + col_j - col) * self.gap_penalty
                score2 = score + self.points.at[index, 'score2']
                if score > 0 and self.points.at[index_ij, 'score2'] < score2:
                    self.points.at[index_ij, 'score2'] = score2
                    self.points.at[index, 'usedtimes2'] += 1
                    self.points.at[index_ij, 'usedtimes2'] += 1
                    self.points.at[index_ij, 'path2'] = self.points.at[index, 'path2'] + [index_ij]
                    gap = min(col_j - col, gap)
                    row_i_old = row_i

    def max_path(self, points):
        """Find the maximum path for the given points."""
        if len(points) == 0:
            self.over_length = 0
            return False
        
        # Initialize path score and index
        self.score, self.path_index = points.loc[points.index[0], ['score', 'path']]
        self.path = points[points.index.isin(self.path_index)]
        self.over_length = len(self.path_index)
        
        # Check if the block overlaps with other blocks
        if self.over_length >= self.over_gap and len(self.path) / self.over_length > self.coverage_ratio:
            points.drop(index=self.path.index, inplace=True)
            [loc1_min, loc2_min], [loc1_max, loc2_max] = self.path[['loc1', 'loc2']].agg(['min', 'max']).to_numpy()

            # Calculate p-value
            gap_init = self.points_init[(loc1_min <= self.points_init['loc1']) & 
                                        (self.points_init['loc1'] <= loc1_max) & 
                                        (loc2_min <= self.points_init['loc2']) & 
                                        (self.points_init['loc2'] <= loc2_max)].copy()
            
            self.p_value = self.p_value_estimated(gap_init, loc1_max - loc1_min + 1, loc2_max - loc2_min + 1)
            self.path = self.path.sort_values(by=['loc1'], ascending=[True])[['loc1', 'loc2']]
            return True
        else:
            points.drop(index=points.index[0], inplace=True)
        return False

    def p_value_estimated(self, gap, L1, L2):
        """Estimate p-value based on the given gap and lengths."""
        N1 = gap['times'].sum()
        N = len(gap)
        self.points_init.loc[gap.index, 'times'] += 1
        m = len(self.path)
        a = (1 - self.score / m / self.grading[0]) * (N1 - m + 1) / N * (L1 - m + 1) * (L2 - m + 1) / L1 / L2
        return round(a, 4)
