import pandas as pd
import numpy as np
def filteredDataframe(csv, hitName):
    '''
    Creates a filtered dataframe for a specific toxin and only includes pairs of two
    Input- Toxin name
    Output- Filtered dataframe
    '''
    df = pd.read_csv(csv)
    domain = df[df['Hit Name'] == hitName]
    contig_counts = domain['Contig'].value_counts()
    valid_contigs = contig_counts[contig_counts == 2].index
    filtered_df = df[(df['Hit Name'] == hitName) & (df['Contig'].isin(valid_contigs))]
    return filtered_df

def uniquePairs(df, g1, g2):
   '''
   Finds and counts the unique pairs or combinations between two columns
   Input- dataframe and two columns
   Output- dataframe with the unique pairs and their counts
   '''
   unique_pairs = df[[g1, g2]].drop_duplicates().reset_index(drop=True)
   pairs = [tuple(x) for x in unique_pairs.values]
   counts = []
   for pair in pairs:
      pair_count = len(df[(df[g1] == pair[0]) & (df[g2] == pair[1])])
      counts.append(pair_count)
   returnFrame = pd.DataFrame({
      g1: [x[0] for x in pairs],
      g2: [x[1] for x in pairs],
      'Counts': counts})
   return returnFrame

def randomizedDF(df, n, col):
   '''
   Randomizes the columns in a dataframe "n" times
   Input- dataframe, number of permutations, column to randomize
   Output- dataframe
   '''
   for _ in range(n):
      df[col] = np.random.permutation(df[col].values)
      print(f"After randomization\n", df[col])
   return df