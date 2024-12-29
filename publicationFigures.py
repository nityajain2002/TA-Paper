import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def createFigure4a(ax, file, xx, y1, y2, title):
   df = pd.read_csv(file)
   df = df.dropna(subset=[xx, y1, y2])
   # Bar width and positions
   x = np.arange(len(df[xx]))  # the label locations
   width = 0.35  # width of the bars

   # Add bars for Value1 and Value2
   bars1 = ax.bar(x - width/2, df[y1], width, label='Observed Counts', color='purple')
   bars2 = ax.bar(x + width/2, df[y2], width, label='Randomized Counts', color='green')
   editedLabels = [
    label.replace("_antitoxin", "").replace("_toxin", "").replace("no domain found", "NDF").replace("No Domain found", "NDF").replace("No domain found", "NDF").replace("_antitox", "")
    for label in df[xx]]
   print(editedLabels)

   # Add x-axis labels, title, and legend
   ax.set_xlabel('Unique Pairs')
   ax.set_ylabel('Counts (log scale)')
   ax.set_title(title)
   ax.set_xticks(x)
   ax.set_xticklabels(editedLabels, rotation=45, ha='right')
   ax.set_yscale('log')
   ax.legend()

def setup4a():
    # Main plotting area for subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))

    # Subplot 1
    createFigure4a(
        axs[0, 0], "final-boss-parE.csv", 'Toxin Pairs', 
        'Observed Counts', 'Randomized Counts x20', 
        'Observed vs Randomized Toxin Counts: ParE'
    )

    # Subplot 2
    createFigure4a(
        axs[0, 1], "final-boss-parE.csv", 'Antitoxin Pairs', 
        'Observed Counts.1', 'Randomized Counts x20.1', 
        'Observed vs Randomized Antitoxin Counts: ParE'
    )

    # Subplot 3
    createFigure4a(
        axs[1, 0], "final-boss-pin-profesh.csv", 'Toxin Pairs', 
        'Observed Counts', 'Randomized Counts', 
        'Observed vs Randomized Toxins Counts: PIN'
    )

    # Subplot 4
    createFigure4a(
        axs[1, 1], "final-boss-pin-profesh.csv", 'Antitoxin Pairs', 
        'Observed Counts.1', 'Randomized Counts.1', 
        'Observed vs Randomized Antitoxins Counts: PIN'
    )

    # Adjust layout and display
    plt.tight_layout()
    plt.show()

def createFigure5a(): 
    '''
    C's code
    '''
    # File paths (update as needed)
    plasmid_data_path = 'Plasmid_Ori_positions_fixed.csv'  # Replace with the uploaded file name
    ta_system_data_path = 'Toxins.csv'  # Replace with the uploaded file name

    # Load datasets
    plasmid_data = pd.read_csv(plasmid_data_path)
    ta_system_data = pd.read_csv(ta_system_data_path)

    # Rename relevant columns for clarity and merging
    plasmid_data = plasmid_data.rename(columns={'contig_name': 'contig'})
    ta_system_data = ta_system_data.rename(columns={'Unnamed: 4': 'contig'})

    # Merge plasmid data with TA system data to check for TA system presence
    plasmid_data_with_ta = plasmid_data.merge(
        ta_system_data[['contig']].drop_duplicates(), 
        on='contig', 
        how='left', 
        indicator=True
    )
    plasmid_data_with_ta['has_ta_system'] = plasmid_data_with_ta['_merge'] == 'both'

    # Summarize data at the plasmid level
    plasmid_summary = plasmid_data_with_ta[['plasmid', 'contig', 'has_ta_system']].drop_duplicates()
    plasmid_level_analysis = plasmid_summary.groupby('plasmid').agg(
        Total_Plasmids=('contig', 'count'),
        Plasmids_with_TA_systems=('has_ta_system', 'sum')
    ).reset_index()

    # Calculate plasmids without TA systems and their percentages
    plasmid_level_analysis['Plasmids_without_TA_systems'] = (
        plasmid_level_analysis['Total_Plasmids'] - plasmid_level_analysis['Plasmids_with_TA_systems']
    )
    plasmid_level_analysis['Percentage_with_TA'] = (
        plasmid_level_analysis['Plasmids_with_TA_systems'] / plasmid_level_analysis['Total_Plasmids'] * 100
    )
    plasmid_level_analysis['Percentage_without_TA'] = (
        plasmid_level_analysis['Plasmids_without_TA_systems'] / plasmid_level_analysis['Total_Plasmids'] * 100
    )

    # Filter plasmid categories with more than 50 plasmids
    filtered_plasmid_analysis = plasmid_level_analysis[plasmid_level_analysis['Total_Plasmids'] > 50]

    # Save the filtered data to CSV
    output_path = 'filtered_plasmid_analysis_from_raw.csv'
    filtered_plasmid_analysis.to_csv(output_path, index=False)

    print(f"Filtered plasmid analysis saved to: {output_path}")
    # Update the analysis to ensure plasmid-level statistics consider the adjusted TA system assignments

    # Recompute the plasmid-level analysis using the corrected `has_ta_system`
    plasmid_summary_updated = plasmid_data_with_ta[['plasmid', 'contig', 'has_ta_system']].drop_duplicates()

    plasmid_level_analysis_updated = plasmid_summary_updated.groupby('plasmid').agg(
        Total_Plasmids=('contig', 'count'),
        Plasmids_with_TA_systems=('has_ta_system', 'sum')
    ).reset_index()

    plasmid_level_analysis_updated['Plasmids_without_TA_systems'] = (
        plasmid_level_analysis_updated['Total_Plasmids'] - plasmid_level_analysis_updated['Plasmids_with_TA_systems']
    )
    plasmid_level_analysis_updated['Percentage_with_TA'] = (
        plasmid_level_analysis_updated['Plasmids_with_TA_systems'] / plasmid_level_analysis_updated['Total_Plasmids'] * 100
    )
    plasmid_level_analysis_updated['Percentage_without_TA'] = (
        plasmid_level_analysis_updated['Plasmids_without_TA_systems'] / plasmid_level_analysis_updated['Total_Plasmids'] * 100
    )

    # Filter plasmids with >50 occurrences
    filtered_plasmid_level_analysis_updated = plasmid_level_analysis_updated[
        plasmid_level_analysis_updated['Total_Plasmids'] > 50
    ]

    # Save updated results to CSV
    output_path_updated = '/mnt/data/filtered_plasmid_analysis_updated.csv'
    filtered_plasmid_level_analysis_updated.to_csv(output_path_updated, index=False)

    # Display updated table to user
    tools.display_dataframe_to_user(name="Filtered Plasmid Analysis with Corrected TA System Presence", dataframe=filtered_plasmid_level_analysis_updated)

    # Recreate the updated graph
    categories_updated = filtered_plasmid_level_analysis_updated['plasmid']
    percentages_with_ta_updated = filtered_plasmid_level_analysis_updated['Percentage_with_TA']
    percentages_without_ta_updated = filtered_plasmid_level_analysis_updated['Percentage_without_TA']

    # Create updated stacked bar chart
    plt.figure(figsize=(12, 8))
    plt.bar(categories_updated, percentages_with_ta_updated, label='With TA Systems', color='#6AA84F')
    plt.bar(categories_updated, percentages_without_ta_updated, bottom=percentages_with_ta_updated, label='Without TA Systems', color='#674EA7')

    # Format the chart
    plt.xlabel('Plasmid Categories', fontsize=14)
    plt.ylabel('Percentage (%)', fontsize=14)
    plt.title('Stacked Percentage of Plasmids with TA Systems by Category (Corrected Data)', fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.legend(title='Legend', fontsize=12)
    plt.grid(False)
    plt.tight_layout()

    # Display the updated chart
    plt.show()

    output_path_updated

def createFigure5b():
    '''
    C's code
    '''
    # Load the data
    plasmid_data_path = 'Plasmid_Ori_positions_fixed - plasmid_data_extracted_with_lengths.csv'  # Update the path if needed
    ta_system_data_path = 'toxins_within_5kb_ori.xlsx'  # Update the path if needed

    plasmid_data = pd.read_csv(plasmid_data_path)
    ta_system_data = pd.read_excel(ta_system_data_path)

    # Clean and merge datasets
    plasmid_data = plasmid_data[['plasmid', 'contig_name', 'positions_in_contig_start', 'positions_in_contig_end']]
    ta_system_data = ta_system_data[['contig_name', 'plasmid', 'within_5kb']]

    # Merge on contig_name and plasmid
    merged_data = pd.merge(plasmid_data, ta_system_data, on=['contig_name', 'plasmid'], how='left')

    # Treat missing within_5kb as False
    merged_data['within_5kb'] = merged_data['within_5kb'].fillna(False).astype(bool)

    # Group by unique ORIs (contig_name, plasmid) and aggregate within_5kb
    collapsed_data = merged_data.groupby(['contig_name', 'plasmid']).agg(
        within_5kb=('within_5kb', 'any')  # True if any TA system is within 5kb
    ).reset_index()

    # Filter plasmid categories with >50 occurrences
    category_counts = collapsed_data['plasmid'].value_counts()
    filtered_categories = category_counts[category_counts > 50].index
    filtered_data = collapsed_data[collapsed_data['plasmid'].isin(filtered_categories)]

    # Calculate metrics for each filtered ORI category
    results = filtered_data.groupby('plasmid').agg(
        Total_ORIs=('contig_name', 'count'),
        ORIs_with_TA_systems=('within_5kb', 'sum')
    ).reset_index()

    results['ORIs_without_TA_systems'] = results['Total_ORIs'] - results['ORIs_with_TA_systems']
    results['Percentage_with_TA'] = (results['ORIs_with_TA_systems'] / results['Total_ORIs']) * 100
    results['Percentage_without_TA'] = (results['ORIs_without_TA_systems'] / results['Total_ORIs']) * 100

    # Prepare data for visualization
    categories = results['plasmid']
    percentages_with_ta = results['Percentage_with_TA']
    percentages_without_ta = results['Percentage_without_TA']

    # Create the stacked bar chart
    plt.figure(figsize=(12, 8))
    plt.bar(categories, percentages_with_ta, label='With TA Systems Nearby', color='#6AA84F')  # Green
    plt.bar(categories, percentages_without_ta, bottom=percentages_with_ta, label='Without TA Systems Nearby', color='#674EA7')  # Purple

    # Formatting the chart
    plt.xlabel('ORI Categories (Plasmid)', fontsize=14)
    plt.ylabel('Percentage (%)', fontsize=14)
    plt.title('Stacked Percentage of ORIs Close to TA Systems by ORI Category (Occurrences > 50)', fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.legend(title='Legend', fontsize=12)
    plt.grid(False)  # Remove gridlines
    plt.tight_layout()

    # Display the chart
    plt.show()

    # Save the resulting table to a CSV file
    output_path = 'filtered_ori_analysis.csv'  # Update the path if needed
    results.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")