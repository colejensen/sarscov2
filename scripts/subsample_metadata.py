# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from epiweeks import Week
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Subsample nextstrain metadata entries keeping only pre-selected samples",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Metadata file from nextstrain")
    parser.add_argument("--keep", required=False, help="List of samples to keep, in all instances")
    parser.add_argument("--remove", required=False, help="List of samples to remove, in all instances")
    parser.add_argument("--scheme", required=True, help="Subsampling scheme")
    parser.add_argument("--output", required=True, help="Selected list of samples")
    parser.add_argument("--report", required=True, help="Report listing samples per category")
    args = parser.parse_args()

    metadata = args.metadata
    keep = args.keep
    remove = args.remove
    scheme = args.scheme
    report = args.report
    output = args.output

    # metadata = path + 'metadata_geo.tsv'
    # keep = path + 'keep.txt'
    # remove = path + 'remove.txt'
    # scheme = path + 'subsampling_scheme.tsv'
    # report = path + 'report.tsv'
    # output = path + 'selected_strains.txt'


    today = time.strftime('%Y-%m-%d', time.gmtime())

    # force genomes to be kept in final dataset
    to_keep = [strain.strip() for strain in open(keep, 'r').readlines() if strain[0] not in ['#', '\n']]

    # subsampling scheme
    dfS = pd.read_csv(scheme, encoding='utf-8', sep='\t', dtype=str)

    results = {}

    ### IGNORING SAMPLES

    # list of rows to be ignored
    ignore = {}
    for idx, val in dfS.loc[dfS['purpose'] == 'ignore', 'name'].to_dict().items():
        key = dfS.iloc[idx]['level']
        if key not in ignore.keys():
            ignore[key] = [val]
        else:
            ignore[key].append(val)
    # print(ignore)

    # nextstrain metadata
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    try:
        dfN = dfN[['strain', 'gisaid_epi_isl', 'genbank_accession', 'date', 'region', 'country', 'division', 'location',
                      'region_exposure', 'country_exposure', 'division_exposure', 'originating_lab', 'submitting_lab']]
    except:
        pass

    dfN.fillna('', inplace=True) # replace empty values by blank

    # drop lines if samples are set to be ignored
    for column, names in ignore.items():
        dfN = dfN[~dfN[column].isin(names)]
    # print(sorted(dfN['country'].unique()))


    # drop rows listed in remove.txt
    to_remove = [strain.strip() for strain in open(remove, 'r').readlines()]
    dfN = dfN[~dfN['strain'].isin(to_remove)]

    # prevent samples already selected in keep.txt from being resampled
    dfN = dfN[~dfN['strain'].isin(to_keep)]

    # drop rows with incomplete dates
    dfN = dfN[dfN['date'].apply(lambda x: len(x.split('-')) == 3)] # accept only full dates
    dfN = dfN[dfN['date'].apply(lambda x: 'X' not in x)] # exclude -XX-XX missing dates

    # convert string dates into date format
    dfN['date'] = pd.to_datetime(dfN['date']) # coverting to datetime format
    dfN = dfN.sort_values(by='date')  # sorting lines by date
    start, end = dfN['date'].min(), today


    # drop columns where 'country' and 'country_exposure' disagree
    dfN['same_division'] = np.where(dfN['division'] == dfN['division_exposure'], 'yes', 'no') # compare values
    dfN.loc[dfN['division_exposure'] == '', 'same_division'] = 'yes'
    dfN = dfN[dfN['same_division'].apply(lambda x: 'yes' in x)] # exclude rows with conflicting place of origin
    # print(dfN[['same_country', 'country', 'country_exposure']])
    dfN.drop(columns=['same_division'])

    # print(dfN[['strain', 'date']].iloc[[0, -1]])

    # get epiweek end date, create column
    dfN['date'] = pd.to_datetime(dfN['date'], errors='coerce')
    dfN['epiweek'] = dfN['date'].apply(lambda x: Week.fromdate(x, system="cdc").enddate())


    ## SAMPLE FOCAL AND CONTEXTUAL SEQUENCES
    purposes = ['focus', 'context']
    subsamplers = [] # list of focal and contextual categories
    for category in purposes:
        query = {} # create a dict for each 'purpose'
        for idx, val in dfS.loc[dfS['purpose'] == category, 'name'].to_dict().items():
            key = dfS.iloc[idx]['level']
            if key not in query.keys():
                query[key] = [val]
            else:
                query[key].append(val)
            subsamplers.append(query)

    # perform subsampling
    for scheme_dict in subsamplers:
        # print(scheme_dict)
        # group by level
        for level, names in scheme_dict.items():
            # print(level, names)
            # create dictionary for level
            if level not in results:
                results[level] = {}

            for name in names:
                # add place name as key in its corresponding level in dict
                if name not in results[level].keys():
                    results[level][name] = []

            glevel = dfN.groupby(level)
            for name, dfLevel in glevel:
                if name in names: # check if name is among focal places
                    min_date, max_date = start, end

                    # define new temporal boundaries, if provided
                    new_start = dfS.loc[dfS['name'] == name, 'start'].values[0]
                    new_end = dfS.loc[dfS['name'] == name, 'end'].values[0]
                    if not pd.isna(new_start):
                        min_date = new_start
                    if not pd.isna(new_end):
                        max_date = new_end

                    # drop any row with dates outside the start/end dates
                    mask = (dfLevel['date'] > min_date) & (dfLevel['date'] <= max_date)
                    dfLevel = dfLevel.loc[mask]  # apply mask

                    gEpiweek = dfLevel.groupby('epiweek')
                    sample_size = int(dfS.loc[dfS['name'] == name, 'size'])
                    total_genomes = dfLevel[level].count()
                    for epiweek, dfEpiweek in gEpiweek:
                        bin_pool = dfEpiweek['epiweek'].count() # genomes in bin
                        sampled = int(np.ceil((bin_pool/total_genomes) * sample_size)) # proportion sampled from bin
                        # print(epiweek, '-', sampled, '/', bin_pool)
                        if sampled > bin_pool: # if # requested samples higher than available genomes, get all
                            sampled = bin_pool

                        # selector
                        random_subset = dfEpiweek.sample(n=sampled)
                        selected = random_subset['strain'].to_list()
                        results[level][name] = results[level][name] + selected

                    # drop pre-selected samples to prevent duplicates
                    dfN = dfN[~dfN[level].isin([name])]


    ### EXPORT RESULTS
    print('\n# Genomes sampled per category in subsampling scheme\n')
    exported = []
    genome_count = 0
    with open(output, 'w') as outfile, open(report, 'w') as outfile2:
        outfile.write('# Genomes selected on ' + today + '\n')
        outfile2.write('sample_size' + '\t' + 'category' + '\n')
        for level, name in results.items():
            for place, entries in name.items():
                if len(entries) > 1:
                    genome_count += len(entries)
                    entry = str(len(entries)) + '\t' + place + ' (' + level + ')'
                    outfile2.write(entry + '\n')
                    print(entry)

                    for strain_name in entries:
                        if strain_name not in [exported + to_keep]:
                            outfile.write(strain_name + '\n')
                            exported.append(strain_name)

        # report selected samples listed in keep.txt
        print('- ' + str(len(to_keep)) + ' genomes added from pre-selected list\n')
        outfile.write('\n# Pre-existing genomes listed in keep.txt\n')
        for strain in to_keep:
            if strain not in exported:
                outfile.write(strain + '\n')
                exported.append(strain)

        print('\n# Genomes matching the criteria below were not found')
        outfile2.write('\n# Categories not found\n')
        for level, name in results.items():
            for place, entries in name.items():
                if len(entries) == 0:
                    entry = str(len(entries)) + '\t' + place + ' (' + level + ')'
                    outfile2.write(entry + '\n')
                    print(entry)

    print('\n' + str(genome_count) + ' genomes exported according to subsampling scheme\n')
