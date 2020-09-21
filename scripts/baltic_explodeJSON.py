import baltic as bt
import pandas as pd
from datetime import datetime
import matplotlib as mpl
import numpy as np


import re, copy, json

import matplotlib as mpl
from matplotlib import pyplot as plt

try:
    from StringIO import StringIO as sio
    from cStringIO import StringIO as csio
except ImportError:
    from io import StringIO as sio
    from io import BytesIO as csio

import requests
from Bio import Phylo

path = 'C:/Users/28Cha/Glab Dropbox/Cole Jensen/sarscov2/auspice'
# treeFile = 'aligned.tree'
# bstree = Phylo.read(path + treeFile, 'newick')

nx_tree = json.load(open(path + 'sarscov2_subsample1.json','r'))

tree = nx_tree['tree'] # tree data
meta = nx_tree['meta'] # metadata

# print(tree)

traitName = 'region' # trait used for colouring the tree
xspan = (0, 0.003)

typeface = 'Helvetica'
mpl.rcParams['font.weight'] = 300
mpl.rcParams['axes.labelweight'] = 300
mpl.rcParams['font.family'] = typeface
mpl.rcParams['font.size'] = 12

json_translation = {'absoluteTime': lambda k: k.traits['node_attrs']['num_date']['value'], 'name': 'name'}  # find correct attributes in JSON, height and name are required at a minimum
# json_translation = {'absoluteTime': lambda k: k.traits['node_attrs']['div'], 'name': 'name'}  # find correct attributes in JSON, height and name are required at a minimum
json_meta = {'file': meta, 'traitName': traitName}  # import the meta file used on nextstrain.org
ll, meta = bt.loadJSON(nx_tree, json_translation, json_meta) # give loadJSON the name of the tree file, the translation dictionary and (optionally) the meta file
ll.sortBranches(descending=True)

print(type(ll))
# print('\n\n\n\n\n')
# print(json_meta['file']['colorings'])
colors = {'ancestor': '#CECECE'}
for metaDict in json_meta['file']['colorings']:
    # print(metaDict)

    if metaDict['key'] == traitName:
        # print(metaDict['scale'])
        for color_pair in metaDict['scale']:
            label, hex = color_pair
            colors[label] = hex

# print(colors)

branchWidth = 2
tipSize = 30

ll.root.parent.traits[traitName] = 'ancestor'  ## add fake trait to root of the tree
subtype_trees = {'ancestor': []}
for element in colors.keys():
    subtype_trees.update({element: []})


results = []
column_list = ['source', 'sink', 'tmrca', 'size', 'tips']
for l in ll.Objects:
    k = l
    kp = l.parent

    ## get current node's and its parent's trait states
    kloc = k.traits[traitName]
    if traitName in k.parent.traits:
        kploc = kp.traits[traitName]
        kpc = kploc
    else:
        kploc = 'ancestor'
        kpc = 'ancestor'

    kc = kloc
    if kc != kpc:
        fields = {column: '' for column in column_list}
        # try:
        traverse_condition = lambda w: w.traits[traitName] == kc
        # print('subtree resulting from %s > %s switch, traversing within %s' % (kpc, kc, kc))
        print(kpc + ' > ' + kc)

        subtree = ll.subtree(k, traverse_condition=traverse_condition)  ## this function returns a new baltic object that contains a trait-traversed subtree, starting from node k, for as long as the traversal stays within the starting trait value state
        print(subtree)
        subtree.traverse_tree()
        # subtree.drawTree(verbose=True)
        subtree.sortBranches()
        # all_dates = [k.absoluteTime for k in subtree.Objects]
        # print(sorted(all_dates))

        # stats
        source = kpc
        sink = kc
        tips = []
        for k in subtree.Objects:
            if k.branchType == 'leaf':
                tips.append(k.name)
            # if k.branchType=='node':
            #     print(k.children)
                # print(k.absoluteTime)
        tmrca = subtree.root.absoluteTime

        if source == 'ancestor':
            tmrca = ll.root.absoluteTime
        print(source, sink, tmrca, len(tips))
        print('')

        date_decimal = float(tmrca)
        date_year = int(date_decimal)
        year_fraction = date_decimal - date_year
        days = int(year_fraction * 365)
        # now convert the year and days into string and then into date (there is probably a better way to do this - without the string step)
        tmrca = datetime.strptime("{}-{}".format(date_year, days),"%Y-%j")


        # fill dict with content to export
        fields['source'] = source
        fields['sink'] = sink
        fields['tmrca'] = tmrca.strftime('%Y-%m-%d')
        fields['size'] = str(len(tips))
        fields['tips'] = ', '.join(tips)
        results.append(fields)

        # tree_strings[kc].append(subtree.toString())  ## remember subtree string, subtree object itself
        subtype_trees[kc].append((kpc, subtree))

        print(fields)

# df = pd.DataFrame(results, columns=list(column_list))
# df = df.sort_values(by=['source'])
# # df = df.sort_values(by=['sink'])
# # df['tmrca'] = df['tmrca'].dt.strftime('%Y-%m-%d')
# # df.sort_values(by=['sink'])
# output = path + 'demes.tsv'
# df.to_csv(output, sep='\t', index=False)
#
#
#
#
#
#
#
# fig, ax = plt.subplots(figsize=(10, 15), facecolor='w')
#
# tipSize = 20
# cumulative_y = 0
#
# x_attr = lambda k: k.absoluteTime
# # c_func=lambda k: ll.cmap[k.traits[json_meta['traitName']]] if traitName in k.traits.keys() else '#CECECE' if k.branchType == 'leaf' else 'grey'  ## colour will be determined by the trait of choice, using the colour map defined by the meta json
# c_func = lambda k: colors[k.traits[traitName]] # colour of branches
#
# s_func = lambda k: tipSize
# z_func = lambda k: 100
#
# su_func = lambda k: tipSize + 30
# cu_func = lambda k: 'k'
# zu_func = lambda k: 99
#
#
#
# # for subtype in ['V', 'Y']:  ## iterate over trait values
# for subtype in subtype_trees.keys(): #['Oceania', 'Eastern Asia', 'Southeastern Asia', 'Northern Europe', 'Western Europe', 'Southern Europe', 'North America']:  ## iterate over trait values
#
#     # print(subtype, tree_strings[subtype])
#     for tr in subtype_trees[subtype]:
#     # for t, tr in enumerate(sorted(subtype_trees[subtype], key=lambda x: (-1 * x[1].root.absoluteTime, len(x[1].Objects)))):  ## iterate over extracted subtrees sorted according to their root height and their size
#         print(subtype, tr[1].root.absoluteTime)
#         origin, loc_tree = tr  ## get origin of subtree, subtree itself
#         y_attr = lambda k: k.y + cumulative_y
#
#         loc_tree.plotTree(ax, x_attr=x_attr, y_attr=y_attr, colour_function=c_func)
#         loc_tree.plotPoints(ax, x_attr=x_attr, y_attr=y_attr, size_function=s_func,
#                             colour_function=c_func, zorder_function=z_func)
#         loc_tree.plotPoints(ax, x_attr=x_attr, y_attr=y_attr, size_function=su_func,
#                             colour_function=cu_func, zorder_function=zu_func)
#
#         oriC = colors[origin] if origin in colors else 'grey'
#         # oriC = ll.cmap[k.traits[json_meta['traitName']]] if traitName in k.traits.keys() else '#CECECE' if k.branchType == 'leaf' else 'grey'
#
#         # oriX = loc_tree.root.absoluteTime
#         # oriY = loc_tree.root.y + cumulative_y
#         # print(oriX, oriY)
#         # # print(cumulative_y)
#         # ax.scatter(oriX, oriY, 100, facecolor=oriC, edgecolor='w', lw=1, zorder=200)  ## add big circle at base of tree to indicate origin
#
#         cumulative_y += loc_tree.ySpan + 5  ## increment y displacement
#
# ax.xaxis.tick_bottom()
# ax.yaxis.tick_left()
#
# [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
#
# ax.tick_params(axis='y', size=0)
# ax.set_yticklabels([])
# ax.set_ylim(-5, cumulative_y)
# # plt.xticks(np.arange(min(xspan), max(xspan), 0.001))
#
# plt.show()
