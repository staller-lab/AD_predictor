import numpy as np
import protfasta
import pandas as pd
import matplotlib.pyplot as plt
import random
from matplotlib.pyplot import figure
import seaborn as sns
sns.set_theme()
sns.set(style = "white", rc={'figure.figsize':(6.5,4)}, font_scale = 1.2)

#Control colors:
LambertTFs_color = "#cccccc"
GSL_color = "#FFD972"
predictions_color = "#BADDE1"



def add_AAcounts_and_allhydros(df):
    AAs=["A","I","L","M","V","F","W","Y","N","C","Q","S","T","D","E","R","H","K","G","P"]
    for aa in AAs:
            df[aa]= df.apply(lambda row: row['ProteinWindowSeq'].count(aa),axis=1)
    df['AllHydros'] = df[["W","Y","F","L"]].sum(axis=1)

def tile_and_plot(column, figsize=(6,4), dpi=80, xtick_interval=10, xtick_labelsize=10, bins_interval=10):
    figure(figsize=figsize, dpi=dpi)
    plt.title(column+" Histogram")
    plt.xlabel(column)
    plt.ylabel('Percent')

    bins=np.arange(min(min(tiled_Lambert[column]),
                             min(tiled_GS_list[column]),
                             min(tiled_predictions[column]))
                   ,max(max(tiled_Lambert[column]),
                             max(tiled_GS_list[column]),
                             max(tiled_predictions[column])),
                   bins_interval)

    plt.xticks(np.arange(min(min(tiled_Lambert[column]),
                             min(tiled_GS_list[column]),
                             min(tiled_predictions[column]))
                   ,max(max(tiled_Lambert[column]),
                             max(tiled_GS_list[column]),
                             max(tiled_predictions[column])),
                   xtick_interval))
    
    plt.tick_params(axis='x', which='major', labelsize=xtick_labelsize)

    tiled_Lambert[column].hist(bins=bins, alpha=0.7, grid=False, weights=np.ones(len(tiled_Lambert.index)) / len(tiled_Lambert.index),label="Lambert TFs", color = LambertTFs_color)
    tiled_GS_list[column].hist(bins=bins, alpha=0.7,grid=False, weights=np.ones(len(tiled_GS_list.index)) / len(tiled_GS_list.index),label="Gold Standard List", color = GSL_color)
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals = 0))

    plt.legend()
    plt.show()
    
# Pass in a dataframe of sequences, will plot the dataframe against human TFs
def plot_heatmap(active_supp_ds_4, table_name = "temp", normalize = True, AAs = ['W', 'F', 'Y', 'L']):
    active_supp_ds_4["GeneID"]=[">sp|"+str(a)+"|"+str(b)+"|"+str(c)+"|"+str(d)
                   for a,b,c,d in zip(active_supp_ds_4["uniprotID"],
                                        active_supp_ds_4["GeneName"],
                                        active_supp_ds_4["Start"],
                                        active_supp_ds_4["End"])]

    active_supp_ds_4_dict=pd.Series(active_supp_ds_4.ProteinRegionSeq.values,index=active_supp_ds_4.GeneID).to_dict()

    protfasta.write_fasta(active_supp_ds_4_dict, '../data/'+table_name+'.fasta')
    
    Prop1Span = np.arange(-39,21,1)# hard coded to match shape of precomputed heatmap
    Prop2Span = np.arange(0,20,1) # hard coded to match shape of precomputed heatmap
    Prop1SpanMesh, Prop2SpanMesh = np.meshgrid(Prop1Span,Prop2Span)

    #Human TFs 
    if AAs == ['W', 'F', 'Y', 'L']:
        Counts_df=pd.read_csv("../data/HumanTFTiles_Counts_v2.csv")
        Counts=Counts_df.values
        Counts=Counts[:,1:]
        Z =np.log10(Counts.transpose())
    else: 
        propset =['Charge','AllHydros']
        tempDF = AD_predictor_tools.makeTilingDF('../data/LambertTFs.fasta', window_size = 39, window_spacing = 1, AAs = AAs)
        tempDF['AllHydros'] = tempDF[AAs].sum(axis=1)
        Prop1SpanMesh,Prop2SpanMesh,Counts_df=AD_predictor_tools.MakeMesh(tempDF,propset, auto_propspan = False)
        Z=np.log10(Counts_df.transpose())

    #Active
    propset =['Charge','AllHydros']
    tempDF = AD_predictor_tools.makeTilingDF('../data/'+table_name+'.fasta',window_size = 39, window_spacing = 1, AAs = AAs)
    tempDF['AllHydros'] = tempDF[AAs].sum(axis=1)
    Prop1SpanMesh2,Prop2SpanMesh2,Counts_df2=AD_predictor_tools.MakeMesh(tempDF,propset, auto_propspan = False)
    Z2=np.log10(Counts_df2.transpose())

    X,Y= Prop1SpanMesh, Prop2SpanMesh,

    plt.figure(figsize=(8,3),dpi=300)
    plt.xlabel('Charge')
    plt.ylabel('Number of '+','.join(AAs))

    if normalize:
        c_cs2 = Z2/Z
    else:
        c_cs2 = Z2
        
    cs1  = plt.scatter(X,Y, c=Z,cmap="Greys", alpha=1,label='All Human TF regions',marker='s',s=21,linewidths=0)
    cs2  = plt.scatter(X,Y, c=c_cs2,cmap="Blues", alpha=1,label='All Human TF regions',marker='s',s=21,linewidths=0)

    LowerCorner= AD_predictor_tools.get_bounds("VP16", AAs)
    UpperCorner= AD_predictor_tools.get_bounds("CITED2", AAs)
    plt.scatter([LowerCorner[0],UpperCorner[0]],[LowerCorner[1],UpperCorner[1]],c='r')

    # Add a color bar which maps values to colors.
    plt.colorbar(cs1, shrink=1, aspect=20,label='log10(Active Tested Abundance)',ticks=[0,1,2,3,4],drawedges=False)
    plt.colorbar(cs2, shrink=1, aspect=20,label='log10(Active Abundance)',ticks=[0,1,2,3,4],drawedges=False)
    #plt.colorbar(cs3, shrink=1, aspect=20,label='log10(Inactive Abundance)',ticks=[0,1,2,3,4],drawedges=False)    
    # 
def tiled_df_to_abundance_df(tiled_Lambert):
    abundance_series = tiled_Lambert.groupby(['Charge', 'AllHydros']).size()
    abundance_df = abundance_series.to_frame()
    abundance_df.reset_index(inplace=True)
    abundance_df = abundance_df.rename(columns = {0: "count"})
    abundance_df["log10count"] = [math.log10(abundance_df["count"][i]) for i in range(len(abundance_df.index))]
    return abundance_df

def plot_marg_hists(tiled_df, label, color, x = "Charge", y = "AllHydros", c = "log10count", cmap = "Greys"):
    abundance_df = tiled_df_to_abundance_df(tiled_df)
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    rect_scatter = [left, bottom, width, height]

    # start with a square Figure
    fig = plt.figure(figsize=(8, 8))
    
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(data = abundance_df, x = x, y = y, c = c, cmap = cmap,
                    alpha=1,label= label,marker='s',s=100,linewidths=0)

    ax_histx.hist(tiled_df[x], bins= np.arange(min(tiled_df[x]), max(tiled_Lambert[x])), color = color)
    ax_histy.hist(tiled_df[y], bins= np.arange(min(tiled_df[y]), max(tiled_Lambert[y])), orientation='horizontal', color = color)

    ax.set_xlabel(x)
    ax.set_ylabel(y)

    plt.show()

def plot_two_marg_hists(tiled_Lambert, tiled_GS_list, cmap1 = "Greys", cmap2 = "Oranges", color2= GSL_color):
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    rect_scatter = [left, bottom, width, height]

    # start with a square Figure
    fig = plt.figure(figsize=(8, 8))
    
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    min_x = min([min(tiled_Lambert["Charge"]), min(tiled_GS_list["Charge"])])
    max_x = max([max(tiled_Lambert["Charge"]), max(tiled_GS_list["Charge"])])

    min_y = min(min(tiled_GS_list["AllHydros"]), min(tiled_Lambert["AllHydros"]))
    max_y = max(max(tiled_GS_list["AllHydros"]), max(tiled_Lambert["AllHydros"]))
    
    LambertTFs_abundanceDF = tiled_df_to_abundance_df(tiled_Lambert)
    GSL_abundanceDF = tiled_df_to_abundance_df(tiled_GS_list)

    ax.scatter(data = LambertTFs_abundanceDF, x = "Charge", y = "AllHydros", c = "log10count", cmap = cmap1,
                    alpha=1,label='All Human TF regions',marker='s',s=100,linewidths=0)
    ax.scatter(data = GSL_abundanceDF, x = "Charge", y = "AllHydros", c = "log10count", cmap = cmap2,
                    alpha=1,label='GSL regions',marker='s',s=100,linewidths=0)

    ax_histx.hist(tiled_Lambert["Charge"], bins= np.arange(min_x, max_x), color = LambertTFs_color, alpha = 0.5, density = True)
    ax_histx.hist(tiled_GS_list["Charge"], bins= np.arange(min_x, max_x), color = color2, alpha = 0.5, density = True)

    ax_histy.hist(tiled_Lambert["AllHydros"], bins= np.arange(min_y, max_y), orientation='horizontal', color = LambertTFs_color, alpha = 0.5, density = True)
    ax_histy.hist(tiled_GS_list["AllHydros"], bins= np.arange(min_y, max_y), orientation='horizontal', color = color2, alpha = 0.5, density = True)

    ax.set_xlabel("Charge")
    ax.set_ylabel("AllHydros")

    plt.show()