import os
import random
import sys
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from kneed import KneeLocator
import random

if __name__ == "__main__":
    start_time = time.time()
    random.seed(1)

    ############ scatter plot showing AS/FGA change in each patient with cutoff  ############
    print('Reading performance data ...')
    data = pd.read_csv('MSK_allInfo.csv')
    data  = data.iloc[:,list(range(21,121))+[2]]
    data_example = data.loc[data['CANCER_TYPE'] == 'Esophagogastric Cancer',]
    sampleNUM = data_example.shape[0]
    cancerType_list = list(set(data['CANCER_TYPE']))
    print(cancerType_list)

    print('Ploting ...')
    output_AUC_fn = '01.AS_FGA_cutoff_response.pdf' # png
    plt.figure(figsize=(3.2, 3.5))
    plt.rcParams['font.size'] = 8
    plt.rcParams["font.family"] = "Arial"
    plt.subplots_adjust(left=0.15, bottom=0.12, right=0.98, top=0.98, wspace=0.45, hspace=0.35)
    color_list = ['0.5', 'pink', 'k', 'r']

    ax1 = plt.subplot(211)
    x_list = np.arange(0.01,0.51,0.01)
    kn_bootstrap = []
    for i in range(sampleNUM):
        y_list = data_example.iloc[i,0:50]
        ax1.plot(x_list,y_list,'k-',alpha=0.5,linewidth=0.3)
    y_list = data_example.iloc[:, 0:50].mean(axis=0)
    ax1.plot(x_list, y_list, 'r-', alpha=1, linewidth=1)
    kn = KneeLocator(x_list, y_list, curve='convex', direction='decreasing')
    for i in range(1000):
        samples_temp = random.choices(range(sampleNUM), k = sampleNUM)
        y_temp = data_example.iloc[samples_temp, 0:50].mean(axis=0)
        kn_temp = KneeLocator(x_list, y_temp, curve='convex', direction='decreasing')
        kn_bootstrap.append(kn_temp.knee)
    kn_mean = np.mean(kn_bootstrap)#np.quantile(kn_bootstrap, 0.5)
    kn_025 = np.quantile(kn_bootstrap, 0.025)
    kn_975 = np.quantile(kn_bootstrap, 0.975)
    ax1.plot([kn_mean, kn_mean], [0, 1], 'g-', alpha=1, linewidth=1)
    ax1.axvspan(kn_025, kn_975, facecolor='g', alpha=0.5)
    print('AS elbow point (mean (95CI) overall) : %.3f (%.3f, %.3f)'%(kn_mean, kn_025, kn_975))
    ax1.set_ylabel("AS", color="k")
    ax1.set_xlabel("Calling cutoff (|log2 copy ratio|)", color="k")
    ax1.set_ylim([0, 1])
    ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax1.set_xticks([0,0.1,0.2,0.3,0.4,0.5])
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    #### AS cutoff of individual cancer types
    for ct in cancerType_list:
        data_ct = data.loc[data['CANCER_TYPE']==ct,]
        y_temp = data_ct.iloc[:, 0:50].mean(axis=0)
        kn = KneeLocator(x_list, y_temp, curve='convex', direction='decreasing')
        sampleNUM_ct = data_ct.shape[0]
        kn_bootstrap_ct = []
        for i in range(1000):
            samples_temp = random.choices(range(sampleNUM_ct), k=sampleNUM_ct)
            y_temp = data_ct.iloc[samples_temp, 0:50].mean(axis=0)
            kn_temp = KneeLocator(x_list, y_temp, curve='convex', direction='decreasing')
            kn_bootstrap_ct.append(kn_temp.knee)
        kn_mean = np.mean(kn_bootstrap_ct)  # np.quantile(kn_bootstrap, 0.5)
        kn_025 = np.quantile(kn_bootstrap_ct, 0.025)
        kn_975 = np.quantile(kn_bootstrap_ct, 0.975)
        print('%50s: %.3f (%.3f, %.3f)' % (ct, kn_mean, kn_025, kn_975))

    ax2 = plt.subplot(212)
    kn_bootstrap = []
    for i in range(sampleNUM):
        y_list = data_example.iloc[i,50:100]
        ax2.plot(x_list,y_list,'k-',alpha=0.5,linewidth=0.3)
    y_list = data_example.iloc[:, 50:100].mean(axis=0)
    ax2.plot(x_list, y_list, 'r-', alpha=1, linewidth=1)
    kn = KneeLocator(x_list, y_list, curve='convex', direction='decreasing')
    for i in range(1000):
        samples_temp = random.choices(range(sampleNUM), k = sampleNUM)
        y_temp = data_example.iloc[samples_temp, 50:100].mean(axis=0)
        kn_temp = KneeLocator(x_list, y_temp, curve='convex', direction='decreasing')
        kn_bootstrap.append(kn_temp.knee)
    kn_mean = np.mean(kn_bootstrap)#np.quantile(kn_bootstrap, 0.5)
    kn_025 = np.quantile(kn_bootstrap, 0.025)
    kn_975 = np.quantile(kn_bootstrap, 0.975)
    ax2.plot([kn_mean, kn_mean], [0, 1], 'b-', alpha=1, linewidth=1)
    ax2.plot([kn_025, kn_025], [0, 1], 'b--', alpha=1, linewidth=1)
    ax2.plot([kn_975, kn_975], [0, 1], 'b--', alpha=1, linewidth=1)
    print('FGA elbow point (mean (95CI) overall) : %.3f (%.3f, %.3f)'%(kn_mean, kn_025, kn_975))
    ax2.set_xlabel("Calling cutoff (|log2 copy ratio|)", color="k")
    ax2.set_ylabel("FGA", color="k")
    ax2.set_ylim([0, 1])
    ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax2.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)
    #### FGA cutoff of individual cancer types
    for ct in cancerType_list:
        data_ct = data.loc[data['CANCER_TYPE']==ct,]
        y_temp = data_ct.iloc[:, 50:100].mean(axis=0)
        kn = KneeLocator(x_list, y_temp, curve='convex', direction='decreasing')
        sampleNUM_ct = data_ct.shape[0]
        kn_bootstrap_ct = []
        for i in range(1000):
            samples_temp = random.choices(range(sampleNUM_ct), k=sampleNUM_ct)
            y_temp = data_ct.iloc[samples_temp, 50:100].mean(axis=0)
            kn_temp = KneeLocator(x_list, y_temp, curve='convex', direction='decreasing')
            kn_bootstrap_ct.append(kn_temp.knee)
        kn_mean = np.mean(kn_bootstrap_ct)  # np.quantile(kn_bootstrap, 0.5)
        kn_025 = np.quantile(kn_bootstrap_ct, 0.025)
        kn_975 = np.quantile(kn_bootstrap_ct, 0.975)
        print('%50s: %.3f (%.3f, %.3f)' % (ct, kn_mean, kn_025, kn_975))

    plt.savefig(output_AUC_fn) # , dpi=400
    plt.close()

    end_time = time.time()
    print('ALL done! Time used: ',end_time-start_time)