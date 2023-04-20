import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator
import random
from scipy.stats import linregress

if __name__ == "__main__":
    start_time = time.time()
    random.seed(1)

    ############ scatter plot showing AS/FGA change in each patient with cutoff  ############
    print('Reading data ...')

    data = pd.read_csv('MSK_allInfo.csv')
    data  = data.iloc[:,list(range(21,121))+[1,2,10]]
    sampleNUM = data.shape[0]
    data['CANCER_TYPE'] = data['CANCER_TYPE'].str.lower().str.capitalize()
    cancerType_list = list(set(data['CANCER_TYPE']))

    print('Ploting ...')
    output_fig_fn = 'tumorPurity_vs_elbowPts.pdf' # png
    plt.figure(figsize=(6, 4))
    plt.rcParams['font.size'] = 12
    plt.rcParams["font.family"] = "Arial"
    plt.subplots_adjust(left=0.12, bottom=0.12, right=0.5, top=0.98, wspace=0.45, hspace=0.35)

    ax1 = plt.subplot(211)
    data_plot = data[['TUMOR_PURITY','CANCER_TYPE']]
    data_plot['TUMOR_PURITY'] = pd.to_numeric(data_plot['TUMOR_PURITY'], errors='coerce')
    data_plot = data_plot.dropna(subset=['TUMOR_PURITY'])

    result = data_plot.groupby('CANCER_TYPE')['TUMOR_PURITY'].agg(['count', 'mean', 'std'])
    result.rename(columns={"count": "sampleNum", "mean": "purity", "std": "purity_std"}, inplace=True)
    result['AS_elbowPt'] = 0
    result['AS_elbowPt_std'] = 0
    result['AS_elbowPt_low'] = 0
    result['AS_elbowPt_up'] = 0
    result['FGA_elbowPt'] = 0
    result['FGA_elbowPt_std'] = 0
    result['FGA_elbowPt_low'] = 0
    result['FGA_elbowPt_up'] = 0

    x_list = np.arange(0.01,0.51,0.01)
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
        kn_std = np.std(kn_bootstrap_ct)
        kn_low = np.quantile(kn_bootstrap_ct, 0.025)
        kn_up = np.quantile(kn_bootstrap_ct, 0.975)
        print('%50s: %.3f (%.3f, %.3f)' % (ct, kn_mean, kn_low, kn_up))
        result.loc[ct, 'AS_elbowPt'] = kn_mean
        result.loc[ct, 'AS_elbowPt_std'] = kn_std
        result.loc[ct, 'AS_elbowPt_low'] = kn_low
        result.loc[ct, 'AS_elbowPt_up'] = kn_up

    data_plot = result
    groups = data_plot.groupby('CANCER_TYPE')
    for name, group in groups:
        ax1.errorbar(group['purity'], group['AS_elbowPt'], xerr=group['purity_std']/np.sqrt(group['sampleNum']),
                     yerr=group['AS_elbowPt_std']/np.sqrt(group['sampleNum']), fmt='o', label=name, markersize=2, linewidth=0.4)

    # Perform linear regression
    slope, intercept, r_value, pvalue, _ = linregress(data_plot['purity'], data_plot['AS_elbowPt'])
    # Create the trend line
    x = np.linspace(min(data_plot['purity']), max(data_plot['purity']), 100)
    y = slope * x + intercept
    ax1.plot(x, y, color='k', zorder=0)
    ax1.text(20, 0.25, 'r = %.2f\np = %.3f'%(r_value, pvalue), fontsize=12)

    ax1.set_xlabel("Tumor purity (%)", color="k")
    ax1.set_ylabel("AS elbow point", color="k")
    ax1.legend(frameon = False, bbox_to_anchor=(1, 1), loc='upper left')
    ax1.set_xlim([10, 90])
    ax1.set_ylim([0.1, 0.3])
    ax1.set_yticks([0.1,0.2,0.3])
    ax1.set_xticks([10, 30,50,70,90])
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)

    ax2 = plt.subplot(212)
    #### FGA cutoff of individual cancer types
    for ct in cancerType_list:
        data_ct = data.loc[data['CANCER_TYPE'] == ct,]
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
        kn_std = np.std(kn_bootstrap_ct)
        kn_low = np.quantile(kn_bootstrap_ct, 0.025)
        kn_up = np.quantile(kn_bootstrap_ct, 0.975)
        print('%50s: %.3f (%.3f, %.3f)' % (ct, kn_mean, kn_low, kn_up))
        result.loc[ct, 'FGA_elbowPt'] = kn_mean
        result.loc[ct, 'FGA_elbowPt_std'] = kn_std
        result.loc[ct, 'FGA_elbowPt_low'] = kn_low
        result.loc[ct, 'FGA_elbowPt_up'] = kn_up

    data_plot = result
    groups = data_plot.groupby('CANCER_TYPE')
    for name, group in groups:
        ax2.errorbar(group['purity'], group['FGA_elbowPt'], xerr=group['purity_std']/np.sqrt(group['sampleNum']),
                     yerr=group['FGA_elbowPt_std']/np.sqrt(group['sampleNum']),fmt='o', label=name, markersize=2, linewidth=0.4)

    # Perform linear regression
    slope, intercept, r_value, pvalue, _ = linregress(data_plot['purity'], data_plot['FGA_elbowPt'])
    # Create the trend line
    x = np.linspace(min(data_plot['purity']), max(data_plot['purity']), 100)
    y = slope * x + intercept
    ax2.plot(x, y, color='k', zorder=3)
    ax2.text(20, 0.25, 'r = %.2f\np = %.3f' % (r_value, pvalue), fontsize=12)

    ax2.set_xlabel("Tumor purity (%)", color="k")
    ax2.set_ylabel("FGA elbow point", color="k")
    ax2.set_xlim([10, 90])
    ax2.set_ylim([0.1, 0.3])
    ax2.set_yticks([0.1, 0.2, 0.3])
    ax2.set_xticks([10, 30, 50, 70, 90])
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)

    plt.savefig(output_fig_fn)  # , dpi=400
    plt.close()

    result_sorted = result.sort_values(by='AS_elbowPt')
    print(result_sorted)
    result_sorted.to_csv('elbow_points.csv', index=True)
    print(np.median(result['AS_elbowPt']), np.mean(result['AS_elbowPt']), np.median(result['FGA_elbowPt']), np.mean(result['FGA_elbowPt']))

    end_time = time.time()
    print('ALL done! Time used: ',end_time-start_time)