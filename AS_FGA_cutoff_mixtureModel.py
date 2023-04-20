import time
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import random

if __name__ == "__main__":
    start_time = time.time()
    random.seed(1)

    ############ scatter plot showing AS/FGA change in each patient with cutoff  ############
    print('Reading data ...')
    data = pd.read_csv('MSK_allInfo.csv')
    data = data.iloc[:,list(range(21,121))+[2]]
    data_example = data.loc[data['CANCER_TYPE'] == 'Esophagogastric Cancer',]
    sampleNUM = data_example.shape[0]
    cancerType_list = list(set(data['CANCER_TYPE']))
    print('sampleNUM: ', data.shape[0], cancerType_list)

    print('AS:')
    ################ AS cutoff of individual cancer types
    x_list = np.arange(0.01,0.51,0.01)
    for ct in cancerType_list:
        data_ct = data.loc[data['CANCER_TYPE']==ct,]
        sampleNUM_ct = data_ct.shape[0]
        kn_bootstrap_ct = []
        for i in range(1000):
            samples_temp = random.choices(range(sampleNUM_ct), k=sampleNUM_ct)
            y_temp = np.array(data_ct.iloc[samples_temp, 0:50].mean(axis=0))
            gmm = GaussianMixture(n_components=2, random_state=0).fit(y_temp.reshape(-1, 1))
            # Predict the probabilities of belonging to each component
            probs = gmm.predict_proba(y_temp.reshape(-1, 1))
            # Find the cutoff point as the maximum probability of belonging to one component
            cutoff = x_list[np.argmax(probs[:, 0])]
            kn_bootstrap_ct.append(cutoff)
        kn_mean = np.mean(kn_bootstrap_ct)  # np.quantile(kn_bootstrap, 0.5)
        kn_std = np.std(kn_bootstrap_ct)
        kn_025 = np.quantile(kn_bootstrap_ct, 0.025)
        kn_975 = np.quantile(kn_bootstrap_ct, 0.975)
        print('%50s: %.2f (%.2f, %.2f)' % (ct, kn_mean, kn_025, kn_975))

    print('FGA:')
    ################ FGA cutoff of individual cancer types
    for ct in cancerType_list:
        data_ct = data.loc[data['CANCER_TYPE']==ct,]
        sampleNUM_ct = data_ct.shape[0]
        kn_bootstrap_ct = []
        for i in range(1000):
            samples_temp = random.choices(range(sampleNUM_ct), k=sampleNUM_ct)
            y_temp = np.array(data_ct.iloc[samples_temp, 50:100].mean(axis=0))
            gmm = GaussianMixture(n_components=2, random_state=0).fit(y_temp.reshape(-1, 1))
            # Predict the probabilities of belonging to each component
            probs = gmm.predict_proba(y_temp.reshape(-1, 1))
            # Find the cutoff point as the maximum probability of belonging to one component
            cutoff = x_list[np.argmax(probs[:, 0])]
            kn_bootstrap_ct.append(cutoff)
        kn_mean = np.mean(kn_bootstrap_ct)  # np.quantile(kn_bootstrap, 0.5)
        kn_std = np.std(kn_bootstrap_ct)
        kn_025 = np.quantile(kn_bootstrap_ct, 0.025)
        kn_975 = np.quantile(kn_bootstrap_ct, 0.975)
        print('%50s: %.2f (%.2f, %.2f)' % (ct, kn_mean, kn_025, kn_975))

    end_time = time.time()
    print('ALL done! Time used: ',end_time-start_time)