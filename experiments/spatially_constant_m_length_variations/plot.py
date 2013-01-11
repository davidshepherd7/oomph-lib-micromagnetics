#!/usr/bin/env python

import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D


# My modules
import scatter_matrix

#def plot():

if __name__ == "__main__":
    data = sp.loadtxt("full_data",  skiprows=1)

    print data.shape

    data_name = sp.array(['eps','renormalise','damping','hk','t',
                          'mx','my','mz', 'r','theta',
                          'phi','exact_t', 'exact_phi','error_t','error_phi'])

    # Remove any lines with extremely large numbers, these are ones that
    # failed to converge properly so dt went to infinity.
    temp = data[(data < 1e5).all(1), :]

    # Remove any lines with mz not near -1, there are cases where the
    # simulation didn't finish.
    clean_data = temp[(data[6] > -0.85),:]

    # Calculate relative errors
    relative_t_error = clean_data[:,13] / clean_data[:,4]
    relative_phi_error = clean_data[:,14] / clean_data[:,10]

    # # 3d scatter plot
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(clean_data[:,0], clean_data[:,2], relative_t_error)
    # plt.xlabel('eps')
    # plt.ylabel('damping')

    # matrix of all scatter plots for the selected columns
    selected_cols = [0,1,2,3,8,13,14]
    fig = scatter_matrix.scatterplot(clean_data[:,selected_cols].transpose(),
                                     data_name[selected_cols])
    plt.show()
