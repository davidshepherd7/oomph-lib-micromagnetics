#!/usr/bin/env python


import scipy as sp
import matplotlib.pyplot as plt
import sys

def relative_error(exact_value, estimate):
    """
    """
    return sp.absolute(exact_value - estimate) / exact_value



def my_read_data(filename):
    """Read in data and create a named array.
    """

    # Read in as a numpy input_array
    input_array = sp.loadtxt(filename)

    # Give all the fields names
    typelist = [('eps', 'd'),
                ('dt', 'd'),
                ('mconstraintmethod', 'd'),
                ('damping', 'd'),
                ('hk', 'd'),
                ('t', 'd'),
                ('mx', 'd'),
                ('my', 'd'),
                ('mz', 'd'),
                ('r', 'd'),
                ('theta', 'd'),
                ('phi', 'd'),
                ('exact_t', 'd'),
                ('exact_phi', 'd'),
                ('error_t', 'd'),
                ('error_phi', 'd'),]

    input_array.dtype = sp.dtype(typelist)

    return input_array


def plot_relerr(ax, data, mconstraintmethod, marker):
    """
    """

    # Get relative error for this constraint method
    d = data[data['mconstraintmethod'] == mconstraintmethod]
    relerr = relative_error(d['exact_t'], d['t'])

    # Plot
    line = ax.scatter(d['dt'], relerr, label='nothing', marker=marker,
                      c='black', s=100)

    return line


def plot_for_one_alpha_value(ax, data, alpha):
    """
    """

    d = data[data['damping'] == alpha]

    # Set up the plot
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.title('alpha = ' + str(alpha))

    l1 = plot_relerr(ax, d, 0, '+')
    l2 = plot_relerr(ax, d, 1, 'x')
    l3 = plot_relerr(ax, d, 2, 'o')

    return (l1, l2, l3)

# def myplot(ax, data, xname, yname, title):
#     """ Given a data set plot a line for each mconstraint value.
#     """
#     l1 = ax.scatter(data[xname], data[yname], c='black', s=100,
#                     marker='+')

#     l2 = ax.scatter(data[xname], data[yname], c='black', s=100,
#                     marker='x')

#     l3 = ax.scatter(data[xname], data[yname], c='black', s=150,
#                     marker='o')
#     plt.title(title)


def main():
    """
    """
    convergencedata = my_read_data("./data_at_time_1")

    fig = plt.figure()
    for (i,alpha) in enumerate([1, 0.5, 0.1, 0.05, 0.01, 0.005]):
        lines = plot_for_one_alpha_value(fig.add_subplot(3,2,i+1), convergencedata, alpha)

    fig.legend(lines, ('nothing', 'renormalised', 'mid-point'),
               scatterpoints=1, frameon=False)

    fig.text(0.1,0.98,'Relative errors vs exact solution at t=1.')

    fig_filename = "plot.pdf"
    fig.savefig(fig_filename, transparent=True, format='pdf')

    plt.show()

    return 0

if __name__ == "__main__":
    sys.exit(main())
