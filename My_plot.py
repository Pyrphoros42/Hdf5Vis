import numpy as np
import matplotlib.pyplot as plt
import h5py
import arviz


"""

Read input from hdf5 files

"""

# Used file for testing
# 
# '/home/shikah/Visualisierung/hdf5Files/MLMCMC_10000_500_200_20_1811.hdf5'


# Searching and plotting up to 12 samples data. The data should not be deeper than two file depths removed from name.
# The corresponding attribute should contain 'samples', '_qois' and / or 'qoi_diff'.

def hdf5_plot(name):
    # creating sample plot
    sample_fig, sample_axes = plt.subplots()
    sample_fig.set_size_inches(19, 15)
    sample_fig.suptitle('Samples of %s' % name, fontsize=20)
    # adjusting layout
    plt.rcParams.update({'figure.autolayout': True})
    plt.axis('off')
    sample_fig.tight_layout(pad=3)

    # creating qoi plot
    qoi_fig, qoi_axes = plt.subplots()
    qoi_fig.set_size_inches(19, 15)
    qoi_fig.suptitle('QOIs of %s' % name, fontsize=20)
    # adjusting layout
    plt.rcParams.update({'figure.autolayout': True})
    plt.axis('off')
    qoi_fig.tight_layout(pad=3)

    # creating qoi_diff / mixing plot
    qoi_diff_fig, qoi_diff_axes = plt.subplots()
    qoi_diff_fig.set_size_inches(19, 15)
    qoi_diff_fig.suptitle('QOI_diff of %s' % name, fontsize=20)
    # adjusting layout
    plt.rcParams.update({'figure.autolayout': True})
    plt.axis('off')
    qoi_diff_fig.tight_layout(pad=3)

    # adjusting layout
    plt.rcParams.update({'figure.autolayout': True})

    qoi_fig.tight_layout(pad=3)
    qoi_diff_fig.tight_layout(pad=3)

    with h5py.File(name, 'r') as f:
        # seperate counter since i is not a number and we keep track of three figures
        sample_count = 0
        qoi_count = 0
        qoi_diff_count = 0

        # stored sample values for later scatter / color plot use
        sample_array = []

        for i in f:

            # If we have samples in the title, add the data to sample_fig

            if 'samples' in i:

                # extract samples
                samples = np.array(f[i]['samples'][()])


                x, y = samples.shape

                # if samples is a 1d array
                if y is None:
                    # plotting data with right side orange
                    ax = sample_fig.add_subplot(3, 4, sample_count + 1)
                    n, bins, patches = ax.hist(samples, 30, density=1, facecolor='orange')

                    # adding samples to sample array
                    sample_array.append(samples)

                    # plotting lines
                    ax.grid()
                    # plt.axvline(-3.0, color='red', linestyle='dashed', linewidth=2, label='ows')
                    ax.axvline(samples.mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                    # setting axis title
                    ax.set_ylabel(f'{len(samples)} iterations')

                    sample_count = sample_count + 1

                # if samples has several dimensions which contain data, loop over all dimensions
                else:
                    for j in range(x):
                        # plotting data
                        ax = sample_fig.add_subplot(3, 4, sample_count + 1)
                        n, bins, patches = ax.hist(samples[j], 30, density=1, facecolor='orange')

                        # adding samples to sample array
                        sample_array.append(samples[j])
                        
                        # drawing sample mean in plot
                        ax.axvline(samples[j].mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                        # setting axis and title
                        ax.set_title(f'%s, number %s' % (i, j))
                        ax.set_ylabel(f'{len(samples[j])} iterations')
                        ax.set_xlabel('ESS = ' + str(round(arviz.ess(samples[j]))))

            # Print various data about sample set

                        print('Variance of sample data set ', sample_count + 1, ' = ', np.var(samples[j]))
                        print('Arithmetic mean of sample data set ', sample_count + 1, ' = ', np.mean(samples[j]))
                        print('ESS of sample data set ', sample_count + 1, ' = ', arviz.ess(samples[j]))
                        plt.legend()

                        sample_count = sample_count + 1

            elif "_qois" in i:

                # extract qois
                qois = np.array(f[i]['samples'][()])

                x, y = qois.shape

                # if samples is a 1d array
                if y is None:
                    # plotting data with right side orange
                    ax = qoi_fig.add_subplot(3, 4, qoi_count + 1)
                    n, bins, patches = ax.hist(qois, 30, density=1, facecolor='orange')

                    # plotting lines
                    ax.grid()
                    
                    # plotting qoi mean
                    ax.axvline(qois.mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                    # adjusting axis and title
                    ax.set_xlabel(f'{len(qois)} iterations')

                    qoi_count = qoi_count + 1

                # if samples has several dimensions which contain data, loop over all dimensions
                else:
                    for j in range(x):
                        # plotting data
                        ax = qoi_fig.add_subplot(3, 4, qoi_count + 1)
                        n, bins, patches = ax.hist(qois[j], 30, density=1, facecolor='orange')

                        # plotting qoi mean
                        ax.grid()
                        ax.axvline(qois[j].mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                        # adjusting axis and title
                        ax.set_title(f'%s, number %s' % (i, j))
                        ax.set_ylabel(f'{len(qois[j])} iterations')
                        ax.set_xlabel('ESS = ' + str(round(arviz.ess(qois[j]))))

                        # printing various data
                        
                        print('Variance of qoi data set ', qoi_count + 1, ' = ', np.var(qois[j]))
                        print('Arithmetic mean of qoi data set ', qoi_count + 1, ' = ', np.mean(qois[j]))
                        print('ESS of qoi data set ', qoi_count + 1, ' = ', arviz.ess(qois[j]))
                        plt.legend()

                        qoi_count = qoi_count + 1

            elif 'qoi_diff' in i:
                # extract qoi_diff / mixing
                qoi_diff = np.array(f[i]['samples'][()])

                x, y = qoi_diff.shape

                # if samples is a 1d array
                if y is None:
                    # plotting data
                    ax = qoi_diff_fig.add_subplot(3, 4, qoi_diff_count + 1)
                    n, bins, patches = ax.hist(qoi_diff, 30, density=1, facecolor='orange')

                    # plotting qoi_diff mean
                    ax.grid()
                    ax.axvline(qoi_diff.mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                    # adjusting axis and title
                    ax.set_xlabel(f'{len(qoi_diff)} iterations')


                    qoi_diff_count = qoi_diff_count + 1

                # if samples has several dimensions which contain data, loop over all dimensions
                else:
                    for j in range(x):
                        # plotting data
                        ax = qoi_diff_fig.add_subplot(3, 4, qoi_diff_count + 1)
                        n, bins, patches = ax.hist(qoi_diff[j], 30, density=1, facecolor='orange')

                        # plotting lines
                        ax.grid()
                        # plt.axvline(-3.0, color='red', linestyle='dashed', linewidth=2, label='ows')
                        ax.axvline(qoi_diff[j].mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                        # adjusting axis and title
                        ax.set_title(f'%s, number %s' % (i, j))
                        ax.set_ylabel(f'{len(qoi_diff[j])} iterations')
                        ax.set_xlabel('ESS = ' + str(round(arviz.ess(qoi_diff[j]))))

                        print('Variance of qoi_diff data set ', qoi_diff_count + 1, ' = ', np.var(qoi_diff[j]))
                        print('Arithmetic mean of qoi_diff data set ', qoi_diff_count + 1, ' = ', np.mean(qoi_diff[j]))
                        print('ESS of qoi_diff data set ', qoi_diff_count + 1, ' = ', arviz.ess(qoi_diff[j]))
                        plt.legend()

                        qoi_diff_count = qoi_diff_count + 1



    sample_fig.savefig("new_samples.png")
    qoi_fig.savefig("new_qoi.png")
    qoi_diff_fig.savefig("new_qoi_diff.png")



    # generate scatter plot of samples if more than one sample chain
    if sample_count > 1:
        scatter_fig, scatter_axes = plt.subplots()
        scatter_fig.set_size_inches(19, 15)
        scatter_fig.suptitle('Scatter plot of %s' % name, fontsize=20)
        scatter_fig.tight_layout(pad=3)
        plt.axis('off')

        # We want to scatter plot at most 4 parameters, less if there are not enough chains.
        param_size = min(sample_count, 4)
        
        if sample_count > 4:
            print(f'Warning: {sample_count} parameters given. Scatter plot will only plot the first four.')

        # adjusting range so that row and column have the values 1 to param_size
        # i as row index, j as column index
        for row in range(1, 1 + param_size):
            for column in range(1, 1 + param_size):

                # Plotting just the samples, no scattering
                if row == column:
                    # plotting data
                    ax = scatter_fig.add_subplot(param_size, param_size, (row - 1) * param_size + column)
                    n, bins, patches = ax.hist(sample_array[row], 30, density=1, facecolor='orange')

                    # plotting lines
                    ax.grid()
                    # plt.axvline(-3.0, color='red', linestyle='dashed', linewidth=2, label='ows')
                    ax.axvline(sample_array[row].mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                    # adjusting axis and title
                    ax.set_xlabel(f'{len(sample_array[row])} iterations')

                else:
                    # Scatter plots need data of the same length, finding smaller length.
                    # Also avoiding plotting more than 1000 data points for visibility.
                    max_len = min(len(sample_array[row]), len(sample_array[column]), 1000)

                    # plot the scatter plots
                    ax = scatter_fig.add_subplot(param_size,param_size, (row - 1) * param_size + column)
                    ax.scatter(sample_array[row][:max_len-1], sample_array[column][:max_len-1])
                    ax.set_xlabel(f'{max_len} iterations')

        scatter_fig.savefig("scatter.png")

    # generate 2d Histogram plot of samples if more than one sample chain

    if sample_count > 1:
         color_fig, color_axes = plt.subplots()
         color_fig.set_size_inches(19, 15)
         color_fig.suptitle('2D Histogram plot of %s' % name, fontsize=20)
         color_fig.tight_layout(pad=3)
         plt.axis('off')

         # We want to scatter plot at most 4 parameters, less if there are not enough chains.
         param_size = min(sample_count, 4)

         if sample_count > 4:
             print(f'Warning: {sample_count} parameters given. Color plot will only plot the first four.')

         # adjusting range so that row and column have the values 1 to param_size
         # i as row index, j as column index
         for row in range(1, 1 + param_size):
             for column in range(1, 1 + param_size):

                 # Plotting just the samples, no scattering
                 if row == column:
                     # plotting data
                     ax = color_fig.add_subplot(param_size, param_size, (row - 1) * param_size + column)
                     n, bins, patches = ax.hist(sample_array[row], 30, density=1, facecolor='orange')

                     # plotting lines
                     ax.grid()
                     # plt.axvline(-3.0, color='red', linestyle='dashed', linewidth=2, label='ows')
                     ax.axvline(sample_array[row].mean(), color='blue', linestyle='dashed', linewidth=2, label='$\mu$')

                     # setting axis title
                     ax.set_xlabel(f'{len(sample_array[row])} iterations')

                 else:
                     # Color plots need data of the same length, finding smaller length.
                     # However no reduction to 1000 data points necessairy
                     max_len = min(len(sample_array[row]), len(sample_array[column]))

                     # plot the scatter plots
                     ax = color_fig.add_subplot(param_size,param_size, (row - 1) * param_size + column)
                     ax.hist2d(sample_array[row][:max_len-1], sample_array[column][:max_len-1])
                     ax.set_xlabel(f'{max_len} iterations')


         color_fig.savefig("color.png")
