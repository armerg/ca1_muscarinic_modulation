import numpy as np
import scipy.stats as sp
import os
import pandas as pd
import h5py

import bokeh.io as bkio
import bokeh.layouts as blay
import bokeh.models as bmod
import bokeh.plotting as bplt

from bokeh.palettes import Category20 as palette
from bokeh.palettes import Category20b as paletteb
import plot_results as plt_res

import frequency_analysis as fan

colrs = palette[20] + paletteb[20] + palette[20] + paletteb[20]


def save_data_to_hdf5(data_folder_path, hdf5_file_path):
    d_paths = [f_file for f_file in os.listdir(data_folder_path) if f_file.endswith('axgt')]

    with pd.HDFStore(hdf5_file_path) as h5file:
        for d_path in d_paths:
            f_path = os.path.join(data_folder_path, d_path)
            d_arr = np.loadtxt(f_path, dtype={'names': ('time', 'Potential', 'Im', 'ACh'),
                                              'formats': ('float', 'float', 'float', 'float')},
                               skiprows=1)



            d_df = pd.DataFrame(d_arr)
            d_name = d_path.split('.')[0].replace(' ', '_').replace('(', '').replace(')', '')
            print(d_name)

            h5file.put('{}/data'.format(d_name), d_df, format='table', data_columns=True)

        print(h5file)


def plot_spike_data(in_h5_file, exclude_list=[]):

    sp_fig = bplt.figure(title='Membrane Potential vs Time')
    sp_fig.xaxis.axis_label = 'time (sec)'
    sp_fig.yaxis.axis_label = 'potential (mV)'
    print('Plotting Potential values from {}'.format(in_h5_file))

    my_lines = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:

        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:
            if 'data' in f_name:

                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    my_lines.append(sp_fig.line(h5_data[f_name]['time'], 1000.0*h5_data[f_name]['Potential'],
                                    line_width=3, color=colrs[f_i])
                                    )
                    legend_items.append((leg_name, [my_lines[-1]]))
                    f_i += 1
    my_legend = bmod.Legend(items=legend_items, location='center')

    sp_fig.add_layout(my_legend, 'right')

    return sp_fig


def plot_selected_ach_data(in_h5_file, select_list):

    sp_fig = bplt.figure(title='Acetylcholine vs Time')
    sp_fig.xaxis.axis_label = 'time (sec)'
    sp_fig.yaxis.axis_label = 'potential (mV)'
    print('Plotting Potential values from {}'.format(in_h5_file))

    my_lines = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:

        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:
            if 'data' in f_name:

                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name in select_list:
                    my_lines.append(sp_fig.line(h5_data[f_name]['time'], h5_data[f_name]['ACh'],
                                    line_width=1, color=colrs[f_i])
                                    )
                    legend_items.append((leg_name, [my_lines[-1]]))
                    f_i += 1
    my_legend = bmod.Legend(items=legend_items, location='center')

    sp_fig.add_layout(my_legend, 'right')

    return sp_fig


def plot_selected_spike_data(in_h5_file, select_list):

    sp_fig = bplt.figure(title='Membrane Potential vs Time')
    sp_fig.xaxis.axis_label = 'time (sec)'
    sp_fig.yaxis.axis_label = 'potential (mV)'
    print('Plotting Potential values from {}'.format(in_h5_file))

    my_lines = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:

        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:
            if 'data' in f_name:

                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name in select_list:
                    my_lines.append(sp_fig.line(h5_data[f_name]['time'], 1000.0*h5_data[f_name]['Potential'],
                                    line_width=1, color=colrs[f_i])
                                    )
                    legend_items.append((leg_name, [my_lines[-1]]))
                    f_i += 1
    my_legend = bmod.Legend(items=legend_items, location='center')

    sp_fig.add_layout(my_legend, 'right')

    return sp_fig


def plot_spike_raster(in_h5_file, exclude_list=[]):
    rast_fig = bplt.figure(title='Spike Raster vs Time')
    rast_fig.xaxis.axis_label = 'time (sec)'
    print('Plotting Spike Raster values from {}'.format(in_h5_file))

    my_circles = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:

        f_i = 1
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:
            if 'spike_times' in f_name:
                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                y_vals = f_i*np.ones(h5_data[f_name].shape)

                if leg_name not in exclude_list:
                    my_circles.append(rast_fig.circle(h5_data[f_name], y_vals,
                                                      line_width=3, color=colrs[f_i-1])
                                      )
                    legend_items.append((leg_name, [my_circles[-1]]))
                    f_i += 1
    my_legend = bmod.Legend(items=legend_items, location='center')

    rast_fig.add_layout(my_legend, 'right')

    return rast_fig


def plot_instantaneous_spike_rate(in_h5_file, exclude_list=[], t_start=0):
    isr_fig = bplt.figure(title='Instantaneous Spike Rate vs Time')
    isr_fig.xaxis.axis_label = 'time (sec)'
    isr_fig.yaxis.axis_label = 'spike rate (Hz)'
    print('Plotting instantaneous spike rate from {}'.format(in_h5_file))

    my_lines = []
    my_circles = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:

        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:
                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    my_lines.append(isr_fig.line(h5_data[f_name]['time'], h5_data[f_name]['ISR'],
                                                line_width=3, color=colrs[f_i])
                                    )
                    my_circles.append(isr_fig.circle(h5_data[f_name]['time'], h5_data[f_name]['ISR'],
                                      size=6, color=colrs[f_i])
                                      )
                    legend_items.append((leg_name, [my_circles[-1], my_lines[-1]]))
                    f_i += 1

    my_legend = bmod.Legend(items=legend_items, location='center')

    isr_fig.add_layout(my_legend, 'right')
    isr_fig.x_range.start = t_start

    return isr_fig


def plot_spike_accel(in_h5_file, exclude_list=[], normalize=False, t_start=0):

    if normalize:
        acc_fig = bplt.figure(title='Normalized Spike Acceleration vs Time')
    else:
        acc_fig = bplt.figure(title='Spike Acceleration vs Time')
    acc_fig.xaxis.axis_label = 'time (sec)'
    acc_fig.yaxis.axis_label = 'spike acceleration (%)'
    print('Plotting spike acceleration from {}'.format(in_h5_file))

    my_lines = []
    my_circles = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:
        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:

                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    if normalize:

                        max_accel = np.max(h5_data[f_name]['Spike_Accel'])
                        my_lines.append(acc_fig.line(h5_data[f_name]['time'], h5_data[f_name]['Spike_Accel']/max_accel,
                                                     line_width=3, color=colrs[f_i])
                                        )
                        my_circles.append(acc_fig.circle(h5_data[f_name]['time'],
                                                         h5_data[f_name]['Spike_Accel']/max_accel,
                                                         size=6, color=colrs[f_i])
                                          )

                    else:
                        my_lines.append(acc_fig.line(h5_data[f_name]['time'], h5_data[f_name]['Spike_Accel'],
                                        line_width=3, color=colrs[f_i])
                                        )
                        my_circles.append(acc_fig.circle(h5_data[f_name]['time'], h5_data[f_name]['Spike_Accel'],
                                          size=6, color=colrs[f_i])
                                          )
                    legend_items.append((leg_name, [my_circles[-1], my_lines[-1]]))
                    f_i += 1
    my_legend = bmod.Legend(items=legend_items, location='center')

    acc_fig.add_layout(my_legend, 'right')
    acc_fig.x_range.start = t_start
    return acc_fig


def plot_spike_accel_aligned(in_h5_file, exclude_list=[], normalize=False):
    if normalize:
        acc_fig = bplt.figure(title='Normalized Spike Acceleration vs Time')
    else:
        acc_fig = bplt.figure(title='Spike Acceleration vs Time')
    acc_fig.xaxis.axis_label = 'time (sec)'
    acc_fig.yaxis.axis_label = 'spike acceleration (%)'
    print('Plotting spike acceleration from {}'.format(in_h5_file))

    my_lines = []
    my_circles = []
    legend_items = []

    with pd.HDFStore(in_h5_file) as h5_data:
        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:
                name = f_name.split('/')[1]
                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    ach_time = h5_data[name + '/ach_times'][0] + 0.5

                    acc_spikes = h5_data[name + '/spike_times'].loc[h5_data[name+'/spike_times'] > ach_time].to_numpy()
                    acc_isr = 1.0 / np.diff(acc_spikes)

                    acc_t = acc_spikes[:-1]
                    sp0 = acc_spikes[0]

                    freq_i = h5_data['frequency_table'].index[h5_data['frequency_table']['Filename'] == name]
                    freq_val = h5_data['frequency_table']['Frequency'][freq_i].values[0]

                    sp_accel = (acc_isr - freq_val)/freq_val*100

                    if normalize:
                        max_accel = np.max(sp_accel)
                        my_lines.append(
                            acc_fig.line(acc_t-sp0, sp_accel / max_accel,
                                         line_width=2, color=colrs[f_i])
                            )
                        my_circles.append(acc_fig.circle(acc_t-sp0, sp_accel / max_accel,
                                                         size=6, color=colrs[f_i])
                                          )

                    else:
                        my_lines.append(acc_fig.line(acc_t-sp0, sp_accel,
                                                     line_width=3, color=colrs[f_i])
                                        )
                        my_circles.append(acc_fig.circle(acc_t-sp0, sp_accel,
                                                         size=6, color=colrs[f_i])
                                          )
                    legend_items.append((leg_name, [my_circles[-1], my_lines[-1]]))
                    f_i += 1
    my_legend = bmod.Legend(items=legend_items, location='center')

    acc_fig.add_layout(my_legend, 'right')

    return acc_fig


def plot_spike_cessation(in_h5_file, exclude_list=[], add_mean=True):

    cess_names = []
    cess_vals = []

    with pd.HDFStore(in_h5_file) as h5_data:

        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:
                name_parts = f_name.split('/')[1].split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    cess_names.append(leg_name)
                    cess_vals.append(1.0/np.min(h5_data[f_name]['ISR']))

    if add_mean:
        mean_cess = np.mean(cess_vals)
        cess_vals.append(mean_cess)
        all_names = cess_names
        mean_name = 'Mean: {0:.2f} sec'.format(mean_cess)
        all_names.append(mean_name)

    else:
        all_names = cess_names
    cess_fig = bplt.figure(x_range=all_names, title='Duration of Spike Cessation after ACh')
    cess_fig.yaxis.axis_label = 'duration (sec)'
    cess_fig.vbar(x=cess_names, top=cess_vals, width=0.9, color=colrs[0])

    if add_mean:
        cess_fig.vbar(x=[mean_name], top=[mean_cess], width=0.9, color='red')

    cess_fig.xaxis.major_label_orientation = np.pi / 2
    cess_fig.y_range.start = 0.0

    return cess_fig


def plot_average_ifr(in_h5_file, exclude_list=[]):

    with pd.HDFStore(in_h5_file) as h5_data:
        h5_df = pd.DataFrame(h5_data['frequency_table'])
        h5_df = h5_df.sort_values(by=['Filename'])

        sel_tab = h5_data['frequency_table'][~h5_data['frequency_table']['Legend'].isin(exclude_list)]

        sel_tab.sort_values('Legend', inplace=True)

        x_names = sel_tab['Legend'].tolist()
        x_names.append('Average')

        cess_fig = bplt.figure(x_range=x_names,
                               title='Average Pre-ACh Frequency and ISR')
        cess_fig.vbar(x=sel_tab['Legend'],
                      top=sel_tab['Frequency'],
                      width=0.9, color='blue', alpha=0.6, legend='Frequency')
        cess_fig.vbar(x=sel_tab['Legend'],
                      top=sel_tab['ISR_Mean'],
                      width=0.6, color='red', alpha=0.6, legend='ISR')

        mean_isr = np.mean(sel_tab['ISR_Mean'])
        mean_freq = np.mean(sel_tab['Frequency'])

        cess_fig.vbar(x=['Average'], top=[mean_freq], width=0.9, color='navy', alpha=0.6)
        cess_fig.vbar(x=['Average'], top=[mean_isr], width=0.6, color='maroon', alpha=0.6)

    cess_fig.xaxis.major_label_orientation = np.pi / 2
    cess_fig.yaxis.axis_label = 'frequency (Hz)'
    cess_fig.y_range.start = 0.0
    cess_fig.legend.location = 'top_right'

    return cess_fig


def plot_average_curve(in_h5_file, time_start=8.5, time_bin_size=0.1, exclude_list=[], spike_acceleration=False,
                       return_curve=False):

    long_time = 0
    with pd.HDFStore(in_h5_file) as h5_data:
        name_sort = list(h5_data.keys())
        name_sort.sort()

        # get longest recorded time
        for f_name in name_sort:
            if 'data' in f_name:
                name = f_name.split('/')[1]
                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    e_time = np.max(h5_data[f_name]['time'])

                    if e_time > long_time:
                        long_time = e_time

        # make array of time bins
        t_bins = np.arange(time_start, long_time+time_bin_size, time_bin_size)
        isr_avg = np.zeros((t_bins.size - 1,))
        acc_avg = np.zeros((t_bins.size - 1,))
        c_count = np.zeros((t_bins.size - 1,))

        for f_name in name_sort:
            if 'spike_times' in f_name:
                name = f_name.split('/')[1]
                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:

                    acc_spikes = h5_data[name + '/spike_times'].loc[h5_data[name + '/spike_times'] > time_start].to_numpy()
                    acc_isrs = 1.0 / np.diff(acc_spikes)

                    acc_t = acc_spikes[:-1]

                    freq_i = h5_data['frequency_table'].index[h5_data['frequency_table']['Filename'] == name]
                    freq_val = h5_data['frequency_table']['Frequency'][freq_i].values[0]

                    sp_accels = (acc_isrs - freq_val) / freq_val * 100

                    sp_is = np.digitize(acc_t, t_bins)

                    for sp_i, sp_acc, sp_isr in zip(sp_is, sp_accels, acc_isrs):
                        isr_avg[sp_i] += sp_isr
                        acc_avg[sp_i] += sp_acc
                        c_count[sp_i] += 1

        isr_avg = np.divide(isr_avg, c_count, where=np.greater(c_count, 0))
        acc_avg = np.divide(acc_avg, c_count, where=np.greater(c_count, 0))

        if spike_acceleration:
            avg_fig = bplt.figure(title='Average Acceleration Versus Time')
            avg_fig.yaxis.axis_label = 'spike acceleration (%)'
            avg_fig.line(t_bins[:-1], acc_avg, line_width=3, color=colrs[0])
            avg_fig.circle(t_bins[:-1], acc_avg, size=12, color=colrs[0])
        else:
            avg_fig = bplt.figure(title='Average Instantaneous Spike Rate Versus Time')
            avg_fig.yaxis.axis_label = 'ISR (Hz)'
            avg_fig.line(t_bins[:-1], isr_avg, line_width=3, color=colrs[0])
            avg_fig.circle(t_bins[:-1], isr_avg, size=12, color=colrs[0])
        avg_fig.xaxis.axis_label = 'time (sec)'

        if return_curve:
            if spike_acceleration:
                return avg_fig, t_bins[:-1], acc_avg
            else:
                return avg_fig, t_bins[:-1], isr_avg
        else:
            return avg_fig


def plot_spike_cessation_vs_isr_variance(in_h5_file, exclude_list=[]):

    cess_names = []

    cess_vals = []
    ifr_vars = []

    with pd.HDFStore(in_h5_file) as h5_data:
        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:
                name = f_name.split('/')[1]
                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    cess_names.append(name)
                    cess_vals.append(1.0 / np.min(h5_data[f_name]['ISR']))

                    c_i = h5_data['frequency_table'].index[h5_data['frequency_table']['Filename'] == name]
                    ifr_vars.append(h5_data['frequency_table']['ISR_Var'][c_i].values[0])

    cess_fig = bplt.figure(title='Spike Cessation vs ISR Variance')
    cess_fig.circle(cess_vals, ifr_vars, size=12, color=colrs[0])
    cess_fig.xaxis.axis_label = 'duration of spike cessation (sec)'
    cess_fig.yaxis.axis_label = 'variance of ISR (Hz)'

    return cess_fig


def plot_peak_acceleration_vs_spike_cessation(in_h5_file, exclude_list=[]):

    fail_acc = []
    fail_cess = []
    fail_names = []

    succ_acc = []
    succ_cess = []
    succ_names = []

    with pd.HDFStore(in_h5_file) as h5_data:
        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:
                name = f_name.split('/')[1]

                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    ach_name = name + '/ach_times'
                    ach_start = h5_data[ach_name][0]

                    cess_val = 1.0 / np.min(h5_data[f_name]['ISR'])
                    acc_i = np.where(h5_data[f_name]['time'] < ach_start)
                    max_acc_pre = np.max(h5_data[f_name].loc[h5_data[f_name]['time'] < ach_start, 'Spike_Accel'].tolist())
                    max_acc = np.max(h5_data[f_name]['Spike_Accel'])
                    if max_acc <= 1.1*max_acc_pre:
                        fail_acc.append(max_acc)
                        fail_cess.append(cess_val)
                        fail_names.append(leg_name)
                    else:
                        succ_acc.append(max_acc)
                        succ_cess.append(cess_val)
                        succ_names.append(leg_name)

    acc_fig = bplt.figure(title='Peak Spike Acceleration vs Duration of Spike Cessation')
    acc_fig.circle(fail_cess, fail_acc, size=12, color='red', legend='no acceleration')
    acc_fig.circle(succ_cess, succ_acc, size=12, color='green', legend='acceleration')
    acc_fig.xaxis.axis_label = 'duration of spike cessation (sec)'
    acc_fig.yaxis.axis_label = 'peak acceleration (%)'

    print('Failed to Demonstrate Spike Acceleration')
    print(fail_names)
    print('Demonstrated at least 10% increase in ISR')
    print(succ_names)

    return acc_fig


def plot_peak_acceleration_vs_isr_variance(in_h5_file, exclude_list=[]):

    acc_vals = []
    var_vals = []
    names = []

    with pd.HDFStore(in_h5_file) as h5_data:
        f_i = 0
        name_sort = list(h5_data.keys())
        name_sort.sort()

        for f_name in name_sort:

            if 'spike_rates' in f_name:
                name = f_name.split('/')[1]

                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                if leg_name not in exclude_list:
                    ach_name = name + '/ach_times'
                    ach_start = h5_data[ach_name][0]

                    c_i = h5_data['frequency_table'].index[h5_data['frequency_table']['Filename'] == name]
                    var_vals.append(h5_data['frequency_table']['ISR_Var'][c_i].values[0])

                    max_acc = np.max(h5_data[f_name]['Spike_Accel'])

                    acc_vals.append(max_acc)
                    names.append(leg_name)

    acc_fig = bplt.figure(title='Peak Spike Acceleration vs ISR Variance')
    acc_fig.circle(var_vals, acc_vals, size=12, color=colrs[0])

    acc_fig.xaxis.axis_label = 'variance of ISR (Hz)'
    acc_fig.yaxis.axis_label = 'peak acceleration (%)'

    return acc_fig


def print_average_table(in_h5_file):

    with pd.HDFStore(in_h5_file) as h5_data:

        h5_df = pd.DataFrame(h5_data['frequency_table'])

        print(h5_df)


def analyze_spike_data_from_hdf5(in_h5_file):

    avg_freqs = []
    avg_isrs = []
    var_isrs = []
    cell_names = []
    legend_names = []

    with pd.HDFStore(in_h5_file) as h5_data:
        for f_i, f_name in enumerate(h5_data.keys()):

            if '/data' in f_name:
                print(f_name)
                name = f_name.split('/')[1]
                sp_name = '{}/spike_times'.format(name)
                isr_name = '{}/spike_rates'.format(name)
                ach_name = '{}/ach_times'.format(name)
                name_parts = name.split('_')
                leg_name = ' '.join(name_parts[:name_parts.index('CA1')])

                # Calculate ACh Times
                ach_i = np.where(h5_data[f_name]['ACh'] > 1e-5)

                ach_times = pd.Series([h5_data[f_name]['time'][ach_i[0][0]], h5_data[f_name]['time'][ach_i[0][-1]]])
                h5_data.put(ach_name, ach_times)

                # Get spike times
                sp_times = fan.get_spike_times(h5_data[f_name]['Potential'], h5_data[f_name]['time'])
                h5_data.put(sp_name, pd.Series(sp_times))

                # Calculate ISR
                isr_vals = np.divide(1.0, np.diff(sp_times))
                isr_ts = (sp_times[:-1] + sp_times[1:]) / 2.0

                # Calculate Average Frequency Before ACh
                pre_spikes = sp_times[np.where(sp_times < ach_times[0])]
                pre_isr = np.divide(1.0, np.diff(pre_spikes))

                avg_isr = np.mean(pre_isr)
                print('Mean of ISR: {}'.format(avg_isr))
                var_isr = np.var(pre_isr)
                print('Variance of ISR: {}'.format(var_isr))
                avg_freq = pre_spikes.size / ach_times[0]
                print('Average Frequency (#/time): {}'.format(avg_freq))
                avg_isrs.append(avg_isr)
                var_isrs.append(var_isr)
                avg_freqs.append(avg_freq)
                cell_names.append(name)
                legend_names.append(leg_name)

                sp_acc = (isr_vals - avg_isr)/avg_isr*100.0
                isr_df = pd.DataFrame(np.vstack((isr_ts, isr_vals, sp_acc)).transpose(),
                                      columns=('time', 'ISR', 'Spike_Accel'))
                h5_data.put(isr_name, isr_df, format='table', data_columns=True)

        freq_dict = {'Filename': cell_names, 'Legend': legend_names, 'ISR_Mean': avg_isrs, 'ISR_Var': var_isrs,
                     'Frequency': avg_freqs}
        h5_data.put('frequency_table', pd.DataFrame(freq_dict), format='table', data_columns=True)


if __name__ == "__main__":

    # Get list of files to convert
    path_parts = os.getcwd().split('/')
    home_i = path_parts.index('ca1_muscarinic_modulation')
    home_path = '/'.join(path_parts[:(home_i+1)])

    mouse_path = os.path.join(home_path,
                              'data_sources/gulledge_lab_recordings/dasari_and_gulledge_2011')
    mouse_h5file = os.path.join(os.path.join(mouse_path, 'dasari_and_gulledge_2011_accel_data.h5'))
    mouse_ex_list = ['Cell 3', 'Cell 11']

    rat_path = os.path.join(home_path,
                            'data_sources/gulledge_lab_recordings/gulledge_and_kawaguchi_2007')
    rat_h5file = os.path.join(os.path.join(rat_path, 'gulledge_and_kawaguchi_2007_accel_data.h5'))
    rat_ex_list = ['Cell 8 Rat', 'Cell 6 trial 4 Rat', 'Cell 6 trial 3 Rat', 'Cell 6 Rat', 'Cell 4 Rat',
                   'Cell 15 Rat', 'Cell 5 trial 2 Rat', 'Cell 5 Rat', 'Cell 14 Rat', 'Cell 13 trial 2 Rat',
                   'Cell 1 trial 2 Rat', 'Cell 13 Rat', 'Cell 4 trial 2 Rat', 'Cell 3 trial 2 Rat',
                   'Cell 25 Rat', 'Cell 2 trial 2 Rat', 'Cell 2 trial 3 Rat', 'Cell 2 trial 4 Rat',
                   'Cell 18 trial 2 Rat', 'Cell 18 trial 3 Rat', 'Cell 22 Rat', 'Cell 19 Rat', 'Cell 15 trial 2 Rat',
                   'Cell 15 trial 4 Rat', 'Cell 15 trial 5 Rat', 'Cell 15 trial 6 Rat', 'Cell 11 Rat', 'Cell 18 Rat',
                   'Cell 6 trial 2 Rat', 'Cell 3 Rat', 'Cell 19 trial 2 Rat', 'Cell 15 trial 3 Rat', 'Cell 16 Rat',
                   'Cell 7 Rat', 'Cell 2 Rat'
                   ]

    h5file = rat_h5file
    ex_list = rat_ex_list

    # save_data_to_hdf5(rat_path, h5file)
    # analyze_spike_data_from_hdf5(h5file)
    # print_average_table(h5file)
    all_figs = []

    all_figs.append(plot_selected_spike_data(h5file, select_list=['Cell 2 Rat', 'Cell 23 Rat']))
    all_figs.append(plot_selected_ach_data(h5file, select_list=['Cell 2 Rat', 'Cell 23 Rat']))
    #all_figs.append(plot_spike_data(h5file, exclude_list=ex_list))
    all_figs.append(plot_spike_raster(h5file, exclude_list=ex_list))
    acc_avg_fig, acc_ts, acc_vals = plot_average_curve(h5file, exclude_list=ex_list,
                                                       spike_acceleration=True, return_curve=True)
    all_figs.append(acc_avg_fig)

    peak_i = np.argmax(acc_vals)
    acc_slope, acc_int, r_val, p_val, acc_err = sp.linregress(acc_ts[peak_i:], acc_vals[peak_i:])
    est_vals = acc_ts*acc_slope + acc_int
    all_figs[-1].line(acc_ts, est_vals, line_width=3, color='green', legend='linear regression')

    lab_source = bmod.ColumnDataSource(dict(x=[11], y=[70],
                                            text=['slope = {0:.2f} (%/sec)\n r = {1:.2f}'.format(acc_slope, r_val)]))
    reg_label = bmod.glyphs.Text(x='x', y='y', text='text', text_color='green')
    all_figs[-1].add_glyph(lab_source, reg_label)

    all_figs.append(plot_spike_accel(h5file, exclude_list=ex_list, t_start=8.5))
    all_figs[-1].line(acc_ts, acc_vals, line_width=4, color='black', legend='Average')
    all_figs[-1].circle(acc_ts, acc_vals, size=16, color='black', legend='Average')

    isr_avg_fig, isr_ts, isr_vals = plot_average_curve(h5file, exclude_list=ex_list,
                                                       spike_acceleration=False, return_curve=True)
    all_figs.append(isr_avg_fig)
    all_figs.append(plot_instantaneous_spike_rate(h5file, exclude_list=ex_list, t_start=8.5))
    all_figs[-1].line(isr_ts, isr_vals, line_width=4, color='black', legend='Average')
    all_figs[-1].circle(isr_ts, isr_vals, size=16, color='black', legend='Average')
    all_figs.append(plot_spike_cessation(h5file, exclude_list=ex_list))
    all_figs.append(plot_average_ifr(h5file, exclude_list=ex_list))

    if all_figs:
        for fig in all_figs:
            plt_res.configure_plot_notebook(fig)

        bkio.show(blay.column(all_figs))
