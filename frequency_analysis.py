import h5py
import numpy as np
import os


def analyze_files_in_dir(in_file_dir, save_to_new=False, force_reanalysis=False, thresh_dict={}, time_int=[]):
    """

    :param in_file_dir:
    :param save_to_new:
    :param force_reanalysis:
    :param thresh_dict:
    :param time_int:
    :return:
    """
    file_list = os.listdir(in_file_dir)
    file_list.sort()

    ana_list = []

    for fil in file_list:

        if fil.endswith('.hdf5') and 'spikes' not in fil:

            if force_reanalysis:
                ana_list.append(fil)
            elif fil.replace('.hdf5', '_spikes.hdf5') not in file_list:
                ana_list.append(fil)

    ana_paths = [os.path.join(in_file_dir, ana_file) for ana_file in ana_list]

    for ana_file in ana_paths:
        print('Analyzing {0}'.format(ana_file))
        analyze_spikes_from_file(ana_file, save_to_new=save_to_new, thresh_dict=thresh_dict,
                                 force_reanalysis=force_reanalysis, time_int=time_int)


def analyze_spikes_from_file(in_file_path, force_reanalysis=False, thresh_dict={}, time_int=[]):
    """

    :param in_file_path:
    :param save_to_new: boolean, if True spike times will be saved to a separate file, if False spike times will be
    saved within the given file as new
    :param force_reanalysis: boolean, sets whether to force a reanalysis of a file even if spike times are already
    present within the given hdf5 file
    :param thresh_dict:
    :param time_int:
    :return: new_file_path
    """
    #
    # if save_to_new:
    #     out_file_path = in_file_path.replace('.hdf5', '_spikes.hdf5')
    #
    # else:
    #     out_file_path = in_file_path

    out_file_path = in_file_path

    res_grps = []
    sp_grps = []
    isi_grps = []
    freq_grps = []

    spike_dict = {}
    isi_dict = {}
    m_isi_dict = {}
    freq_dict = {}
    m_freq_dict = {}
    time_dict = {}

    m_cell_dict = {}

    with h5py.File(in_file_path, 'r') as in_file:

        for group_name in in_file.keys():
            if 'result' in group_name:
                res_grps.append(group_name)
            elif 'spikes' in group_name:
                sp_grps.append(group_name)
            elif 'isi' in group_name:
                isi_grps.append(group_name)
            elif 'freq' in group_name:
                freq_grps.append(group_name)

        for sp_grp in sp_grps:
            print('Analyzing spike times for group: {0}'.format(sp_grp))
            cur_group = in_file[sp_grp]
            r_grp = sp_grp.replace('spikes', 'results')
            isi_dict[sp_grp] = {}
            m_isi_dict[sp_grp] = {}
            freq_dict[sp_grp] = {}
            m_freq_dict[sp_grp] = {}
            m_cell_dict[sp_grp] = {}

            for dataset in in_file[r_grp].keys():

                if 'time' in dataset:
                    if time_int:
                        sim_dur = float(time_int[1] - time_int[0]) / 1000  # duration of simulation in seconds
                    else:
                        sim_dur = (in_file[r_grp][dataset][-1] - in_file[r_grp][dataset][0]) / 1000

            for dataset in cur_group.keys():
                sp_data = np.array(cur_group[dataset])

                if 'rand_int' in dataset or 'set_times' in dataset or 'set_freq' in dataset:
                    if 'rand_int' in dataset:
                        key_base = dataset.replace('_rand_int', '')
                    elif 'set_freq' in dataset:
                        key_base = dataset.replace('_set_freq', '')
                    else:
                        key_base = dataset.replace('_set_times', '')
                else:
                    key_base = dataset

                cell_name = key_base.replace('_spikes', '')
                cell_base = cell_name.rstrip('0123456789_')

                isi_key = key_base.replace('spikes', 'isi')
                freq_key = key_base.replace('spikes', 'freq')
                cfreq_key = key_base.replace('spikes', 'cfreq')
                m_isi_key = key_base.replace('spikes', 'mean_isi')
                m_freq_key = key_base.replace('spikes', 'mean_freq')
                s_freq_key = key_base.replace('spikes', 'std_freq')

                if time_int:
                    int_spike = np.intersect1d(np.where(sp_data > time_int[0]),
                                               np.where(sp_data < time_int[1]))
                    ana_spikes = cur_group[dataset][int_spike]
                else:
                    ana_spikes = sp_data

                n_spikes = ana_spikes.size
                freq_dict[sp_grp][cfreq_key] = float(n_spikes) / sim_dur

                if len(ana_spikes) > 1:
                    isi_dict[sp_grp][isi_key] = np.diff(ana_spikes)
                    freq_dict[sp_grp][freq_key] = np.array(
                        [1000 / isi for isi in isi_dict[sp_grp][isi_key] if isi != 0])

                    if cell_base not in m_cell_dict[sp_grp]:
                        m_cell_dict[sp_grp][cell_base] = [cell_name]
                    else:
                        m_cell_dict[sp_grp][cell_base].append(cell_name)

                    m_isi_dict[sp_grp][m_isi_key] = np.mean(isi_dict[sp_grp][isi_key])
                    m_freq_dict[sp_grp][m_freq_key] = np.mean(freq_dict[sp_grp][freq_key])
                    m_freq_dict[sp_grp][s_freq_key] = np.std(freq_dict[sp_grp][freq_key])
                else:
                    m_isi_dict[sp_grp][m_isi_key] = sim_dur
                    m_freq_dict[sp_grp][m_freq_key] = 0
                    m_freq_dict[sp_grp][s_freq_key] = 0

    for sp_group, cell_dict in m_cell_dict.iteritems():

        for cell_type, cell_list in cell_dict.iteritems():
            m_isi = 0
            m_freq = 0
            m_cfreq = 0

            cfreq_list = []

            freq_list = []

            for cell_name in cell_list:
                freq_key = cell_name + '_freq'
                m_isi_key = cell_name + '_mean_isi'
                m_freq_key = cell_name + '_mean_freq'
                cfreq_key = cell_name + '_cfreq'

                m_isi += m_isi_dict[sp_group][m_isi_key]
                m_freq += m_freq_dict[sp_group][m_freq_key]
                m_cfreq += freq_dict[sp_group][cfreq_key]

                cfreq_list.append(freq_dict[sp_group][cfreq_key])
                freq_list += list(freq_dict[sp_group][freq_key])

            isi_tot_k = cell_type + '_total_mean_isi'
            freq_tot_k = cell_type + '_total_mean_freq'
            freq_std_k = cell_type + '_total_std_freq'
            cfreq_std_k = cell_type + '_std_cfreq'
            cfreq_tot_k = cell_type + '_total_cfreq'

            m_isi_dict[sp_group][isi_tot_k] = np.array(m_isi/len(cell_list))
            m_freq_dict[sp_group][freq_tot_k] = np.array(m_freq/len(cell_list))
            m_freq_dict[sp_group][cfreq_tot_k] = np.array(m_cfreq/len(cell_list))
            m_freq_dict[sp_group][cfreq_std_k] = np.array([np.std(cfreq_list)])
            m_freq_dict[sp_group][freq_std_k] = np.array([np.std(freq_list)])

    print('Spike analysis complete. Saving to {0}'.format(out_file_path))

    with h5py.File(out_file_path, 'a') as s_file:

        for grp_name, grp_dict in freq_dict.iteritems():

            freq_name = grp_name.replace('spikes', 'frequencies')

            if freq_name in s_file:
                del s_file[freq_name]

            freq_grp = s_file.create_group(freq_name)

            for src_name, freq_arr in grp_dict.iteritems():
                freq_grp.create_dataset(src_name, data=freq_arr)

        for grp_name, grp_dict in isi_dict.iteritems():

            isi_name = grp_name.replace('spikes', 'intervals')

            if isi_name in s_file:
                del s_file[isi_name]

            isi_grp = s_file.create_group(isi_name)

            for src_name, isi_arr in grp_dict.iteritems():
                isi_grp.create_dataset(src_name, data=isi_arr)

        for grp_name, grp_dict in m_isi_dict.iteritems():

            isi_name = grp_name.replace('spikes', 'mean_intervals')

            if isi_name in s_file:
                del s_file[isi_name]

            isi_grp = s_file.create_group(isi_name)

            for src_name, isi_arr in grp_dict.iteritems():
                isi_grp.create_dataset(src_name, data=isi_arr)

        for grp_name, grp_dict in m_freq_dict.iteritems():

            freq_name = grp_name.replace('spikes', 'mean_frequencies')

            if freq_name in s_file:
                del s_file[freq_name]

            freq_grp = s_file.create_group(freq_name)

            for src_name, freq_arr in grp_dict.iteritems():
                freq_grp.create_dataset(src_name, data=freq_arr)

    print('Saving to file complete.')

    return out_file_path


def analyze_spikes_from_file_old(in_file_path, force_reanalysis=False, save_to_new=False, thresh_dict={}, time_int=[]):
    """
    analyzes files created by network_neuron_parallel.py before 2017-08-23
    :param in_file_path:
    :param save_to_new: boolean, if True spike times will be saved to a separate file, if False spike times will be
    saved within the given file as new
    :param force_reanalysis: boolean, sets whether to force a reanalysis of a file even if spike times are already
    present within the given hdf5 file
    :param thresh_dict:
    :param time_int:
    :return: new_file_path
    """

    if save_to_new:
        out_file_path = in_file_path.replace('.hdf5', '_spikes.hdf5')

    else:
        out_file_path = in_file_path

    res_grps = []
    sp_grps = []
    isi_grps = []
    freq_grps = []

    spike_dict = {}
    isi_dict = {}
    m_isi_dict = {}
    freq_dict = {}
    m_freq_dict = {}
    time_dict = {}

    m_cell_dict = {}

    with h5py.File(in_file_path, 'r') as in_file:

        for group_name in in_file.keys():
            if 'result' in group_name:
                res_grps.append(group_name)
            elif 'spikes' in group_name:
                sp_grps.append(group_name)
            elif 'isi' in group_name:
                isi_grps.append(group_name)
            elif 'freq' in group_name:
                freq_grps.append(group_name)

        if len(res_grps) != len(sp_grps) or force_reanalysis:
            # Delete spike data groups to remove chance of conflict
            # for del_grp in sp_grps:
            #     del in_file[del_grp]
            # for del_grp in isi_grps:
            #     del in_file[del_grp]
            # for del_grp in freq_grps:
            #     del in_file[del_grp]

            for res_grp in res_grps:  # for each simulation run using a different seed
                print('Analyzing spike times for group: {0}'.format(res_grp))
                cur_group = in_file[res_grp]
                spike_dict[res_grp] = {}
                time_dict[res_grp] = {}
                isi_dict[res_grp] = {}
                m_isi_dict[res_grp] = {}
                freq_dict[res_grp] = {}
                m_freq_dict[res_grp] = {}
                m_cell_dict[res_grp] = {}

                time_names = []
                time_lengths = []

                # get time arrays and note the lengths of the arrays in case the lengths differ due to
                # differences between node calculations

                for dataset in cur_group.keys():
                    if 'time' in dataset:
                        time_dict[res_grp][dataset] = np.array(cur_group[dataset])
                        time_names.append(dataset)
                        time_lengths.append(time_dict[res_grp][dataset].size)

                        if time_int:
                            sim_dur = float(time_int[1] - time_int[0])/1000  # duration of simulation in seconds
                        else:
                            sim_dur = (time_dict[res_grp][dataset][-1] - time_dict[res_grp][dataset][0])/1000

                for dataset in cur_group.keys():

                    # print('Analyzing {0}'.format(dataset))

                    if 'soma_v' in dataset:
                        cur_data = cur_group[dataset]
                        d_size = cur_data.size
                        key_base = dataset.replace('soma_v_', '')
                        cell_base = key_base.rstrip('0123456789_')
                        spike_key = key_base + '_spikes'

                        if 'CA1P' in dataset:
                            t = 0

                        isi_key = key_base + '_isi'
                        freq_key = key_base + '_freq'
                        cfreq_key = key_base + '_cfreq'

                        # finds the appropriate time array for this data array
                        d_time = time_names[time_lengths.index(d_size)]

                        if not thresh_dict:
                            spike_dict[res_grp][spike_key] = get_spike_times(np.array(cur_data), time_dict[res_grp][d_time])

                        else:
                            id_base = dataset.replace('soma_v_', '').rstrip('0123456789')

                            if id_base in thresh_dict:
                                spike_dict[res_grp][spike_key] = get_spike_times(np.array(cur_data), time_dict[res_grp][d_time],
                                                                                 threshold=thresh_dict[id_base])
                            else:
                                spike_dict[res_grp][spike_key] = get_spike_times(np.array(cur_data), time_dict[res_grp][d_time])

                        if time_int:
                            int_spike = np.intersect1d(np.where(spike_dict[res_grp][spike_key] > time_int[0]),
                                                       np.where(spike_dict[res_grp][spike_key] < time_int[1]))
                            ana_spikes = spike_dict[res_grp][spike_key][int_spike]
                        else:
                            ana_spikes = spike_dict[res_grp][spike_key]

                        n_spikes = ana_spikes.size
                        freq_dict[res_grp][cfreq_key] = n_spikes / sim_dur

                        m_isi_key = key_base + '_mean_isi'
                        m_freq_key = key_base + '_mean_freq'
                        s_freq_key = key_base + '_std_freq'

                        if len(spike_dict[res_grp][spike_key]) > 1:
                            isi_dict[res_grp][isi_key] = np.diff(ana_spikes)
                            freq_dict[res_grp][freq_key] = np.array(
                                [1000/isi for isi in isi_dict[res_grp][isi_key] if isi != 0])

                            if cell_base not in m_cell_dict[res_grp]:
                                m_cell_dict[res_grp][cell_base] = [key_base]
                            else:
                                m_cell_dict[res_grp][cell_base].append(key_base)

                            m_isi_dict[res_grp][m_isi_key] = np.mean(isi_dict[res_grp][isi_key])
                            m_freq_dict[res_grp][m_freq_key] = np.mean(freq_dict[res_grp][freq_key])
                            m_freq_dict[res_grp][s_freq_key] = np.std(freq_dict[res_grp][freq_key])
                        else:
                            m_isi_dict[res_grp][m_isi_key] = sim_dur
                            m_freq_dict[res_grp][m_freq_key] = 0
                            m_freq_dict[res_grp][s_freq_key] = 0

                    elif 'rand_int' in dataset or 'set_times' in dataset or 'set_freq' in dataset:
                        if 'rand_int' in dataset:
                            key_base = dataset.replace('_rand_int', '')
                        elif 'set_freq' in dataset:
                            key_base = dataset.replace('_set_freq', '')
                        else:
                            key_base = dataset.replace('_set_times', '')

                        cell_base = key_base.rstrip('0123456789')

                        spike_key = dataset

                        if save_to_new:
                            spike_dict[res_grp][dataset] = np.array(cur_group[dataset])

                        isi_key = key_base + '_isi'
                        freq_key = key_base + '_freq'
                        cfreq_key = key_base + '_cfreq'

                        if time_int:
                            #  look
                            int_spike = np.intersect1d(np.where(spike_dict[res_grp][spike_key] > time_int[0]),
                                                       np.where(spike_dict[res_grp][spike_key] < time_int[1]))
                            # spikes that occur within requested time interval
                            ana_spikes = spike_dict[res_grp][spike_key][int_spike]
                        else:
                            # spikes to be analyzed
                            ana_spikes = spike_dict[res_grp][spike_key]

                        n_spikes = ana_spikes.size
                        freq_dict[res_grp][cfreq_key] = n_spikes/sim_dur

                        if len(spike_dict[res_grp][spike_key]) > 1:
                            isi_dict[res_grp][isi_key] = np.diff(ana_spikes)
                            freq_dict[res_grp][freq_key] = np.array([1000/isi for isi in isi_dict[res_grp][isi_key]
                                                                     if isi != 0])

                            m_isi_key = key_base + '_mean_isi'
                            m_freq_key = key_base + '_mean_freq'

                            if cell_base not in m_cell_dict[res_grp]:
                                m_cell_dict[res_grp][cell_base] = [key_base]
                            else:
                                m_cell_dict[res_grp][cell_base].append(key_base)

                            m_isi_dict[res_grp][m_isi_key] = np.mean(isi_dict[res_grp][isi_key])
                            m_freq_dict[res_grp][m_freq_key] = np.mean(freq_dict[res_grp][freq_key])

        else:

            print('File {0} has already been analyzed for spike times'.format(in_file_path))
            exit(0)

    for res_grp_n, cell_dict in m_cell_dict.iteritems():

        for cell_type, cell_list in cell_dict.iteritems():
            m_isi = 0
            m_freq = 0
            m_cfreq = 0

            cfreq_list = []

            freq_list = []

            for cell_name in cell_list:
                freq_key = cell_name + '_freq'
                m_isi_key = cell_name + '_mean_isi'
                m_freq_key = cell_name + '_mean_freq'
                cfreq_key = cell_name + '_cfreq'

                m_isi += m_isi_dict[res_grp][m_isi_key]
                m_freq += m_freq_dict[res_grp][m_freq_key]
                m_cfreq += freq_dict[res_grp][cfreq_key]

                cfreq_list.append(freq_dict[res_grp][cfreq_key])
                freq_list += list(freq_dict[res_grp][freq_key])

            isi_tot_k = cell_type + '_total_mean_isi'
            freq_tot_k = cell_type + '_total_mean_freq'
            freq_std_k = cell_type + '_total_std_freq'
            cfreq_std_k = cell_type + '_std_cfreq'
            cfreq_tot_k = cell_type + '_total_cfreq'

            m_isi_dict[res_grp][isi_tot_k] = np.array(m_isi/len(cell_list))
            m_freq_dict[res_grp][freq_tot_k] = np.array(m_freq/len(cell_list))
            m_freq_dict[res_grp][cfreq_tot_k] = np.array(m_cfreq/len(cell_list))
            m_freq_dict[res_grp][cfreq_std_k] = np.array([np.std(cfreq_list)])
            m_freq_dict[res_grp][freq_std_k] = np.array([np.std(freq_list)])

    print('Spike analysis complete. Saving to {0}'.format(out_file_path))

    with h5py.File(out_file_path, 'w') as s_file:

        for grp_name, grp_dict in spike_dict.iteritems():

            if grp_name == 'results':
                sp_name = 'spikes'
            elif 'results_' in grp_name:
                sp_name = grp_name.replace('results_', 'spikes_')

            sp_grp = s_file.create_group(sp_name)

            for src_name, sp_arr in grp_dict.iteritems():
                sp_grp.create_dataset(src_name, data=sp_arr)

        for grp_name, grp_dict in freq_dict.iteritems():

            if grp_name == 'results':
                freq_name = 'frequencies'
            elif 'results_' in grp_name:
                freq_name = grp_name.replace('results_', 'frequencies_')

            freq_grp = s_file.create_group(freq_name)

            for src_name, freq_arr in grp_dict.iteritems():
                freq_grp.create_dataset(src_name, data=freq_arr)

        for grp_name, grp_dict in isi_dict.iteritems():

            if grp_name == 'results':
                isi_name = 'intervals'
            elif 'results_' in grp_name:
                isi_name = grp_name.replace('results_', 'intervals_')

            isi_grp = s_file.create_group(isi_name)

            for src_name, isi_arr in grp_dict.iteritems():
                isi_grp.create_dataset(src_name, data=isi_arr)

        for grp_name, grp_dict in m_isi_dict.iteritems():

            if grp_name == 'results':
                isi_name = 'mean_intervals'
            elif 'results_' in grp_name:
                isi_name = grp_name.replace('results_', 'mean_intervals_')

            isi_grp = s_file.create_group(isi_name)

            for src_name, isi_arr in grp_dict.iteritems():
                isi_grp.create_dataset(src_name, data=isi_arr)

        for grp_name, grp_dict in m_freq_dict.iteritems():

            if grp_name == 'results':
                freq_name = 'mean_frequencies'
            elif 'results_' in grp_name:
                freq_name = grp_name.replace('results_', 'mean_frequencies_')

            freq_grp = s_file.create_group(freq_name)

            for src_name, freq_arr in grp_dict.iteritems():
                freq_grp.create_dataset(src_name, data=freq_arr)

    print('Saving to file complete.')

    return out_file_path


def analyze_spikes_from_file_parallel(in_file_path, save_to_new=True, force_reanalysis=False, thresh_dict={}, time_int=[]):
    """

    :param in_file_path:
    :param save_to_new: boolean, if True spike times will be saved to a separate file, if False spike times will be
    saved within the given file as new
    :param force_reanalysis: boolean, sets whether to force a reanalysis of a file even if spike times are already
    present within the given hdf5 file
    :param thresh_dict:
    :param time_int:
    :return:
    """

    num_procs = comm.size()
    my_rank = comm.rank()

    if save_to_new:
        out_file_path = in_file_path.replace('.hdf5', '_spikes.hdf5')

    else:
        out_file_path = in_file_path

    res_grps = []
    sp_grps = []
    isi_grps = []
    freq_grps = []

    spike_dict = {}
    isi_dict = {}
    freq_dict = {}
    time_dict = {}

    with h5py.File(in_file_path, 'r', driver='mpio', comm=MPI.COMM_WORLD) as in_file:



        for group_name in in_file.keys():
            if 'result' in group_name:
                res_grps.append(group_name)
            elif 'spikes' in group_name:
                sp_grps.append(group_name)
            elif 'isi' in group_name:
                isi_grps.append(group_name)
            elif 'freq' in group_name:
                freq_grps.append(group_name)

        if len(res_grps) != len(sp_grps) or force_reanalysis or save_to_new:
            # Delete spike data groups to remove chance of conflict
            # for del_grp in sp_grps:
            #     del in_file[del_grp]
            # for del_grp in isi_grps:
            #     del in_file[del_grp]
            # for del_grp in freq_grps:
            #     del in_file[del_grp]

            for res_grp in res_grps:  # for each simulation run using a different seed
                print('Analyzing spike times for group: {0}'.format(res_grp))
                cur_group = in_file[res_grp]
                spike_dict[res_grp] = {}
                time_dict[res_grp] = {}
                isi_dict[res_grp] = {}
                freq_dict[res_grp] = {}

                time_names = []
                time_lengths = []

                # get time arrays and note the lengths of the arrays in case the lengths differ due to
                # differences between node calculations

                all_dsets = cur_group.keys()
                all_dsets.sort()

                my_dsets = [all_dsets[i] for i in range(my_rank, len(all_dsets), num_procs)]

                for dataset in my_dsets:

                    if 'time' in dataset:
                        time_dict[res_grp][dataset] = np.array(cur_group[dataset])
                        time_names.append(dataset)
                        time_lengths.append(time_dict[res_grp][dataset].size)

                for dataset in my_dsets:

                    print('Analyzing {0}'.format(dataset))

                    if 'soma_v' in dataset:
                        cur_data = cur_group[dataset]
                        d_size = cur_data.size
                        key_base = dataset.replace('soma_v_', '')
                        spike_key = key_base + '_spikes'

                        isi_key = key_base + '_isi'
                        freq_key = key_base + '_freq'

                        # finds the appropriate time array for this data array
                        d_time = time_names[time_lengths.index(d_size)]

                        if not thresh_dict:
                            spike_vals = get_spike_times(np.array(cur_data), time_dict[res_grp][d_time])
                            if spike_vals:
                                spike_dict[res_grp][spike_key] = spike_vals
                        else:
                            id_base = dataset.replace('soma_v_', '').rstrip('0123456789')

                            if id_base in thresh_dict:
                                spike_vals = get_spike_times(cur_data, time_dict[res_grp][d_time], threshold=thresh_dict[id_base])
                                if spike_vals:
                                    spike_dict[res_grp][spike_key] = spike_vals
                            else:
                                spike_vals = get_spike_times(cur_data, time_dict[res_grp][d_time])
                                if spike_vals:
                                    spike_dict[res_grp][spike_key] = spike_vals

                        if len(spike_dict[res_grp][spike_key]) > 1:
                            isi_dict[res_grp][isi_key] = np.diff(spike_dict[res_grp][spike_key])
                            freq_dict[res_grp][freq_key] = np.array(
                                [1000/isi for isi in isi_dict[res_grp][isi_key] if isi != 0])

                    elif 'rand_int' in dataset or 'set_times' in dataset or 'set_freq' in dataset:
                        if 'rand_int' in dataset:
                            key_base = dataset.replace('rand_int', '')
                        elif 'set_freq' in dataset:
                            key_base = dataset.replace('set_freq', '')
                        else:
                            key_base = dataset.replace('set_times', '')

                        spike_key = dataset

                        if save_to_new:
                            spike_dict[res_grp][dataset] = np.array(cur_group[dataset])

                        isi_key = key_base + 'isi'
                        freq_key = key_base + 'freq'

                        if len(spike_dict[res_grp][spike_key]) > 1:
                            isi_dict[res_grp][isi_key] = np.diff(spike_dict[res_grp][spike_key])
                            freq_dict[res_grp][freq_key] = np.array([1000/isi for isi in isi_dict[res_grp][isi_key]
                                                                     if isi != 0])

        else:

            print('File {0} has already been analyzed for spike times'.format(in_file_path))
            exit(0)

    print('Spike analysis complete. Saving to {0}'.format(out_file_path))

    with h5py.File(out_file_path, 'w') as s_file:

        for grp_name, grp_dict in spike_dict.iteritems():

            sp_name = grp_name.replace('results_', 'spikes_')

            sp_grp = s_file.create_group(sp_name)

            for src_name, sp_arr in grp_dict.iteritems():
                sp_grp.create_dataset(src_name, data=sp_arr)

        for grp_name, grp_dict in freq_dict.iteritems():
            freq_name = grp_name.replace('results_', 'frequencies_')

            freq_grp = s_file.create_group(freq_name)

            for src_name, freq_arr in grp_dict.iteritems():
                freq_grp.create_dataset(src_name, data=freq_arr)

        for grp_name, grp_dict in isi_dict.iteritems():
            isi_name = grp_name.replace('results_', 'intervals_')

            isi_grp = s_file.create_group(isi_name)

            for src_name, isi_arr in grp_dict.iteritems():
                isi_grp.create_dataset(src_name, data=isi_arr)

    print('Saving to file complete.')


def calc_freq_from_interval(interval, unit='msec'):
    """
    given an interval in the provided unit, returns the frequency
    :param interval:
    :param unit: unit of time interval given, valid values include ('msec', 'sec')
    :return: frequency in hertz
    """
    if interval <= 0:
        print('Invalid interval: {0} provided for calculating frequency'.format(interval))
        return -1

    if unit == 'msec':
        return 1000.0/interval
    elif unit == 'sec':
        return 1.0/interval
    else:
        print('Invalid unit: {0} provided for calculating frequency of interval'.format(unit))
        return -1


def calc_mean_frequency(spike_train, unit='second', time_limits=[]):
    """

    :param spike_train:
    :param unit: ['second', 'msec']
    :param time_limits: limits of time to perform analysis on
    :return:
    """
    if time_limits:
        first_time = time_limits[0]
        last_time = time_limits[1]
    else:
        first_time = spike_train[0]
        last_time = spike_train[-1]

    dt = last_time-first_time

    if unit =='second':
        multi = 1
    elif unit == 'msec':
        multi = 1000

    return len(spike_train)/dt*multi


def calc_mean_isi(spike_train, unit='second', time_limits=[]):
    """

    :param spike_train:
    :param unit: ['second', 'msec']
    :param time_limits: limits of time to perform analysis on
    :return:
    """
    train_intvs = np.diff(spike_train)

    return np.mean(train_intvs)


def calc_number_of_spikes(in_array, threshold=0.0):
    """

    :param in_array: array of values
    :param threshold: value that if in_array passes, counts as a spike
    :return: num_spikes: integer value of the number of spikes
    """
    num_spikes = 0

    for in_ind, in_val in enumerate(in_array[0:-2]):
        if in_val < threshold < in_array[in_ind+1]:
            num_spikes += 1

    return num_spikes


def gauss_func(x, mu=0, sigma=1.0):
    """
    Simple implementation of gaussian function
    :param mu:
    :param sigma:
    :return: returns value of function given input parameters
    """

    my_exp = -0.5*((x-mu)/sigma)*((x-mu)/sigma)
    my_mult = 1/(sigma*np.sqrt(2*np.pi))
    return my_mult*np.exp(my_exp)


def get_all_spike_times(spike_dict):
    """

    :param spike_dict:
    :return: spike_list
    """

    spike_list = []

    for k, vals in spike_dict.iteritems():
        spike_list += vals

    return spike_list


def get_spike_times(in_array, time_array, threshold=0.0):
    """

    :param in_array:
    :param time_array:
    :param threshold:
    :return:
    """

    spike_list = []

    high_inds = np.where(in_array > threshold)

    if high_inds:
        for in_ind in high_inds[0]:
            if in_ind > 0 and in_array[(in_ind-1)] < threshold:
                spike_list.append(time_array[in_ind])

    return np.array(spike_list)


def measure_spike_trains_correlation_coefficient(spike_train_1, spike_train_2, std_val=5.0, sample_period=0.01):
    """

    :param spike_train_1: numpy array holding the spike times. times are assumed to be in msec
    :param spike_train_2: numpy array holding the spike times. times are assumed to be in msec
    :param std_val: standard deviation used to create the gaussian function which will be convolved with spike trains
    :param sample_period: period between samples. This is used for creating a time series for each spike train. Value is
        assumed to be in msec
    :return: coe: coefficient of correlation for the given spike trains calculated with the given parameters
             arr_dict: dictionary of arrays that could be used for plotting the time series created to calculate coe
                'spike1': array holding spike time series 1,
                'conv1': array holding spike time series 1 convolved with gaussian function,
                'spike2': array holding spike time series 2,
                'conv2': array holding spike time series 2 convolved with gaussian function,
                't': array holding time steps
    """

    first_sp = min((np.min(spike_train_1), np.min(spike_train_2)))
    last_sp = max((np.max(spike_train_1), np.max(spike_train_2)))

    t_begin = first_sp - 5*std_val
    t_end = last_sp + 5*std_val

    t_vect = np.arange(t_begin, t_end, sample_period)

    sp_arr1 = np.zeros(t_vect.shape)

    for sp_time in spike_train_1:
        t_val = sp_time-t_begin
        sp_ind = np.floor_divide(t_val, sample_period)
        sp_arr1[sp_ind] = 1

    sp_arr2 = np.zeros(t_vect.shape)

    for sp_time in spike_train_2:
        t_val = sp_time-t_begin
        sp_ind = np.floor_divide(t_val, sample_period)
        sp_arr2[sp_ind] = 1

    gaus_t = np.arange(-5*std_val, 5*std_val+sample_period, sample_period)
    gaus_func = np.array(map(lambda x: gauss_func(x, sigma=std_val), gaus_t))

    sp_conv1 = np.convolve(sp_arr1, gaus_func, mode='same')
    sp_conv2 = np.convolve(sp_arr2, gaus_func, mode='same')

    sp1sp1 = sp_conv1.dot(sp_conv1)
    sp2sp2 = sp_conv2.dot(sp_conv2)
    sp1sp2 = sp_conv1.dot(sp_conv2)

    coe = sp1sp2/np.sqrt(sp1sp1*sp2sp2)

    arr_dict = {'spike1': sp_arr1,
                'conv1': sp_conv1,
                'spike2': sp_arr2,
                'conv2': sp_conv2,
                't': t_vect}

    return coe, arr_dict


def measure_spike_trains_correlation_coefficient_file(input_file, desired_ctypes, std_val=5.0, sample_period=0.01):
    """
    hdf5 file
    :param input_file:
    :param desired_ctypes: list of cell types desired for comparison
    :param std_val:
    :param sample_period:
    :return:
    """

    cor_dict = {}
    # cor_dict[cell_type] = correlation coefficient dictionary as returned by
    #   measure_spike_trains_correlation_coefficient_multi() for that cell type
    #

    with h5py.File(input_file, "r") as in_fobj:
        sp_groups = [sp_grp for sp_grp in in_fobj.keys() if 'spikes_' in sp_grp]
        sp_g = sp_groups[0]

        c_grp = in_fobj[sp_g]

        for cell_i, cell_type in enumerate(desired_ctypes):
            sp_dict = {}

            sp_vars = [sp_name for sp_name in c_grp.keys() if cell_type in sp_name]

            if sp_vars:

                print("Performing correlation coefficient analysis on {0} spike trains".format(len(sp_vars)))

                for sp_var in sp_vars:
                    if c_grp[sp_var].size >= 2:
                        sp_dict[sp_var] = np.array(c_grp[sp_var])

                if len(sp_dict.keys()) > 1:
                    cor_dict[cell_type] = measure_spike_trains_correlation_coefficient_multi(sp_dict,
                                                                                             std_val=std_val,
                                                                                             sample_period=sample_period)
                else:
                    print("Not enough spike trains found for {0}".format(cell_type))
                    cor_dict[cell_type] = 0
            else:
                print("No spikes found for {0}".format(cell_type))
                cor_dict[cell_type] = 0

    return cor_dict


def measure_spike_trains_correlation_coefficient_multi(spike_train_dict, std_val=5.0, sample_period=0.01):
    """

    :param spike_train_dict: dictionary where every key is the name of a cell and every value is a numpy array with the
        spike times in msec
    :param std_val: standard deviation used to create the gaussian function which will be convolved with spike trains
    :param sample_period: period between samples. This is used for creating a time series for each spike train. Value is
        assumed to be in msec
    :return: coe_dict: dictionary of the correlation coefficients for every combination of the spike trains
    """

    spike_names = spike_train_dict.keys()
    n_trains = len(spike_names)

    coe_dict = {}

    avg = 0
    n_comps = 0

    for i in range(n_trains-1):
        for j in range(i+1, n_trains):

            temp_coe, _ = measure_spike_trains_correlation_coefficient(spike_train_dict[spike_names[i]],
                                                                       spike_train_dict[spike_names[j]],
                                                                       std_val=std_val,
                                                                       sample_period=sample_period)

            comb_k = '{0} & {1}'.format(spike_names[i], spike_names[j])
            coe_dict[comb_k] = temp_coe
            n_comps += 1

    for k, v in coe_dict.iteritems():
        avg += v
    avg /= n_comps

    coe_dict['average'] = avg
    print(avg)

    coe_dict['sigma'] = std_val

    return coe_dict


def measure_spike_trains_isi():
    pass


if __name__ == "__main__":
    pass

    res_file = '/home/adam/Documents/Python/MEMORY_NETWORK/trunk/results/neuron_network/2017-12-05/full_dendrite_test/0p0001ACh/ca1_test_17_43_00.hdf5'
    c_list = ['CA1P', 'BCCCK', 'BCPV']

    rpath = '/home/adam/Documents/Python/MEMORY_NETWORK/trunk/results/ca3_sinusoidal_test/2018-03-10/'

    # res_file1 = os.path.join(rpath, 'a_sr_distr_lognormal_sinusoidal_10_so_alpha_n7p0_20_59_50.hdf5')
    # res_file2 = os.path.join(rpath, 'a_sr_distr_lognormal_sinusoidal_5_so_alpha_n7p0_20_28_38.hdf5')
    res_file3 = os.path.join(rpath, 'a_sr_distr_lognormal_nonsinusoidal_so_alpha_n7p0_14_46_14.hdf5')

    th_d = {'CA1P': 0, 'BCPV': -10.0, 'BCCK': -10.0, 'BIS': -10.0, 'OLM': -10.0}
    # analyze_spikes_from_file(res_file1, thresh_dict=th_d, time_int=[])
    # analyze_spikes_from_file(res_file2, thresh_dict=th_d, time_int=[])
    analyze_spikes_from_file(res_file3, thresh_dict=th_d, time_int=[])

    # measure_spike_trains_correlation_coefficient_file(res_file, c_list, std_val=5.0, sample_period=0.01)

    #
    # c_dir = os.getcwd()
    # m_dir = os.path.dirname(c_dir)
    #
    # th_d = {'CA1P': 0, 'BCPV': -10.0, 'BCCK': -10.0, 'BIS': -10.0, 'OLM': -10.0}
    #
    # # ta = np.array([-70.0, -70.0, 10.0, 50.0, -50, -70.0])
    # # tim_a = np.array([0, 10, 20, 30, 40, 50])
    # # print get_spike_times(ta, tim_a)
    #
    # r_string = 'results/epsp_comparison/2017-08-15'
    # # r_string = 'ca1_test_14_08_58.hdf5'
    # rpath = '/home/adam/Documents/Python/MEMORY_NETWORK/trunk/results/ach_variation/2017-06-25/full_nocck_redo/ach_full_0p000000_10_26_51.hdf5'
    # r_path = os.path.join(m_dir, r_string)
    #
    # r_path = '/home/adam/Documents/Python/MEMORY_NETWORK/trunk/ca1_test_17_07_41.hdf5'
    # analyze_spikes_from_file(r_path)
    # new_path = analyze_spikes_from_file(rpath, save_to_new=True, force_reanalysis=True)
    # print new_path
    # analyze_files_in_dir(r_path, save_to_new=True, force_reanalysis=True, thresh_dict=th_d, time_int=[50, 1000])

