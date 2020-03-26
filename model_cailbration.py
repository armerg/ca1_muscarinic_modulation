import datetime as dat
import numpy as np

import os
import pandas as pd
import scipy as sp
import scipy.stats as scat
import sys

print(os.getcwd())

c_dir = os.path.join(os.getcwd(), 'model_sources/migliore_etal_2018/MiglioreEtAl2018PLOSCompBiol2018')
res_dir = os.path.join(os.getcwd(), 'results/Tuning Calcium Model/RXD')

os.chdir(c_dir)

from mpl_toolkits.mplot3d import Axes3D
from itertools import cycle
from neuron import h
from os.path import join
from scipy import integrate as spint
from scipy import optimize as spopt

import matplotlib.pyplot as plt

import bokeh.io as bkio
import bokeh.layouts as blay
import bokeh.models as bmod
import bokeh.plotting as bplt
import various_scripts.frequency_analysis as fan
from bokeh.palettes import Category20 as palette
from bokeh.palettes import Category20b as paletteb
from selenium import webdriver

colrs = palette[20] + paletteb[20]

import plot_results as plt_res
#
# module_path = os.getcwd() # os.path.abspath(os.path.join('../..'))
# if module_path not in sys.path:
#     sys.path.append(module_path)

import migliore_python as mig_py

mu = u'\u03BC'
delta = u'\u0394'
tau = u'\u03C4'

# chrome_path = r'/usr/local/bin/chromedriver'
# my_opts = webdriver.chrome.options.Options()
# my_opts.add_argument('start-maximized')
# my_opts.add_argument('disable-infobars')
# my_opts.add_argument('--disable-extensions')
# my_opts.add_argument('window-size=1200,1000')
# my_opts.add_argument('--headless')


def single_exponential(x_data, tau):
    x_data = np.array(x_data)

    return np.exp(-x_data / tau)


def single_exponential_rise(x_data, tau, asym, tshift):
    return asym * (1.0 - np.exp(-(x_data - tshift) / tau))


def estimate_decay_constant(t_vals, f_vals, change_tol=0.001):
    max_i = np.argmax(f_vals)
    plus_i = range(max_i, f_vals.size)
    pos_i = np.where(f_vals > 0)

    targ_i = np.intersect1d(plus_i, pos_i)

    min_val = 0.999 * np.min(f_vals[targ_i])

    dec_t = t_vals[targ_i] - t_vals[max_i]

    dec_norm = (f_vals[targ_i] - min_val) / (f_vals[max_i] - min_val)

    opt_params, pcov = sp.optimize.curve_fit(single_exponential, dec_t, dec_norm)

    return opt_params


def run_recharge_simulation(param_dict, sim_dur=180000):

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')
    t_steps = np.arange(0, sim_dur, 100)

    cell = mig_py.MyCell(t_path, True, param_dict)

    # Calreticulin Parameter Values
    total_car = param_dict['car_total']
    KD_car = param_dict['car_KD']
    car_kon = param_dict['car_kon']
    car_koff = KD_car * car_kon

    carca_init = total_car * 0.001 / (KD_car + 0.01)
    car_init = total_car - carca_init

    ap_secs = [cell.apical[i] for i in range(10)]

    t_secs = cell.somatic + ap_secs

    h.CVode().active(True)
    h.finitialize(-65.0)

    for sec in t_secs:
        cell.ca[cell.er].nodes(sec).concentration = 0.001
        cell.carca.nodes(sec).concentration = carca_init
        cell.car.nodes(sec).concentration = car_init

    h.CVode().re_init()

    # Create Numpy Arrays For Storing Data

    ca_cyt_arr = np.zeros((t_steps.shape[0], 3))
    ca_er_arr = np.zeros((t_steps.shape[0], 3))
    carca_arr = np.zeros(t_steps.shape)

    cbdhca_arr = np.zeros(t_steps.shape)
    cbdlca_arr = np.zeros(t_steps.shape)
    ogb1ca_arr = np.zeros(t_steps.shape)

    # dyeca_arr = np.zeros(t_steps.shape)
    #
    # e1_arr = np.zeros(t_steps.shape)
    # e1_2ca_arr = np.zeros(t_steps.shape)
    # e1_2ca_p_arr = np.zeros(t_steps.shape)
    # e2_2ca_p_arr = np.zeros(t_steps.shape)
    # e2_p_arr = np.zeros(t_steps.shape)
    # e2_arr = np.zeros(t_steps.shape)
    #
    # c1_arr = np.zeros(t_steps.shape)
    # c2_arr = np.zeros(t_steps.shape)
    # c3_arr = np.zeros(t_steps.shape)
    # c4_arr = np.zeros(t_steps.shape)
    # o5_arr = np.zeros(t_steps.shape)
    # o6_arr = np.zeros(t_steps.shape)

    for t_i, t_step in enumerate(t_steps):
        h.continuerun(t_step)

        ca_cyt_arr[t_i, 0] = cell.ca[cell.cyt].nodes(cell.somatic[0])(0.5)[0].concentration
        ca_cyt_arr[t_i, 1] = cell.ca[cell.cyt].nodes(cell.apical[0])(0.5)[0].concentration
        ca_cyt_arr[t_i, 2] = cell.ca[cell.cyt].nodes(cell.apical[9])(0.5)[0].concentration
        ca_er_arr[t_i, 0] = cell.ca[cell.er].nodes(cell.somatic[0])(0.5)[0].concentration
        ca_er_arr[t_i, 1] = cell.ca[cell.er].nodes(cell.apical[0])(0.5)[0].concentration
        ca_er_arr[t_i, 2] = cell.ca[cell.er].nodes(cell.apical[9])(0.5)[0].concentration

        cbdhca_arr[t_i] = cell.cbdhca[cell.cyt].nodes(cell.somatic[0])(0.5)[0].concentration
        cbdlca_arr[t_i] = cell.cbdlca[cell.cyt].nodes(cell.somatic[0])(0.5)[0].concentration
        ogb1ca_arr[t_i] = cell.dyeca[cell.cyt].nodes(cell.somatic[0])(0.5)[0].concentration

        carca_arr[t_i] = cell.carca[cell.er].nodes(cell.somatic[0])(0.5)[0].concentration
        # dyeca_arr[t_i] = dyeca.nodes(soma)(0.5)[0].concentration
        #
        # e1_arr[t_i] = e1_serca.nodes(soma)(0.5)[0].concentration
        # e1_2ca_arr[t_i] = e1_2ca_serca.nodes(soma)(0.5)[0].concentration
        # e1_2ca_p_arr[t_i] = e1_2ca_p_serca.nodes(soma)(0.5)[0].concentration
        # e2_2ca_p_arr[t_i] = e2_2ca_p_serca.nodes(soma)(0.5)[0].concentration
        # e2_p_arr[t_i] = e2_p_serca.nodes(soma)(0.5)[0].concentration
        # e2_arr[t_i] = e2_serca.nodes(soma)(0.5)[0].concentration
        #
        # c1_arr[t_i] = c1_ip3r.nodes(soma)(0.5)[0].concentration
        # c2_arr[t_i] = c2_ip3r.nodes(soma)(0.5)[0].concentration
        # c3_arr[t_i] = c3_ip3r.nodes(soma)(0.5)[0].concentration
        # c4_arr[t_i] = c4_ip3r.nodes(soma)(0.5)[0].concentration
        # o5_arr[t_i] = o5_ip3r.nodes(soma)(0.5)[0].concentration
        # o6_arr[t_i] = o6_ip3r.nodes(soma)(0.5)[0].concentration

    print('Final SERCA States')
    # e1 = cell.e1_serca.nodes(cell.somatic[0])(0.5)[0].concentration
    # e1_2ca = cell.e1_2ca_serca.nodes(cell.somatic[0])(0.5)[0].concentration
    # e1_2ca_p = cell.e1_2ca_p_serca.nodes(cell.somatic[0])(0.5)[0].concentration
    # e2_2ca_p = cell.e2_2ca_p_serca.nodes(cell.somatic[0])(0.5)[0].concentration
    # e2_p = cell.e2_p_serca.nodes(cell.somatic[0])(0.5)[0].concentration
    # e2 = cell.e2_serca.nodes(cell.somatic[0])(0.5)[0].concentration

    # total = e1 + e1_2ca + e1_2ca_p + e2_2ca_p + e2_p + e2
    #
    # print('e1: {}'.format(e1/total))
    # print('e1-2ca: {}'.format(e1_2ca / total))
    # print('e1-2ca-p: {}'.format(e1_2ca_p / total))
    # print('e2-2ca-p: {}'.format(e2_2ca_p / total))
    # print('e2-p: {}'.format(e2_p / total))
    # print('e2: {}'.format(e2 / total))

    result_d = {'t': t_steps,
                'sec names': ['Soma', 'Proximal Apical Trunk', 'Distal Apical Trunk'],
                'cyt ca': ca_cyt_arr,
                'er ca': ca_er_arr,
                'carca': carca_arr,
                'cbdhca': cbdhca_arr,
                'cbdlca': cbdlca_arr,
                'ogb1ca': ogb1ca_arr,
                }

    return result_d


def plot_recharge(result_d):
    t_steps = result_d['t'] / 1000.0

    ret_figs = []

    cyt_ca_fig = bplt.figure(title='Cytosol Calcium vs Time')
    ret_figs.append(cyt_ca_fig)
    cyt_ca_fig.xaxis.axis_label = 'time (seconds)'
    cyt_ca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)

    er_ca_fig = bplt.figure(title='Endoplasmic Reticulum Calcium vs Time')
    ret_figs.append(er_ca_fig)
    er_ca_fig.xaxis.axis_label = 'time (seconds)'
    er_ca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)

    carca_fig = bplt.figure(title='Bound Calreticulin vs Time')
    ret_figs.append(carca_fig)
    carca_fig.xaxis.axis_label = 'time (seconds)'
    carca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    carca_fig.line(t_steps, 1000.0*result_d['carca'][:], line_width=3, color=colrs[0], legend='Somatic ER')

    cbd_fig = bplt.figure(title='Bound Calbindin-D28k vs Time')
    ret_figs.append(cbd_fig)
    cbd_fig.xaxis.axis_label = 'time (seconds)'
    cbd_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    cbd_fig.line(t_steps, 1000.0 * result_d['cbdhca'][:], line_width=3, color=colrs[0], legend='High')
    cbd_fig.line(t_steps, 1000.0 * result_d['cbdlca'][:], line_width=3, color=colrs[2], legend='Low')

    ogb1_fig = bplt.figure(title='Bound OGB-1 vs Time')
    ret_figs.append(ogb1_fig)
    ogb1_fig.xaxis.axis_label = 'time (seconds)'
    ogb1_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ogb1_fig.line(t_steps, 1000.0 * result_d['ogb1ca'][:], line_width=3, color=colrs[0], legend='High')

    # for l_i, loc_name in enumerate(result_d['sec names']):
    #     if 'soma' in loc_name:
    #         cyt_ca_fig.line(t_steps, 1000.0*result_d['cyt ca'][:, l_i], line_width=4, color='black', legend='model result')
    #         er_ca_fig.line(t_steps, 1000.0*result_d['er ca'][:, l_i], line_width=4, color='black', legend='model result')

    for l_i, loc_name in enumerate(result_d['sec names']):

        cyt_ca_fig.line(t_steps, 1000.0*result_d['cyt ca'][:, l_i], line_width=4, color=colrs[l_i], legend=loc_name)
        er_ca_fig.line(t_steps, 1000.0*result_d['er ca'][:, l_i], line_width=4, color=colrs[l_i], legend=loc_name)

    cyt_ca_fig.legend.location = 'bottom_right'

    return ret_figs    


def run_current_injection(rxd_sim, param_dict={}, sim_dur=500, c_int=[50, 100], c_amp=1):
    """

    :param cell_type: type of neuron cell to test\n
    :param conf_obj: configuration dictionary which holds all of the parameters for the neuron cell to be tested\n
    :param swc_path: location of swc_file, required for CA1PyramidalCell class\n
    :param sim_dur: simulation duration (msec)\n
    :param c_int: interval of time current injection pulse is active\n
    :param c_amp: amplitude of current injection\n
    :return: t_array, v_array, i_array: t_array\n
        time array = numpy array of simulation time steps\n
        v_array = numpy array potentials for soma[0](0.5) of cell\n
        i_array = numpy array of values for the injected current\n
    """

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = mig_py.MyCell(t_path, rxd_sim, param_dict)

    t_curr = h.IClamp(cell.somatic[0](0.5))

    t_curr.delay = c_int[0]
    t_curr.amp = c_amp
    t_curr.dur = c_int[1] - c_int[0]

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values
    i_rec = h.Vector().record(t_curr._ref_i)
    t = h.Vector().record(h._ref_t)
    soma_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    apic_v = h.Vector().record(cell.apical[9](0.5)._ref_v)
    axon_v = h.Vector().record(cell.axonal[0](0.5)._ref_v)

    s_ica_l = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cal)
    a_ica_l = h.Vector().record(cell.apical[9](0.5)._ref_ica_cal)
    s_ica_n = h.Vector().record(cell.somatic[0](0.5)._ref_ica_can)
    a_ica_n = h.Vector().record(cell.apical[9](0.5)._ref_ica_can)
    s_ica_t = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cat)
    a_ica_t = h.Vector().record(cell.apical[9](0.5)._ref_ica_cat)

    if not rxd_sim:
        s_ca_cyt = h.Vector().record(cell.somatic[0](0.5)._ref_cai)
        a_ca_cyt = h.Vector().record(cell.apical[9](0.5)._ref_cai)
    else:
        s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
        a_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)
        ax_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.axonal[0])[0]._ref_concentration)

        s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
        a_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[9])[0]._ref_concentration)

        s_dyeca = h.Vector().record(cell.dyeca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_dyeca = h.Vector().record(cell.dyeca.nodes(cell.apical[9])[0]._ref_concentration)

        s_cbdhca = h.Vector().record(cell.cbdhca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_cbdhca = h.Vector().record(cell.cbdhca.nodes(cell.apical[9])[0]._ref_concentration)

        s_cbdlca = h.Vector().record(cell.cbdlca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_cbdlca = h.Vector().record(cell.cbdlca.nodes(cell.apical[9])[0]._ref_concentration)

        s_carca = h.Vector().record(cell.carca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_carca = h.Vector().record(cell.carca.nodes(cell.apical[9])[0]._ref_concentration)

        s_ip3r = h.Vector().record(cell.ro_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
        a_ip3r = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[9])[0]._ref_concentration)

    h.cvode.active(1)

    h.v_init = -69.4

    h.tstop = sim_dur
    h.celsius = 34.0

    # h.load_file('negative_init.hoc')
    # h.init()

    h.stdinit()

    print('Running current injection, amplitude = {0} nA'.format(c_amp))

    h.continuerun(sim_dur)

    print('Final IP3R Values')
    print('r: {0}'.format(cell.r_ip3r.nodes(cell.somatic[0]).concentration))
    print('ri: {0}'.format(cell.ri_ip3r.nodes(cell.somatic[0]).concentration))
    print('ro: {0}'.format(cell.ro_ip3r.nodes(cell.somatic[0]).concentration))
    print('rc: {0}'.format(cell.rc_ip3r.nodes(cell.somatic[0]).concentration))
    print('rc2: {0}'.format(cell.rc2_ip3r.nodes(cell.somatic[0]).concentration))
    print('rc3: {0}'.format(cell.rc3_ip3r.nodes(cell.somatic[0]).concentration))
    print('rc4: {0}'.format(cell.rc4_ip3r.nodes(cell.somatic[0]).concentration))

    print('Final IP3 Production Numbers')
    print('ip3: {0}'.format(cell.ip3.nodes(cell.somatic[0]).concentration))
    print('plc: {0}'.format(cell.plc_m1.nodes(cell.somatic[0]).concentration))
    print('ga_gtp: {0}'.format(cell.ga_gtp_m1.nodes(cell.somatic[0]).concentration))
    print('ga_gtp_plc: {0}'.format(cell.ga_gtp_plc_m1.nodes(cell.somatic[0]).concentration))
    print('ip5p: {0}'.format(cell.ip5p.nodes(cell.somatic[0]).concentration))
    print('ip5p_ip3: {0}'.format(cell.ip5p_ip3.nodes(cell.somatic[0]).concentration))
    print('ip3k: {0}'.format(cell.ip3k.nodes(cell.somatic[0]).concentration))
    print('ip3k_2ca: {0}'.format(cell.ip3k_2ca.nodes(cell.somatic[0]).concentration))
    print('ip3k_2ca_ip3: {0}'.format(cell.ip3k_2ca_ip3.nodes(cell.somatic[0]).concentration))

    res_dict = {'t': np.array(t.as_numpy()),
                'i_stim': np.array(i_rec.as_numpy()),
                'soma_v': np.array(soma_v),
                'soma_cyt': np.array(s_ca_cyt),
                'axon_v': np.array(axon_v),
                'axon_cyt': np.array(ax_ca_cyt),
                'apic9_v': np.array(apic_v),
                'apic9_cyt': np.array(a_ca_cyt),
                'soma_ica_l': np.array(s_ica_l),
                'apic9_ica_l': np.array(a_ica_l),
                'soma_ica_n': np.array(s_ica_n),
                'apic9_ica_n': np.array(a_ica_n),
                'soma_ica_t': np.array(s_ica_t),
                'apic9_ica_t': np.array(a_ica_t),
    }

    if rxd_sim:
        res_dict['soma_er'] = np.array(s_ca_er)
        res_dict['apic9_er'] = np.array(a_ca_er)
        res_dict['soma_ip3r_open'] = np.array(s_ip3r)
        res_dict['apic9_ip3r_open'] = np.array(a_ip3r)
        res_dict['soma_dyeca'] = np.array(s_dyeca)
        res_dict['apic9_dyeca'] = np.array(a_dyeca)
        res_dict['soma_cbdhca'] = np.array(s_cbdhca)
        res_dict['apic9_cbdhca'] = np.array(s_cbdhca)
        res_dict['soma_cbdlca'] = np.array(s_cbdlca)
        res_dict['apic9_cbdlca'] = np.array(s_cbdlca)
        res_dict['soma_carca'] = np.array(s_carca)
        res_dict['apic9_carca'] = np.array(a_carca)
    return res_dict

    
def run_current_injection_series(rxd_sim, param_dict={}, sim_dur=500, pulse_times=[50], pulse_amps=[1.0], pulse_length=10):

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = mig_py.MyCell(t_path, rxd_sim, param_dict)

    t_curr = h.IClamp(cell.somatic[0](0.5))

    amp_list = []
    amp_time = []
    for p_time, p_amp in zip(pulse_times, pulse_amps):
        amp_list += [p_amp, 0.0]
        amp_time += [p_time, p_time+pulse_length]

    c_vec = h.Vector().from_python(amp_list)
    c_time = h.Vector().from_python(amp_time)

    t_curr.delay = 0
    t_curr.dur = 1e9
    c_vec.play(t_curr._ref_amp, c_time)

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values
    i_rec = h.Vector().record(t_curr._ref_i)
    t = h.Vector().record(h._ref_t)
    soma_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    apic_v = h.Vector().record(cell.apical[9](0.5)._ref_v)
    axon_v = h.Vector().record(cell.axonal[0](0.5)._ref_v)

    s_ica_l = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cal)
    a_ica_l = h.Vector().record(cell.apical[9](0.5)._ref_ica_cal)
    s_ica_n = h.Vector().record(cell.somatic[0](0.5)._ref_ica_can)
    a_ica_n = h.Vector().record(cell.apical[9](0.5)._ref_ica_can)
    s_ica_t = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cat)
    a_ica_t = h.Vector().record(cell.apical[9](0.5)._ref_ica_cat)

    if not rxd_sim:
        s_ca_cyt = h.Vector().record(cell.somatic[0](0.5)._ref_cai)
        a_ca_cyt = h.Vector().record(cell.apical[9](0.5)._ref_cai)

    else:
        s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
        a_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)
        ax_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.axonal[0])[0]._ref_concentration)
        s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
        a_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[9])[0]._ref_concentration)

        s_dyeca = h.Vector().record(cell.dyeca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_dyeca = h.Vector().record(cell.dyeca.nodes(cell.apical[9])[0]._ref_concentration)

        s_cbdhca = h.Vector().record(cell.cbdhca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_cbdhca = h.Vector().record(cell.cbdhca.nodes(cell.apical[9])[0]._ref_concentration)

        s_cbdlca = h.Vector().record(cell.cbdlca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_cbdlca = h.Vector().record(cell.cbdlca.nodes(cell.apical[9])[0]._ref_concentration)

        s_carca = h.Vector().record(cell.carca.nodes(cell.somatic[0])[0]._ref_concentration)
        a_carca = h.Vector().record(cell.carca.nodes(cell.apical[9])[0]._ref_concentration)

    h.cvode.active(1)

    h.v_init = -69.4

    h.tstop = sim_dur
    h.celsius = 34.0

    # h.load_file('negative_init.hoc')
    # h.init()

    h.stdinit()

    print('Running current injection')

    h.continuerun(sim_dur)

    res_dict = {'t': np.array(t.as_numpy()),
                'i_stim': np.array(i_rec.as_numpy()),
                'soma_v': np.array(soma_v),
                'soma_cyt': np.array(s_ca_cyt),
                'apic9_v': np.array(apic_v),
                'apic9_cyt': np.array(a_ca_cyt),
                'axon_v': np.array(axon_v),
                'axon_cyt': np.array(ax_ca_cyt),
                'soma_ica_l': np.array(s_ica_l),
                'apic9_ica_l': np.array(a_ica_l),
                'soma_ica_n': np.array(s_ica_n),
                'apic9_ica_n': np.array(a_ica_n),
                'soma_ica_t': np.array(s_ica_t),
                'apic9_ica_t': np.array(a_ica_t)
                }

    if rxd_sim:
        res_dict['soma_er'] = np.array(s_ca_er)
        res_dict['apic9_er'] = np.array(a_ca_er)
        res_dict['soma_dyeca'] = np.array(s_dyeca)
        res_dict['apic9_dyeca'] = np.array(a_dyeca)
        res_dict['soma_cbdhca'] = np.array(s_cbdhca)
        res_dict['apic9_cbdhca'] = np.array(a_cbdhca)
        res_dict['soma_cbdlca'] = np.array(s_cbdlca)
        res_dict['apic9_cbdlca'] = np.array(a_cbdlca)
        res_dict['soma_carca'] = np.array(s_carca)
        res_dict['apic9_carca'] = np.array(a_carca)

    return res_dict


def plot_current_injection(result_d, rxd_bool, t_inj, t_ignore=0):
    result_figs = []
    t_ig = np.squeeze(np.where(result_d['t'] > t_ignore))
    i_inj = np.squeeze(np.where(result_d['t'] >= t_inj))[0]
    t_arr = result_d['t'][t_ig] - t_ignore

    v_fig = bplt.figure(title='Membrane Potential vs Time')
    result_figs.append(v_fig)
    v_fig.line(t_arr, result_d['soma_v'][t_ig], line_width=3, color='blue', legend='Soma')
    v_fig.line(t_arr, result_d['apic9_v'][t_ig], line_width=3, color='green', legend='Apical Trunk',
               line_dash='dashed')

    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'

    i_fig = bplt.figure(title='Current Injection vs Time')
    result_figs.append(i_fig)
    i_fig.line(t_arr, result_d['i_stim'][t_ig], line_width=3, color='blue', legend='Soma')

    i_fig.xaxis.axis_label = 'time (msec)'
    i_fig.yaxis.axis_label = 'current (nA)'

    cai_fig = bplt.figure(title='Intracellular Calcium vs Time')
    result_figs.append(cai_fig)
    cai_fig.line(t_arr, result_d['soma_cyt'][t_ig] * 1000.0, line_width=3, color='blue', legend='Soma')
    cai_fig.line(t_arr, result_d['apic9_cyt'][t_ig] * 1000.0, line_width=3, color='green', legend='Apical Trunk',
                 line_dash='dashed')
    cai_fig.xaxis.axis_label = 'time (msec)'
    cai_fig.yaxis.axis_label = '[Ca] ({}M)'.format(mu)

    ica_fig = bplt.figure(title='Calcium Currents vs Time')
    result_figs.append(ica_fig)
    ica_fig.line(t_arr, result_d['soma_ica_l'][t_ig], line_width=3, color=colrs[0], legend='Soma - L')
    ica_fig.line(t_arr, result_d['apic9_ica_l'][t_ig], line_width=3, color=colrs[1], legend='Apical Trunk - L')
    ica_fig.line(t_arr, result_d['soma_ica_n'][t_ig], line_width=3, color=colrs[2], legend='Soma - N')
    ica_fig.line(t_arr, result_d['apic9_ica_n'][t_ig], line_width=3, color=colrs[3], legend='Apical Trunk - N')
    ica_fig.line(t_arr, result_d['soma_ica_t'][t_ig], line_width=3, color=colrs[4], legend='Soma - T')
    ica_fig.line(t_arr, result_d['apic9_ica_t'][t_ig], line_width=3, color=colrs[5], legend='Apical Trunk - T')

    ica_fig.xaxis.axis_label = 'time (msec)'
    ica_fig.yaxis.axis_label = 'current (mA/cm^2)'

    if rxd_bool:
        caer_fig = bplt.figure(title='Endoplasmic Reticulum Calcium vs Time')
        result_figs.append(caer_fig)
        caer_fig.line(t_arr, result_d['soma_er'][t_ig] * 1000.0, line_width=3, color='blue', legend='Soma')
        caer_fig.line(t_arr, result_d['apic9_er'][t_ig] * 1000.0, line_width=3, color='green', legend='Apical Trunk')
        caer_fig.xaxis.axis_label = 'time (msec)'
        caer_fig.yaxis.axis_label = '[Ca]_ER ({}M)'.format(mu)

        cbd_fig = bplt.figure(title='Bound Calbindin D28k vs Time')
        result_figs.append(cbd_fig)
        cbd_fig.line(t_arr, result_d['soma_cbdhca'][t_ig] * 1000.0, line_width=3, color='blue', legend='Soma - HA')
        cbd_fig.line(t_arr, result_d['apic9_cbdhca'][t_ig] * 1000.0, line_width=3, color='green',
                     legend='Apical Trunk - HA')
        cbd_fig.line(t_arr, result_d['soma_cbdlca'][t_ig] * 1000.0, line_width=3, color='blue', legend='Soma - LA',
                     line_dash='dashed')
        cbd_fig.line(t_arr, result_d['apic9_cbdlca'][t_ig] * 1000.0, line_width=3, color='green',
                     legend='Apical Trunk- LA', line_dash='dashed')
        cbd_fig.xaxis.axis_label = 'time (msec)'
        cbd_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)

        dye_fig = bplt.figure(title='Change in Fluorescence vs Time')
        result_figs.append(dye_fig)

        dF_s = 100 * (result_d['soma_dyeca'][t_ig] - result_d['soma_dyeca'][i_inj]) / result_d['soma_dyeca'][i_inj]
        dF_a = 100 * (result_d['apic9_dyeca'][t_ig] - result_d['apic9_dyeca'][i_inj]) / result_d['apic9_dyeca'][i_inj]

        dye_fig.line(t_arr, dF_s, line_width=3, color='blue', legend='Soma')
        dye_fig.line(t_arr, dF_a, line_width=3, color='green', legend='Apical Trunk')
        dye_fig.xaxis.axis_label = 'time (msec)'
        dye_fig.yaxis.axis_label = '{}F (%)'.format(delta)

        carca_fig = bplt.figure(title='Bound Calreticulin vs Time')
        result_figs.append(carca_fig)
        carca_fig.line(t_arr, result_d['soma_carca'][t_ig] * 1000.0, line_width=3, color='blue', legend='Soma')
        carca_fig.line(t_arr, result_d['apic9_carca'][t_ig] * 1000.0, line_width=3, color='green', legend='Distal Apical Trunk')
        carca_fig.xaxis.axis_label = 'time (msec)'
        carca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)

        ip3r_fig = bplt.figure(title='Open IP3R vs Time')
        result_figs.append(ip3r_fig)
        ip3r_fig.line(t_arr, result_d['soma_ip3r_open'][t_ig]*100, line_width=3, color='blue', legend='Soma')
        ip3r_fig.line(t_arr, result_d['apic9_ip3r_open'][t_ig]*100, line_width=3, color='green',
                      legend='Distal Apical Trunk')
        ip3r_fig.xaxis.axis_label = 'time (msec)'
        ip3r_fig.yaxis.axis_label = 'Percent Open'

    return result_figs


def run_ip3_pulse(t_steps, param_dict={}, pulse_times=[50], pulse_amps=[0.001], pulse_amps_high=[0.001],
                  pulse_length=10, current_dict={'amp': 0.0, 'dur': 100.0, 'start': 0.0},
                  im_inhib_dict={'perc': 1.0, 'start': 0.0, 'end': 0.0}):
    p_times = []
    p_amps_high = []
    p_amps_f = []
    p_amps_b = []

    for p_t, p_a, p_h in zip(pulse_times, pulse_amps, pulse_amps_high):
        p_times += [p_t, p_t+pulse_length]
        p_amps_high += [p_h, 0]
        p_amps_f += [p_a, 0]
        p_amps_b += [0, p_a]

    print('Pulse Times: {}'.format(p_times))
    print('Pulse Amps: {}'.format(p_amps_f))

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = mig_py.MyCell(t_path, True, param_dict)

    if current_dict['amp']:
        t_curr = h.IClamp(cell.somatic[0](0.5))
        t_curr.amp = current_dict['amp']
        t_curr.delay = current_dict['start']
        t_curr.dur = current_dict['dur']

    n_nodes = len(cell.ca[cell.cyt].nodes)
    apical_trunk_inds = [0, 8, 9, 11, 13, 19]
    sec_list = [sec for sec in cell.somatic] + [cell.apical[i] for i in apical_trunk_inds]
    apic_names = ['apical_{0}'.format(num) for num in apical_trunk_inds]
    sec_names = ['soma_0'] + apic_names

    h.distance(0, cell.somatic[0](0.5))

    node_dists = []

    for sec in sec_list:
        for seg in sec:
            node_dists.append(h.distance(sec(seg.x)))

    node_locs = []

    for sec_name in sec_names:
        sec = cell.get_sec_by_name(sec_name)

        for node in cell.ca[cell.cyt].nodes:
            if node.satisfies(sec):
                node_locs.append('{0}({1:.3})'.format(sec_name, node.x))

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values
    t_vec = h.Vector().record(h._ref_t)

    s_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    a0_v = h.Vector().record(cell.apical[0](0.5)._ref_v)
    a9_v = h.Vector().record(cell.apical[9](0.5)._ref_v)

    s_isk = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kca)
    a0_isk = h.Vector().record(cell.apical[0](0.5)._ref_ik_kca)
    a9_isk = h.Vector().record(cell.apical[9](0.5)._ref_ik_kca)

    s_im = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kmb_inh)
    ax_im = h.Vector().record(cell.axonal[0](0.5)._ref_ik_kmb_inh)

    s_ica_l = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cal)
    a_ica_l = h.Vector().record(cell.apical[9](0.5)._ref_ica_cal)
    s_ica_n = h.Vector().record(cell.somatic[0](0.5)._ref_ica_can)
    a_ica_n = h.Vector().record(cell.apical[9](0.5)._ref_ica_can)
    s_ica_t = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cat)
    a_ica_t = h.Vector().record(cell.apical[9](0.5)._ref_ica_cat)

    s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[0])[0]._ref_concentration)
    a9_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)

    s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[0])[0]._ref_concentration)
    a9_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[9])[0]._ref_concentration)

    s_ip3 = h.Vector().record(cell.ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[0])[0]._ref_concentration)
    a9_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[9])[0]._ref_concentration)

    s_po = h.Vector().record(cell.ro_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    a9_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[9])[0]._ref_concentration)

    cyt_vals = np.zeros((t_steps.size, n_nodes))
    er_vals = np.zeros((t_steps.size, n_nodes))
    ip3_vals = np.zeros((t_steps.size, n_nodes))
    ip3r_open_vals = np.zeros((t_steps.size, n_nodes))

    cv_act = 1
    print('CVode: {0}'.format(cv_act))
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    h.stdinit()

    print('Running IP3 pulse')

    for t_i, t_step in enumerate(t_steps):
        h.continuerun(t_step)

        if t_step in p_times:
            p_i = p_times.index(t_step)
            print('time: {0}, ip3 rate: {1}'.format(t_step, p_amps_f[p_i]))

            cell.ip3.nodes(cell.cyt).concentration = p_amps_f[p_i]

            cell.ip3.nodes(cell.apical[apical_trunk_inds[2]]).concentration = p_amps_high[p_i]

            # ip3_f = p_amps_f[p_i]
            # ip3_b = p_amps_b[p_i]
            #
            # cell.ip3_prod.b_rate = ip3_f
            # cell.ip3_prod.b_rate = ip3_b

            h.CVode().re_init()

            # cell.ip3.concentration = ip3_amp

        if im_inhib_dict['perc'] < 1.0:
            if t_step == im_inhib_dict['start']:
                for soma_sec in cell.somatic:
                    for seg in soma_sec:
                        seg.perc_test_kmb_inh = im_inhib_dict['perc']
            elif t_step == im_inhib_dict['end']:
                for soma_sec in cell.somatic:
                    for seg in soma_sec:
                        seg.perc_test_kmb_inh = 1.0

        cyt_vals[t_i, :] = cell.ca[cell.cyt].nodes.concentration
        er_vals[t_i, :] = cell.ca[cell.er].nodes.concentration
        ip3_vals[t_i, :] = cell.ip3.nodes.concentration
        # ip3r_open_vals[t_i, :] = np.array(cell.o5_ip3r.nodes.concentration) + np.array(cell.o6_ip3r.nodes.concentration)
        ip3r_open_vals[t_i, :] = np.array(cell.ro_ip3r.nodes.concentration)

    res_dict = {'node names': node_locs,
                'node distances': node_dists,
                'ip3 vals': ip3_vals,
                'cyt vals': cyt_vals,
                'er vals': er_vals,
                'soma v': np.array(s_v),
                'apical0 v': np.array(a0_v),
                'apical9 v': np.array(a9_v),
                'ip3r open vals': ip3r_open_vals,
                'soma cyt time': np.array(s_ca_cyt),
                'apical0 cyt time': np.array(a0_ca_cyt),
                'apical9 cyt time': np.array(a9_ca_cyt),
                'soma er time': np.array(s_ca_er),
                'apical0 er time': np.array(a0_ca_er),
                'apical9 er time': np.array(a9_ca_er),
                'soma ip3 time': np.array(s_ip3),
                'apical0 ip3 time': np.array(a0_ip3),
                'apical9 ip3 time': np.array(a9_ip3),
                'soma im': np.array(s_im),
                'axon im': np.array(ax_im),
                'soma isk': np.array(s_isk),
                'apical0 isk': np.array(a0_isk),
                'apical9 isk': np.array(a9_isk),
                'soma ica_l': np.array(s_ica_l),
                'apic9 ica_l': np.array(a_ica_l),
                'soma ica_n': np.array(s_ica_n),
                'apic9 ica_n': np.array(a_ica_n),
                'soma ica_t': np.array(s_ica_t),
                'apic9 ica_t': np.array(a_ica_t),
                'soma open probability': np.array(s_po),
                'apical0 open probability': np.array(a0_po),
                'apical9 open probability': np.array(a9_po),
                't': np.array(t_vec),
                }

    return res_dict


def plot_ip3_pulse(result_d, t_is=[0, -1]):
    ret_figs = []
    
    cyt_t_fig = bplt.figure(title='Cytosol Calcium vs Time')
    cyt_t_fig.xaxis.axis_label = 'time (msec)'
    cyt_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(cyt_t_fig)
    er_t_fig = bplt.figure(title='ER Calcium vs Time')
    er_t_fig.xaxis.axis_label = 'time (msec)'
    er_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(er_t_fig)
    ip3_t_fig = bplt.figure(title='IP3 vs Time')
    ip3_t_fig.xaxis.axis_label = 'time (msec)'
    ip3_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ip3_t_fig)
    ip3r_open_t_fig = bplt.figure(title='Open IP3R vs Time')
    ip3r_open_t_fig.xaxis.axis_label = 'time (msec)'
    ip3r_open_t_fig.yaxis.axis_label = 'probability (%)'
    ret_figs.append(ip3r_open_t_fig)
    isk_t_fig = bplt.figure(title='SK Current vs Time')
    isk_t_fig.xaxis.axis_label = 'time (msec)'
    isk_t_fig.yaxis.axis_label = 'current (mA/cm^2)'.format(mu)
    ret_figs.append(isk_t_fig)
    im_t_fig = bplt.figure(title='M Current vs Time')
    im_t_fig.xaxis.axis_label = 'time (msec)'
    im_t_fig.yaxis.axis_label = 'current (mA/cm^2)'.format(mu)
    ret_figs.append(im_t_fig)
    ica_t_fig = bplt.figure(title='Ca Current vs Time')
    ica_t_fig.xaxis.axis_label = 'time (msec)'
    ica_t_fig.yaxis.axis_label = 'current (mA/cm^2)'.format(mu)
    ret_figs.append(ica_t_fig)
    v_fig = bplt.figure(title='Membrane Potential vs Time')
    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(v_fig)

    cyt_t_fig.line(result_d['t'], result_d['soma cyt time'] * 1000.0, line_width=3, color=colrs[0], legend='somatic[0]')
    cyt_t_fig.line(result_d['t'], result_d['apical0 cyt time'] * 1000.0, line_width=3, color=colrs[1], legend='apical[0]')
    cyt_t_fig.line(result_d['t'], result_d['apical9 cyt time'] * 1000.0, line_width=3, color=colrs[2], legend='apical[9]',
                   line_dash='dashed')

    er_t_fig.line(result_d['t'], result_d['soma er time'] * 1000.0, line_width=3, color=colrs[0], legend='somatic[0]')
    er_t_fig.line(result_d['t'], result_d['apical0 er time'] * 1000.0, line_width=3, color=colrs[1], legend='apical[0]')
    er_t_fig.line(result_d['t'], result_d['apical9 er time'] * 1000.0, line_width=3, color=colrs[2], legend='apical[9]',
                  line_dash='dashed')

    ip3_t_fig.line(result_d['t'], result_d['soma ip3 time'] * 1000.0, line_width=3, color=colrs[0], legend='somatic[0]')
    ip3_t_fig.line(result_d['t'], result_d['apical0 ip3 time'] * 1000.0, line_width=3, color=colrs[1], legend='apical[0]')
    ip3_t_fig.line(result_d['t'], result_d['apical9 ip3 time'] * 1000.0, line_width=3, color=colrs[2], legend='apical[9]',
                   line_dash='dashed')

    ip3r_open_t_fig.line(result_d['t'], result_d['soma open probability'], line_width=3, color=colrs[0], legend='somatic[0]')
    ip3r_open_t_fig.line(result_d['t'], result_d['apical0 open probability'], line_width=3, color=colrs[1],
                   legend='apical[0]')
    ip3r_open_t_fig.line(result_d['t'], result_d['apical9 open probability'], line_width=3, color=colrs[2],
                   legend='apical[9]',
                   line_dash='dashed')

    im_t_fig.line(result_d['t'], result_d['soma im'], line_width=3, color=colrs[0], legend='somatic[0]')
    im_t_fig.line(result_d['t'], result_d['axon im'], line_width=3, color=colrs[0], legend='axon[0]')

    isk_t_fig.line(result_d['t'], result_d['soma isk'], line_width=3, color=colrs[0], legend='somatic[0]')
    isk_t_fig.line(result_d['t'], result_d['apical0 isk'], line_width=3, color=colrs[1],
                  legend='apical[0]')
    isk_t_fig.line(result_d['t'], result_d['apical9 isk'], line_width=3, color=colrs[2],
                  legend='apical[9]', line_dash='dashed')

    ica_t_fig.line(result_d['t'], result_d['soma ica_l'], line_width=3, color=colrs[0], legend='somatic[0], i_cal')
    ica_t_fig.line(result_d['t'], result_d['apic9 ica_l'], line_width=3, color=colrs[1],
                  legend='apical[9], i_cal')
    ica_t_fig.line(result_d['t'], result_d['soma ica_n'], line_width=3, color=colrs[2], legend='somatic[0], i_can')
    ica_t_fig.line(result_d['t'], result_d['apic9 ica_n'], line_width=3, color=colrs[3],
                   legend='apical[9], i_can')
    ica_t_fig.line(result_d['t'], result_d['soma ica_t'], line_width=3, color=colrs[4], legend='somatic[0], i_cat')
    ica_t_fig.line(result_d['t'], result_d['apic9 ica_t'], line_width=3, color=colrs[5],
                   legend='apical[9], i_cat')

    v_fig.line(result_d['t'], result_d['soma v'], line_width=3, color=colrs[0], legend='soma[0](0.5)')
    v_fig.line(result_d['t'], result_d['apical0 v'], line_width=3, color=colrs[1], legend='apical[0](0.5)')
    v_fig.line(result_d['t'], result_d['apical9 v'], line_width=3, color=colrs[2], legend='apical[9](0.5)',
               line_dash='dashed')

    # cal_fig = bplt.figure(title='Cytosol Calcium vs Location')
    # cal_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(cal_fig)
    # er_fig = bplt.figure(title='ER Calcium vs Location')
    # er_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(er_fig)
    # ip3_fig = bplt.figure(title='IP3 vs Location')
    # ip3_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(ip3_fig)
    # ip3r_open_fig = bplt.figure(title='Open IP3R vs Location')
    # ip3r_open_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(ip3r_open_fig)
    #
    # for s_i, t_i in enumerate(t_is):
    #     cal_fig.circle(result_d['node distances'], 1000.0*result_d['cyt vals'][t_i, :], size=20-s_i, color=colrs[s_i],
    #                    legend='t={0} msec'.format(my_steps[t_i]))
    #     cal_fig.line(result_d['node distances'], 1000.0 * result_d['cyt vals'][t_i, :], line_width=3,
    #                  color=colrs[s_i])
    #     er_fig.circle(result_d['node distances'], 1000.0*result_d['er vals'][t_i, :], size=20-s_i, color=colrs[s_i],
    #                   legend='t={0} msec'.format(my_steps[t_i]))
    #     er_fig.line(result_d['node distances'], 1000.0 * result_d['er vals'][t_i, :], line_width=3,
    #                 color=colrs[s_i])
    #     ip3_fig.circle(result_d['node distances'], 1000.0*result_d['ip3 vals'][t_i, :], size=20-s_i, color=colrs[s_i],
    #                    legend='t={0} msec'.format(my_steps[t_i]))
    #     ip3_fig.line(result_d['node distances'], 1000.0 * result_d['ip3 vals'][t_i, :], line_width=3,
    #                  color=colrs[s_i])
    #     ip3r_open_fig.circle(result_d['node distances'], 1e6*result_d['ip3r open vals'][t_i, :], size=20-s_i,
    #                          color=colrs[s_i], legend='t={0} msec'.format(my_steps[t_i]))
    #     ip3r_open_fig.line(result_d['node distances'], 1e6*result_d['ip3r open vals'][t_i, :],
    #                        line_width=3, color=colrs[s_i])

    return ret_figs


def run_ach_pulse(t_steps, param_dict={}, pulse_times=[50], pulse_amps=[0.001],
                  pulse_length=10, current_dict={'amp': 0.0, 'start': 0.0, 'dur': 0.0}):
    p_times = []
    p_amps_high = []
    p_amps_f = []
    p_amps_b = []

    for p_t, p_a in zip(pulse_times, pulse_amps):
        p_times += [p_t, p_t + pulse_length]
        p_amps_f += [p_a, 0]
        p_amps_b += [0, p_a]

    print('Pulse Times: {}'.format(p_times))
    print('Pulse Amps: {}'.format(p_amps_f))

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = mig_py.MyCell(t_path, True, param_dict)

    if current_dict['amp']:
        t_curr = h.IClamp(cell.somatic[0](0.5))
        t_curr.amp = current_dict['amp']
        t_curr.delay = current_dict['start']
        t_curr.dur = current_dict['dur']
        print('Current Injection at t = {0} msec, amplitude = {1} nA'.format(current_dict['start'], current_dict['amp']))

    n_nodes = len(cell.ca[cell.cyt].nodes)
    apical_trunk_inds = [0, 8, 9, 11, 13, 19]
    sec_list = [sec for sec in cell.somatic] + [cell.apical[i] for i in apical_trunk_inds]
    apic_names = ['apical_{0}'.format(num) for num in apical_trunk_inds]
    sec_names = ['soma_0'] + apic_names

    h.distance(0, cell.somatic[0](0.5))

    node_dists = []

    for sec in sec_list:
        for seg in sec:
            node_dists.append(h.distance(sec(seg.x)))

    node_locs = []

    for sec_name in sec_names:
        sec = cell.get_sec_by_name(sec_name)

        for node in cell.ca[cell.cyt].nodes:
            if node.satisfies(sec):
                node_locs.append('{0}({1:.3})'.format(sec_name, node.x))

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)

    s_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    a0_v = h.Vector().record(cell.apical[0](0.5)._ref_v)
    a9_v = h.Vector().record(cell.apical[9](0.5)._ref_v)
    ax_v = h.Vector().record(cell.axonal[0](0.5)._ref_v)

    s_isk = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kca)
    a0_isk = h.Vector().record(cell.apical[0](0.5)._ref_ik_kca)
    a9_isk = h.Vector().record(cell.apical[9](0.5)._ref_ik_kca)

    s_im = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kmb_inh)
    ax_im = h.Vector().record(cell.axonal[0](0.5)._ref_ik_kmb_inh)

    s_ach = h.Vector().record(cell.ach.nodes(cell.somatic[0])[0]._ref_concentration)

    s_dag = h.Vector().record(cell.dag.nodes(cell.somatic[0])[0]._ref_concentration)

    s_pip2 = h.Vector().record(cell.pip2.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_pip2 = h.Vector().record(cell.pip2.nodes(cell.apical[0])[0]._ref_concentration)
    a9_pip2 = h.Vector().record(cell.pip2.nodes(cell.apical[9])[0]._ref_concentration)
    ax_pip2 = h.Vector().record(cell.pip2.nodes(cell.axonal[0])[0]._ref_concentration)
    s_pip2_b = h.Vector().record(cell.pip2_bound.nodes(cell.somatic[0])[0]._ref_concentration)
    ax_pip2_b = h.Vector().record(cell.pip2_bound.nodes(cell.axonal[0])[0]._ref_concentration)

    s_pi4p = h.Vector().record(cell.pi4p.nodes(cell.somatic[0])[0]._ref_concentration)

    s_pi = h.Vector().record(cell.pi.nodes(cell.somatic[0])[0]._ref_concentration)

    s_pip2_kcnq = h.Vector().record(cell.pip2_kcnq.nodes(cell.somatic[0])[0]._ref_concentration)
    s_perc_i = h.Vector().record(cell.somatic[0](0.5)._ref_perc_i_kmb_inh)
    ax_perc_i = h.Vector().record(cell.axonal[0](0.5)._ref_perc_i_kmb_inh)
    s_ga_gtp = h.Vector().record(cell.ga_gtp_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_plc = h.Vector().record(cell.plc_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_active_plc = h.Vector().record(cell.ga_gtp_plc_m1.nodes(cell.somatic[0])[0]._ref_concentration)

    s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[0])[0]._ref_concentration)
    a9_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)
    s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[0])[0]._ref_concentration)
    a9_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[9])[0]._ref_concentration)
    s_ip3 = h.Vector().record(cell.ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[0])[0]._ref_concentration)
    a9_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[9])[0]._ref_concentration)
    s_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    a9_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.apical[9])[0]._ref_concentration)

    # s_po = h.Vector().record(cell.ro_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    # a9_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[9])[0]._ref_concentration)

    s_po = h.Vector().record(cell.x10_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_po = h.Vector().record(cell.x10_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    a9_po = h.Vector().record(cell.x10_ip3r.nodes(cell.apical[9])[0]._ref_concentration)

    s_ip5p = h.Vector().record(cell.ip5p.nodes(cell.somatic[0])[0]._ref_concentration)
    s_ip5p_ip3 = h.Vector().record(cell.ip5p_ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    s_ip3k = h.Vector().record(cell.ip3k.nodes(cell.somatic[0])[0]._ref_concentration)
    s_ip3k_2ca = h.Vector().record(cell.ip3k_2ca_ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    s_ip3k_2ca_ip3 = h.Vector().record(cell.ip3k_2ca_ip3.nodes(cell.somatic[0])[0]._ref_concentration)

    s_r = h.Vector().record(cell.r_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_rl = h.Vector().record(cell.rl_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_rlg = h.Vector().record(cell.rlg_m1.nodes(cell.somatic[0])[0]._ref_concentration)

    s_g = h.Vector().record(cell.g_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_gbg = h.Vector().record(cell.gbg_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_rg = h.Vector().record(cell.rg_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    s_rgbg = h.Vector().record(cell.rgbg_m1.nodes(cell.somatic[0])[0]._ref_concentration)

    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    h.stdinit()

    print('Running ACh pulse')

    for t_i, t_step in enumerate(t_steps):
        h.continuerun(t_step)

        if t_step in p_times:
            p_i = p_times.index(t_step)
            print('time: {0}, [ACh]: {1} {2}M'.format(t_step, p_amps_f[p_i]*1000.0, mu))

            cell.ach.nodes(cell.cyt).concentration = p_amps_f[p_i]

            h.CVode().re_init()

    res_dict = {'node names': node_locs,
                'node distances': node_dists,
                # 'ip3 vals': ip3_vals,
                # 'cyt vals': cyt_vals,
                # 'er vals': er_vals,
                # 'ip3r open vals': ip3r_open_vals,
                'soma_v': np.array(s_v),
                'apical0_v': np.array(a0_v),
                'apical9_v': np.array(a9_v),
                'axon_v': np.array(ax_v),
                'soma_cyt': np.array(s_ca_cyt),
                'apical0_cyt': np.array(a0_ca_cyt),
                'apical9_cyt': np.array(a9_ca_cyt),
                'soma_er': np.array(s_ca_er),
                'apical0_er': np.array(a0_ca_er),
                'apical9_er': np.array(a9_ca_er),
                'soma_ip3': np.array(s_ip3),
                'apical0_ip3': np.array(a0_ip3),
                'apical9_ip3': np.array(a9_ip3),
                'soma_ip3r_ri': np.array(s_ri),
                'apical0_ip3r_ri': np.array(a0_ri),
                'apical9_ip3r_ri': np.array(a9_ri),
                'soma_ip3r_open_probability': np.array(s_po),
                'apical0_ip3r_open_probability': np.array(a0_po),
                'apical9_ip3r_open_probability': np.array(a9_po),
                'soma_im': np.array(s_im),
                'axon_im': np.array(ax_im),
                'soma_isk': np.array(s_isk),
                'apical0_isk': np.array(a0_isk),
                'apical9_isk': np.array(a9_isk),
                'soma_dag': np.array(s_dag),
                'soma_pip2': np.array(s_pip2),
                'soma_pi': np.array(s_pi),
                'soma_pi4p': np.array(s_pi4p),
                'apical0_pip2': np.array(a0_pip2),
                'apical9_pip2': np.array(a9_pip2),
                'axon_pip2': np.array(ax_pip2),
                'axon_pip2_bound': np.array(ax_pip2_b),
                'soma_ach': np.array(s_ach),
                'soma_pip2_bound': np.array(s_pip2_b),
                'soma_pip2_kcnq': np.array(s_pip2_kcnq),
                'soma_perc_i': np.array(s_perc_i),
                'axon_perc_i': np.array(ax_perc_i),
                'soma_ga_gtp': np.array(s_ga_gtp),
                'soma_plc': np.array(s_plc),
                'soma_active_plc': np.array(s_active_plc),
                'soma_r': np.array(s_r),
                'soma_g': np.array(s_g),
                'soma_gbg': np.array(s_gbg),
                'soma_rg': np.array(s_rg),
                'soma_rgbg': np.array(s_rgbg),
                'soma_rl': np.array(s_rl),
                'soma_rlg': np.array(s_rlg),
                'soma_ip5p': np.array(s_ip5p),
                'soma_ip5p_ip3': np.array(s_ip5p_ip3),
                'soma_ip3k': np.array(s_ip3k),
                'soma_ip3k_2ca': np.array(s_ip3k_2ca),
                'soma_ip3k_2ca_ip3': np.array(s_ip3k_2ca_ip3),
                # 'soma dumb1': np.array(s_d1),
                # 'soma dumb2': np.array(s_d2),
                't': np.array(t_vec),
                }
    return res_dict


def plot_ach_pulse(result_dict, t_is=[0, -1], ach_times=[100, 150], t_ignore=0):
    ret_figs = []

    cyt_t_fig = bplt.figure(title='Cytosol Calcium vs Time')
    cyt_t_fig.xaxis.axis_label = 'time (msec)'
    cyt_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(cyt_t_fig)
    er_t_fig = bplt.figure(title='ER Calcium vs Time')
    er_t_fig.xaxis.axis_label = 'time (msec)'
    er_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(er_t_fig)
    ach_t_fig = bplt.figure(title='ACh vs Time')
    ach_t_fig.xaxis.axis_label = 'time (msec)'
    ach_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ach_t_fig)
    m1_t_fig = bplt.figure(title='Soma M1 States vs Time')
    m1_t_fig.xaxis.axis_label = 'time (msec)'
    m1_t_fig.yaxis.axis_label = 'density ({}m^2)'.format(mu)
    ret_figs.append(m1_t_fig)
    act_plc_t_fig = bplt.figure(title='Activated PLC vs Time')
    act_plc_t_fig.xaxis.axis_label = 'time (msec)'
    act_plc_t_fig.yaxis.axis_label = 'density ({}m^-2)'.format(mu)
    ret_figs.append(act_plc_t_fig)
    pip2_t_fig = bplt.figure(title='PIP2 vs Time')
    pip2_t_fig.xaxis.axis_label = 'time (msec)'
    pip2_t_fig.yaxis.axis_label = 'density ({}m^-2)'.format(mu)
    ret_figs.append(pip2_t_fig)
    pi_flux_fig = bplt.figure(title='PI fluxes vs Time')
    pi_flux_fig.xaxis.axis_label = 'time (msec)'
    pi_flux_fig.yaxis.axis_label = 'flux rate (um^-2 msec^-1)'
    ret_figs.append(pi_flux_fig)
    ip3_t_fig = bplt.figure(title='IP3 vs Time')
    ip3_t_fig.xaxis.axis_label = 'time (msec)'
    ip3_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ip3_t_fig)
    ip3r_open_t_fig = bplt.figure(title='Open IP3R vs Time')
    ip3r_open_t_fig.xaxis.axis_label = 'time (msec)'
    ip3r_open_t_fig.yaxis.axis_label = 'probability (%)'
    ret_figs.append(ip3r_open_t_fig)
    ip3_kinase_fig = bplt.figure(title='IP3 Breakdown vs Time')
    ip3_kinase_fig.xaxis.axis_label = 'time (msec)'
    ip3_kinase_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ip3_kinase_fig)
    pip2_kcnq_t_fig = bplt.figure(title='Percent I_M vs Time')
    pip2_kcnq_t_fig.xaxis.axis_label = 'time (msec)'
    pip2_kcnq_t_fig.yaxis.axis_label = 'percentage (%)'.format(mu)
    ret_figs.append(pip2_kcnq_t_fig)
    ik_t_fig = bplt.figure(title='Potassium Currents vs Time')
    ik_t_fig.xaxis.axis_label = 'time (msec)'
    ik_t_fig.yaxis.axis_label = 'current (mA/cm^2)'.format(mu)
    ret_figs.append(ik_t_fig)
    v_fig = bplt.figure(title='Somatic Membrane Potential vs Time')
    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(v_fig)

    sp_times = fan.get_spike_times(result_dict['soma_v'], result_dict['t'])

    t_bool = result_dict['t'] > t_ignore
    t_ig = result_dict['t'][t_bool] - t_ignore

    if sp_times.size > 3:

        isr_fig = bplt.figure(title='Instantaneous Spiking Rate vs Time', x_range=[t_ig[0], t_ig[-1]])
        isr_fig.xaxis.axis_label = 'time (msec)'
        isr_fig.yaxis.axis_label = 'ISR (Hz)'
        ret_figs.append(isr_fig)

        isr = np.diff(sp_times)
        isr_vals = np.divide(1000.0, isr)
        isr_ts = (sp_times[:-1] + sp_times[1:])/2.0

        isr_is = np.where(isr_ts > t_ignore)

        isr_fig.circle(isr_ts[isr_is]-t_ignore, isr_vals[isr_is], size=12, color=colrs[0], legend='Simulation Result')
        isr_fig.line(isr_ts[isr_is]-t_ignore, isr_vals[isr_is], line_width=3, color=colrs[0])

        label_source = bmod.ColumnDataSource(dict(x=[ach_times[0] + 2000], y=[8],
                                                  text=['cessation = {0:.2f} (msec)'.format(np.max(isr))]))
        reg_label = bmod.glyphs.Text(x='x', y='y', text='text', text_color='black')
        isr_fig.add_glyph(label_source, reg_label)

        acc_fig = bplt.figure(title='Spike Acceleration vs Time', x_range=[t_ig[0], t_ig[-1]])
        acc_fig.xaxis.axis_label = 'time (msec)'
        acc_fig.yaxis.axis_label = 'spike acceleration (%)'
        ret_figs.append(acc_fig)

        pre_isr = isr_vals[np.intersect1d(np.where(isr_ts < ach_times[0]), np.where(isr_ts > t_ignore))]
        pre_avg_isr = np.mean(pre_isr)
        post_isr = isr_vals[np.where(isr_ts > ach_times[0])]

        acc_vals = (post_isr - pre_avg_isr)/pre_avg_isr*100.0
        acc_ts = isr_ts[np.where(isr_ts > ach_times[0])]

        peak_i = np.argmax(acc_vals)
        print('Peak Spike Acceleration = {0}'.format(acc_vals[peak_i]))
        acc_slope, acc_int, r_val, p_val, acc_err = scat.linregress(acc_ts[peak_i:], acc_vals[peak_i:])
        est_vals = acc_ts * acc_slope + acc_int

        acc_fig.circle(acc_ts - t_ignore, acc_vals, size=12, color=colrs[0], legend='Simulation Result')
        acc_fig.line(acc_ts - t_ignore, acc_vals, line_width=3, color=colrs[0])

        acc_fig.line(acc_ts - t_ignore, est_vals, line_width=3, color='green', legend='Linear Regression')

        lable_source = bmod.ColumnDataSource(dict(x=[ach_times[0]+3000], y=[0],
                                                text=[
                                                    'slope = {0:.2f} (%/sec)\n r = {1:.2f}'.format(acc_slope*1000.0, r_val)]))
        reg_label = bmod.glyphs.Text(x='x', y='y', text='text', text_color='green')
        acc_fig.add_glyph(lable_source, reg_label)

    # dumb_fig = bplt.figure(title='Dummy Variables to Test Density/Concentration')
    # dumb_fig.xaxis.axis_label = 'time ('msec)'
    # ret_figs.append(dumb_fig)

    ach_spans = []

    for ach_time in ach_times:
        ach_spans.append(bmod.Span(location=ach_time-t_ignore, dimension='height',
                                   line_color='green', line_dash='dashed', line_width=3)
                         )

    cyt_t_fig.line(t_ig, result_dict['apical9_cyt'][t_bool] * 1000.0, line_width=3, color=colrs[2],
                   legend='apical[9]')
    cyt_t_fig.line(t_ig, result_dict['apical0_cyt'][t_bool] * 1000.0, line_width=3, color=colrs[0],
                   legend='apical[0]')
    cyt_t_fig.line(t_ig, result_dict['soma_cyt'][t_bool] * 1000.0, line_width=3, color='black')

    for ach_sp in ach_spans:
        cyt_t_fig.add_layout(ach_sp)

    er_t_fig.line(t_ig, result_dict['soma_er'][t_bool] * 1000.0, line_width=3, color='black', legend='somatic[0]')
    er_t_fig.line(t_ig, result_dict['apical0_er'][t_bool] * 1000.0, line_width=3, color=colrs[0], legend='apical[0]')
    er_t_fig.line(t_ig, result_dict['apical9_er'][t_bool] * 1000.0, line_width=3, color=colrs[2], legend='apical[9]')

    for ach_sp in ach_spans:
        er_t_fig.add_layout(ach_sp)

    ach_t_fig.line(t_ig, result_dict['soma_ach'][t_bool]*1000.0, line_width=3, color='green', legend='somatic[0](0.5)')

    m1_t_fig.line(t_ig, result_dict['soma_r'][t_bool], line_width=3, color=colrs[0], legend='R')
    m1_t_fig.line(t_ig, result_dict['soma_g'][t_bool], line_width=3, color=colrs[2], legend='G')
    m1_t_fig.line(t_ig, result_dict['soma_gbg'][t_bool], line_width=3, color=colrs[3], legend='G-bg')
    m1_t_fig.line(t_ig, result_dict['soma_rg'][t_bool], line_width=3, color=colrs[4], legend='RG')
    m1_t_fig.line(t_ig, result_dict['soma_rgbg'][t_bool], line_width=3, color=colrs[5], legend='RG-bg')
    m1_t_fig.line(t_ig, result_dict['soma_rl'][t_bool], line_width=3, color=colrs[6], legend='RL')
    m1_t_fig.line(t_ig, result_dict['soma_rlg'][t_bool], line_width=3, color=colrs[7], legend='RLG')

    for ach_sp in ach_spans:
        m1_t_fig.add_layout(ach_sp)

    act_plc_t_fig.line(t_ig, result_dict['soma_rlg'][t_bool], line_width=3, color=colrs[4], legend='RLG')
    act_plc_t_fig.line(t_ig, result_dict['soma_ga_gtp'][t_bool], line_width=3, color=colrs[0],
                       legend='ga-gtp')
    # act_plc_t_fig.line(t_ig, result_d['soma plc'], line_width=3, color=colrs[1], legend='plc')
    act_plc_t_fig.line(t_ig, result_dict['soma_active_plc'][t_bool], line_width=3, color=colrs[2], legend='ga-gtp-plc')

    for ach_sp in ach_spans:
        act_plc_t_fig.add_layout(ach_sp)

    pip2_t_fig.line(t_ig, result_dict['soma_pip2'][t_bool], line_width=3, color=colrs[0], legend='Somatic PIP2')
    pip2_t_fig.line(t_ig, result_dict['soma_pip2_bound'][t_bool], line_width=3, color=colrs[5], legend='Somatic PIP2 bound')
    # pip2_t_fig.line(t_ig, result_dict['apical0_pip2'], line_width=3, color=colrs[1], legend='Proximal Apical PIP2')
    # pip2_t_fig.line(t_ig, result_dict['apical9_pip2'], line_width=3, color=colrs[2], legend='Distal Apical PIP2')
    # pip2_t_fig.line(t_ig, result_dict['axon_pip2'], line_width=3, color=colrs[3], legend='Axonal PIP2',
    #                 line_dash='dashed')
    # pip2_t_fig.line(t_ig, result_dict['axon_pip2_bound'], line_width=3, color=colrs[6], legend='Axonal PIP2 bound',
    #                 line_dash='dashed')
    pip2_t_fig.line(t_ig, result_dict['soma_pi4p'][t_bool], line_width=3, color=colrs[4], legend='somatic[0](0.5) PI4P')

    for ach_sp in ach_spans:
        pip2_t_fig.add_layout(ach_sp)

    k_PLC = 0.3  # 0.0003  # (um^2 msec^-1) dropped by a fact
    k_4K = 0.0008 * 2  # 0.0000008  # (msec^-1)
    k_4P = 0.12  # 0.00012  # (msec^-1)
    k_5K = 0.02  # 0.00002  # (msec^-1)
    k_5P = 0.028  # 0.000028  # (msec^-1)

    pi_flux_fig.line(t_ig, result_dict['soma_pi'][t_bool]*k_4K,line_width=3, color=colrs[0], legend='somatic[0](0.5) PI->PI4P')
    pi_flux_fig.line(t_ig, result_dict['soma_pi4p'][t_bool] * k_4P, line_width=3, color=colrs[1],
                     legend='somatic[0](0.5) PI4P->PI')
    pi_flux_fig.line(t_ig, result_dict['soma_pi4p'][t_bool] * k_5K, line_width=3, color=colrs[2],
                     legend='somatic[0](0.5) PI4P->PIP2')
    pi_flux_fig.line(t_ig, result_dict['soma_pip2'][t_bool] * k_5P, line_width=3, color=colrs[3],
                     legend='somatic[0](0.5) PIP2->PI4P')
    pi_flux_fig.line(t_ig, result_dict['soma_active_plc'][t_bool]*result_dict['soma_pip2'][t_bool]*k_PLC*3, line_width=3,
                     color=colrs[4], legend='somatic[0](0.5) PIP2->DAG+IP3')

    # ip3_t_fig.line(t_ig, dip3, line_width=3, color=colrs[0], legend='somatic[0] flux')
    ip3_t_fig.line(t_ig, result_dict['soma_ip3'][t_bool] * 1000.0, line_width=3, color='black', legend='somatic[0]')
    ip3_t_fig.line(t_ig, result_dict['apical0_ip3'][t_bool] * 1000.0, line_width=3, color=colrs[0],
                   legend='apical[0]')
    ip3_t_fig.line(t_ig, result_dict['apical9_ip3'][t_bool] * 1000.0, line_width=3, color=colrs[2],
                   legend='apical[9]')
    for ach_sp in ach_spans:
        ip3_t_fig.add_layout(ach_sp)

    ip3r_open_t_fig.line(t_ig, 100*result_dict['soma_ip3r_open_probability'][t_bool], line_width=3, color='black',
                         legend='somatic[0]')
    # ip3r_open_t_fig.line(t_ig, 100*result_dict['soma_ip3r_ri'], line_width=3, color='black',
    #                      legend='somatic[0] ri', line_dash='dashed')
    ip3r_open_t_fig.line(t_ig, 100*result_dict['apical0_ip3r_open_probability'][t_bool], line_width=3, color=colrs[0],
                         legend='apical[0]')
    # ip3r_open_t_fig.line(t_ig, 100 * result_dict['apical0_ip3r_ri'], line_width=3, color=colrs[0],
    #                      legend='apical[0] ri', line_dash='dashed')
    ip3r_open_t_fig.line(t_ig, 100*result_dict['apical9_ip3r_open_probability'][t_bool], line_width=3, color=colrs[1],
                         legend='apical[9]')
    # ip3r_open_t_fig.line(t_ig, 100 * result_dict['apical9_ip3r_ri'], line_width=3, color=colrs[1],
    #                      legend='apical[9] ri', line_dash='dashed')
    for ach_sp in ach_spans:
        ip3r_open_t_fig.add_layout(ach_sp)

    ip3_kinase_fig.line(t_ig, 1000.0 * result_dict['soma_ip5p'][t_bool],
                        line_width=3, color=colrs[0], legend='ip5p')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_dict['soma_ip5p_ip3'][t_bool],
                        line_width=3, line_dash='dashed', color=colrs[1], legend='ip5p_ip3')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_dict['soma_ip3k'][t_bool],
                        line_width=3, color=colrs[2], legend='ip3k')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_dict['soma_ip3k_2ca'][t_bool],
                        line_width=3, line_dash='dashed',color=colrs[3], legend='ip3k_2ca')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_dict['soma_ip3k_2ca_ip3'][t_bool],
                        line_width=3, line_dash='dotted', color=colrs[4], legend='ip3k_2ca_ip3')

    pip2_kcnq_t_fig.line(t_ig, 100.0*result_dict['soma_perc_i'][t_bool], line_width=3, color=colrs[0], legend='somatic[0](0.5)')
    pip2_kcnq_t_fig.line(t_ig, 100.0 * result_dict['axon_perc_i'][t_bool], line_width=3, color=colrs[2],
                         legend='axon[0](0.5)')
    for ach_sp in ach_spans:
        pip2_kcnq_t_fig.add_layout(ach_sp)

    ik_t_fig.line(t_ig, result_dict['soma_im'][t_bool], line_width=3, color='black', legend='Soma Im', line_dash='dashed')
    # ik_t_fig.line(t_ig, result_d['axon im'][i_ig], line_width=3, color=colrs[2], legend='axonal[0]')

    ik_t_fig.line(t_ig, result_dict['apical9_isk'][t_bool], line_width=3, color=colrs[2],
                   legend='apical[9] Isk')
    ik_t_fig.line(t_ig, result_dict['apical0_isk'][t_bool], line_width=3, color=colrs[0],
                   legend='apical[0] Isk')
    ik_t_fig.line(t_ig, result_dict['soma_isk'][t_bool], line_width=3, color='black', legend='Soma Isk')
    for ach_sp in ach_spans:
        ik_t_fig.add_layout(ach_sp)

    v_fig.line(t_ig, result_dict['soma_v'][t_bool], line_width=3, color='black')
    # v_fig.line(t_ig, result_d['apical0 v'][i_ig], line_width=3, color=colrs[0], legend='apical[0](0.5)',
    #            line_dash='dashed')
    # v_fig.line(t_ig, result_d['apical9 v'][i_ig], line_width=3, color=colrs[2], legend='apical[9](0.5)',
    #            line_dash='dotted')
    # v_fig.line(t_ig, result_d['axon v'][i_ig], line_width=3, color=colrs[3], legend='axonal[0](0.5)',
    #            line_dash='dotted')
    for ach_sp in ach_spans:
        v_fig.add_layout(ach_sp)
    v_fig.legend.location = 'bottom_right'

    # dumb_fig.line(result_d['t'], result_d['soma dumb1'], line_width=3, color=colrs[0], legend='density')
    # dumb_fig.line(result_d['t'], result_d['soma dumb2'], line_width=3, color=colrs[2], legend='concentration')

    # cal_fig = bplt.figure(title='Cytosol Calcium vs Location')
    # cal_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(cal_fig)
    # er_fig = bplt.figure(title='ER Calcium vs Location')
    # er_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(er_fig)
    # ip3_fig = bplt.figure(title='IP3 vs Location')
    # ip3_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(ip3_fig)
    # ip3r_open_fig = bplt.figure(title='Open IP3R vs Location')
    # ip3r_open_fig.xaxis.axis_label = 'distance from soma ({}m)'.format(mu)
    # ret_figs.append(ip3r_open_fig)
    #
    # for s_i, t_i in enumerate(t_is):
    #     cal_fig.circle(result_d['node distances'], 1000.0 * result_d['cyt vals'][t_i, :], size=20 - s_i,
    #                    color=colrs[s_i],
    #                    legend='t={0} msec'.format(my_steps[t_i]))
    #     cal_fig.line(result_d['node distances'], 1000.0 * result_d['cyt vals'][t_i, :], line_width=3,
    #                  color=colrs[s_i])
    #     er_fig.circle(result_d['node distances'], 1000.0 * result_d['er vals'][t_i, :], size=20 - s_i, color=colrs[s_i],
    #                   legend='t={0} msec'.format(my_steps[t_i]))
    #     er_fig.line(result_d['node distances'], 1000.0 * result_d['er vals'][t_i, :], line_width=3,
    #                 color=colrs[s_i])
    #     ip3_fig.circle(result_d['node distances'], 1000.0 * result_d['ip3 vals'][t_i, :], size=20 - s_i,
    #                    color=colrs[s_i],
    #                    legend='t={0} msec'.format(my_steps[t_i]))
    #     ip3_fig.line(result_d['node distances'], 1000.0 * result_d['ip3 vals'][t_i, :], line_width=3,
    #                  color=colrs[s_i])
    #     ip3r_open_fig.circle(result_d['node distances'], 1e6 * result_d['ip3r open vals'][t_i, :], size=20 - s_i,
    #                          color=colrs[s_i], legend='t={0} msec'.format(my_steps[t_i]))
    #     ip3r_open_fig.line(result_d['node distances'], 1e6 * result_d['ip3r open vals'][t_i, :],
    #                        line_width=3, color=colrs[s_i])

    return ret_figs


def add_dasari_isr_data(in_fig):
    das_spikes = np.array([0.0423397381, 0.3421186873, 0.499234824, 0.684917531, 0.8563169529, 1.0705662302,
                           1.2275123278, 1.3989117497, 1.6274443122, 1.8700901207, 2.0130930114, 2.0695459956,
                           2.9981295698, 3.0693759565, 3.1265090971, 3.1836422377, 3.2550586635, 3.297908519,
                           3.3693249447, 3.4121748002, 3.4834211869, 3.5259309641, 3.5971773508, 3.6687638157,
                           3.7258969563, 3.7680666553, 3.825199796, 3.8966162217, 3.9537493624, 4.010882503,
                           4.0678456045, 4.1533752763, 4.2394150655, 4.2962081279, 4.396191124, 4.4533242646,
                           4.5104574052, 4.581873831, 4.6675735419, 4.7389899677, 4.8245196395, 4.881482741,
                           4.9530692059, 5.0244856317, 5.0957320184, 5.1385818738, 5.2099982996, 5.2671314402,
                           5.3528311512, 5.4387009012, 5.524230573, 5.5956469988, 5.7240265261, 5.809726237,
                           5.895425948, 5.9668423737, 6.0525420847, 6.1523550417, 6.2525080769, 6.3380377487,
                           6.4381907839, 6.5521169869, 6.6379867369, 6.7665363033, 6.8522360143, 6.9522190104,
                           7.0663152525, 7.1807515729, 7.294847815, 7.4376806666, 7.5517769087, 7.6374766196,
                           7.7517429009, 7.8660091821, 7.9659921782, 8.0945417446, 8.2229212719, 8.3657541234,
                           8.4800204047, 8.6087400102, 8.7087230063, 8.8661792212, 8.9806155416, 9.0948818228,
                           9.2091481041, 9.308961061, 9.4517939126, 9.566230233, 9.7233463697, 9.8949158306,
                           9.9948988267, 10.1379017174, 10.3093011393])

    das_t = das_spikes*1000.0
    das_isr = np.divide(1000.0, np.diff(das_t))
    das_ints = (das_t[:-1] + das_t[1:])/2.0

    in_fig.circle(das_ints, das_isr, size=12, color=colrs[2], legend='Dasari Data')
    in_fig.line(das_ints, das_isr, line_width=3, color=colrs[2])


def plot_dasari_and_gulledge_data():
    dg_figs = []

    mV_span = bmod.Span(location=-70.0, dimension='width', line_color='grey',
                        line_dash='dashed', line_width=3)

    dg_fig1a2_rmp = np.array([[-3.69, 5.16],[-3.55, 4.82],[-3.39, 5.64],[-3.35, 4.71],[-3.16, 5.04],[-3.12, 5.62],
                              [-3.06, 4.69],[-3.01, 5.16],[-2.97, 4.98],[-2.91, 5.53],[-2.85, 5.37],[-2.79, 5.60],
                              [-2.76, 5.23],[-2.67, 5.99],[-2.59, 8.63],[-2.52, 6.93],[-2.43, 5.90],[-2.39, 4.77],
                              [-2.36, 3.87],[-2.31, 2.89],[-2.28, 1.84],[-2.25, 1.34],[-2.20, 1.43],[-2.13, 1.30],
                              [-2.07, 1.28],[-2.00, 1.06],[-1.95, 1.96],[-1.88, 1.82],[-1.79, 2.35],[-1.71, 2.82],
                              [-1.63, 3.01],[-1.59, 3.48],[-1.57, 3.25],[-1.53, 3.97],[-1.48, 4.18],[-1.43, 4.51],
                              [-1.39, 5.00],[-1.34, 5.47],[-1.31, 5.92],[-1.28, 5.37],[-1.25, 6.25],[-1.19, 6.42],
                              [-1.14, 6.60],[-1.12, 6.21],[-1.07, 7.67],[-1.04, 7.46],[-0.987, 8.02],[-0.943, 8.26],
                              [-0.911, 8.06],[-0.871, 8.31],[-0.813, 7.98],[-0.776, 7.24],[-0.723, 7.92],[-0.665, 8.57],
                              [-0.598, 8.16],[-0.546, 7.73],[-0.477, 7.63],[-0.445, 8.22],[-0.378, 8.31],[-0.318, 7.14],
                              [-0.256, 7.90],[-0.189, 8.39],[-0.126, 8.24],[-0.0856, 7.63],[0.00142, 7.69],
                              [0.0826, 7.77],[0.135, 7.83],[0.170, 8.45],[0.210, 7.67],[0.282, 7.24],[0.353, 8.31],
                              [0.460, 7.94],[0.492, 8.78],[0.555, 7.98],[0.608, 8.53],[0.662, 7.90],[0.749, 7.65],
                              [0.790, 8.26],[0.845, 7.71],[0.912, 7.69],[0.958, 8.43],[1.01, 8.12],[1.07, 7.36],
                              [1.11, 8.06],[1.20, 7.48],[1.27, 8.26],[1.34, 7.75],[1.40, 7.73],[1.44, 7.24],
                              [1.51, 7.79],[1.58, 7.98],[1.62, 8.24],[1.67, 8.02],[1.73, 7.85],[1.79, 7.69],
                              [1.83, 7.67]])

    dg_fig1a2_t = (dg_fig1a2_rmp[:, 0] - dg_fig1a2_rmp[0, 0]) * 1000.0
    dg_fig1a2_v = dg_fig1a2_rmp[:, 1] - dg_fig1a2_rmp[0, 1] - 70.0

    dg_fig1a2_ach = (np.array([-2.70, -2.66]) - dg_fig1a2_rmp[0,0])*1000.0

    ach_span_1 = bmod.Span(location=dg_fig1a2_ach[0], dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_2 = bmod.Span(location=dg_fig1a2_ach[1], dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    dg_figs.append(bplt.figure(title='Dasari and Gulledge Fig 1A2'))
    dg_figs[-1].line(dg_fig1a2_t, dg_fig1a2_v, line_width=3, color='blue', legend='soma')
    dg_figs[-1].xaxis.axis_label = 'time (msec)'
    dg_figs[-1].yaxis.axis_label = 'membrane potential (mV)'
    dg_figs[-1].add_layout(mV_span)
    dg_figs[-1].add_layout(ach_span_1)
    dg_figs[-1].add_layout(ach_span_2)

    dg_fig2a_rmp = np.array([[-5.43, 1.73],[-5.33, 1.73],[-5.24, 1.60],[-5.09, 1.93],[-4.95, 1.60],[-4.82, 1.93],[-4.80, 2.40],
                    [-4.65, 1.80],[-4.49, 2.00],[-4.46, 4.07],[-4.41, 1.73],[-4.28, 3.00],[-4.21, 2.20],[-4.00, 2.27],
                    [-3.92, 1.73],[-3.84, 2.00],[-3.82, 2.80],[-3.63, -1.47],[-3.53, -2.13],[-3.40, -2.20],
                    [-3.37, -1.40],[-3.33, -2.40],[-3.26, -2.27],[-3.21, -2.40],[-3.04, -2.07],[-2.88, -2.07],
                    [-2.85, -1.40],[-2.81, -1.87],[-2.65, -1.53],[-2.54, -1.27],[-2.42, -1.13],[-2.29, -0.600],
                    [-2.20, -0.0667],[-2.08, 0.600],[-1.97, 1.27],[-1.85, 2.00],[-1.72, 2.93],[-1.62, 3.13],
                    [-1.51, 3.27],[-1.45, 3.60],[-1.29, 3.53],[-1.12, 3.33],[-1.11, 4.73],[-1.02, 3.40],[-0.789, 3.60],
                    [-0.716, 3.80],[-0.600, 3.80],[-0.326, 3.93],[-0.221, 3.73],[-0.0632, 3.47],[0.0316, 3.87],
                    [0.0737, 4.73],[0.105, 4.00],[0.189, 3.53],[0.295, 3.87],[0.379, 4.20],[0.400, 5.00],[0.453, 4.20],
                    [0.547, 4.13]])
    dg_fig2a_t = (dg_fig2a_rmp[:,0] - dg_fig2a_rmp[0,0])*1000.0
    dg_fig2a_v = dg_fig2a_rmp[:,1] - dg_fig2a_rmp[0, 1] - 70.0

    dg_fig2a_ach = (np.array([-4.35, -4.31]) - dg_fig2a_rmp[0,0])*1000.0

    ach_span_3 = bmod.Span(location=dg_fig2a_ach[0], dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_4 = bmod.Span(location=dg_fig2a_ach[1], dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    dg_figs.append(bplt.figure(title='Dasari and Gulledge Fig 2A '))
    dg_figs[-1].line(dg_fig2a_t, dg_fig2a_v, line_width=3, color='blue', legend='soma')
    dg_figs[-1].xaxis.axis_label='time (msec)'
    dg_figs[-1].yaxis.axis_label='membrane potential (mV)'
    dg_figs[-1].add_layout(ach_span_3)
    dg_figs[-1].add_layout(ach_span_4)
    dg_figs[-1].add_layout(mV_span)

    return dg_figs


def ode_m1(species, t, k_f, k_r):
    """
    species = ['R', 'L', 'RL']
    """

    r1 = species[0]*species[1] * k_f - species[2] * k_r

    return [-r1, 0, r1]


def sig_func(x, k, x0):
    return 1.0/(1.0 + np.exp(-k*(x - x0)))


def run_m1_activation_oxom_vs_ach(conc_vals, t_sim_end=5000):
    """
    conc_vals: numpy array with the concentrations of agonists to be simulated
    t_sim_end: length of simulations used to calculate RL values

    :return: r_dict: dictionary of results for simulations
        r_dict['ACh'][concentration index] =

    """

    # Run with agonist at some concentration for some amount of time.
    # Check to see if it's at a steady state. If so use that as value for plotting concentration
    # response.

    r_dict = {}

    init_vals = [15.87, 0.0, 0.0]
    t_vals = np.linspace(0, t_sim_end, 100)  # Simulation time steps in msec

    ach_vals = np.zeros(conc_vals.shape)
    oxom_vals = np.zeros(conc_vals.shape)

    for ag in ['ACh', 'Oxom']:
        print(ag)
        if ag == 'ACh':
            k_f = 0.00278  # (uM^-1 msec^-1)
            k_r = 0.00215  # (msec^-1)

        else:
            k_f = 0.00278  # (uM^-1 msec^-1)
            k_r = 0.0055  # (msec^-1)

        r_dict[ag] = {}

        for c_i, c_val in enumerate(conc_vals):
            print('{0}, {1} {2}M'.format(c_i, c_val, mu))
            init_vals[1] = c_val

            r_dict[ag][c_i] = spint.odeint(ode_m1, init_vals, t_vals, args=(k_f, k_r))

            if ag == 'ACh':
                ach_vals[c_i] = r_dict[ag][c_i][-1, 2]
            else:
                oxom_vals[c_i] = r_dict[ag][c_i][-1, 2]

    r_dict['concentrations'] = conc_vals
    r_dict['t'] = t_vals
    r_dict['ACh final values'] = ach_vals
    r_dict['Oxom final values'] = oxom_vals
    return r_dict


def plot_m1_activation_oxom_vs_ach(res_dict):

    conc_vals = res_dict['concentrations']
    ach_vals = res_dict['ACh final values']
    ach_norm = ach_vals/np.max(ach_vals)
    oxom_vals = res_dict['Oxom final values']
    oxom_norm = oxom_vals/np.max(oxom_vals)

    # Estimate sigmoidal curves to fit response

    # Calculate estimated curve values
    l_concs = np.logspace(np.log10(conc_vals[0]), np.log10(conc_vals[-1]), 39)

    conc_log = np.log10(conc_vals)
    ach_opt, _ = spopt.curve_fit(sig_func, conc_log, ach_norm)
    print('ACh\nk = {0}\nE50 = {1}'.format(ach_opt[0], ach_opt[1]))
    oxom_opt, _ = spopt.curve_fit(sig_func, conc_log, oxom_norm)
    print('Oxom\nk = {0}\nE50 = {1}'.format(oxom_opt[0], oxom_opt[1]))

    ach_ec50 = 10**ach_opt[1]
    oxom_ec50 = 10**oxom_opt[1]

    # Plot values vs time
    ach_est = sig_func(np.log10(l_concs), ach_opt[0], ach_opt[1])
    oxom_est = sig_func(np.log10(l_concs), oxom_opt[0], oxom_opt[1])

    figs = []

    for ag in ['ACh', 'Oxom']:
        figs.append(bplt.figure(title='RL versus time, {}'.format(ag)))

        for c_i, c_val in enumerate(conc_vals):
            figs[-1].line(res_dict['t'], res_dict[ag][c_i][:, 2], line_width=3, color=colrs[c_i],
                          legend='{0}, {1} {2}M'.format(ag, c_val, mu))

    # Plot concentration response along with estimated curves

    conc_fig = bplt.figure(title='M1 Activation (RL) vs Agonist', x_axis_type='log')
    figs.append(conc_fig)
    conc_fig.xaxis.axis_label = 'Agonist Concentration ({}M)'.format(mu)
    conc_fig.yaxis.axis_label = 'Percent Bound (%)'

    conc_fig.circle(conc_vals, oxom_norm*100.0, size=12, color='blue', legend='Oxo-M Parameters')
    conc_fig.circle(conc_vals, ach_norm*100.0, size=12, color='green', legend='ACh Parameters')
    conc_fig.legend.location = 'bottom_right'

    conc_fig.line(l_concs, oxom_est*100.0, line_width=3, color='blue', legend='Oxo-M Fit, EC50={0:.3f}'.format(oxom_ec50))
    conc_fig.line(l_concs, ach_est*100.0, line_width=3, color='green', legend='ACh Fit, EC50={0:.3f}'.format(ach_ec50))
    return figs


def plot_model_vs_gulledge(in_trace, bio_trace):

    in_trace = 0

    v_fig = bplt.figure()
    acc_fig = bplt.figure()

    return [v_fig, acc_fig]


if __name__ == "__main__":

    # Ensure directory exists for saving figures to

    all_figs = []
    all_names = []

    ext_to_er_ratio = 1.0

    ncx_tot = 0.0
    g_ext_leak = 3.2e-5
    pmca_tot = 2.56e-5

    # Increasing er_ratio increases input/output so higher values leads to faster recharge of er
    g_er_leak = 2100.0
    g_serca = 63525.0

    print('g_ext_leak: {0}'.format(g_ext_leak))
    print('pmca_tot: {0}'.format(pmca_tot))
    print('ncx tot: {0}'.format(ncx_tot))
    print('g_er_leak: {0}'.format(g_er_leak))
    print('g_serca: {0}'.format(g_serca))

    # Simulation Options
    # 'recharge_er':  baseline recharge endoplasmic reticulum
    # 'single_ap': test calcium transient after a single action potential
    # 'four_aps': test calcium transient after four action potentials
    # 'ach_pulse': test entire reaction after a pulse of acetylcholine
    # 'test_refill_spiking': test repeated challenges with ACh while current injection induces spikes to see if recharge
    #   is sufficient to cause another inhibition
    # 'test_refill_resting': test repeated challenges with ACh while at resting membrane potential to see if recharge is
    #   sufficient to cause another inhibition

    run_sim = 'test_refill_resting'
    save_result = False
    use_file_to_plot = False
    save_figs = True

    plot_dasari = False
    rxd_bool = True

    sim_dict = {'frac_cyt': 0.9,
                'frac_er': 0.1,
                'ca ext val': 2.0,
                'g_serca': g_serca,
                'g ext leak': g_ext_leak,
                'ca_cyt_val': 0.1*0.001,
                'ca_er_val': 175.0*0.001,
                'apical calcium mult': 300.0,
                'soma calcium mult': 360.0,
                'soma kca mult': 5.5, # 4.5,
                'apical kca mult': 5.5, # 4.5,
                'soma im mult': 1.0,
                'axon im mult': 1.0,
                'ca_diff': 0.03,
                'ip3_diff': 0.3,
                'ip3_init': 0.000003,
                'ip3k_total': 0.05,  # 0.007,
                'ip5p_total': 0.01,  # 0.0002,
                'ip3r count': 1e-3,
                'cbd_total': 0.045,
                'cbdh_KD': 0.000237,
                'cbdh_kon': 11.0,
                'cbdl_KD': 0.000411,
                'cbdl_kon': 87.0,
                'car_total': 86.0,
                'car_KD': 2.0,
                'car_kon': 0.01,
                'dye identity': 'ogb1',
                'ogb1_total': 0.0,
                'ogb1_KD': 0.00043,  # 0.0004
                'ogb1_kon': 100.0,  # 10.0
                'ogb5_total': 0.0,
                'ogb5_KD': 0.01,
                'ogb5_kon': 10.0,
                'g_ip3r': 7.5e5,  # 0.9e6,  # 1e5,
                'g_er_leak': g_er_leak,  # 27.5, # 2800.0,  # 1.8
                'total pmca': pmca_tot,
                'total ncx': ncx_tot,
                'ip3 lock': False,
                'tau_leak_soma': 1.0,
                'tau_leak_apical': 1.0,
                'm1 init': 15.87*5.0,
                'g init': 40.0,
                'plc init': 3.1198*5.0,
                'pip2 init': 3232,  # 3232,
                'pi4p init': 4540,  # 4540,
                'pi init': 226975.0,  # 226975
                }

    if run_sim == 'm1_agonist':
        ag_concs = np.logspace(-3, 3, 13)
        t_end = 1500
        r_dict = run_m1_activation_oxom_vs_ach(ag_concs, t_sim_end=t_end)

        m1_figs = plot_m1_activation_oxom_vs_ach(r_dict)

        all_figs += m1_figs
        all_names += ['rl_vs_t_ach', 'rl_vs_t_oxom', 'rl_vs_concentration']

    elif run_sim == 'recharge_er':
        #######################
        # Recharge Simulation
        #######################

        recharge_dict = run_recharge_simulation(sim_dict, sim_dur=180000)
        recharge_figs = plot_recharge(recharge_dict)

        t_est = np.arange(0, 180001, 20000)
        v_est = 177.0 * (1.0 - np.exp(-t_est / 59000))

        recharge_figs[1].line(t_est/1000.0, v_est, line_width=3, color='black', line_dash='dashed', legend='target, {} = 59 sec'.format(tau))
        recharge_figs[1].legend.location = 'bottom_right'
        all_figs += recharge_figs
        all_names += ['recharge_cytosol', 'recharge_er']

    elif run_sim == 'single_ap':

        ###########################
        # Single Action Potential
        ###########################
        # Target Single Action Potential
        # Apical Dendrite dF = 22%, decay time constant = 276 msec
        # Soma dF = 9%, decay time constant = 391 msec

        sim_dict['cbd_total'] = 0.045*0.2
        sim_dict['ogb1_total'] = 0.050

        # Stimulate Action Potential with Current Injection
        charge_int = [2000, 2010]
        r_d = run_current_injection(rxd_bool, param_dict=sim_dict, sim_dur=5200, c_int=charge_int, c_amp=0.5)

        time_ig = 1000
        reg_figs = plot_current_injection(r_d, rxd_bool, charge_int[0], t_ignore=time_ig)

        targ_max_s = 9.0  # (%)
        targ_tau_s = 391.0
        targ_max_a = 22.0  # (%)
        targ_tau_a = 276.0

        # Data taken from Dendrite Figure 5a from Power and Sah, 2002
        fig5a_dend = np.array([[-2.0408163265306136, -8.925244010647745], [-1.7346938775510221, -8.958518189884654],
                               [-1.7142857142857153, 2.995785270629984], [-1.6938775510204085, 10.6022626441881],
                               [-1.6530612244897966, 13.858695652173902],
                               [-1.612244897959183, 9.07165039929014], [-1.571428571428573, 6.241126885536815],
                               [-1.4897959183673475, 3.6235581188997266],
                               [-1.408163265306122, 1.0059893522626382], [-1.2653061224489797, -1.8356255545696598],
                               [-1.1224489795918373, -3.8076752440106496],
                               [-1.0408163265306136, -6.207852706299914], [-0.8163265306122458, -6.6670363797693],
                               [-0.6938775510204085, -7.33251996450754],
                               [-0.5102040816326543, -7.569875776397517], [-0.3061224489795915, -8.244232475598931],
                               [-0.20408163265306278, -8.255323868677902],
                               [-0.1020408163265305, -8.266415261756872], [0.02040816326530681, -9.14929015084294],
                               [0.22448979591836782, -9.606255545696534],
                               [0.26530612244898144, -9.610692102928123], [0.387755102040817, -9.189219165927234],
                               [0.4897959183673457, -9.852484472049685],
                               [0.6938775510204085, -10.961623779946756], [0.8571428571428577, -10.544587400177456],
                               [1.0612244897959204, -9.479813664596268]])

        fig5a_dend[:, 0] = (fig5a_dend[:, 0] - fig5a_dend[0, 0])*1000 + 700.0
        i_dend = np.squeeze(np.where(fig5a_dend[:, 0] < 2000))
        fig5a_dend[:, 1] = fig5a_dend[:, 1] - fig5a_dend[0, 1]

        # Data taken from Soma Figure 5a from Power and Sah, 2002
        fig5a_soma = np.array([[-2.0921658986174307, -0.03599999999999967], [-1.7952442396312627, -0.03714285714285687],
                               [-1.7679631336404782, 0.05142857142857172], [-1.742525345622047, 0.054285714285714576],
                               [-1.6913548387096053, 0.057714285714285996], [-1.6505806451612202, 0.04171428571428599],
                               [-1.6505806451612202, 0.04171428571428599], [-1.5726451612902554, 0.03771428571428598],
                               [-1.5186728110598402, 0.019428571428571677], [-1.4399999999999338, 0.009714285714285939],
                               [-1.2969585253455609, 0.001142857142857362],
                               [-1.2448294930874972, -0.0028571428571426555],
                               [-1.1263410138248275, -0.021142857142856963],
                               [-1.0632995391704512, -0.009714285714285523],
                               [-1.0100645161289776, -0.022285714285714117],
                               [-0.8286082949308256, -0.02857142857142843],
                               [-0.7909308755759881, -0.02057142857142842], [-0.7909308755759881, -0.02057142857142842],
                               [-0.7637972350229933, -0.030857142857142722],
                               [-0.7387281105990304, -0.02514285714285701],
                               [-0.6479631336405074, -0.028571428571428456],
                               [-0.6479631336405074, -0.028571428571428456],
                               [-0.5957603686635515, -0.03314285714285703], [-0.5320552995391292, -0.02685714285714276],
                               [-0.4280184331796839, -0.03314285714285706],
                               [-0.37611059907830224, -0.03542857142857135],
                               [-0.3247188940091803, -0.03371428571428564],
                               [-0.25931797235019616, -0.04057142857142851],
                               [-0.20962211981563428, -0.025714285714285648],
                               [-0.20962211981563428, -0.025714285714285648],
                               [-0.1820460829492756, -0.03942857142857137], [-0.1436313364054964, -0.0371428571428571],
                               [-0.09172350230411652, -0.03942857142857138],
                               [-0.09172350230411652, -0.03942857142857138],
                               [-0.013714285714257812, -0.04399999999999997],
                               [-0.013714285714257812, -0.04399999999999997],
                               [-0.0024331797234733443, -0.0314285714285714],
                               [0.07572350230417335, -0.03714285714285713],
                               [0.07572350230417335, -0.03714285714285713], [0.1399447004608554, -0.03485714285714285],
                               [0.24265437788020705, -0.030857142857142875], [0.2957419354838908, -0.0422857142857143],
                               [0.3467649769585446, -0.03771428571428574], [0.3730875576037054, -0.04171428571428573],
                               [0.47631336405531677, -0.041714285714285745],
                               [0.47631336405531677, -0.041714285714285745],
                               [0.4884055299539334, -0.035428571428571476], [0.5271889400921808, -0.036000000000000046],
                               [0.5796129032258204, -0.04228571428571434], [0.5796129032258204, -0.04228571428571434],
                               [0.6175115207373398, -0.036000000000000046], [0.6571797235023169, -0.0434285714285715],
                               [0.6571797235023169, -0.0434285714285715], [0.7079815668202869, -0.037142857142857214],
                               [0.8121658986175202, -0.044571428571428665], [0.8121658986175202, -0.044571428571428665],
                               [0.8497695852534637, -0.03600000000000009], [0.953290322580651, -0.038285714285714395],
                               [1.0046082949308772, -0.0360000000000001], [1.0443502304147483, -0.04400000000000012],
                               [1.095078341013826, -0.03714285714285727]])

        fig5a_soma[:, 0] = (fig5a_soma[:, 0] - fig5a_soma[0, 0])*1000 + 700
        i_soma = np.squeeze(np.where(fig5a_soma[:, 0] < 2000))
        fig5a_soma[:, 1] = (fig5a_soma[:, 1] - fig5a_soma[0, 1])*100

        i_ig = np.squeeze(np.where(r_d['t'] > time_ig))
        i_max_s = np.argmax(r_d['apic9_dyeca'][i_ig])
        t_max_a = r_d['t'][i_ig[i_max_s]] - time_ig
        t_est = np.arange(r_d['t'][i_ig[i_max_s]], r_d['t'][-1], 20) - time_ig
        est_val_s = targ_max_s * np.exp(-(t_est - t_max_a) / targ_tau_s)
        est_val_a = targ_max_a * np.exp(-(t_est - t_max_a) / targ_tau_a)

        # reg_figs[6].line(t_est, est_val_s, line_width=3, color='darkorange', legend='Target-Soma', line_dash='dashed')
        # reg_figs[6].line(t_est, est_val_a, line_width=3, color='red', legend='Target-Apical Trunk', line_dash='dashed')
        reg_figs[6].line(fig5a_soma[:, 0][i_soma], fig5a_soma[:, 1][i_soma], line_width=3,
                         color='darkorange', legend='Target-Soma', line_dash='dashed')
        reg_figs[6].line(fig5a_dend[:, 0][i_dend], fig5a_dend[:, 1][i_dend], line_width=3,
                         color='red', legend='Target-Apical Trunk', line_dash='dashed')

        all_figs += reg_figs
        if rxd_bool:
            plot_names = ['one_ap_v_vs_t', 'one_ap_i_ing_vs_t', 'one_ap_cyt_vs_t', 'one_ap_calcium_currents_vs_t',
                          'one_ap_er_vs_t', 'one_ap_cbd_ca_vs_t', 'one_ap_dF_vs_t', 'one_ap_car_ca_vs_t',
                          'one_ap_rO_ip3r_vs_t']
        else:
            plot_names = ['one_ap_v_vs_t', 'one_ap_i_ing_vs_t', 'one_ap_cyt_vs_t', 'one_ap_calcium_currents_vs_t']
        all_names += plot_names

    elif run_sim == 'four_aps':

        ###############################
        # Four 20 Hz Action Potentials
        ###############################
        # Target Four 20 Hz Action Potentials
        # Apical Dendrite dF = 100%, decay time constant = 211 msec
        # Soma dF = 22%, decay time constant = 343 msec

        r_d2 = run_current_injection_series(rxd_bool, param_dict=sim_dict, sim_dur=500, pulse_times=[50, 100, 150, 200], pulse_amps=[0.7, 0.7, 0.7, 0.7])
        reg_figs2 = plot_current_injection(r_d2, rxd_bool)

        targ_max_s = 22.0  # (%)
        targ_tau_s = 343.0
        targ_max_a = 100.0  # (%)
        targ_tau_a = 211.0
        i_max_s = np.argmax(r_d2['soma_dyeca'])
        t_est = np.arange(r_d2['t'][i_max_s], r_d2['t'][-1], 10)
        est_val_s = targ_max_s * np.exp(-(t_est - r_d2['t'][i_max_s]) / targ_tau_s)
        est_val_a = targ_max_a * np.exp(-(t_est - r_d2['t'][i_max_s]) / targ_tau_a)

        reg_figs2[6].line(t_est, est_val_s, line_width=3, color='red', legend='target soma', line_dash='dashed')
        reg_figs2[6].line(t_est, est_val_a, line_width=3, color='darkorange', legend='target apical', line_dash='dashed')

        all_figs += reg_figs2
        if rxd_bool:
            plot_names = ['four_ap_v_vs_t', 'four_ap_i_ing_vs_t', 'four_ap_cyt_vs_t', 'four_ap_calcium_currents_vs_t',
                          'four_ap_er_vs_t', 'four_ap_cbd_ca_vs_t', 'four_ap_dF_vs_t', 'four_ap_car_ca_vs_t']
        else:
            plot_names = ['four_ap_v_vs_t', 'four_ap_i_ing_vs_t', 'four_ap_cyt_vs_t', 'four_ap_calcium_currents_vs_t']
        all_names += plot_names

    elif run_sim == 'ach_pulse':
        ############
        # ACh Pulse
        ############

        sim_dur = 5000
        t_ignore = 0
        h5_path = os.path.join(res_dir, 'ach_pulse_spikes.h5')

        my_steps = np.arange(0, sim_dur+1, 50)

        c_dict = {'amp': 0.0, # 0.525,
                  'start': 1.0,
                  'dur': sim_dur,
                  }

        if not use_file_to_plot:
            rd_ach = run_ach_pulse(my_steps, sim_dict, pulse_times=[200], pulse_amps=[0.1],
                                   pulse_length=50, current_dict=c_dict)

            if save_result:

                save_dict = {}
                # Only save numpy arrays holding time series
                for k, v in rd_ach.items():
                    if type(v) == np.ndarray:
                        save_dict[k] = v

                save_df = pd.DataFrame.from_dict(save_dict)

                with pd.HDFStore(h5_path) as h5file:
                    h5file.put('ach_pulse', save_df, format='table', data_columns=True)
                print('Results saved to {}'.format(h5_path))

            ach_figs = plot_ach_pulse(rd_ach, t_is=[0, 1], ach_times=[200, 250], t_ignore=t_ignore)

        else:
            with pd.HDFStore(h5_path) as h5file:
                ton_dict = {'soma_v_dict': {},
                            'ca_cyt_dict': {},
                            'ca_er_dict': {},
                            'soma_im_dict': {},
                            'soma_isk_dict': {},
                            'soma_pip2_dict': {},
                            't': {}}

                for k in h5file.keys():
                    if 'soma_v' in k:
                        ton_dict['soma_v_dict'] = h5file[k]
                    elif 'soma_im' in k:
                        ai = int(k.split('_')[-1])
                        ton_dict['soma_im_dict'][ai] = h5file[k]
                    elif 'soma_isk' in k:
                        ai = int(k.split('_')[-1])
                        ton_dict['soma_isk_dict'][ai] = h5file[k]
                    elif 'ca_cyt' in k:
                        ai = int(k.split('_')[-1])
                        ton_dict['ca_cyt_dict'][ai] = h5file[k]
                    elif 'ca_er' in k:
                        ai = int(k.split('_')[-1])
                        ton_dict['ca_er_dict'][ai] = h5file[k]
                    elif 't_dict' in k:
                        ai = int(k.split('_')[-1])
                        ton_dict['t'][ai] = h5file[k]

                rd_ach = h5file['ach_pulse']

                ach_figs = plot_ach_pulse(rd_ach, t_is=[0, 1], ach_times=[2200, 2250], t_ignore=t_ignore)

        ach_figs[-1].legend.visible = False

        ach_names = ['ach_pulse_cyt_vs_t', 'ach_pulse_er_vs_t', 'ach_pulse_ach_vs_t', 'ach_pulse_rl_vs_t',
                     'ach_pulse_active_plc_vs_t', 'ach_pulse_pip2_vs_t', 'ach_pulse_pi_flux_vs_t', 'ach_pulse_ip3_vs_t',
                     'ach_pulse_open_ip3r_vs_t', 'ach_pulse_ip3_breakdown_vs_t', 'ach_pulse_perc_im_vs_t',
                     'ach_pulse_ik_vs_t', 'ach_pulse_v_vs_t']

        # if len(ach_figs) > len(ach_names):
        #     add_dasari_isr_data(ach_figs[-1])
        #
        #     ach_names.append('ach_pulse_ifr_vs_t')
        #     ach_names.append('ach_pulse_spike_accel_vs_t')

        ach_names_dose = [a_name + '_100uM_50msec_ach' for a_name in ach_names]

        all_figs += ach_figs
        all_names += ach_names_dose

    elif run_sim == 'test_refill_spiking':
        #################
        # Test ER refill
        #################

        sim_dur = 24200

        h5_path = os.path.join(res_dir, 'ach_pulse_refill_test_spikes.h5')

        my_steps = np.arange(0, sim_dur + 1, 50)

        c_dict = {'amp': 0.525,
                  'start': 1.0,
                  'dur': sim_dur,
                  }

        if not use_file_to_plot:
            rd_ach = run_ach_pulse(my_steps, sim_dict, pulse_times=[2200, 10200, 18200], pulse_amps=[0.1, 0.1, 0.1],
                                   pulse_length=50, current_dict=c_dict)

            if save_result:

                save_dict = {}
                # Only save numpy arrays holding time series
                for k, v in rd_ach.items():
                    if type(v) == np.ndarray:
                        save_dict[k] = v

                save_df = pd.DataFrame.from_dict(save_dict)

                with pd.HDFStore(h5_path) as h5file:
                    h5file.put('ach_pulse', save_df, format='table', data_columns=True)
                print('Results saved to {}'.format(h5_path))

            ach_figs = plot_ach_pulse(rd_ach, t_is=[0, 1], ach_times=[2200, 2250, 10200, 10250, 18200, 18250], t_ignore=200)

        else:
            with pd.HDFStore(h5_path) as h5file:

                rd_ach = h5file['ach_pulse']

                ach_figs = plot_ach_pulse(rd_ach, t_is=[0, 1], ach_times=[2200, 2250, 10200, 10250, 18200, 18250], t_ignore=200)

        ach_figs[-1].legend.visible = False

        ach_names = ['ach_refill_test_cyt_vs_t', 'ach_refill_test_er_vs_t', 'ach_refill_test_ach_vs_t', 'ach_refill_test_rl_vs_t',
                     'ach_refill_test_active_plc_vs_t', 'ach_refill_test_pip2_vs_t', 'ach_refill_test_pi_flux_vs_t', 'ach_refill_test_ip3_vs_t',
                     'ach_refill_test_open_ip3r_vs_t', 'ach_refill_test_ip3_breakdown_vs_t', 'ach_refill_test_perc_im_vs_t',
                     'ach_refill_test_ik_vs_t', 'ach_refill_test_v_vs_t']

        # if len(ach_figs) > len(ach_names):
        #     add_dasari_isr_data(ach_figs[-1])
        #
        #     ach_names.append('ach_pulse_ifr_vs_t')
        #     ach_names.append('ach_pulse_spike_accel_vs_t')

        ach_names_dose = [a_name + '_refill_100uM_50msec_0p125Hz_ach_spikes' for a_name in ach_names]

        all_figs += ach_figs
        all_names += ach_names_dose

    elif run_sim == 'test_refill_resting':
        #################
        # Test ER refill
        #################

        sim_dur = 24200

        h5_path = os.path.join(res_dir, 'ach_pulse_refill_test_resting.h5')

        my_steps = np.arange(0, sim_dur + 1, 50)

        c_dict = {'amp': 0.0,
                  'start': 1.0,
                  'dur': sim_dur,
                  }

        if not use_file_to_plot:
            rd_ach = run_ach_pulse(my_steps, sim_dict, pulse_times=[2200, 10200, 18200], pulse_amps=[0.1, 0.1, 0.1],
                                   pulse_length=50, current_dict=c_dict)

            if save_result:

                save_dict = {}
                # Only save numpy arrays holding time series
                for k, v in rd_ach.items():
                    if type(v) == np.ndarray:
                        save_dict[k] = v

                save_df = pd.DataFrame.from_dict(save_dict)

                with pd.HDFStore(h5_path) as h5file:
                    h5file.put('ach_pulse', save_df, format='table', data_columns=True)
                print('Results saved to {}'.format(h5_path))

            ach_figs = plot_ach_pulse(rd_ach, t_is=[0, 1], ach_times=[2200, 2250, 10200, 10250, 18200, 18250], t_ignore=200)

        else:
            with pd.HDFStore(h5_path) as h5file:

                rd_ach = h5file['ach_pulse']

                ach_figs = plot_ach_pulse(rd_ach, t_is=[0, 1], ach_times=[2200, 2250, 10200, 10250, 18200, 18250], t_ignore=200)

        ach_figs[-1].legend.visible = False

        ach_names = ['ach_refill_test__cyt_vs_t', 'ach_refill_test_er_vs_t', 'ach_refill_test_ach_vs_t', 'ach_refill_test_rl_vs_t',
                     'ach_refill_test_active_plc_vs_t', 'ach_refill_test_pip2_vs_t', 'ach_refill_test_pi_flux_vs_t', 'ach_refill_test_ip3_vs_t',
                     'ach_refill_test_open_ip3r_vs_t', 'ach_refill_test_ip3_breakdown_vs_t', 'ach_refill_test_perc_im_vs_t',
                     'ach_refill_test_ik_vs_t', 'ach_refill_test_v_vs_t']

        # if len(ach_figs) > len(ach_names):
        #     add_dasari_isr_data(ach_figs[-1])
        #
        #     ach_names.append('ach_pulse_ifr_vs_t')
        #     ach_names.append('ach_pulse_spike_accel_vs_t')

        ach_names_dose = [a_name + '_refill_100uM_50msec_0p125Hz_ach_resting' for a_name in ach_names]

        all_figs += ach_figs
        all_names += ach_names_dose
    
    if plot_dasari:
        #################################
        # Plot Dasari and Gulledge Data
        #################################
        ifr_fig = bplt.figure(title='Instantaneous Firing Rate vs Time')
        add_dasari_isr_data(ifr_fig)

        ifr_fig.xaxis.axis_label = 'time (msec)'
        ifr_fig.yaxis.axis_label = 'IFR (Hz)'

        all_figs += [ifr_fig]
        all_figs += plot_dasari_and_gulledge_data()
        all_names += ['dasari_and_gulledge_ifr', 'dasari_and_gulledge_1a2_rmp', 'dasari_and_gulledge_2a_rmp']

    for fig in all_figs:
        plt_res.configure_plot(fig)

    # plt.show()
    bkio.show(blay.column(all_figs))

    if save_figs:

        driver_path = r'/usr/local/bin/geckodriver'
        my_opts = webdriver.firefox.options.Options()
        # my_opts.add_argument('start-maximized')
        # my_opts.add_argument('disable-infobars')
        my_opts.add_argument('--disable-extensions')
        my_opts.add_argument('window-size=1200,1000')
        my_opts.headless = True

        my_driver = webdriver.Firefox(options=my_opts, executable_path=driver_path)

        for fig, name in zip(all_figs, all_names):
            fig_path = os.path.join(res_dir, '{}.png'.format(name))
            bkio.export_png(fig, filename=fig_path, webdriver=my_driver)

