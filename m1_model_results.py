import datetime as dat
import frequency_analysis as fan
import ca1_pyramidal as ca1p
import config_utils as conf
import numpy as np
import os
import pandas as pd
import pickle
import plot_results as plt_res
import time

from scipy import integrate as scin
from scipy import optimize as spopt
from neuron import h
from os.path import join

import bokeh.io as bkio
import bokeh.layouts as blay
import bokeh.models as bmod
import bokeh.plotting as bplt
import bokeh.palettes as bpal

from bokeh.palettes import Category20 as palette
from bokeh.palettes import Category20b as paletteb
# from selenium import webdriver

colrs = palette[20] + paletteb[20]

mu = u'\u03BC'
delta = u'\u0394'
ohm = u'\u03A9'


def run_ach_pulse(t_steps, param_dict={}, pulse_times=[50], pulse_amps=[0.001], pulse_amps_high=[0.001],
                  pulse_length=10, current_dict={'amp': 0.0, 'start': 0.0, 'dur': 0.0}):
    p_times = []
    p_amps_high = []
    p_amps_f = []
    p_amps_b = []

    for p_t, p_a, p_h in zip(pulse_times, pulse_amps, pulse_amps_high):
        p_times += [p_t, p_t + pulse_length]
        p_amps_high += [p_h, 0]
        p_amps_f += [p_a, 0]
        p_amps_b += [0, p_a]

    print('Pulse Times: {}'.format(p_times))
    print('Pulse Amps: {}'.format(p_amps_f))

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = ca1p.MyCell(t_path, param_dict)
    cell.insert_rxd()

    if current_dict['amp']:
        t_curr = h.IClamp(cell.somatic[0](0.5))
        t_curr.amp = current_dict['amp']
        t_curr.delay = current_dict['start']
        t_curr.dur = current_dict['dur']
        print(
            'Current Injection at t = {0} msec, amplitude = {1} nA'.format(current_dict['start'], current_dict['amp']))

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

    # s_dag = h.Vector().record(cell.somatic[0](0.5)._ref_DAG_m1_kmb)
    # s_pip2 = h.Vector().record(cell.somatic[0](0.5)._ref_PIP2_m1_kmb)
    # s_pip2_kcnq = h.Vector().record(cell.somatic[0](0.5)._ref_PIP2_KCNQ_m1_kmb)
    # s_ga_gtp = h.Vector().record(cell.somatic[0](0.5)._ref_Ga_GTP_m1_kmb)
    # s_plc = h.Vector().record(cell.somatic[0](0.5)._ref_PLC_m1_kmb)
    # s_active_plc = h.Vector().record(cell.somatic[0](0.5)._ref_Ga_GTP_PLC_m1_kmb)

    s_dag = h.Vector().record(cell.dag.nodes(cell.somatic[0])[0]._ref_concentration)

    s_pip2 = h.Vector().record(cell.pip2.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_pip2 = h.Vector().record(cell.pip2.nodes(cell.apical[0])[0]._ref_concentration)
    a9_pip2 = h.Vector().record(cell.pip2.nodes(cell.apical[9])[0]._ref_concentration)
    ax_pip2 = h.Vector().record(cell.pip2.nodes(cell.axonal[0])[0]._ref_concentration)

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
    # s_ip3 = h.Vector().record(cell.somatic[0](0.5)._ref_IP3i)
    a0_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[0])[0]._ref_concentration)
    a9_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[9])[0]._ref_concentration)
    s_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    a9_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.apical[9])[0]._ref_concentration)
    s_po = h.Vector().record(cell.ro_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    a0_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    a9_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[9])[0]._ref_concentration)

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
            print('time: {0}, [ACh]: {1} {2}M'.format(t_step, p_amps_f[p_i] * 1000.0, mu))

            cell.ach.nodes(cell.cyt).concentration = p_amps_f[p_i]

            h.CVode().re_init()

    res_dict = {'node names': node_locs,
                'node distances': node_dists,
                'soma v': np.array(s_v),
                'apical0 v': np.array(a0_v),
                'apical9 v': np.array(a9_v),
                'axon v': np.array(ax_v),
                'soma cyt time': np.array(s_ca_cyt),
                'apical0 cyt time': np.array(a0_ca_cyt),
                'apical9 cyt time': np.array(a9_ca_cyt),
                'soma er time': np.array(s_ca_er),
                'apical0 er time': np.array(a0_ca_er),
                'apical9 er time': np.array(a9_ca_er),
                'soma ip3 time': np.array(s_ip3),
                'apical0 ip3 time': np.array(a0_ip3),
                'apical9 ip3 time': np.array(a9_ip3),
                'soma ri': np.array(s_ri),
                'apical0 ri': np.array(a0_ri),
                'apical9 ri': np.array(a9_ri),
                'soma open probability': np.array(s_po),
                'apical0 open probability': np.array(a0_po),
                'apical9 open probability': np.array(a9_po),
                'soma im': np.array(s_im),
                'axon im': np.array(ax_im),
                'soma isk': np.array(s_isk),
                'apical0 isk': np.array(a0_isk),
                'apical9 isk': np.array(a9_isk),
                'soma dag': np.array(s_dag),
                'soma pip2': np.array(s_pip2),
                'soma pi': np.array(s_pi),
                'soma pi4p': np.array(s_pi4p),
                'apical0 pip2': np.array(a0_pip2),
                'apical9 pip2': np.array(a9_pip2),
                'axon pip2': np.array(ax_pip2),
                'soma ach': np.array(s_ach),
                'soma pip2-kcnq': np.array(s_pip2_kcnq),
                'soma perc_i': np.array(s_perc_i),
                'axon perc_i': np.array(ax_perc_i),
                'soma ga-gtp': np.array(s_ga_gtp),
                'soma plc': np.array(s_plc),
                'soma active plc': np.array(s_active_plc),
                'soma r': np.array(s_r),
                'soma g': np.array(s_g),
                'soma gbg': np.array(s_gbg),
                'soma rg': np.array(s_rg),
                'soma rgbg': np.array(s_rgbg),
                'soma rl': np.array(s_rl),
                'soma rlg': np.array(s_rlg),
                'soma ip5p': np.array(s_ip5p),
                'soma ip5p_ip3': np.array(s_ip5p_ip3),
                'soma ip3k': np.array(s_ip3k),
                'soma ip3k_2ca': np.array(s_ip3k_2ca),
                'soma ip3k_2ca_ip3': np.array(s_ip3k_2ca_ip3),
                't': np.array(t_vec),
                }
    return res_dict


def plot_ach_pulse(result_d, ach_times=[100, 150], t_ignore=0):
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
    im_t_fig = bplt.figure(title='M Current vs Time')
    im_t_fig.xaxis.axis_label = 'time (msec)'
    im_t_fig.yaxis.axis_label = 'current (mA/cm^2)'.format(mu)
    ret_figs.append(im_t_fig)
    isk_t_fig = bplt.figure(title='SK Current vs Time')
    isk_t_fig.xaxis.axis_label = 'time (msec)'
    isk_t_fig.yaxis.axis_label = 'current (mA/cm^2)'.format(mu)
    ret_figs.append(isk_t_fig)
    v_fig = bplt.figure(title='Membrane Potential vs Time')
    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(v_fig)

    sp_times = fan.get_spike_times(result_d['soma v'], result_d['t'])

    t_arr = result_d['t']

    i_ig = np.where(t_arr > t_ignore)
    t_ig = t_arr[i_ig] - t_ignore

    if sp_times.size > 3:
        isr_fig = bplt.figure(title='Instantaneous Spiking Frequency vs Time', x_range=[t_ig[0], t_ig[-1]])
        isr_fig.xaxis.axis_label = 'time (msec)'
        isr_fig.yaxis.axis_label = 'ISF (Hz)'
        ret_figs.append(isr_fig)

        isr_vals = np.divide(1000.0, np.diff(sp_times))
        isr_ts = (sp_times[:-1] + sp_times[1:]) / 2.0

        isr_is = np.where(isr_ts > t_ignore)

        isr_fig.circle(isr_ts[isr_is] - t_ignore, isr_vals[isr_is], size=12, color=colrs[0], legend='Simulation Result')
        isr_fig.line(isr_ts[isr_is] - t_ignore, isr_vals[isr_is], line_width=3, color=colrs[0])

    # dumb_fig = bplt.figure(title='Dumby Variables to Test Density/Concentration')
    # dumb_fig.xaxis.axis_label = 'time (msec)'
    # ret_figs.append(dumb_fig)

    ach_span_1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    cyt_t_fig.line(t_ig, result_d['soma cyt time'][i_ig] * 1000.0, line_width=3, color=colrs[0], legend='somatic[0]')
    cyt_t_fig.line(t_ig, result_d['apical0 cyt time'][i_ig] * 1000.0, line_width=3, color=colrs[1],
                   legend='apical[0]')
    cyt_t_fig.line(t_ig, result_d['apical9 cyt time'][i_ig] * 1000.0, line_width=3, color=colrs[2],
                   legend='apical[9]',
                   line_dash='dashed')
    cyt_t_fig.add_layout(ach_span_1)
    cyt_t_fig.add_layout(ach_span_2)

    er_t_fig.line(t_ig, result_d['soma er time'][i_ig] * 1000.0, line_width=3, color=colrs[0], legend='somatic[0]')
    er_t_fig.line(t_ig, result_d['apical0 er time'][i_ig] * 1000.0, line_width=3, color=colrs[1], legend='apical[0]')
    er_t_fig.line(t_ig, result_d['apical9 er time'][i_ig] * 1000.0, line_width=3, color=colrs[2], legend='apical[9]',
                  line_dash='dashed')
    er_t_fig.add_layout(ach_span_1)
    er_t_fig.add_layout(ach_span_2)

    ach_t_fig.line(t_ig, result_d['soma ach'][i_ig] * 1000.0, line_width=3, color='green', legend='somatic[0](0.5)')

    m1_t_fig.line(t_ig, result_d['soma r'][i_ig], line_width=3, color=colrs[0], legend='R')
    m1_t_fig.line(t_ig, result_d['soma g'][i_ig], line_width=3, color=colrs[2], legend='G')
    m1_t_fig.line(t_ig, result_d['soma gbg'][i_ig], line_width=3, color=colrs[3], legend='G-bg')
    m1_t_fig.line(t_ig, result_d['soma rg'][i_ig], line_width=3, color=colrs[4], legend='RG')
    m1_t_fig.line(t_ig, result_d['soma rgbg'][i_ig], line_width=3, color=colrs[5], legend='RG-bg')
    m1_t_fig.line(t_ig, result_d['soma rl'][i_ig], line_width=3, color=colrs[6], legend='RL')
    m1_t_fig.line(t_ig, result_d['soma rlg'][i_ig], line_width=3, color=colrs[7], legend='RLG')
    # m1_t_fig.line(t_ig, result_d['soma dag'], line_width=3, color=colrs[8], legend='DAG')
    m1_t_fig.add_layout(ach_span_1)
    m1_t_fig.add_layout(ach_span_2)

    act_plc_t_fig.line(t_ig, result_d['soma rlg'][i_ig], line_width=3, color=colrs[4], legend='RLG')
    act_plc_t_fig.line(t_ig, result_d['soma ga-gtp'][i_ig], line_width=3, color=colrs[0],
                       legend='ga-gtp')
    # act_plc_t_fig.line(t_ig, result_d['soma plc'], line_width=3, color=colrs[1], legend='plc')
    act_plc_t_fig.line(t_ig, result_d['soma active plc'][i_ig], line_width=3, color=colrs[2], legend='ga-gtp-plc')
    act_plc_t_fig.add_layout(ach_span_1)
    act_plc_t_fig.add_layout(ach_span_2)

    pip2_t_fig.line(t_ig, result_d['soma pip2'][i_ig], line_width=3, color=colrs[0], legend='somatic[0](0.5) PIP2')
    pip2_t_fig.line(t_ig, result_d['apical0 pip2'][i_ig], line_width=3, color=colrs[1], legend='apical[0](0.5) PIP2')
    pip2_t_fig.line(t_ig, result_d['apical9 pip2'][i_ig], line_width=3, color=colrs[2], legend='apical[9](0.5) PIP2')
    pip2_t_fig.line(t_ig, result_d['axon pip2'][i_ig], line_width=3, color=colrs[3], legend='axonal[0](0.5) PIP2',
                    line_dash='dashed')
    pip2_t_fig.line(t_ig, result_d['soma pi4p'][i_ig], line_width=3, color=colrs[4], legend='somatic[0](0.5) PI4P')

    pip2_t_fig.add_layout(ach_span_1)
    pip2_t_fig.add_layout(ach_span_2)

    k_PLC = 0.3  # 0.0003  # (um^2 msec^-1) dropped by a fact
    k_4K = 0.0008 * 2  # 0.0000008  # (msec^-1)
    k_4P = 0.12  # 0.00012  # (msec^-1)
    k_5K = 0.02  # 0.00002  # (msec^-1)
    k_5P = 0.028  # 0.000028  # (msec^-1)

    pi_flux_fig.line(t_ig, result_d['soma pi'][i_ig] * k_4K, line_width=3, color=colrs[0],
                     legend='somatic[0](0.5) PI->PI4P')
    pi_flux_fig.line(t_ig, result_d['soma pi4p'][i_ig] * k_4P, line_width=3, color=colrs[1],
                     legend='somatic[0](0.5) PI4P->PI')
    pi_flux_fig.line(t_ig, result_d['soma pi4p'][i_ig] * k_5K, line_width=3, color=colrs[2],
                     legend='somatic[0](0.5) PI4P->PIP2')
    pi_flux_fig.line(t_ig, result_d['soma pip2'][i_ig] * k_5P, line_width=3, color=colrs[3],
                     legend='somatic[0](0.5) PIP2->PI4P')
    pi_flux_fig.line(t_ig, result_d['soma active plc'][i_ig] * result_d['soma pip2'][i_ig] * k_PLC * 3, line_width=3,
                     color=colrs[4], legend='somatic[0](0.5) PIP2->DAG+IP3')

    s_ip3_flux = result_d['soma active plc'] * result_d['soma pip2'] * k_PLC * 3 / 602214.129

    dip3 = np.diff(result_d['soma ip3 time'][i_ig]) / np.diff(t_ig)

    # ip3_t_fig.line(t_ig, dip3, line_width=3, color=colrs[0], legend='somatic[0] flux')
    ip3_t_fig.line(t_ig, result_d['soma ip3 time'][i_ig] * 1000.0, line_width=3, color=colrs[0], legend='somatic[0]')
    ip3_t_fig.line(t_ig, result_d['apical0 ip3 time'][i_ig] * 1000.0, line_width=3, color=colrs[1],
                   legend='apical[0]')
    ip3_t_fig.line(t_ig, result_d['apical9 ip3 time'][i_ig] * 1000.0, line_width=3, color=colrs[2],
                   legend='apical[9]', line_dash='dashed')
    ip3_t_fig.add_layout(ach_span_1)
    ip3_t_fig.add_layout(ach_span_2)

    ip3r_open_t_fig.line(t_ig, 100 * result_d['soma open probability'][i_ig], line_width=3, color=colrs[0],
                         legend='somatic[0]')
    ip3r_open_t_fig.line(t_ig, 100 * result_d['soma ri'][i_ig], line_width=3, color=colrs[0],
                         legend='somatic[0] ri', line_dash='dashed')
    ip3r_open_t_fig.line(t_ig, 100 * result_d['apical0 open probability'][i_ig], line_width=3, color=colrs[1],
                         legend='apical[0]')
    ip3r_open_t_fig.line(t_ig, 100 * result_d['apical0 ri'][i_ig], line_width=3, color=colrs[1],
                         legend='apical[0] ri', line_dash='dashed')
    ip3r_open_t_fig.line(t_ig, 100 * result_d['apical9 open probability'][i_ig], line_width=3, color=colrs[2],
                         legend='apical[9]')
    ip3r_open_t_fig.line(t_ig, 100 * result_d['apical9 ri'][i_ig], line_width=3, color=colrs[2],
                         legend='apical[9] ri', line_dash='dashed')
    ip3r_open_t_fig.add_layout(ach_span_1)
    ip3r_open_t_fig.add_layout(ach_span_2)

    ip3_kinase_fig.line(t_ig, 1000.0 * result_d['soma ip5p'][i_ig],
                        line_width=3, color=colrs[0], legend='ip5p')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_d['soma ip5p_ip3'][i_ig],
                        line_width=3, line_dash='dashed', color=colrs[1], legend='ip5p_ip3')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_d['soma ip3k'][i_ig],
                        line_width=3, color=colrs[2], legend='ip3k')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_d['soma ip3k_2ca'][i_ig],
                        line_width=3, line_dash='dashed', color=colrs[3], legend='ip3k_2ca')
    ip3_kinase_fig.line(t_ig, 1000.0 * result_d['soma ip3k_2ca_ip3'][i_ig],
                        line_width=3, line_dash='dotted', color=colrs[4], legend='ip3k_2ca_ip3')

    pip2_kcnq_t_fig.line(t_ig, 100.0 * result_d['soma perc_i'][i_ig], line_width=3, color=colrs[0],
                         legend='somatic[0](0.5)')
    pip2_kcnq_t_fig.line(t_ig, 100.0 * result_d['axon perc_i'][i_ig], line_width=3, color=colrs[2],
                         legend='axon[0](0.5)')
    pip2_kcnq_t_fig.add_layout(ach_span_1)
    pip2_kcnq_t_fig.add_layout(ach_span_2)

    im_t_fig.line(t_ig, result_d['soma im'][i_ig], line_width=3, color=colrs[0], legend='somatic[0]')
    im_t_fig.line(t_ig, result_d['axon im'][i_ig], line_width=3, color=colrs[2], legend='axonal[0]')
    im_t_fig.add_layout(ach_span_1)
    im_t_fig.add_layout(ach_span_2)

    isk_t_fig.line(t_ig, result_d['soma isk'][i_ig], line_width=3, color=colrs[0], legend='somatic[0]')
    isk_t_fig.line(t_ig, result_d['apical0 isk'][i_ig], line_width=3, color=colrs[1],
                   legend='apical[0]')
    isk_t_fig.line(t_ig, result_d['apical9 isk'][i_ig], line_width=3, color=colrs[2],
                   legend='apical[9]', line_dash='dashed')
    isk_t_fig.add_layout(ach_span_1)
    isk_t_fig.add_layout(ach_span_2)

    v_fig.line(t_ig, result_d['soma v'][i_ig], line_width=3, color=colrs[0], legend='soma[0](0.5)')
    v_fig.line(t_ig, result_d['apical0 v'][i_ig], line_width=3, color=colrs[1], legend='apical[0](0.5)',
               line_dash='dashed')
    v_fig.line(t_ig, result_d['apical9 v'][i_ig], line_width=3, color=colrs[2], legend='apical[9](0.5)',
               line_dash='dotted')
    v_fig.line(t_ig, result_d['axon v'][i_ig], line_width=3, color=colrs[3], legend='axonal[0](0.5)',
               line_dash='dotted')
    v_fig.add_layout(ach_span_1)
    v_fig.add_layout(ach_span_2)

    return ret_figs


def run_ahp_test(param_dict={}, ach_levels=[0.0, 100.0], run_time=2000,
                 current_dict={'amp': 0.0, 'start': 0.0, 'dur': 0.0}):

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = ca1p.MyCell(t_path, param_dict)
    cell.insert_rxd()

    if current_dict['amp']:
        t_curr = h.IClamp(cell.somatic[0](0.5))
        t_curr.amp = current_dict['amp']
        t_curr.delay = current_dict['start']
        t_curr.dur = current_dict['dur']
        print(
            'Current Injection at t = {0} msec, amplitude = {1} nA'.format(current_dict['start'], current_dict['amp']))

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)

    s_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    ax_v = h.Vector().record(cell.axonal[0](0.5)._ref_v)
    ap_v = h.Vector().record(cell.apical[9](0.5)._ref_v)

    s_ina = h.Vector().record(cell.somatic[0](0.5)._ref_ina_nax)
    s_isk = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kca)
    s_ikdr = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kdr)
    s_ibk = h.Vector().record(cell.somatic[0](0.5)._ref_ik_cagk)
    s_ih = h.Vector().record(cell.somatic[0](0.5)._ref_i_hd)
    s_im = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kmb_inh)
    s_ipas = h.Vector().record(cell.somatic[0](0.5)._ref_i_pas)
    s_ikap = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kap)
    s_ical = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cal)
    s_ican = h.Vector().record(cell.somatic[0](0.5)._ref_ica_can)
    s_icat = h.Vector().record(cell.somatic[0](0.5)._ref_ica_cat)
    ax_ikap = h.Vector().record(cell.axonal[0](0.5)._ref_ik_kap)
    ax_im = h.Vector().record(cell.axonal[0](0.5)._ref_ik_kmb_inh)
    ax_ikdr = h.Vector().record(cell.axonal[0](0.5)._ref_ik_kdr)

    # s_ach = h.Vector().record(cell.ach.nodes(cell.somatic[0])[0]._ref_concentration)

    # s_dag = h.Vector().record(cell.somatic[0](0.5)._ref_DAG_m1_kmb)
    # s_pip2 = h.Vector().record(cell.somatic[0](0.5)._ref_PIP2_m1_kmb)
    # s_pip2_kcnq = h.Vector().record(cell.somatic[0](0.5)._ref_PIP2_KCNQ_m1_kmb)
    # s_ga_gtp = h.Vector().record(cell.somatic[0](0.5)._ref_Ga_GTP_m1_kmb)
    # s_plc = h.Vector().record(cell.somatic[0](0.5)._ref_PLC_m1_kmb)
    # s_active_plc = h.Vector().record(cell.somatic[0](0.5)._ref_Ga_GTP_PLC_m1_kmb)

    # s_dag = h.Vector().record(cell.dag.nodes(cell.somatic[0])[0]._ref_concentration)
    #
    # s_pip2 = h.Vector().record(cell.pip2.nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_pip2 = h.Vector().record(cell.pip2.nodes(cell.apical[0])[0]._ref_concentration)
    # a9_pip2 = h.Vector().record(cell.pip2.nodes(cell.apical[9])[0]._ref_concentration)
    # ax_pip2 = h.Vector().record(cell.pip2.nodes(cell.axonal[0])[0]._ref_concentration)
    #
    # s_pi4p = h.Vector().record(cell.pi4p.nodes(cell.somatic[0])[0]._ref_concentration)
    #
    # s_pi = h.Vector().record(cell.pi.nodes(cell.somatic[0])[0]._ref_concentration)
    #
    # s_pip2_kcnq = h.Vector().record(cell.pip2_kcnq.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_perc_i = h.Vector().record(cell.somatic[0](0.5)._ref_perc_i_kmb_inh)
    # ax_perc_i = h.Vector().record(cell.axonal[0](0.5)._ref_perc_i_kmb_inh)
    # s_ga_gtp = h.Vector().record(cell.ga_gtp_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_plc = h.Vector().record(cell.plc_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_active_plc = h.Vector().record(cell.ga_gtp_plc_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    #
    s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[0])[0]._ref_concentration)
    # a9_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)
    # s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[0])[0]._ref_concentration)
    # a9_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[9])[0]._ref_concentration)
    # s_ip3 = h.Vector().record(cell.ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[0])[0]._ref_concentration)
    # a9_ip3 = h.Vector().record(cell.ip3.nodes(cell.apical[9])[0]._ref_concentration)
    # s_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    # a9_ri = h.Vector().record(cell.ri_ip3r.nodes(cell.apical[9])[0]._ref_concentration)
    # s_po = h.Vector().record(cell.ro_ip3r.nodes(cell.somatic[0])[0]._ref_concentration)
    # a0_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[0])[0]._ref_concentration)
    # a9_po = h.Vector().record(cell.ro_ip3r.nodes(cell.apical[9])[0]._ref_concentration)
    #
    # s_ip5p = h.Vector().record(cell.ip5p.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_ip5p_ip3 = h.Vector().record(cell.ip5p_ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_ip3k = h.Vector().record(cell.ip3k.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_ip3k_2ca = h.Vector().record(cell.ip3k_2ca_ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_ip3k_2ca_ip3 = h.Vector().record(cell.ip3k_2ca_ip3.nodes(cell.somatic[0])[0]._ref_concentration)
    #
    # s_r = h.Vector().record(cell.r_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_rl = h.Vector().record(cell.rl_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_rlg = h.Vector().record(cell.rlg_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    #
    # s_g = h.Vector().record(cell.g_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_gbg = h.Vector().record(cell.gbg_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_rg = h.Vector().record(cell.rg_m1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_rgbg = h.Vector().record(cell.rgbg_m1.nodes(cell.somatic[0])[0]._ref_concentration)

    # s_r = h.Vector().record(cell.somatic[0](0.5)._ref_R_m1_kmb)
    # s_rl = h.Vector().record(cell.somatic[0](0.5)._ref_RLA_m1_kmb)
    # s_rlg = h.Vector().record(cell.somatic[0](0.5)._ref_RLGA_m1_kmb)
    #
    # s_g = h.Vector().record(cell.somatic[0](0.5)._ref_G_m1_kmb)
    # s_gbg = h.Vector().record(cell.somatic[0](0.5)._ref_Gbg_m1_kmb)
    # s_rg = h.Vector().record(cell.somatic[0](0.5)._ref_RG_m1_kmb)
    # s_rgbg = h.Vector().record(cell.somatic[0](0.5)._ref_RGbg_m1_kmb)

    # s_d1 = h.Vector().record(cell.dumb1.nodes(cell.somatic[0])[0]._ref_concentration)
    # s_d2 = h.Vector().record(cell.dumb2.nodes(cell.somatic[0])[0]._ref_concentration)

    # cyt_vals = np.zeros((t_steps.size, n_nodes))
    # er_vals = np.zeros((t_steps.size, n_nodes))
    # ip3_vals = np.zeros((t_steps.size, n_nodes))
    # ip3r_open_vals = np.zeros((t_steps.size, n_nodes))
    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    ax_v_dict = {}
    ap_v_dict = {}
    t_dict = {}
    ca_dict = {}
    isk_dict = {}
    im_dict = {}
    ibk_dict = {}
    ikdr_dict = {}
    ikap_dict = {}
    ina_dict = {}
    ih_dict = {}
    ipas_dict = {}
    ican_dict = {}
    ical_dict = {}
    icat_dict = {}
    ax_im_dict = {}
    ax_ikdr_dict = {}
    ax_ikap_dict = {}

    for a_i, ach_val in enumerate(ach_levels):

        h.stdinit()

        print('Running ACh pulse')

        cell.ach.nodes(cell.cyt).concentration = ach_val
        h.CVode().re_init()

        h.continuerun(run_time)

        s_v_dict[a_i] = np.array(s_v)
        ax_v_dict[a_i] = np.array(ax_v)
        ap_v_dict[a_i] = np.array(ap_v)
        t_dict[a_i] = np.array(t_vec)
        ca_dict[a_i] = np.array(s_ca_cyt)
        isk_dict[a_i] = np.array(s_isk)
        ax_im_dict[a_i] = np.array(ax_im)
        ax_ikdr_dict[a_i] = np.array(ax_ikdr)
        ax_ikap_dict[a_i] = np.array(ax_ikap)
        im_dict[a_i] = np.array(s_im)
        ibk_dict[a_i] = np.array(s_ibk)
        ikdr_dict[a_i] = np.array(s_ikdr)
        ikap_dict[a_i] = np.array(s_ikap)
        ipas_dict[a_i] = np.array(s_ipas)
        ina_dict[a_i] = np.array(s_ina)
        ih_dict[a_i] = np.array(s_ih)
        ical_dict[a_i] = np.array(s_ical)
        ican_dict[a_i] = np.array(s_ican)
        icat_dict[a_i] = np.array(s_icat)

    res_dict = {
                # 'ip3 vals': ip3_vals,
                # 'cyt vals': cyt_vals,
                # 'er vals': er_vals,
                # 'ip3r open vals': ip3r_open_vals,
                'soma v dict': s_v_dict,
                'axon v dict': ax_v_dict,
                'apical v dict': ap_v_dict,
                'ach levels': ach_levels,
                'ca cyt dict': ca_dict,
                # 'apical0 v': np.array(a0_v),
                # 'apical9 v': np.array(a9_v),
                # 'axon v': np.array(ax_v),
                # 'soma cyt time': np.array(s_ca_cyt),
                # 'apical0 cyt time': np.array(a0_ca_cyt),
                # 'apical9 cyt time': np.array(a9_ca_cyt),
                # 'soma er time': np.array(s_ca_er),
                # 'apical0 er time': np.array(a0_ca_er),
                # 'apical9 er time': np.array(a9_ca_er),
                # 'soma ip3 time': np.array(s_ip3),
                # 'apical0 ip3 time': np.array(a0_ip3),
                # 'apical9 ip3 time': np.array(a9_ip3),
                # 'soma ri': np.array(s_ri),
                # 'apical0 ri': np.array(a0_ri),
                # 'apical9 ri': np.array(a9_ri),
                # 'soma open probability': np.array(s_po),
                # 'apical0 open probability': np.array(a0_po),
                # 'apical9 open probability': np.array(a9_po),
                'soma ina dict': ina_dict,
                'soma ih dict': ih_dict,
                'soma im dict': im_dict,
                'soma ipas dict': ipas_dict,
                # 'axon im': np.array(ax_im),
                'soma isk dict': isk_dict,
                'soma ibk dict': ibk_dict,
                'soma ikdr dict': ikdr_dict,
                'soma ikap dict': ikap_dict,
                'axon im dict': ax_im_dict,
                'axon ikdr dict': ax_ikdr_dict,
                'axon ikap dict': ax_ikap_dict,
                'soma ical dict': ical_dict,
                'soma ican dict': ican_dict,
                'soma icat dict': icat_dict,

                # 'apical0 isk': np.array(a0_isk),
                # 'apical9 isk': np.array(a9_isk),
                # 'soma dag': np.array(s_dag),
                # 'soma pip2': np.array(s_pip2),
                # 'soma pi': np.array(s_pi),
                # 'soma pi4p': np.array(s_pi4p),
                # 'apical0 pip2': np.array(a0_pip2),
                # 'apical9 pip2': np.array(a9_pip2),
                # 'axon pip2': np.array(ax_pip2),
                # 'soma ach': np.array(s_ach),
                # 'soma pip2-kcnq': np.array(s_pip2_kcnq),
                # 'soma perc_i': np.array(s_perc_i),
                # 'axon perc_i': np.array(ax_perc_i),
                # 'soma ga-gtp': np.array(s_ga_gtp),
                # 'soma plc': np.array(s_plc),
                # 'soma active plc': np.array(s_active_plc),
                # 'soma r': np.array(s_r),
                # 'soma g': np.array(s_g),
                # 'soma gbg': np.array(s_gbg),
                # 'soma rg': np.array(s_rg),
                # 'soma rgbg': np.array(s_rgbg),
                # 'soma rl': np.array(s_rl),
                # 'soma rlg': np.array(s_rlg),
                # 'soma ip5p': np.array(s_ip5p),
                # 'soma ip5p_ip3': np.array(s_ip5p_ip3),
                # 'soma ip3k': np.array(s_ip3k),
                # 'soma ip3k_2ca': np.array(s_ip3k_2ca),
                # 'soma ip3k_2ca_ip3': np.array(s_ip3k_2ca_ip3),
                # 'soma dumb1': np.array(s_d1),
                # 'soma dumb2': np.array(s_d2),
                't': t_dict,
                }
    return res_dict


def plot_ahp_test(result_d, z_times=[1250, 1350], t_ignore=0):

    i_igs = {}
    t_igs = {}

    for a_i, ach_conc in enumerate(result_d['ach levels']):
        t_arr = result_d['t'][a_i]

        i_igs[a_i] = np.where(t_arr > t_ignore)
        t_igs[a_i] = t_arr[i_igs[a_i]] - t_ignore

    ret_figs = []

    v_fig = bplt.figure(title='Somatic Membrane Potential')
    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(v_fig)

    z_fig = bplt.figure(title='Zoom on AHP/ADP', x_range=z_times)
    z_fig.xaxis.axis_label = 'time (msec)'
    z_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(z_fig)

    ca_fig = bplt.figure(title='Cytosol Calcium')
    ca_fig.xaxis.axis_label = 'time (msec)'
    ca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ca_fig)

    ik_fig = bplt.figure(title='Potassium Currents', x_range=z_times, y_range=[0, 0.1])
    ik_fig.xaxis.axis_label = 'time (msec)'
    ik_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(ik_fig)

    ina_fig = bplt.figure(title='Sodium Currents', x_range=z_times)
    ina_fig.xaxis.axis_label = 'time (msec)'
    ina_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(ina_fig)

    ica_fig = bplt.figure(title='Calcium Currents', x_range=z_times)
    ica_fig.xaxis.axis_label = 'time (msec)'
    ica_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(ica_fig)

    iax_fig = bplt.figure(title='Currents in Axon', x_range=z_times)
    iax_fig.xaxis.axis_label = 'time (msec)'
    iax_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(iax_fig)

    itot_fig = bplt.figure(title='Total Soma Current', x_range=z_times, y_range=[-0.1, 0.1])
    itot_fig.xaxis.axis_label = 'time (msec)'
    itot_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(itot_fig)

    huh = ['solid', 'dashed', 'dotted']

    tot_dict = {}

    for a_i, ach_conc in enumerate(result_d['ach levels']):

        tot_dict[a_i] = result_d['soma isk dict'][a_i] + result_d['soma ibk dict'][a_i] + result_d['soma im dict'][a_i]\
                        + result_d['soma ikdr dict'][a_i] + result_d['soma ikap dict'][a_i] + \
                        result_d['soma ipas dict'][a_i] + result_d['soma ina dict'][a_i] + \
                        result_d['soma ih dict'][a_i] + result_d['soma ical dict'][a_i] + \
                        result_d['soma ican dict'][a_i] + result_d['soma icat dict'][a_i]

        v_fig.line(t_igs[a_i], result_d['soma v dict'][a_i][i_igs[a_i]], color=colrs[3*a_i],
                   legend='Soma ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        v_fig.line(t_igs[a_i], result_d['axon v dict'][a_i][i_igs[a_i]], color=colrs[3*a_i+1],
                   legend='Axon ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        v_fig.line(t_igs[a_i], result_d['apical v dict'][a_i][i_igs[a_i]], color=colrs[3*a_i+2],
                   legend='Apical ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])

        z_fig.line(t_igs[a_i], result_d['soma v dict'][a_i][i_igs[a_i]], color=colrs[2*a_i],
                   legend='ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        z_fig.line(t_igs[a_i], result_d['axon v dict'][a_i][i_igs[a_i]], color=colrs[2*a_i+1],
                   legend='ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])

        ca_fig.line(result_d['t'][a_i], 1000.0*result_d['ca cyt dict'][a_i], color=colrs[a_i],
                    legend='ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3)

        # ik_fig.line(result_d['t'][a_i], result_d['soma isk dict'][a_i],
        #             color=colrs[6 * a_i],
        #             legend='I_SK ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])

        ik_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma isk dict'][a_i][i_igs[a_i]], color=colrs[6*a_i],
                    legend='I_SK ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ik_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma ibk dict'][a_i][i_igs[a_i]], color=colrs[6 * a_i+1],
                    legend='I_BK ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ik_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma im dict'][a_i][i_igs[a_i]], color=colrs[6 * a_i+2],
                    legend='I_M ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ik_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma ikdr dict'][a_i][i_igs[a_i]], color=colrs[6 * a_i + 3],
                    legend='I_KDR ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ik_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma ikap dict'][a_i][i_igs[a_i]], color=colrs[6 * a_i + 4],
                    legend='I_KA ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ik_fig.line(result_d['t'][a_i][i_igs[a_i]] - t_ignore, result_d['soma ipas dict'][a_i][i_igs[a_i]],
                    color=colrs[6 * a_i + 5], legend='I_Passive ACh: {0:.2f} {1}M'.format(ach_conc, mu),
                    line_width=3, line_dash=huh[a_i])

        ina_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma ina dict'][a_i][i_igs[a_i]], color=colrs[2*a_i],
                     legend='I_Na ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ina_fig.line(result_d['t'][a_i][i_igs[a_i]] - t_ignore, result_d['soma ih dict'][a_i][i_igs[a_i]],
                     color=colrs[2*a_i + 1], legend='I_H ACh: {0:.2f} {1}M'.format(ach_conc, mu),
                     line_width=3, line_dash=huh[a_i])

        ica_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['soma ical dict'][a_i][i_igs[a_i]], color=colrs[3*a_i],
                     legend='I_CaL ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        ica_fig.line(result_d['t'][a_i][i_igs[a_i]] - t_ignore, result_d['soma ican dict'][a_i][i_igs[a_i]],
                     color=colrs[3*a_i + 1], legend='I_CaN ACh: {0:.2f} {1}M'.format(ach_conc, mu),
                     line_width=3, line_dash=huh[a_i])
        ica_fig.line(result_d['t'][a_i][i_igs[a_i]] - t_ignore, result_d['soma icat dict'][a_i][i_igs[a_i]],
                     color=colrs[3 * a_i + 2], legend='I_CaT ACh: {0:.2f} {1}M'.format(ach_conc, mu),
                     line_width=3, line_dash=huh[a_i])

        iax_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, result_d['axon ikdr dict'][a_i][i_igs[a_i]], color=colrs[3*a_i],
                     legend='I_KDR ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        iax_fig.line(result_d['t'][a_i][i_igs[a_i]] - t_ignore, result_d['axon im dict'][a_i][i_igs[a_i]],
                     color=colrs[3 * a_i + 1],
                     legend='I_KM ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])
        iax_fig.line(result_d['t'][a_i][i_igs[a_i]] - t_ignore, result_d['axon ikap dict'][a_i][i_igs[a_i]],
                     color=colrs[3 * a_i + 2],
                     legend='I_KA ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3, line_dash=huh[a_i])

        itot_fig.line(result_d['t'][a_i][i_igs[a_i]]-t_ignore, tot_dict[a_i][i_igs[a_i]], color=colrs[a_i],
                      legend='ACh: {0:.2f} {1}M'.format(ach_conc, mu), line_width=3)

    v_fig.legend.location = 'bottom_right'

    return ret_figs


def run_buffer_comparison(param_dict, calbindin_concs, ach_times, ach_concs=[100.0], run_time=4000):
    """
    param: param_dict: dictionary of values used to define the cell model
    param: calbindin_concs: array of concentrations of calbindin to simulate for
    param: ach_times: array of times that acetylcholine pulse will start and stop
    param: ach_conc: concentration of acetylcholine to be used in simulations
    param: run_time: length of simulated time in milliseconds

    return: res_dict: dictionary holding results for simulations

    """
    morph_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = ca1p.MyCell(morph_path, param_dict)
    cell.insert_rxd()

    s_v_vec = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    ca_cyt_vec = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    ca_er_vec = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)

    # Record Values
    t_vec = h.Vector().record(h._ref_t)

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    ca_cyt_dict = {}
    ca_er_dict = {}
    t_dict = {}

    ca_cyt_val = param_dict['ca_cyt_val']

    for a_i, ach_conc in enumerate(ach_concs):
        s_v_dict[a_i] = {}
        ca_cyt_dict[a_i] = {}
        ca_er_dict[a_i] = {}
        t_dict[a_i] = {}

        for c_i, cal_val in enumerate(calbindin_concs):
            print('Running for Calbindin = {0} mM'.format(cal_val))
            total_cbd = cal_val
            total_cbdh = total_cbd * 2.0  # 2 High Affinity Binding Sites
            total_cbdl = total_cbd * 2.0  # 2 Low Affinity Binding Sites

            h.stdinit()

            KD_cbdh = 0.000237
            kf_cbdh = 11.0
            cbdhca_init = total_cbdh * ca_cyt_val / (KD_cbdh + ca_cyt_val)
            cbdh_init = total_cbdh - cbdhca_init

            KD_cbdl = 0.000411
            kf_cbdl = 87.0
            cbdlca_init = total_cbdl * ca_cyt_val / (KD_cbdl + ca_cyt_val)
            cbdl_init = total_cbdl - cbdlca_init

            cell.cbdh.nodes(cell.cyt).concentration = cbdh_init
            cell.cbdhca.nodes(cell.cyt).concentration = cbdhca_init
            cell.cbdl.nodes(cell.cyt).concentration = cbdl_init
            cell.cbdlca.nodes(cell.cyt).concentration = cbdlca_init

            h.CVode().re_init()

            if ach_times[0]:
                h.continuerun(ach_times[0])

            cell.ach.nodes(cell.cyt).concentration = ach_conc * 0.001

            h.CVode().re_init()

            if ach_times[1] < run_time:
                h.continuerun(ach_times[1])

            cell.ach.nodes(cell.cyt).concentration = 0.0

            h.CVode().re_init()

            h.continuerun(run_time)

            s_v_dict[a_i][c_i] = np.array(s_v_vec)
            ca_cyt_dict[a_i][c_i] = np.array(ca_cyt_vec)
            ca_er_dict[a_i][c_i] = np.array(ca_er_vec)
            t_dict[a_i][c_i] = np.array(t_vec)

    res_dict = {
        'soma v dict': s_v_dict,
        'cytosol calcium': ca_cyt_dict,
        'endoplasmic calcium': ca_er_dict,
        'calbindin values': np.array(calbindin_concs),
        't': t_dict,
        'ach values': ach_concs,
        'ach times': ach_times,
    }

    return res_dict


def plot_buffer_comparison(result_d, t_ignore=0):

    # Make figure showing voltage for every synaptic delay value

    ret_figs = []

    ach_concs = result_d['ach values']

    ach_span_1 = bmod.Span(location=result_d['ach times'][0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_2 = bmod.Span(location=result_d['ach times'][1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    for a_i, ach_conc in enumerate(ach_concs):

        v_fig = bplt.figure(title='Somatic Potential, ACh={0} {1}M'.format(ach_conc, mu))
        v_fig.xaxis.axis_label = 'time (msec)'
        v_fig.yaxis.axis_label = 'potential (mV)'
        ret_figs.append(v_fig)
        v_fig.add_layout(ach_span_1)
        v_fig.add_layout(ach_span_2)

        ca_cyt_fig = bplt.figure(title='Cytosol Calcium vs Time, ACh={0} {1}M'.format(ach_conc, mu))
        ca_cyt_fig.xaxis.axis_label = 'time (msec)'
        ca_cyt_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
        ret_figs.append(ca_cyt_fig)
        ca_cyt_fig.add_layout(ach_span_1)
        ca_cyt_fig.add_layout(ach_span_2)

        ca_er_fig = bplt.figure(title='Endoplasmic Reticulum Calcium vs Time, ACh={0} {1}M'.format(ach_conc, mu))
        ca_er_fig.xaxis.axis_label = 'time (msec)'
        ca_er_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
        ret_figs.append(ca_er_fig)
        ca_er_fig.add_layout(ach_span_1)
        ca_er_fig.add_layout(ach_span_2)

        peak_fig = bplt.figure(title='Time to Peak Calcium vs Buffer, ACh={0} {1}M'.format(ach_conc, mu))
        peak_fig.xaxis.axis_label = 'concentration ({}M)'.format(mu)
        peak_fig.yaxis.axis_label = 'time to peak (msec)'
        ret_figs.append(peak_fig)

        int_fig = bplt.figure(title='Integrated Cytosol Calcium vs Buffer, ACh={0} {1}M'.format(ach_conc, mu))
        int_fig.xaxis.axis_label = 'concentration ({}M)'.format(mu)
        int_fig.yaxis.axis_label = 'total calcium'.format(mu)
        ret_figs.append(int_fig)

        peak_ts = np.zeros(result_d['calbindin values'].shape)

        for cbd_i, cbd_val in enumerate(result_d['calbindin values']):
            i_ig = np.squeeze(np.where(result_d['t'][a_i][cbd_i] > t_ignore))

            v_fig.line(result_d['t'][a_i][cbd_i][i_ig]-t_ignore, result_d['soma v dict'][a_i][cbd_i][i_ig],
                       line_width=3, color=colrs[2*cbd_i], legend='Calbindin = {0} {1}M'.format(cbd_val, mu))

            ca_cyt_fig.line(result_d['t'][a_i][cbd_i][i_ig] - t_ignore, 1000.0*result_d['cytosol calcium'][a_i][cbd_i][i_ig],
                            line_width=3, color=colrs[2*cbd_i], legend='Calbindin = {0} {1}M'.format(cbd_val, mu))

            ca_er_fig.line(result_d['t'][a_i][cbd_i][i_ig] - t_ignore, 1000.0*result_d['endoplasmic calcium'][a_i][cbd_i][i_ig],
                           line_width=3, color=colrs[2*cbd_i], legend='Calbindin = {0} {1}M'.format(cbd_val, mu))

            peak_ts[cbd_i] = result_d['t'][a_i][cbd_i][np.argmax(result_d['cytosol calcium'][a_i][cbd_i])]

        peak_fig.circle(result_d['calbindin values'], peak_ts - result_d['ach times'][0], size=12, color=colrs[0])
        peak_fig.line(result_d['calbindin values'], peak_ts - result_d['ach times'][0], line_width=3, color=colrs[0])

        # Integrate area under curves of cytosol calcium
        ca_int = np.zeros(result_d['calbindin values'].shape)

        for cbd_i, cbd_val in enumerate(result_d['calbindin values']):

            i_int = np.intersect1d(np.where(result_d['t'][a_i][cbd_i] > result_d['ach times'][0]),
                                   np.where(result_d['t'][a_i][cbd_i] < (result_d['ach times'][0] + 3000.0))
                                   )

            test = scin.trapz(result_d['cytosol calcium'][a_i][cbd_i][i_int],
                              result_d['t'][a_i][cbd_i][i_int], axis=0)
            ca_int[cbd_i] = test

        int_fig.circle(result_d['calbindin values'], ca_int, size=12, color=colrs[0])
        int_fig.line(result_d['calbindin values'], ca_int, line_width=3, color=colrs[0])

    return ret_figs


def run_calcium_wave(param_dict, time_values, ach_times, ach_conc, run_time=3000):
    all_t_vals = np.sort(np.array(time_values + ach_times))

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')
    cell = ca1p.MyCell(t_path, param_dict)
    cell.insert_rxd()

    line_dict = get_line_segs_calcium(cell)

    r_dict = {}
    r_dict['time values'] = time_values
    r_dict['ach conc'] = ach_conc

    t_vec = h.Vector().record(h._ref_t)
    s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[0])[0]._ref_concentration)
    a9_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)
    s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
    a0_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[0])[0]._ref_concentration)
    a9_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.apical[9])[0]._ref_concentration)

    calcium_arr = np.zeros((line_dict['number total segs'], len(time_values)))

    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    h.stdinit()

    print('Running ACh pulse')

    ach_high = False

    for t_i, t_step in enumerate(all_t_vals):
        h.continuerun(t_step)

        if t_step in ach_times:
            if ach_high:
                cell.ach.nodes(cell.cyt).concentration = 0.0
            else:
                cell.ach.nodes(cell.cyt).concentration = ach_conc
                ach_high = True

            h.CVode().re_init()
        if t_step in time_values:
            rec_i = time_values.index(t_step)
            cur_i = 0

            for s_i, sec in enumerate(cell.cal_list):
                s_n = sec.name()
                if 'axon' not in s_n:
                    for seg_i in range(line_dict[s_n]['nseg']):
                        calcium_arr[cur_i, rec_i] = cell.ca[cell.cyt].nodes(sec)[seg_i].concentration
                        cur_i += 1

    h.continuerun(run_time)

    r_dict['soma ca cyt'] = np.array(s_ca_cyt)
    r_dict['apical0 ca cyt'] = np.array(a0_ca_cyt)
    r_dict['apical9 ca cyt'] = np.array(a9_ca_cyt)
    r_dict['soma ca er'] = np.array(s_ca_er)
    r_dict['apical0 ca er'] = np.array(a0_ca_er)
    r_dict['apical9 ca er'] = np.array(a9_ca_er)
    r_dict['calcium array'] = calcium_arr
    r_dict['t'] = np.array(t_vec)
    return line_dict, r_dict


def get_line_segs_calcium(cell):

    lseg_dict = {}

    current_peak = 0

    n_total_segs = 0

    for s_name, sec in zip(cell.cal_names, cell.cal_list):
        nseg = sec.nseg
        n_total_segs += nseg
        length = sec.L
        npts = int(h.n3d(sec=sec))

        if not npts:
            print('{} has no 3d pts'.format(s_name))
            continue

        xs = np.zeros((nseg + 1,))
        ys = np.arange(0, length + 1, length / float(nseg)) + current_peak

        pts = np.vstack((xs, ys)).transpose()

        lseg_dict[s_name] = {}
        lseg_dict[s_name]['pts'] = {}

        for i in range(nseg):
            lseg_dict[s_name]['pts'][i] = pts[i:i+2, :]
        lseg_dict[s_name]['nseg'] = nseg
        lseg_dict[s_name]['diam'] = sec.diam
        lseg_dict[s_name]['seg values'] = np.arange(0.5 / nseg, 0.99, 1.0 / nseg)

        current_peak += length
        # temp_pts = []
        # for i in range(npts):
        #     temp_pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])
        # pt_arr = np.array(temp_pts)

        # if nseg == 1:
        #     lseg_dict[s_name]['pts'][0] = pt_arr
        #     lseg_dict[s_name]['seg values'] = np.array([0.5])
        # else:
        #     print('{0} has {1} segments'.format(s_name, nseg))
        #
        #     pt_dists = np.linalg.norm(np.diff(pt_arr, axis=0), axis=1)
        #     prop_dists = np.cumsum(pt_dists)/np.sum(pt_dists)
        #
        #     last_min = 0
        #
        #     lseg_dict[s_name]['seg values'] = np.arange(0.5/nseg, 0.99, 1.0/nseg)
        #
        #     new_pts = np.array(pt_arr)
        #
        #     for trans_i, trans_x in enumerate(np.arange(1.0/nseg, 0.99, 1.0/nseg)):
        #         t_i = np.searchsorted(prop_dists, trans_x)
        #         pt0 = pt_arr[t_i-1]
        #         pt1 = pt_arr[t_i]
        #
        #         p0 = prop_dists[t_i-1]
        #         p1 = prop_dists[t_i]
        #
        #         prop_seg = (trans_x - p0)/(p1 - p0)
        #         seg_vec = pt1 - pt0
        #         new_pt = pt0 + prop_seg*seg_vec
        #         new_pts = np.insert(new_pts, t_i+1, new_pt, axis=0)
        #         prop_dists = np.insert(prop_dists, t_i, trans_x)
        #
        #         lseg_dict[s_name]['pts'][trans_i] = np.array(new_pts[last_min:(t_i+1), :])
        #         last_min = t_i
        #         last_pt = new_pt
        #
        #     lseg_dict[s_name]['pts'][trans_i+1] = np.array(new_pts[last_min:, :])

    lseg_dict['number total segs'] = n_total_segs

    return lseg_dict


def plot_calcium_wave(result_d, line_dict, t_ignore=0):

    r_figs = []

    # for sec_name in line_dict:
        # if 'soma[0]' in sec_name:
        #     x0 = line_dict[sec_name]['pts'][0][0, 0]
        #     y0 = line_dict[sec_name]['pts'][0][0, 1]
        #     z0 = line_dict[sec_name]['pts'][0][0, 2]
        #     break

    cyt_t_fig = bplt.figure(title='Cytsosol Calcium vs Time')
    cyt_t_fig.xaxis.axis_label = 'time (msec)'
    cyt_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    r_figs.append(cyt_t_fig)

    i_ig = np.where(result_d['t'] > t_ignore)

    cyt_t_fig.line(result_d['t'][i_ig]-t_ignore, 1000.0*result_d['soma ca cyt'][i_ig], line_width=3, color=colrs[0],
                   legend='Soma')
    cyt_t_fig.line(result_d['t'][i_ig] - t_ignore, 1000.0*result_d['apical0 ca cyt'][i_ig], line_width=3, color=colrs[2],
                   legend='Proximal Dendrite')
    cyt_t_fig.line(result_d['t'][i_ig] - t_ignore, 1000.0*result_d['apical9 ca cyt'][i_ig], line_width=3, color=colrs[4],
                   legend='Distal Dendrite')


    er_t_fig = bplt.figure(title='ER Calcium vs Time')
    er_t_fig.xaxis.axis_label = 'time (msec)'
    er_t_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    r_figs.append(er_t_fig)

    er_t_fig.line(result_d['t'][i_ig] - t_ignore, 1000.0*result_d['soma ca er'][i_ig], line_width=3, color=colrs[0],
                   legend='Soma')
    er_t_fig.line(result_d['t'][i_ig] - t_ignore, 1000.0*result_d['apical0 ca er'][i_ig], line_width=3, color=colrs[2],
                   legend='Proximal Dendrite')
    er_t_fig.line(result_d['t'][i_ig] - t_ignore, 1000.0*result_d['apical9 ca er'][i_ig], line_width=3, color=colrs[4],
                   legend='Distal Dendrite')

    n_colors = 256
    my_cols = bpal.viridis(n_colors)
    calcium_vals = np.linspace(0.0, 2.0, n_colors - 1)

    color_map = bmod.LinearColorMapper(palette='Viridis256', low=0.0, high=2.0)
    color_bar = bmod.ColorBar(color_mapper=color_map, label_standoff=20, major_label_text_font_size='18pt',
                              border_line_color=None, location=(0, 0))

    for t_i, t_step in enumerate(result_d['time values']):

        cal_plt = bplt.figure(title='time = {0} msec'.format(t_step-t_ignore),
                              x_range=[-10, 10])
        cal_plt.plot_width = 400
        cal_plt.plot_height = 600
        cal_plt.xaxis.axis_label_text_font = 'arial'
        cal_plt.xaxis.axis_label_text_font_style = 'bold'
        cal_plt.xaxis.axis_label_text_font_size = '18pt'
        cal_plt.xaxis.major_label_text_font = 'arial'
        cal_plt.xaxis.major_label_text_font_size = "18pt"
        cal_plt.xaxis.major_label_text_font_style = 'bold'
        cal_plt.xaxis.axis_label = 'distance ({}m)'.format(mu)

        cal_plt.yaxis.axis_label_text_font = 'arial'
        cal_plt.yaxis.axis_label_text_font_style = 'bold'
        cal_plt.yaxis.axis_label_text_font_size = '18pt'
        cal_plt.yaxis.major_label_text_font = 'arial'
        cal_plt.yaxis.major_label_text_font_size = "18pt"
        cal_plt.yaxis.major_label_text_font_style = 'bold'
        cal_plt.yaxis.axis_label = 'distance ({}m)'.format(mu)

        cal_plt.title.text_font = 'arial'
        cal_plt.title.text_font_size = '18pt'
        cal_plt.title.text_font_style = 'bold'

        cal_plt.toolbar.logo = None
        cal_plt.toolbar_location = None
        cal_plt.background_fill_color = None
        cal_plt.border_fill_color = None
        cal_plt.outline_line_color = None

        r_figs.append(cal_plt)

        c_i = 0

        for s_i, sec_name in enumerate(line_dict):
            if sec_name != 'number total segs':
                for x_i, x_val in enumerate(line_dict[sec_name]['seg values']):

                    h = line_dict[sec_name]['pts'][x_i][1, 1] - line_dict[sec_name]['pts'][x_i][0, 1]
                    c_y = (line_dict[sec_name]['pts'][x_i][1, 1] - line_dict[sec_name]['pts'][x_i][0, 1])*0.5 + line_dict[sec_name]['pts'][x_i][0, 1]

                    cal_val = result_d['calcium array'][c_i, t_i]*1000.0

                    cal_plt.rect(x=line_dict[sec_name]['pts'][x_i][0,0],
                                 y=c_y,
                                 height=h,
                                 width=line_dict[sec_name]['diam'],
                                 color=my_cols[np.digitize(cal_val, calcium_vals)])

                    c_i += 1

                    # cal_plt.line(line_dict[sec_name]['pts'][x_i][:, 0] - x0, line_dict[sec_name]['pts'][x_i][:, 1] - y0,
                    #              line_width=line_dict[sec_name]['diam'], color=colrs[x_i])

    cal_plt.add_layout(color_bar, 'right')

    return r_figs


def run_tonic_test(param_dict={}, ach_levels=[0.0, 100.0], run_time=2000, ach_times=[],
                   current_dict={'amp': 0.0, 'start': 0.0, 'dur': 0.0}):

    if not ach_times:
        ach_times = [0, run_time]

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = ca1p.MyCell(t_path, param_dict)
    cell.insert_rxd()

    if current_dict['amp']:
        t_curr = h.IClamp(cell.somatic[0](0.5))
        t_curr.amp = current_dict['amp']
        t_curr.delay = current_dict['start']
        t_curr.dur = current_dict['dur']
        print(
            'Current Injection at t = {0} msec, amplitude = {1} nA'.format(current_dict['start'], current_dict['amp']))

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)

    s_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    # ax_v = h.Vector().record(cell.axonal[0](0.5)._ref_v)
    # ap_v = h.Vector().record(cell.apical[9](0.5)._ref_v)

    s_im = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kmb_inh)
    s_isk = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kca)

    s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)
    # s_pip2 = h.Vector().record(cell.pip2.nodes(cell.somatic[0])[0]._ref_concentration)
    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    # ax_v_dict = {}
    # ap_v_dict = {}
    t_dict = {}
    ca_dict = {}
    er_dict = {}
    isk_dict = {}
    im_dict = {}
    # pip2_dict = {}

    for a_i, ach_val in enumerate(ach_levels):

        h.stdinit()

        print('Running ACh pulse')

        if ach_times[0]:
            h.continuerun(ach_times[0])

        cell.ach.nodes(cell.cyt).concentration = ach_val*0.001
        h.CVode().re_init()

        if ach_times[1] < run_time:
            h.continuerun(ach_times[1])

            cell.ach.nodes(cell.cyt).concentration = 0.0
            h.CVode().re_init()

        h.continuerun(run_time)

        s_v_dict[a_i] = np.array(s_v)
        # ax_v_dict[a_i] = np.array(ax_v)
        # ap_v_dict[a_i] = np.array(ap_v)
        t_dict[a_i] = np.array(t_vec)
        ca_dict[a_i] = np.array(s_ca_cyt)
        er_dict[a_i] = np.array(s_ca_er)
        isk_dict[a_i] = np.array(s_isk)
        im_dict[a_i] = np.array(s_im)
        # pip2_dict[a_i] = np.array(s_pip2)

    res_dict = {
                'soma_v_dict': s_v_dict,
                # 'axon_v_dict': ax_v_dict,
                # 'apical_v_dict': ap_v_dict,
                'ach_levels': ach_levels,
                'ca_cyt_dict': ca_dict,
                'ca_er_dict': er_dict,
                'soma_im_dict': im_dict,
                'soma_isk_dict': isk_dict,
                # 'soma_pip2_dict': pip2_dict,
                't': t_dict,
                }

    return res_dict


def run_input_resistance_test(param_dict, ach_times, ach_levels=[0.0, 100.0], run_time=2000,
                              current_dict={'amps': [0.0], 'start': 0.0, 'dur': 0.0}):

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = ca1p.MyCell(t_path, param_dict)
    cell.insert_rxd()

    t_curr = h.IClamp(cell.somatic[0](0.5))
    t_curr.delay = current_dict['start']
    t_curr.dur = current_dict['dur']

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)

    s_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    s_im = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kmb_inh)
    s_isk = h.Vector().record(cell.somatic[0](0.5)._ref_ik_kca)
    s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
    s_ca_er = h.Vector().record(cell.ca[cell.er].nodes(cell.somatic[0])[0]._ref_concentration)

    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    t_dict = {}
    ca_dict = {}
    er_dict = {}
    im_dict = {}
    isk_dict = {}

    for a_i, ach_val in enumerate(ach_levels):
        print('Running ACh pulse')

        s_v_dict[a_i] = {}
        t_dict[a_i] = {}
        ca_dict[a_i] = {}
        er_dict[a_i] = {}
        im_dict[a_i] = {}
        isk_dict[a_i] = {}

        for c_i, current_amp in enumerate(current_dict['amps']):

            t_curr.amp = current_amp

            h.stdinit()

            if ach_times[0]:
                h.continuerun(ach_times[0])

            cell.ach.nodes(cell.cyt).concentration = ach_val*0.001
            h.CVode().re_init()

            if ach_times[1] < run_time:
                h.continuerun(ach_times[1])

                cell.ach.nodes(cell.cyt).concentration = 0.0
                h.CVode().re_init()

            h.continuerun(run_time)

            s_v_dict[a_i][c_i] = np.array(s_v)
            t_dict[a_i][c_i] = np.array(t_vec)
            ca_dict[a_i][c_i] = np.array(s_ca_cyt)
            er_dict[a_i][c_i] = np.array(s_ca_er)
            im_dict[a_i][c_i] = np.array(s_im)
            isk_dict[a_i][c_i] = np.array(s_isk)

    res_dict = {
                'soma_v_dict': s_v_dict,
                'ach_levels': np.array(ach_levels),
                'current_amplitudes': np.array(current_dict['amps']),
                'ca_cyt_dict': ca_dict,
                'ca_er_dict': er_dict,
                'soma_im_dict': im_dict,
                'soma_isk_dict': isk_dict,
                't': t_dict,
                }

    return res_dict


def plot_input_resistance_test(result_d, curr_times, ach_times, t_ignore=0):
    """

    return: ret_figs: list of figures
            slopes: numpy array containing the calculated input resistance values
    """
    dv_dict = {}
    slopes = np.zeros(len(result_d['ach_levels']))
    intercepts = np.zeros(len(result_d['ach_levels']))

    # Make figure showing voltage for every ach concentration

    for a_i, ach_conc in enumerate(result_d['ach_levels']):
        dv_dict[a_i] = np.zeros((len(result_d['current_amplitudes']),))

        for c_i, curr_amp in enumerate(result_d['current_amplitudes']):
            p_is = np.squeeze(np.intersect1d(np.where(result_d['t'][a_i][c_i] > curr_times[0]),
                                             np.where(result_d['t'][a_i][c_i] < curr_times[1])))


            dv = result_d['soma_v_dict'][a_i][c_i][p_is[-1]] - result_d['soma_v_dict'][a_i][c_i][p_is[0] - 1]
            dv_dict[a_i][c_i] = dv

        dv_matrix = np.vstack((np.array(result_d['current_amplitudes']), np.ones(len(result_d['current_amplitudes'])))).T

        print('ACh: {0}'.format(ach_conc))
        print(dv_matrix)
        print(dv_dict[a_i])

        slopes[a_i], intercepts[a_i] = np.linalg.lstsq(dv_matrix, dv_dict[a_i])[0]

    ret_figs = []

    spans1 = []
    spans2 = []

    ispans1 = []
    ispans2 = []

    for a_i, ach_conc in enumerate(result_d['ach_levels']):

        v_fig = bplt.figure(title='Somatic Membrane Potential, ACh={0:.3f} {1}M'.format(ach_conc, mu))
        v_fig.xaxis.axis_label = 'time (msec)'
        v_fig.yaxis.axis_label = 'potential (mV)'
        ret_figs.append(v_fig)

        spans1.append(bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3))
        spans2.append(bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3))
        v_fig.add_layout(spans1[-1])
        v_fig.add_layout(spans2[-1])

        i_fig = bplt.figure(title='Currents vs Time, ACh={0:.3f} {1}M'.format(ach_conc, mu))
        i_fig.xaxis.axis_label = 'time (msec)'
        i_fig.yaxis.axis_label = 'current (mA/cm^2)'
        ret_figs.append(i_fig)

        ispans1.append(bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3))
        ispans2.append(bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3))
        i_fig.add_layout(ispans1[-1])
        i_fig.add_layout(ispans2[-1])

        for c_i, c_amp in enumerate(result_d['current_amplitudes']):
            v_fig.line(result_d['t'][a_i][c_i], result_d['soma_v_dict'][a_i][c_i], color=colrs[c_i],
                       legend='I_inj: {0:.2f} nA'.format(c_amp), line_width=3)

            i_fig.line(result_d['t'][a_i][c_i], result_d['soma_isk_dict'][a_i][c_i], line_width=3, color=colrs[2 * c_i],
                       legend='Isk {0:.2f} nA'.format(c_amp))
            i_fig.line(result_d['t'][a_i][c_i], result_d['soma_im_dict'][a_i][c_i], line_width=3, color=colrs[2 * c_i + 1],
                       legend='Im {0:.2f} nA'.format(c_amp))

        v_fig.legend.location = 'top_right'

    dv_fig = bplt.figure(title='Depolarization vs Current Amplitude')
    dv_fig.xaxis.axis_label = 'current (nA)'
    dv_fig.yaxis.axis_label = 'depolarization (mV)'
    ret_figs.append(dv_fig)

    for a_i, ach_conc in enumerate(result_d['ach_levels']):
        dv_fig.circle(result_d['current_amplitudes'], dv_dict[a_i], size=12, color=colrs[a_i],
                      legend='ACh: {0:.3f} {1}M'.format(ach_conc, mu))

        lin_approx = slopes[a_i]*np.array(result_d['current_amplitudes']) + intercepts[a_i]

        dv_fig.line(result_d['current_amplitudes'], lin_approx, line_width=3, color=colrs[a_i])

    dv_fig.legend.location = 'bottom_right'

    rin_fig = bplt.figure(title='Input Resistance vs Tonic Acetylcholine'.format(slopes[0], ohm, mu),
                          x_axis_type='log')
    rin_fig.xaxis.axis_label = 'Acetylcholine Concentration ({}M)'.format(mu)
    rin_fig.yaxis.axis_label = 'Input Resistance (M{})'.format(ohm)
    ret_figs.append(rin_fig)

    # inp_opts = get_ach_curve_params(result_d['ach_levels'][1:], slopes[1:], positive=False)
    # print('Input Resistance EC50: ', 10 ** inp_opts[1])

    rin_fig.line(result_d['ach_levels'][1:], slopes[1:], color=colrs[0], line_width=3)
    rin_fig.circle(result_d['ach_levels'][1:], slopes[1:], color=colrs[0], size=12)

    return ret_figs, slopes


def run_epsp_sweep(in_cell, ach_times, stim_times, synapse_location, ach_conc=100.0, run_time=4000,
                   synapse_dict={'amplitude': 0.0, 'e': 0.0, 'tau1': 0.0, 'tau2': 10.0, 'delay': 1.0}):
    """
    param: in_cell: ca1 pyramidal model to be used in simulations
    param: ach_times: array of times that mark the beginning and end of the acetylcholine pulse
    param: stim_times: array of times that synapse will be activated
    param: synapse_location: dictionary of this structure
        synapse_location['section name'] = location
        where 'section name' is a string with the
    param: ach_conc: concentration of acetylcholine to be used in simulations
    param: run_time: length of simulated time in milliseconds
    param: synapse_dict: dictionary containing parameters for synapse

    return: res_dict: dictionary holding results for simulations

    """

    netcon_dict = {}
    syn_dict = {}

    vect = h.Vector()
    stim = h.VecStim()

    for sec_name in synapse_location:
        netcon_dict[sec_name] = {}
        syn_dict[sec_name] = {}

        if 'apic' in sec_name:
            sec_list = in_cell.apical
        elif 'dend' in sec_name:
            sec_list = in_cell.basal
        elif 'soma' in sec_name:
            sec_list = in_cell.somatic
        else:
            print('Poop: {0}'.format(sec_name))

        sec_num = int(sec_name.rstrip(']').split('[')[-1])

        netcon_dict[sec_name] = {}
        syn_dict[sec_name] = {}

        for loc in synapse_location[sec_name]:

            syn_dict[sec_name][loc] = h.Exp2Syn(loc, sec=sec_list[sec_num])
            syn_dict[sec_name][loc].tau1 = synapse_dict['tau1']
            syn_dict[sec_name][loc].tau2 = synapse_dict['tau2']
            syn_dict[sec_name][loc].e = synapse_dict['e']

            netcon_dict[sec_name][loc] = h.NetCon(stim,
                                                  syn_dict[sec_name][loc],
                                                  1,
                                                  synapse_dict['delay'],
                                                  synapse_dict['amplitude'],
                                                  sec=sec_list[sec_num]
                                                  )

            last_syn = [sec_name, loc]

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)
    s_v = h.Vector().record(in_cell.somatic[0](0.5)._ref_v)

    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    t_dict = {}
    dep_list = []
    dep_dict = {}

    n_delays = len(stim_times)

    for d_i, delay_val in enumerate(stim_times):

        vect.from_python([delay_val])
        stim.play(vect)

        print('{0} of {1} delay values'.format(d_i+1, n_delays))

        s_v_dict[d_i] = {}
        t_dict[d_i] = {}
        dep_dict[d_i] = {}

        h.stdinit()

        if ach_times[0]:
            h.continuerun(ach_times[0])

        in_cell.ach.nodes(in_cell.cyt).concentration = ach_conc * 0.001
        h.CVode().re_init()

        if ach_times[1] < run_time:
            h.continuerun(ach_times[1])

            in_cell.ach.nodes(in_cell.cyt).concentration = 0.0
            h.CVode().re_init()

        h.continuerun(run_time)

        s_v_arr = np.array(s_v)
        t_arr = np.array(t_vec)

        s_v_dict[d_i] = s_v_arr
        t_dict[d_i] = t_arr

        dep_i = np.intersect1d(np.where(t_arr >= delay_val),
                               np.where(t_arr <= (delay_val+50.0))
                               )

        dep_list.append(np.max(s_v_arr[dep_i]) - s_v_arr[dep_i[0]-1])

    res_dict = {
        'soma v dict': s_v_dict,
        't': t_dict,
        'ach val': ach_conc,
        'ach times': ach_times,
        'delay vals': stim_times,
        'depolarization values': np.array(dep_list),
    }

    return res_dict


def plot_epsp_sweep(result_d, t_ignore=0):

    # Make figure showing voltage for every synaptic delay value

    ret_figs = []

    ach_span_1 = bmod.Span(location=result_d['ach times'][0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_2 = bmod.Span(location=result_d['ach times'][1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    dv_vals = result_d['depolarization values']

    v_fig = bplt.figure(title='Somatic Potential, ACh = {0} {1}M'.format(result_d['ach val'], mu))
    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(v_fig)

    for d_i, delay_val in enumerate(result_d['delay vals']):
        i_ig = np.where(result_d['t'][d_i] > t_ignore)
        v_fig.line(result_d['t'][d_i][i_ig]-t_ignore, result_d['soma v dict'][d_i][i_ig], line_width=3, color=colrs[d_i])

    delays = np.array(result_d['delay vals']) - result_d['ach times'][0]
    dv_fig = bplt.figure(title='Depolarization vs ACh Pulse Time')
    dv_fig.xaxis.axis_label = 'time of ACh'
    dv_fig.yaxis.axis_label = 'depolarization (mV)'

    dv_fig.line(delays, dv_vals, line_width=3, color=colrs[0])
    dv_fig.circle(delays, dv_vals, size=12, color=colrs[0])

    ret_figs.append(dv_fig)

    p_diff = (dv_vals - dv_vals[0])/dv_vals[0]*100.0

    perc_fig = bplt.figure(title='Percent Change in Depolarization vs ACh Pulse Time')
    perc_fig.xaxis.axis_label = 'time of ACh'
    perc_fig.yaxis.axis_label = 'percent change (%)'

    perc_fig.line(delays, p_diff, line_width=3, color=colrs[0])
    perc_fig.circle(delays, p_diff, size=12, color=colrs[0])

    ret_figs.append(perc_fig)

    return ret_figs


def run_epsp_test(in_cell, ach_times, synapse_locations, ach_levels=[0.0, 100.0], run_time=4000, syn_delay =3000,
                  synapse_dict={'amplitude': 0.0, 'e': 0.0, 'tau1': 0.0, 'tau2': 10.0, 'delay': 1.0}):

    vect = h.Vector()
    t_list = [syn_delay]
    vect.from_python(t_list)

    stim = h.VecStim()
    stim.play(vect)

    netcon_dict = {}
    syn_dict = {}

    for sec_name in synapse_locations:
        netcon_dict[sec_name] = {}
        syn_dict[sec_name] = {}

        if 'apic' in sec_name:
            sec_list = in_cell.apical
        elif 'dend' in sec_name:
            sec_list = in_cell.basal
        elif 'soma' in sec_name:
            sec_list = in_cell.somatic
        else:
            print('Poop: {0}'.format(sec_name))

        sec_num = int(sec_name.rstrip(']').split('[')[-1])

        netcon_dict[sec_name] = {}
        syn_dict[sec_name] = {}

        for loc in synapse_locations[sec_name]:

            syn_dict[sec_name][loc] = h.Exp2Syn(loc, sec=sec_list[sec_num])
            syn_dict[sec_name][loc].tau1 = synapse_dict['tau1']
            syn_dict[sec_name][loc].tau2 = synapse_dict['tau2']
            syn_dict[sec_name][loc].e = synapse_dict['e']

            netcon_dict[sec_name][loc] = h.NetCon(stim,
                                                  syn_dict[sec_name][loc],
                                                  1,
                                                  synapse_dict['delay'],
                                                  0.0,
                                                  sec=sec_list[sec_num]
                                                  )

            last_syn = [sec_name, loc]

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)
    s_v = h.Vector().record(in_cell.somatic[0](0.5)._ref_v)

    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    t_dict = {}
    dep_dict = {}

    min_dv = 100.0
    max_dv = 0.0

    n_levels = len(ach_levels)

    for a_i, ach_val in enumerate(ach_levels):
        print('{0} of {1} ACh Levels'.format(a_i+1, n_levels))

        s_v_dict[a_i] = {}
        t_dict[a_i] = {}
        dep_dict[a_i] = {}

        for sec_name in synapse_locations:

            s_v_dict[a_i][sec_name] = {}
            t_dict[a_i][sec_name] = {}
            dep_dict[a_i][sec_name] = {}

            s_v_dict[a_i][sec_name] = {}
            t_dict[a_i][sec_name] = {}
            dep_dict[a_i][sec_name] = {}

            for loc in synapse_locations[sec_name]:

                netcon_dict[last_syn[0]][last_syn[1]].weight[0] = 0.0
                netcon_dict[sec_name][loc].weight[0] = synapse_dict['amplitude']
                last_syn = [sec_name, loc]

                h.stdinit()

                if ach_times[0]:
                    h.continuerun(ach_times[0])

                in_cell.ach.nodes(in_cell.cyt).concentration = ach_val * 0.001
                h.CVode().re_init()

                if ach_times[1] < run_time:
                    h.continuerun(ach_times[1])

                    in_cell.ach.nodes(in_cell.cyt).concentration = 0.0
                    h.CVode().re_init()

                h.continuerun(run_time)

                # Turn off last synapse
                # Activate Appropriate synapse (set amp to nonzero value)

                # Run Simulation
                # Measure depolarization at synapse
                s_v_arr = np.array(s_v)
                t_arr = np.array(t_vec)

                s_v_dict[a_i][sec_name][loc] = s_v_arr
                t_dict[a_i][sec_name][loc] = t_arr

                dep_i = np.intersect1d(np.where(t_arr >= syn_delay),
                                       np.where(t_arr <= (syn_delay+500.0))
                                       )

                dep_dict[a_i][sec_name][loc] = np.max(s_v_arr[dep_i]) - s_v_arr[dep_i[0]-1]

                if dep_dict[a_i][sec_name][loc] > max_dv:
                    max_dv = dep_dict[a_i][sec_name][loc]
                if dep_dict[a_i][sec_name][loc] < min_dv:
                    min_dv = dep_dict[a_i][sec_name][loc]

    res_dict = {
        'soma v dict': s_v_dict,
        't': t_dict,
        'ach levels': ach_levels,
        'depolarization values': dep_dict,
        'min dV': min_dv,
        'max dV': max_dv,
    }

    return res_dict


def plot_epsp_test(result_d, line_dict, ach_times, t_ignore=0):

    # Make figure showing voltage for every ach concentration

    ret_figs = []

    ach_span_1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    dv_dict = result_d['depolarization values']

    n_colors = 31
    my_cols = bpal.viridis(n_colors)
    dv_vals = np.linspace(result_d['min dV'], result_d['max dV'], n_colors-1)

    dV_list = []

    color_map = bmod.LinearColorMapper(palette='Viridis256', low=result_d['min dV'], high=result_d['max dV'])
    color_bar = bmod.ColorBar(color_mapper=color_map, label_standoff=20, major_label_text_font_size='18pt',
                              border_line_color=None, location=(0, 0))

    ticker = bmod.SingleIntervalTicker(interval=100)

    for a_i, ach_conc in enumerate(result_d['ach levels']):
        temp_list = []

        dv_fig = bplt.figure(title='EPSP at Soma'.format(ach_conc, mu))
        ret_figs.append(dv_fig)
        dv_fig.xaxis.ticker = ticker
        dv_fig.xaxis.axis_label = 'time (msec)'
        dv_fig.plot_width = 400
        dv_fig.plot_height = 600
        dv_fig.xaxis.axis_label_text_font = 'arial'
        dv_fig.xaxis.axis_label_text_font_style = 'bold'
        dv_fig.xaxis.axis_label_text_font_size = '22pt'
        dv_fig.xaxis.major_label_text_font = 'arial'
        dv_fig.xaxis.major_label_text_font_size = "18pt"
        dv_fig.xaxis.major_label_text_font_style = 'bold'
        dv_fig.xaxis.axis_label = 'distance ({}m)'.format(mu)
        dv_fig.add_layout(color_bar, 'right')

        dv_fig.yaxis.axis_label_text_font = 'arial'
        dv_fig.yaxis.axis_label_text_font_style = 'bold'
        dv_fig.yaxis.axis_label_text_font_size = '22pt'
        dv_fig.yaxis.major_label_text_font = 'arial'
        dv_fig.yaxis.major_label_text_font_size = "18pt"
        dv_fig.yaxis.major_label_text_font_style = 'bold'
        dv_fig.yaxis.axis_label = 'distance ({}m)'.format(mu)

        dv_fig.title.text_font = 'arial'
        dv_fig.title.text_font_size = '24pt'
        dv_fig.title.text_font_style = 'bold'

        dv_fig.toolbar.logo = None
        dv_fig.toolbar_location = None

        print(len(line_dict.keys()))

        for sec_name in line_dict:
            if 'soma[0]' in sec_name:
                x0 = line_dict[sec_name]['pts'][0, 0]
                y0 = line_dict[sec_name]['pts'][0, 1]
                break

        if not x0:
            x0, y0 = 0.0, 0.0

        for sec_name in dv_dict[a_i]:
            for loc_val in dv_dict[a_i][sec_name]:
                dv_val = dv_dict[a_i][sec_name][loc_val]
                temp_list.append(dv_val)
                break

            dv_fig.line(line_dict[sec_name]['pts'][:, 0] - x0, line_dict[sec_name]['pts'][:, 1] - y0,
                        line_width=line_dict[sec_name]['diam'],
                        color=my_cols[np.digitize(dv_val, dv_vals)])

        dV_list.append(np.array(temp_list))

    dv_dot_fig = bplt.figure(title='EPSP Amplitude vs [ACh]', x_axis_type='log')
    dv_dot_fig.xaxis.axis_label = 'concentration ({}M)'.format(mu)
    dv_dot_fig.yaxis.axis_label = 'depolarization (mV)'
    ret_figs.append(dv_dot_fig)

    dv_soma_v_fig = bplt.figure(title='Soma EPSP vs Time')
    dv_soma_v_fig.xaxis.axis_label = 'time (msec)'
    dv_soma_v_fig.yaxis.axis_label = 'membrane potential (mV)'
    ret_figs.append(dv_soma_v_fig)

    sec_names = list(result_d['soma v dict'][0].keys())
    locs = list(result_d['soma v dict'][0][sec_names[0]].keys())

    for a_i, a_val in enumerate(result_d['ach levels']):

        dv_soma_v_fig.line(result_d['t'][a_i][sec_names[0]][locs[0]],
                           result_d['soma v dict'][a_i][sec_names[0]][locs[0]],
                           line_width=3, color=colrs[a_i], legend='[ACh]: {0:.2f} ({1}M)'.format(a_val, mu))

        if a_i:
            dv_dot_fig.circle(a_val*np.ones((len(dV_list[a_i]),)), dV_list[a_i], color=colrs[a_i-1], size=12)

    return ret_figs


def run_rheobase_test(param_dict, ach_times, ach_levels=[0.0, 100.0], run_time=2000, starting_int=[0.4, 0.5],
                       desired_int=0.01, current_dict={'amp': [0.0], 'start': 0.0, 'dur': 0.0}):
    def calc_half(x, y):
        return (x + y) / 2.0

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = ca1p.MyCell(t_path, param_dict)
    cell.insert_rxd()

    t_curr = h.IClamp(cell.somatic[0](0.5))
    t_curr.delay = current_dict['start']
    t_curr.dur = current_dict['dur']

    # Record Values
    print('Creating Vectors for recording variables')
    t_vec = h.Vector().record(h._ref_t)
    s_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)

    print('Vectors created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    t_dict = {}
    current_dict = {}
    thresh_vals = []

    for a_i, ach_val in enumerate(ach_levels):
        s_v_dict[a_i] = {}
        t_dict[a_i] = {}
        current_dict[a_i] = []

        print('Searching for rheobase ACh pulse: {0} {1}M'.format(ach_val, mu))

        ### Test if High edge of current range causes spikes

        t_curr.amp = starting_int[1]

        h.stdinit()

        if ach_times[0]:
            h.continuerun(ach_times[0])

        cell.ach.nodes(cell.cyt).concentration = ach_val * 0.001
        h.CVode().re_init()

        if ach_times[1] < run_time:
            h.continuerun(ach_times[1])

            cell.ach.nodes(cell.cyt).concentration = 0.0
            h.CVode().re_init()

        h.continuerun(run_time)

        s_v_dict[a_i]['h'] = np.array(s_v)
        t_dict[a_i]['h'] = np.array(t_vec)

        current_dict[a_i].append(starting_int[1])
        h_spikes = fan.get_spike_times(s_v_dict[a_i]['h'], t_dict[a_i]['h'])

        h_b = len(h_spikes) > 0

        # Test if Low edge of current range does not cause spikes

        t_curr.amp = starting_int[0]

        h.stdinit()

        if ach_times[0]:
            h.continuerun(ach_times[0])

        cell.ach.nodes(cell.cyt).concentration = ach_val * 0.001
        h.CVode().re_init()

        if ach_times[1] < run_time:
            h.continuerun(ach_times[1])

            cell.ach.nodes(cell.cyt).concentration = 0.0
            h.CVode().re_init()

        h.continuerun(run_time)

        s_v_dict[a_i]['l'] = np.array(s_v)
        t_dict[a_i]['l'] = np.array(t_vec)

        current_dict[a_i].append(starting_int[0])
        l_spikes = fan.get_spike_times(s_v_dict[a_i]['l'], t_dict[a_i]['l'])

        l_b = len(l_spikes) > 0

        if h_b == l_b:
            print('Interval does not contain rheobase.')
            if l_b:
                print('Both edges caused spikes')
            else:
                print('Neither edge caused spikes')

            thresh_vals.append(0.0)
            continue

        b_int = [l_b, h_b]
        work_int = list(starting_int)
        int_size = work_int[1] - work_int[0]

        c_var = 2

        while int_size > desired_int:
            new_val = calc_half(work_int[0], work_int[1])

            current_dict[a_i].append(new_val)

            t_curr.amp = new_val

            h.stdinit()

            if ach_times[0]:
                h.continuerun(ach_times[0])

            cell.ach.nodes(cell.cyt).concentration = ach_val * 0.001
            h.CVode().re_init()

            if ach_times[1] < run_time:
                h.continuerun(ach_times[1])

                cell.ach.nodes(cell.cyt).concentration = 0.0
                h.CVode().re_init()

            h.continuerun(run_time)

            s_v_dict[a_i][c_var] = np.array(s_v)
            t_dict[a_i][c_var] = np.array(t_vec)

            m_spikes = fan.get_spike_times(s_v_dict[a_i][c_var], t_dict[a_i][c_var])

            h_m = len(m_spikes) > 0
            work_int[b_int.index(h_m)] = new_val
            int_size = work_int[1] - work_int[0]
            c_var += 1

        thresh_vals.append(work_int[1])
        print('Rheobase Found: {0} nA'.format(work_int[1]))

    res_dict = {
        'soma_v_dict': s_v_dict,
        't': t_dict,
        'current_dict': current_dict,
        'ach_levels': np.array(ach_levels),
        'rheobase_values': np.array(thresh_vals),
    }

    return res_dict


def plot_cell_bokeh(line_dict, color_dict={}, title='bokeh test'):

    dv_fig = bplt.figure(title=title, match_aspect=True)
    dv_fig.plot_width = 400
    dv_fig.plot_height = 600
    dv_fig.xaxis.axis_label_text_font = 'arial'
    dv_fig.xaxis.axis_label_text_font_style = 'bold'
    dv_fig.xaxis.axis_label_text_font_size = '22pt'
    dv_fig.xaxis.major_label_text_font = 'arial'
    dv_fig.xaxis.major_label_text_font_size = "18pt"
    dv_fig.xaxis.major_label_text_font_style = 'bold'
    dv_fig.xaxis.axis_label = 'distance ({}m)'.format(mu)

    dv_fig.yaxis.axis_label_text_font = 'arial'
    dv_fig.yaxis.axis_label_text_font_style = 'bold'
    dv_fig.yaxis.axis_label_text_font_size = '22pt'
    dv_fig.yaxis.major_label_text_font = 'arial'
    dv_fig.yaxis.major_label_text_font_size = "18pt"
    dv_fig.yaxis.major_label_text_font_style = 'bold'
    dv_fig.yaxis.axis_label = 'distance ({}m)'.format(mu)

    dv_fig.title.text_font = 'arial'
    dv_fig.title.text_font_size = '24pt'
    dv_fig.title.text_font_style = 'bold'

    dv_fig.toolbar.logo = None
    dv_fig.toolbar_location = None
    dv_fig.background_fill_color = None
    dv_fig.border_fill_color = None
    dv_fig.outline_line_color = None

    print(len(line_dict.keys()))

    for sec_name in line_dict:

        if 'soma[0]' in sec_name:
            x0 = line_dict[sec_name]['pts'][0, 0]
            y0 = line_dict[sec_name]['pts'][0, 1]
            break

    if not x0:
        x0, y0 = 0.0, 0.0

    for sec_name in line_dict:

        if color_dict:
            my_col = color_dict[sec_name]
        else:
            my_col=colrs[0]

        dv_fig.line(line_dict[sec_name]['pts'][:, 0] - x0, line_dict[sec_name]['pts'][:, 1] - y0,
                    line_width=line_dict[sec_name]['diam'], color=my_col)

    return dv_fig


def plot_rheobase_test(result_d, ach_times, curr_times, t_ignore=0):

    slopes = np.zeros(len(result_d['ach_levels']))
    intercepts = np.zeros(len(result_d['ach_levels']))

    # Make figure showing voltage for every ach concentration

    ret_figs = []

    ach_span_1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    z_times = [curr_times[0]-100, curr_times[1]+100]

    for a_i, ach_conc in enumerate(result_d['ach_levels']):

        v_fig = bplt.figure(title='Somatic Membrane Potential, ACh: {0:.3f} {1}M'.format(ach_conc, mu),
                            x_range=z_times)
        v_fig.xaxis.axis_label = 'time (msec)'
        v_fig.yaxis.axis_label = 'potential (mV)'
        ret_figs.append(v_fig)
        v_fig.legend.location = 'bottom_right'
        # v_fig.add_layout(ach_span_1)
        # v_fig.add_layout(ach_span_2)

        arg_sort = np.argsort(result_d['current_dict'][a_i])

        for c_i in arg_sort:
            c_amp = result_d['current_dict'][a_i][c_i]

            if c_i == 0:
                c_val = 'h'
            elif c_i == 1:
                c_val = 'l'
            else:
                c_val = c_i

            v_fig.line(result_d['t'][a_i][c_val], result_d['soma_v_dict'][a_i][c_val], color=colrs[c_i],
                       legend='{0:.3f} nA'.format(c_amp), line_width=3)

    thresh_fig = bplt.figure(title='Action Potential Rheobase vs [ACh]', x_axis_type='log')
    thresh_fig.xaxis.axis_label = '[ACh] ({}M)'.format(mu)
    thresh_fig.yaxis.axis_label = 'rheobase (nA)'

    rheo_opts = get_ach_curve_params(result_d['ach_levels'][1:], result_d['rheobase_values'][1:], positive=False)
    print('Rheobase EC50: ', 10 ** rheo_opts[1])

    ret_figs.append(thresh_fig)

    thresh_fig.line(result_d['ach_levels'][1:], result_d['rheobase_values'][1:], color=colrs[0], line_width=3)
    thresh_fig.circle(result_d['ach_levels'][1:], result_d['rheobase_values'][1:], color=colrs[0], size=12)

    return ret_figs


def plot_spike_acc_test(result_d, ach_times=[0, 50], t_ignore=0, plot_inds=[]):

    spike_dict = {}
    ifr_dict = {}
    acc_dict = {}
    t_dict = {}
    peak_acc = []
    peak_ca_vals = []
    max_hault = []

    for a_i, ach_conc in enumerate(result_d['ach_levels']):
        if a_i > 0:
            peak_ca_vals.append(np.max(result_d['ca_cyt_dict'][a_i]))

        spike_dict[a_i] = fan.get_spike_times(result_d['soma_v_dict'][a_i], result_d['t'][a_i])
        i_ig = np.squeeze(np.where(spike_dict[a_i] > t_ignore))

        if spike_dict[a_i][i_ig].size > 2:

            ifr_dict[a_i] = np.divide(1000.0, np.diff(spike_dict[a_i][i_ig]))

            t_dict[a_i] = (spike_dict[a_i][i_ig[1:]] + spike_dict[a_i][i_ig[:-1]])/2.0 - t_ignore

            if a_i == 0:
                normal_ifr = np.mean(ifr_dict[a_i])
            else:
                acc_dict[a_i] = 100.0*(ifr_dict[a_i] - normal_ifr)/normal_ifr
                peak_acc.append(np.max(acc_dict[a_i]))
                max_hault.append(1.0 / np.min(ifr_dict[a_i]))

        else:
            print('Not Enough Spikes in ACh = {0} {1}M'.format(ach_conc, mu))
            ifr_dict[a_i] = np.array([0.0])
            t_dict[a_i] = np.array([0.0])

    peak_acc = np.array(peak_acc)
    max_hault = np.array(max_hault)

    ret_figs = []

    ach_span_v1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_v2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    freq_fig = bplt.figure(title='Instantaneous Firing Rate')
    freq_fig.xaxis.axis_label = 'time (msec)'
    freq_fig.yaxis.axis_label = 'IFR (Hz)'
    ret_figs.append(freq_fig)
    freq_fig.legend.location = 'bottom_right'
    freq_fig.add_layout(ach_span_v1)
    freq_fig.add_layout(ach_span_v2)

    ach_span_a1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)
    ach_span_a2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)

    acc_fig = bplt.figure(title='Spike Acceleration')
    acc_fig.xaxis.axis_label = 'time (msec)'
    acc_fig.yaxis.axis_label = 'acceleration (%)'
    acc_fig.legend.location = 'bottom_right'
    ret_figs.append(acc_fig)
    acc_fig.add_layout(ach_span_a1)
    acc_fig.add_layout(ach_span_a2)

    acc_ach_fig = bplt.figure(title='Peak Spike Acceleration vs [ACh]', x_axis_type='log')
    acc_ach_fig.xaxis.axis_label = '[ACh] ({}M)'.format(mu)
    acc_ach_fig.yaxis.axis_label = 'acceleration (%)'
    ret_figs.append(acc_ach_fig)

    hault_ach_fig = bplt.figure(title='Duration of Spike Inhibition vs [ACh]', x_axis_type='log')
    hault_ach_fig.xaxis.axis_label = '[ACh] ({}M)'.format(mu)
    hault_ach_fig.yaxis.axis_label = '(seconds)'
    ret_figs.append(hault_ach_fig)

    ach_span_ca1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)
    ach_span_ca2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)

    ca_fig = bplt.figure(title='Somatic Cytosol Calcium')
    ca_fig.xaxis.axis_label = 'time (msec)'
    ca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ca_fig)
    ca_fig.add_layout(ach_span_ca1)
    ca_fig.add_layout(ach_span_ca2)

    ach_span_er1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)
    ach_span_er2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)

    er_fig = bplt.figure(title='Somatic Endoplasmic Reticulum Calcium')
    er_fig.xaxis.axis_label = 'time (msec)'
    er_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(er_fig)
    er_fig.add_layout(ach_span_er1)
    er_fig.add_layout(ach_span_er2)

    ach_span_ik1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)
    ach_span_ik2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)

    ik_fig = bplt.figure(title='Potassium Currents')
    ik_fig.xaxis.axis_label = 'time (msec)'
    ik_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(ik_fig)
    ik_fig.add_layout(ach_span_ik1)
    ik_fig.add_layout(ach_span_ik2)

    pca_fig = bplt.figure(title='Peak Cytosol Calcium vs ACh', x_axis_type='log')
    pca_fig.xaxis.axis_label = 'ACh ({}M)'.format(mu)
    pca_fig.yaxis.axis_label = '[Ca] ({}M)'.format(mu)
    ret_figs.append(pca_fig)

    freq_list = []
    acc_list = []
    cyt_list = []
    er_list = []
    ik_list = []

    if plot_inds:
        my_is = plot_inds
    else:
        my_is = range(len(result_d['ach_levels']))

    spans1 = []
    spans2 = []

    dashes = ['solid', 'solid', 'dashed', 'dashed', 'dotted', 'dotted', 'dashdot', 'dashdot']
    markers = ['circle', 'triangle', 'diamond', 'hex', 'inverted_triangle', 'square', 'asterisk']

    for i_i, a_i in enumerate(my_is):
        ach_conc = result_d['ach_levels'][a_i]

        v_fig = bplt.figure(title='Somatic Membrane Potential, ACh = {0:.3f} {1}M'.format(ach_conc, mu))

        v_fig.plot_width = 950
        v_fig.yaxis.ticker = [-70.0, 0.0, 40.0]
        v_fig.yaxis.axis_label = '(mV)'
        v_fig.xgrid.grid_line_color = None
        v_fig.ygrid.grid_line_color = None
        ret_figs.append(v_fig)
        v_fig.legend.location = 'bottom_right'

        if ach_conc > 0.0:
            spans1.append(bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                                    line_dash='dashed', line_width=3))
            spans2.append(bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                                    line_dash='dashed', line_width=3))

            v_fig.add_layout(spans1[-1])
            v_fig.add_layout(spans2[-1])

        i_ig = np.squeeze(np.where(result_d['t'][a_i] > t_ignore))
        t_arr = np.array(result_d['t'][a_i][i_ig]) - t_ignore

        v_fig.line(t_arr, result_d['soma_v_dict'][a_i][i_ig], color='black', line_width=1)

        freq_fig.line(t_dict[a_i], ifr_dict[a_i], color=colrs[i_i], line_width=1)
        f_scat = freq_fig.scatter(t_dict[a_i], ifr_dict[a_i], marker=markers[i_i], color=colrs[i_i], size=12, alpha=0.5)
        freq_list.append(('{0:.3f} {1}M'.format(ach_conc, mu), [f_scat]))

        cyt_line = ca_fig.line(t_arr, 1000.0*result_d['ca_cyt_dict'][a_i][i_ig], color=colrs[i_i], line_width=1,
                               line_dash=dashes[i_i])
        cyt_list.append(('{0:.3f} {1}M'.format(ach_conc, mu), [cyt_line]))

        er_line = er_fig.line(t_arr, 1000.0*result_d['ca_er_dict'][a_i][i_ig], color=colrs[i_i], line_width=1,
                               line_dash=dashes[i_i])
        er_list.append(('{0:.3f} {1}M'.format(ach_conc, mu), [er_line]))

        ik_line = ik_fig.line(t_arr, result_d['soma_isk_dict'][a_i][i_ig], color=colrs[2*i_i], line_width=3)
        ik_list.append(('Isk {0:.3f} {1}M'.format(ach_conc, mu), [ik_line]))
        ik_line2 = ik_fig.line(t_arr, result_d['soma_im_dict'][a_i][i_ig], color=colrs[2*i_i+1], line_width=3)
        ik_list.append(('Im {0:.3f} {1}M'.format(ach_conc, mu), [ik_line2]))

        if a_i > 0.0:
            acc_line = acc_fig.line(t_dict[a_i], acc_dict[a_i], color=colrs[i_i], line_width=3, line_dash=dashes[i_i])
            acc_list.append(('ACh: {0:.3f} {1}M'.format(ach_conc, mu), [acc_line]))
            acc_fig.scatter(t_dict[a_i], acc_dict[a_i], color=colrs[i_i], marker=markers[i_i], size=12, alpha=0.5)

    f_leg = bmod.Legend(items=freq_list, location='center')
    freq_fig.add_layout(f_leg, 'right')

    cyt_leg = bmod.Legend(items=cyt_list, location='center')
    ca_fig.add_layout(cyt_leg, 'right')

    er_leg = bmod.Legend(items=er_list, location='center')
    er_fig.add_layout(er_leg, 'right')

    ik_leg = bmod.Legend(items=ik_list, location='center')
    ik_fig.add_layout(ik_leg, 'right')

    # acc_opts = get_ach_curve_params(result_d['ach_levels'][1:], peak_acc)
    # print('Acceleration EC50: ', 10 ** acc_opts[1])

    nz_inds = np.nonzero(result_d['ach_levels'])
    ach_arr = np.array(result_d['ach_levels'])
    acc_ach_fig.line(ach_arr[nz_inds], peak_acc, color=colrs[0], line_width=3)
    acc_ach_fig.circle(ach_arr[nz_inds], peak_acc, color=colrs[0], size=12)

    hault_ach_fig.line(ach_arr[nz_inds], max_hault, color=colrs[0], line_width=3)
    hault_ach_fig.circle(ach_arr[nz_inds], max_hault, color=colrs[0], size=12)

    peak_ca_vals = np.array(peak_ca_vals)
    pca_fig.line(ach_arr[nz_inds], peak_ca_vals * 1000.0, line_width=3, color='green')
    pca_fig.circle(ach_arr[nz_inds], peak_ca_vals * 1000.0, size=12, color='green')

    return ret_figs


def plot_ahp_adp_test(result_d, phasic_or_tonic, ach_times=[0, 50], t_ignore=0, plot_indices=[]):

    dep_vals = []
    hyp_vals = []
    peak_ca_vals = []

    v_list = []
    cyt_list = []
    er_list = []

    for a_i, ach_conc in enumerate(result_d['ach_levels']):

        i_3k = np.squeeze(np.where(result_d['t'][a_i] > 3000.0))

        if a_i == 0:
            rmp = result_d['soma_v_dict'][a_i][i_3k[0]]

        else:
            i_hyp = np.intersect1d(np.where(result_d['t'][a_i] > t_ignore), np.where(result_d['t'][a_i] < (t_ignore + 3000.0)))
            hyp_val = np.min(result_d['soma_v_dict'][a_i][i_hyp]) - rmp

            hyp_vals.append(hyp_val)

            dep_val = np.max(result_d['soma_v_dict'][a_i]) - rmp
            dep_vals.append(dep_val)

            peak_ca_vals.append(np.max(result_d['ca_cyt_dict'][a_i]))

    print('Resting Membrane Potential: {0} mV'.format(rmp))
    print('Acetylcholine Values: ', result_d['ach_levels'][1:])
    dep_arr = np.array(dep_vals)
    print('Depolarization Values: ', dep_arr)
    if dep_arr.size > 3:
        dep_opts = get_ach_curve_params(result_d['ach_levels'][1:], dep_arr)
        print('Depolarization EC50: ', 10**dep_opts[1])
    hyp_arr = np.array(hyp_vals)
    print('Hyperpolarization Values: ', hyp_arr)
    if hyp_arr.size > 3:
        peak_ca_vals = np.array(peak_ca_vals)
        hyp_opts = get_ach_curve_params(result_d['ach_levels'][1:], hyp_arr, positive=False)
        print('Hyperpolarization EC50: ', 10**hyp_opts[1])

    ret_figs = []

    ach_span_v1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)
    ach_span_v2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                            line_dash='dashed', line_width=3)

    v_fig = bplt.figure(title='Somatic Membrane Potential')
    v_fig.xaxis.axis_label = 'time (msec)'
    v_fig.yaxis.axis_label = 'potential (mV)'
    ret_figs.append(v_fig)
    v_fig.add_layout(ach_span_v1)
    v_fig.add_layout(ach_span_v2)

    ach_span_ca1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_ca2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    ca_fig = bplt.figure(title='Cytosol Calcium')
    ca_fig.xaxis.axis_label = 'time (msec)'
    ca_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(ca_fig)
    ca_fig.add_layout(ach_span_ca1)
    ca_fig.add_layout(ach_span_ca2)

    ach_span_er1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_er2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    er_fig = bplt.figure(title='Endoplasmic Reticulum Calcium')
    er_fig.xaxis.axis_label = 'time (msec)'
    er_fig.yaxis.axis_label = 'concentration ({}M)'.format(mu)
    ret_figs.append(er_fig)
    er_fig.add_layout(ach_span_er1)
    er_fig.add_layout(ach_span_er2)

    ach_span_ik1 = bmod.Span(location=ach_times[0] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)
    ach_span_ik2 = bmod.Span(location=ach_times[1] - t_ignore, dimension='height', line_color='green',
                           line_dash='dashed', line_width=3)

    ik_fig = bplt.figure(title='Potassium Currents')
    ik_fig.xaxis.axis_label = 'time (msec)'
    ik_fig.yaxis.axis_label = 'current density (mA/cm^2)'
    ret_figs.append(ik_fig)
    ik_fig.add_layout(ach_span_ik1)
    ik_fig.add_layout(ach_span_ik2)

    dep_fig = bplt.figure(title='Depolarization/Hyperpolarization vs {} ACh'.format(phasic_or_tonic), x_axis_type='log')
    dep_fig.xaxis.axis_label = 'ACh ({}M)'.format(mu)
    dep_fig.yaxis.axis_label = 'change in potential (mV)'

    ret_figs.append(dep_fig)

    pca_fig = bplt.figure(title='Peak Cytosol Calcium vs {} ACh'.format(phasic_or_tonic), x_axis_type='log')
    pca_fig.xaxis.axis_label = 'ACh ({}M)'.format(mu)
    pca_fig.yaxis.axis_label = '[Ca] ({}M)'.format(mu)
    ret_figs.append(pca_fig)

    dashes = ['solid', 'solid', 'solid', 'solid','solid', 'solid','solid', 'solid','solid', 'solid',]

    if plot_indices:
        my_is = plot_indices
    else:
        my_is = range(len(result_d['ach_levels']))

    for i_i, a_i in enumerate(my_is):
        ach_conc = result_d['ach_levels'][a_i]
        i_ig = np.squeeze(np.where(result_d['t'][a_i] > t_ignore))

        if not plot_indices or a_i in plot_indices:
            vline = v_fig.line(result_d['t'][a_i][i_ig] - t_ignore, result_d['soma_v_dict'][a_i][i_ig],
                               color=colrs[i_i], line_width=3, line_dash=dashes[i_i])
            v_list.append(('{0:.3f} {1}M'.format(ach_conc, mu), [vline,]))

            cytline = ca_fig.line(result_d['t'][a_i][i_ig] - t_ignore, 1000.0*result_d['ca_cyt_dict'][a_i][i_ig],
                                  color=colrs[i_i], line_width=3, line_dash=dashes[i_i])
            cyt_list.append(('{0:.3f} {1}M'.format(ach_conc, mu), [cytline,]))

            erline = er_fig.line(result_d['t'][a_i][i_ig] - t_ignore, 1000.0 * result_d['ca_er_dict'][a_i][i_ig],
                                 color=colrs[i_i], line_width=3, line_dash=dashes[i_i])
            er_list.append(('{0:.3f} {1}M'.format(ach_conc, mu), [erline, ]))

            ik_fig.line(result_d['t'][a_i][i_ig] - t_ignore, result_d['soma_isk_dict'][a_i][i_ig],
                        color=colrs[2*i_i], line_dash=dashes[i_i],
                        legend='I_SK ACh: {0:.3g} {1}M'.format(ach_conc, mu), line_width=3)
            ik_fig.line(result_d['t'][a_i][i_ig] - t_ignore, result_d['soma_isk_dict'][a_i][i_ig],
                        color=colrs[2*i_i+1], line_dash=dashes[i_i],
                        legend='I_M ACh: {0:.3g} {1}M'.format(ach_conc, mu), line_width=3)

    v_leg = bmod.Legend(items=v_list, location='center')
    v_fig.add_layout(v_leg, 'right')

    cyt_leg = bmod.Legend(items=cyt_list, location='center')
    ca_fig.add_layout(cyt_leg, 'right')

    er_leg = bmod.Legend(items=er_list, location='center')
    er_fig.add_layout(er_leg, 'right')

    dep_items = []

    dep_line = dep_fig.line(result_d['ach_levels'][1:], dep_arr, line_width=3, color='blue')
    # dep_items.append(('Depolarization', [dep_line,]))
    dep_fig.square(result_d['ach_levels'][1:], dep_arr, size=12, color='blue', legend='Depolarization')
    hyp_line = dep_fig.line(result_d['ach_levels'][1:], hyp_arr, line_width=3, color='red')
    # dep_items.append(('Hyperpolarization', [hyp_line, ]))
    dep_fig.circle(result_d['ach_levels'][1:], hyp_arr, size=12, color='red', legend='Hyperpolarization')
    dep_fig.legend.location = 'center_right'
    # dep_leg = bmod.Legend(items=dep_items, location='center')
    # dep_fig.add_layout(dep_leg, 'right')

    if len(peak_ca_vals) > 3:
        pca_fig.line(result_d['ach_levels'][1:], peak_ca_vals*1000.0, line_width=3, color='green')
        pca_fig.diamond(result_d['ach_levels'][1:], peak_ca_vals*1000.0, size=12, color='green')

    v_fig.legend.location = 'bottom_right'

    return ret_figs


def sig_func(x, k, x0):
    return 1.0/(1.0 + np.exp(-k*(x - x0)))


def neg_sig_func(x, k, x0, y_max, y_min):
    return -1.0/(1.0 + np.exp(-k*(x - x0)))


def get_ach_curve_params(ach_levels, y_vals, positive=True):

    conc_log = np.log10(ach_levels)

    if positive:
        val_shift = np.array(y_vals) - np.min(y_vals)
        val_norm = val_shift/np.max(val_shift)

        ach_opt, _ = spopt.curve_fit(sig_func, conc_log, val_norm)
    else:
        val_shift = np.array(y_vals) - np.max(y_vals)
        val_norm = -val_shift / np.min(val_shift)

        ach_opt, _ = spopt.curve_fit(neg_sig_func, conc_log, val_norm)

    return ach_opt


if __name__ == '__main__':

    all_figs = []
    all_names = []

    # mpg1412009_A_idA is the file for the parameters that the M1 model was tuned for
    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    # t_path = join(os.getcwd(), 'morphologies/mpg141208_B_idA.asc')

    path_parts = os.getcwd().split('/')
    home_i = path_parts.index('ca1_muscarinic_modulation')
    home_path = '/'.join(path_parts[:(home_i + 1)])

    config_name = 'config_mpg141209_A_idA_test_doi' # Original Model
    # config_name = 'config_mpg141208_B_idA'
    spec_name = config_name + '_spec'
    conf_path = os.path.join(home_path, 'config_files', config_name)
    spec_path = conf_path + '_spec'
    conf_path += '.txt'
    spec_path += '.txt'

    config_obj = conf.load_config(conf_path, spec_path)

    res_dir = os.path.join(home_path,
                           'results/M1 Model Results')
    res_h5file = os.path.join(os.path.join(res_dir, 'm1_model_simulations.h5'))

    # Simulation Options
    # 'ahp test': Limited AHP Test
    # 'tonic depolarization': Tonic ACh Hyperpolarization and Depolarization
    # 'phasic depolarization': Phasic ACh Hyperpolarization and Depolarization
    # 'phasic acceleration': Phasic ACh Spike Acceleration
    # 'tonic acceleration': Tonic ACh Spike Acceleration
    # 'tonic input resistance': Tonic Input Resistance
    # 'tonic rheobase': Tonic ACh rheobase
    # 'epsp modulation': EPSP Modulation
    # 'epsp sweep': EPSP Amplitude vs Time of ACh Pulse
    # 'buffer test': Test Calcium Dynamics vs Presence of Intracellular Buffer

    simulation_choice = 'tonic acceleration'
    save_figs = False
    save_result = True
    use_file_to_plot = True

    if simulation_choice == 'ahp test':
        #############
        ## Limited AHP Test
        #############

        ahp_dict = run_ahp_test(config_obj, run_time=5000, current_dict={'amp': 0.0, 'start': 10100, 'dur': 250.0})
        # ahp_dict = run_ahp_test(config_obj, run_time=5000, current_dict={'amp': 1.2, 'start': 10100, 'dur': 250.0})

        ahp_figs = plot_ahp_test(ahp_dict, z_times=[350, 450], t_ignore=0)
        ahp_names = ['ahp_test_v', 'ahp_test_zoom', 'ahp_test_ca_cyt']

        all_figs += ahp_figs
        all_names += ahp_names

    elif simulation_choice == 'tonic depolarization':
        ###############################################
        ## Tonic Hyperpolarization and Depolarization
        ###############################################

        ach_levs = np.logspace(-3, 2, 11)
        v_plot_inds = [0, 1, 3, 5, 7, 9, 11]
        h5_path = os.path.join(res_dir, 'tonic_depolarization_1min.h5')

        if not use_file_to_plot:
            ton_dict = run_tonic_test(config_obj, run_time=61500, ach_levels=[0.0] + list(ach_levs),
                                      ach_times=[1500, 61500])

            if save_result:

                with pd.HDFStore(h5_path) as h5file:
                    for k, v in ton_dict['soma_v_dict'].items():
                        v_series = pd.Series(v)
                        h5file.put('soma_v_dict_{0}'.format(k), v_series, format='table', data_columns=True)
                    for k, v in ton_dict['ca_cyt_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_cyt_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in ton_dict['ca_er_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_er_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in ton_dict['soma_im_dict'].items():
                        im_series = pd.Series(v)
                        h5file.put('soma_im_dict_{0}'.format(k), im_series, format='table', data_columns=True)
                    for k, v in ton_dict['soma_isk_dict'].items():
                        isk_series = pd.Series(v)
                        h5file.put('soma_isk_dict_{0}'.format(k), isk_series, format='table', data_columns=True)
                    # for k, v in ton_dict['soma_pip2_dict'].items():
                    #     pip2_series = pd.Series(v)
                    #     h5file.put('soma_pip2_dict_{0}'.format(k), pip2_series, format='table', data_columns=True)
                    for k, v in ton_dict['t'].items():
                        t_series = pd.Series(v)
                        h5file.put('t_dict_{0}'.format(k), t_series, format='table', data_columns=True)

                    ach_series = pd.Series(ton_dict['ach_levels'])
                    h5file.put('ach_levels', ach_series, format='table', data_columns=True)

                print('Results saved to {}'.format(h5_path))

                ton_figs = plot_ahp_adp_test(ton_dict, ach_times=[1500, 61500], t_ignore=500, phasic_or_tonic='Tonic',
                                             plot_indices=v_plot_inds)
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
                    if 'ach_levels' in k:
                        ton_dict['ach_levels'] = h5file[k].to_numpy()
                    elif 'soma_v' in k:
                        ai = int(k.split('_')[-1])
                        ton_dict['soma_v_dict'][ai] = h5file[k]
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

                ton_figs = plot_ahp_adp_test(ton_dict, ach_times=[1500, 61500], t_ignore=500, phasic_or_tonic='Tonic',
                                             plot_indices=v_plot_inds)

        ton_names = ['tonic_test_v', 'tonic_test_cyt', 'tonic_test_er', 'tonic_test_ik', 'tonic_test_depol',
                     'tonic_test_peak_ca']
        all_figs += ton_figs
        all_names += ton_names

    elif simulation_choice == 'phasic depolarization':
        ################################################
        ## Phasic Hyperpolarization and Depolarization
        ################################################

        ach_levs = np.logspace(-3, 2, 11)
        v_plot_inds = [0, 1, 3, 5, 7, 9, 11]

        h5_path = os.path.join(res_dir, 'phasic_depolarization.h5')

        if not use_file_to_plot:

            phasic_dict = run_tonic_test(config_obj, run_time=10500, ach_levels=[0.0] + list(ach_levs),
                                         ach_times=[1500, 1550])

            if save_result:
                save_dict = {}
                # Only save numpy arrays holding time series
                with pd.HDFStore(h5_path) as h5file:
                    for k, v in phasic_dict['soma_v_dict'].items():
                        v_series = pd.Series(v)
                        h5file.put('soma_v_dict_{0}'.format(k), v_series, format='table', data_columns=True)
                    for k, v in phasic_dict['ca_cyt_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_cyt_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in phasic_dict['ca_er_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_er_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in phasic_dict['soma_im_dict'].items():
                        im_series = pd.Series(v)
                        h5file.put('soma_im_dict_{0}'.format(k), im_series, format='table', data_columns=True)
                    for k, v in phasic_dict['soma_isk_dict'].items():
                        isk_series = pd.Series(v)
                        h5file.put('soma_isk_dict_{0}'.format(k), isk_series, format='table', data_columns=True)
                    # for k, v in phasic_dict['soma_pip2_dict'].items():
                    #     pip2_series = pd.Series(v)
                    #     h5file.put('soma_pip2_dict_{0}'.format(k), pip2_series, format='table', data_columns=True)
                    for k, v in phasic_dict['t'].items():
                        t_series = pd.Series(v)
                        h5file.put('t_dict_{0}'.format(k), t_series, format='table', data_columns=True)

                    ach_series = pd.Series(phasic_dict['ach_levels'])
                    h5file.put('ach_levels', ach_series, format='table', data_columns=True)

                print('Results saved to {}'.format(h5_path))

            phasic_figs = plot_ahp_adp_test(phasic_dict, ach_times=[1500, 1550], t_ignore=500, phasic_or_tonic='Phasic',
                                            plot_indices=v_plot_inds)
        else:
            with pd.HDFStore(h5_path) as h5file:
                phasic_dict = {'soma_v_dict': {},
                               'ca_cyt_dict': {},
                               'ca_er_dict': {},
                               'soma_im_dict': {},
                               'soma_isk_dict': {},
                               'soma_pip2_dict': {},
                               't': {}}

                for k in h5file.keys():
                    if 'ach_levels' in k:
                        phasic_dict['ach_levels'] = h5file[k].to_numpy()
                    elif 'soma_v' in k:
                        ai = int(k.split('_')[-1])
                        phasic_dict['soma_v_dict'][ai] = h5file[k]
                    elif 'soma_im' in k:
                        ai = int(k.split('_')[-1])
                        phasic_dict['soma_im_dict'][ai] = h5file[k]
                    elif 'soma_isk' in k:
                        ai = int(k.split('_')[-1])
                        phasic_dict['soma_isk_dict'][ai] = h5file[k]
                    elif 'ca_cyt' in k:
                        ai = int(k.split('_')[-1])
                        phasic_dict['ca_cyt_dict'][ai] = h5file[k]
                    elif 'ca_er' in k:
                        ai = int(k.split('_')[-1])
                        phasic_dict['ca_er_dict'][ai] = h5file[k]
                    elif 't_dict' in k:
                        ai = int(k.split('_')[-1])
                        phasic_dict['t'][ai] = h5file[k]

                phasic_figs = plot_ahp_adp_test(phasic_dict, ach_times=[1500, 1550], t_ignore=500,
                                                phasic_or_tonic='Phasic', plot_indices=v_plot_inds)
        phasic_names = ['phasic_ahpadp_v', 'phasic_ahpadp_cyt', 'phasic_ahpadp_er', 'phasic_ahpadp_ik', 'phasic_ahpadp_depol',
                        'phasic_ahpdadp_peak_ca']

        all_figs += phasic_figs
        all_names += phasic_names

    elif simulation_choice == 'phasic acceleration':
        ##############################
        ## Phasic Spike Acceleration
        ##############################

        ach_levs = np.logspace(-3, 2, 11)
        v_plot_inds = [0,1,3,5,7,9,11]
        h5_path = os.path.join(res_dir, 'phasic_acceleration.h5')
        sim_dur = 10500
        # Time Allowed for model to achieve steady state. Time series before this value are not plotted
        time_ignore = 500
        c_dict = {'amp': 0.42, 'start': 0.0, 'dur': sim_dur}

        ach_levs = [0.0] + list(ach_levs)

        if not use_file_to_plot:

            acc_dict = run_tonic_test(config_obj, run_time=sim_dur, ach_times=[1500, 1550],
                                      current_dict=c_dict, ach_levels=ach_levs)
            if save_result:

                with pd.HDFStore(h5_path) as h5file:
                    for k, v in acc_dict['soma_v_dict'].items():
                        v_series = pd.Series(v)
                        h5file.put('soma_v_dict_{0}'.format(k), v_series, format='table', data_columns=True)
                    for k, v in acc_dict['ca_cyt_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_cyt_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in acc_dict['ca_er_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_er_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in acc_dict['soma_im_dict'].items():
                        im_series = pd.Series(v)
                        h5file.put('soma_im_dict_{0}'.format(k), im_series, format='table', data_columns=True)
                    for k, v in acc_dict['soma_isk_dict'].items():
                        isk_series = pd.Series(v)
                        h5file.put('soma_isk_dict_{0}'.format(k), isk_series, format='table', data_columns=True)
                    # for k, v in acc_dict['soma_pip2_dict'].items():
                    #     pip2_series = pd.Series(v)
                    #     h5file.put('soma_pip2_dict_{0}'.format(k), pip2_series, format='table', data_columns=True)
                    for k, v in acc_dict['t'].items():
                        t_series = pd.Series(v)
                        h5file.put('t_dict_{0}'.format(k), t_series, format='table', data_columns=True)

                    ach_series = pd.Series(acc_dict['ach_levels'])
                    h5file.put('ach_levels', ach_series, format='table', data_columns=True)

                print('Results saved to {}'.format(h5_path))

            acc_figs = plot_spike_acc_test(acc_dict, ach_times=[1500, 1550], t_ignore=time_ignore,
                                           plot_inds=v_plot_inds)
        else:
            with pd.HDFStore(h5_path) as h5file:

                acc_dict = {'soma_v_dict': {},
                            'ca_cyt_dict': {},
                            'ca_er_dict': {},
                            'soma_im_dict': {},
                            'soma_isk_dict': {},
                            'soma_pip2_dict': {},
                            't': {}}

                for k in h5file.keys():
                    if 'ach_levels' in k:
                        acc_dict['ach_levels'] = h5file[k].to_numpy()
                    elif 'soma_v' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['soma_v_dict'][ai] = np.array(h5file[k])
                    elif 'soma_im' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['soma_im_dict'][ai] = np.array(h5file[k])
                    elif 'soma_isk' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['soma_isk_dict'][ai] = np.array(h5file[k])
                    elif 'ca_cyt' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['ca_cyt_dict'][ai] = np.array(h5file[k])
                    elif 'ca_er' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['ca_er_dict'][ai] = np.array(h5file[k])
                    elif 't_dict' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['t'][ai] = np.array(h5file[k])

                acc_figs = plot_spike_acc_test(acc_dict, ach_times=[1500, 1550], t_ignore=time_ignore,
                                               plot_inds=v_plot_inds)
        acc_names = ['phasic_accel_test_ifr', 'phasic_accel_test_accel',
                     'phasic_accel_test_accel_vs_ach', 'phasic_accel_test_cessation_vs_ach', 'phasic_accel_test_cyt',
                     'phasic_accel_test_er', 'phasic_accel_test_ik', 'phasic_peak_cyt_ca']

        if v_plot_inds:
            for v_ind in v_plot_inds:
                ach_level = ach_levs[v_ind]
                lev_str = str(ach_level).replace('.', 'p')
                acc_names.append('phasic_accel_test_potential_{0}'.format(lev_str))
        else:
            for ach_level in ach_levs:
                lev_str = str(ach_level).replace('.', 'p')
                acc_names.append('phasic_accel_test_potential_{0}'.format(lev_str))

        all_figs += acc_figs
        all_names += acc_names

    elif simulation_choice == 'tonic acceleration':
        #############################
        ## Tonic Spike Acceleration
        #############################

        ach_levs = np.logspace(-3, 2, 11)
        h5_path = os.path.join(res_dir, 'tonic_acceleration.h5')
        v_plot_inds = [0, 1, 3, 5, 7, 9, 11]
        sim_dur = 10500
        c_dict = {'amp': 0.420, 'start': 0.0, 'dur': sim_dur}

        ach_levs = [0.0] + list(ach_levs)

        if not use_file_to_plot:

            acc_dict = run_tonic_test(config_obj, run_time=sim_dur, ach_levels=ach_levs,
                                      ach_times=[1500, sim_dur], current_dict=c_dict)

            if save_result:

                with pd.HDFStore(h5_path) as h5file:
                    for k, v in acc_dict['soma_v_dict'].items():
                        v_series = pd.Series(v)
                        h5file.put('soma_v_dict_{0}'.format(k), v_series, format='table', data_columns=True)
                    for k, v in acc_dict['ca_cyt_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_cyt_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in acc_dict['ca_er_dict'].items():
                        ca_series = pd.Series(v)
                        h5file.put('ca_er_dict_{0}'.format(k), ca_series, format='table', data_columns=True)
                    for k, v in acc_dict['soma_im_dict'].items():
                        im_series = pd.Series(v)
                        h5file.put('soma_im_dict_{0}'.format(k), im_series, format='table', data_columns=True)
                    for k, v in acc_dict['soma_isk_dict'].items():
                        isk_series = pd.Series(v)
                        h5file.put('soma_isk_dict_{0}'.format(k), isk_series, format='table', data_columns=True)
                    # for k, v in acc_dict['soma_pip2_dict'].items():
                    #     pip2_series = pd.Series(v)
                    #     h5file.put('soma_pip2_dict_{0}'.format(k), pip2_series, format='table', data_columns=True)
                    for k, v in acc_dict['t'].items():
                        t_series = pd.Series(v)
                        h5file.put('t_dict_{0}'.format(k), t_series, format='table', data_columns=True)

                    ach_series = pd.Series(acc_dict['ach_levels'])
                    h5file.put('ach_levels', ach_series, format='table', data_columns=True)

                print('Results saved to {}'.format(h5_path))
                acc_figs = plot_spike_acc_test(acc_dict, ach_times=[1500, 20500], t_ignore=500, plot_inds=v_plot_inds)
        else:
            with pd.HDFStore(h5_path) as h5file:
                acc_dict = {'soma_v_dict': {},
                            'ca_cyt_dict': {},
                            'ca_er_dict': {},
                            'soma_im_dict': {},
                            'soma_isk_dict': {},
                            'soma_pip2_dict': {},
                            't': {}}

                for k in h5file.keys():
                    if 'ach_levels' in k:
                        acc_dict['ach_levels'] = h5file[k].to_numpy()
                    elif 'soma_v' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['soma_v_dict'][ai] = np.array(h5file[k])
                    elif 'soma_im' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['soma_im_dict'][ai] = np.array(h5file[k])
                    elif 'soma_isk' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['soma_isk_dict'][ai] = np.array(h5file[k])
                    elif 'ca_cyt' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['ca_cyt_dict'][ai] = np.array(h5file[k])
                    elif 'ca_er' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['ca_er_dict'][ai] = np.array(h5file[k])
                    elif 't_dict' in k:
                        ai = int(k.split('_')[-1])
                        acc_dict['t'][ai] = np.array(h5file[k])

                acc_figs = plot_spike_acc_test(acc_dict, ach_times=[1500, 20500], t_ignore=500, plot_inds=v_plot_inds)

        acc_names = ['tonic_accel_test_ifr', 'tonic_accel_test_accel',
                     'tonic_accel_test_accel_vs_ach', 'tonic_accel_test_cessation_vs_ach', 'tonic_accel_test_cyt',
                     'tonic_accel_test_er', 'tonic_accel_test_ik', 'tonic_peak_cyt_ca']

        for ach_level in ach_levs:
            lev_str = str(ach_level).replace('.', 'p')
            acc_names.append('tonic_accel_test_v_{0}'.format(lev_str))

        all_figs += acc_figs
        all_names += acc_names

    elif simulation_choice == 'tonic input resistance':
        ###########################
        ## Tonic Input Resistance
        ###########################

        ach_levs = np.logspace(-3, 2, 11)
        h5_path = os.path.join(res_dir, 'tonic_input_resistance.h5')
        sim_dur = 10100
        c_dict = {'amps': [-0.1, 0.0, 0.1], 'start': 6100, 'dur': 400}

        ach_levels = [0.0] + list(ach_levs)

        if not use_file_to_plot:
            inr_dict = run_input_resistance_test(config_obj, run_time=sim_dur, ach_levels=ach_levels,
                                                 ach_times=[200, sim_dur - 200], current_dict=c_dict)

            if save_result:

                with pd.HDFStore(h5_path) as h5file:
                    for k, v in inr_dict['soma_v_dict'].items():
                        for k2, v2 in v.items():
                            v_series = pd.Series(v2)
                            h5file.put('soma_v_dict_{0}_{1}'.format(k, k2), v_series, format='table', data_columns=True)
                    for k, v in inr_dict['ca_cyt_dict'].items():
                        for k2, v2 in v.items():
                            ca_series = pd.Series(v2)
                            h5file.put('ca_cyt_dict_{0}_{1}'.format(k, k2), ca_series, format='table', data_columns=True)
                    for k, v in inr_dict['ca_er_dict'].items():
                        for k2, v2 in v.items():
                            ca_series = pd.Series(v2)
                            h5file.put('ca_er_dict_{0}_{1}'.format(k, k2), ca_series, format='table', data_columns=True)
                    for k, v in inr_dict['soma_im_dict'].items():
                        for k2, v2 in v.items():
                            im_series = pd.Series(v2)
                            h5file.put('soma_im_dict_{0}_{1}'.format(k, k2), im_series, format='table', data_columns=True)
                    for k, v in inr_dict['soma_isk_dict'].items():
                        for k2, v2 in v.items():
                            isk_series = pd.Series(v2)
                            h5file.put('soma_isk_dict_{0}_{1}'.format(k, k2), isk_series, format='table', data_columns=True)
                    for k, v in inr_dict['t'].items():
                        for k2, v2 in v.items():
                            t_series = pd.Series(v2)
                            h5file.put('t_dict_{0}_{1}'.format(k, k2), t_series, format='table', data_columns=True)

                    cur_series = pd.Series(inr_dict['current_amplitudes'])
                    h5file.put('current_amplitudes', cur_series, format='table', data_columns=True)

                    ach_series = pd.Series(inr_dict['ach_levels'])
                    h5file.put('ach_levels', ach_series, format='table', data_columns=True)

                print('Results saved to {}'.format(h5_path))

            inr_plots, inr_vals = plot_input_resistance_test(inr_dict, curr_times=[6100, 6500],
                                                             ach_times=[200, sim_dur - 200], t_ignore=100)
        else:
            with pd.HDFStore(h5_path) as h5file:
                inr_dict = {'soma_v_dict': {}, 'soma_im_dict': {}, 'soma_isk_dict': {},
                            'ca_cyt_dict': {}, 'ca_er_dict': {}, 't': {}}
                for k in h5file.keys():
                    if 'soma_v_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)
                        if a_num not in inr_dict['soma_v_dict']:
                            inr_dict['soma_v_dict'][a_num] = {}
                        inr_dict['soma_v_dict'][a_num][c_num] = h5file[k]
                    elif 'ca_cyt_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)
                        if a_num not in inr_dict['ca_cyt_dict']:
                            inr_dict['ca_cyt_dict'][a_num] = {}
                        inr_dict['ca_cyt_dict'][a_num][c_num] = h5file[k]
                    elif 'ca_er_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)
                        if a_num not in inr_dict['ca_er_dict']:
                            inr_dict['ca_er_dict'][a_num] = {}
                        inr_dict['ca_er_dict'][a_num][c_num] = h5file[k]
                    elif 'soma_im_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)
                        if a_num not in inr_dict['soma_im_dict']:
                            inr_dict['soma_im_dict'][a_num] = {}
                        inr_dict['soma_im_dict'][a_num][c_num] = h5file[k]
                    elif 'soma_isk_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)
                        if a_num not in inr_dict['soma_isk_dict']:
                            inr_dict['soma_isk_dict'][a_num] = {}
                        inr_dict['soma_isk_dict'][a_num][c_num] = h5file[k]
                    elif 'current_amplitudes' in k:
                        inr_dict['current_amplitudes'] = h5file[k].to_numpy()
                    elif 'ach_levels' in k:
                        inr_dict['ach_levels'] = h5file[k].to_numpy()
                    elif 't_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)
                        if a_num not in inr_dict['t']:
                            inr_dict['t'][a_num] = {}
                        inr_dict['t'][a_num][c_num] = h5file[k]

                inr_plots, inr_vals = plot_input_resistance_test(inr_dict, curr_times=[6100, 6500],
                                                         ach_times=[200, sim_dur-200], t_ignore=100)

        print('Minimum Input Resistance: {0} M{1}'.format(inr_vals[0], ohm))
        print('Maximum Input Resistance: {0} M{1}'.format(inr_vals[-1], ohm))
        perc_inc = 100*(inr_vals[-1] - inr_vals[0])/inr_vals[0]
        print('Percent Increase: {0} %'.format(perc_inc))

        all_figs += inr_plots

        inr_names = []

        for ach_level in ach_levels:
            lev_str = str(ach_level).replace('.', 'p')
            inr_names.append('input_resistance_test_v_{0}'.format(lev_str))
            inr_names.append('input_resistance_test_i_{0}'.format(lev_str))

        inr_names += ['input_resistance_test_depolarization', 'input_resistance_test_rin']
        all_names += inr_names

    elif simulation_choice == 'tonic rheobase':
        ####################
        ## Tonic rheobase
        ####################

        h5_path = os.path.join(res_dir, 'tonic_rheobase_test.h5')

        ach_levs = np.logspace(-3, 2, 11)

        sim_dur = 5300
        c_dict = {'amps': [0.0], 'start': 5000, 'dur': 200}

        ach_levels = [0.0] + list(ach_levs)

        if not use_file_to_plot:

            rheo_dict = run_rheobase_test(config_obj, run_time=sim_dur, ach_levels=ach_levels, starting_int=[0.0, 0.4],
                                          desired_int=0.01, ach_times=[100, sim_dur - 100], current_dict=c_dict)

            if save_result:
                with pd.HDFStore(h5_path) as h5file:
                    for k, v in rheo_dict['soma_v_dict'].items():
                        for k2, v2 in v.items():
                            v_series = pd.Series(v2)
                            h5file.put('rheobase/soma_v_dict_{0}_{1}'.format(k, k2), v_series, format='table', data_columns=True)

                    for k, v in rheo_dict['t'].items():
                        for k2, v2 in v.items():
                            t_series = pd.Series(v2)
                            h5file.put('rheobase/t_dict_{0}_{1}'.format(k, k2), t_series, format='table', data_columns=True)

                    for k, v in rheo_dict['current_dict'].items():
                        save_df = pd.Series(v)
                        h5file.put('rheobase/current_dict_{}'.format(k), save_df, format='table', data_columns=True)

                    ach_vals = pd.Series(rheo_dict['ach_levels'])
                    h5file.put('rheobase/ach_levels'.format(k), ach_vals, format='table', data_columns=True)
                    rheo_vals = pd.Series(rheo_dict['rheobase_values'])
                    h5file.put('rheobase/rheobase_values'.format(k), rheo_vals, format='table', data_columns=True)

                print('Results saved to {}'.format(h5_path))

        else:
            with pd.HDFStore(h5_path) as h5file:

                rheo_dict = {'soma_v_dict': {}, 'current_dict': {}, 't': {}}

                for k in h5file.keys():

                    if 'soma_v_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)

                        if a_num not in rheo_dict['soma_v_dict']:
                            rheo_dict['soma_v_dict'][a_num] = {}
                        rheo_dict['soma_v_dict'][a_num][c_num] = h5file[k]
                    elif 'current_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-1])
                        rheo_dict['current_dict'][a_num] = h5file[k]

                    elif 't_dict' in k:
                        k_parts = k.split('_')
                        a_num = int(k_parts[-2])
                        c_num = k_parts[-1]
                        if c_num.isdigit():
                            c_num = int(c_num)

                        if a_num not in rheo_dict['t']:
                            rheo_dict['t'][a_num] = {}
                        rheo_dict['t'][a_num][c_num] = h5file[k]
                    elif 'ach_levels' in k:
                        rheo_dict['ach_levels'] = h5file[k].to_numpy()
                    elif 'rheobase_values' in k:
                        rheo_dict['rheobase_values'] = h5file[k].to_numpy()

        rheo_plots = plot_rheobase_test(rheo_dict, ach_times=[100, sim_dur - 100], curr_times=[5000, 5200])

        all_figs += rheo_plots

        thresh_names = []

        print('Minimum Rheobase: {0} nA'.format(rheo_dict['rheobase_values'][-1]))
        print('Maximum Rheobase: {0} nA'.format(rheo_dict['rheobase_values'][0]))
        perc_inc = 100*(rheo_dict['rheobase_values'][-1]-rheo_dict['rheobase_values'][0])/rheo_dict['rheobase_values'][0]
        print('Percent Increase: {0} %'.format(perc_inc))

        for ach_level in ach_levels:
            lev_str = str(ach_level).replace('.', 'p')
            thresh_names.append('tonic_rheobase_test_v_{0}'.format(lev_str))

        thresh_names += ['tonic_rheobase_test_rheobase_vs_ach']
        all_names += thresh_names

    elif simulation_choice == 'epsp modulation':
        ####################
        ## EPSP modulation
        ####################

        ach_levs = np.logspace(-2, 2, 5)

        sim_dur = 5200

        ach_levels = [0.0] + list(ach_levs)

        t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')
        dis_cell = ca1p.MyCell(t_path, config_obj)
        dis_cell.insert_rxd()

        syn_locs = {}

        for sec_list in [dis_cell.somatic, dis_cell.apical, dis_cell.basal]:
            for sec in sec_list:
                syn_locs[sec.name()] = [0.5,]

        dv_dict = run_epsp_test(dis_cell, [100, sim_dur - 100], syn_locs,
                                run_time=sim_dur,
                                syn_delay=4000,
                                ach_levels=ach_levels,
                                synapse_dict={'amplitude': 0.00115, 'tau1': 0.1, 'tau2': 5.0, 'e': 0.0, 'delay': 1.0})

        min_dv_dict = {}
        min_dv_dict['min dV'] = dv_dict['min dV']
        min_dv_dict['max dV'] = dv_dict['max dV']
        min_dv_dict['ach levels'] = dv_dict['ach levels']
        min_dv_dict['depolarization values'] = dv_dict['depolarization values']
        min_dv_dict['soma v dict'] = {}
        min_dv_dict['t'] = {}

        for a_i, a_val in enumerate(ach_levels):
            min_dv_dict['soma v dict'][a_i] = {}
            min_dv_dict['t'][a_i] = {}
            for sec_name in dv_dict['soma v dict'][a_i]:
                if 'soma' in sec_name:
                    min_dv_dict['soma v dict'][a_i][sec_name] = dv_dict['soma v dict'][a_i][sec_name]
                    min_dv_dict['t'][a_i][sec_name] = dv_dict['t'][a_i][sec_name]

        # my_dict = min_dv_dict
        now = dat.datetime.now()
        date_str = now.strftime('%Y_%m_%d')
        file_name = 'dv_results_{}.p'.format(date_str)
        print(file_name)

        with open(file_name, 'wb') as p_obj:
            pickle.dump(min_dv_dict, p_obj)

        l_file = 'dv_results_2019_10_03.p'

        with open(l_file, 'rb') as p_obj:
            my_dict = pickle.load(p_obj)

        line_dict = dis_cell.get_line_segs()
        dv_plots = plot_epsp_test(my_dict, line_dict, ach_times=[100, sim_dur - 100])

        all_figs += dv_plots

        dv_names = []

        for ach_level in ach_levels:
            lev_str = str(ach_level).replace('.', 'p')
            dv_names.append('dV_vs_location_{0}'.format(lev_str))

        dv_names += ['dV_vs_ACh']
        all_names += dv_names

    elif simulation_choice == 'epsp sweep':
        ########################################
        # EPSP amplitude vs Time of ACh Pulse
        ########################################

        t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')
        dis_cell = ca1p.MyCell(t_path, config_obj)
        dis_cell.insert_rxd()

        sim_dur = 5000
        ach_level = 100
        ach_times = [1200, 1250]
        stim_times = [1000, 1250, 1400, 1700, 2200, 2700, 3200]

        targ_sec = {'somatic[0]': [0.5]}

        sw_dict = run_epsp_sweep(dis_cell, ach_times, stim_times, targ_sec, run_time=sim_dur,
                                 synapse_dict={'amplitude': 0.00115, 'tau1': 0.1, 'tau2': 5.0, 'e': 0.0, 'delay': 1.0})
        sweep_figs = plot_epsp_sweep(sw_dict, t_ignore=200)
        all_figs += sweep_figs

        sweep_names = ['epsp_sweep_soma_v', 'epsp_sweep_dv', 'epsp_sweep_perc']
        all_names += sweep_names

    elif simulation_choice == 'buffer test':
        ##############################################################
        # Test Calcium Dynamics vs Presence of Intracellular Buffer
        ##############################################################

        sim_dur = 5200
        ach_levels = [100.0]  # (uM)
        ach_times = [1200, 1250]

        calbindin_concs = [0.0, 0.045, 0.090]  # (mM)

        buff_dict = run_buffer_comparison(config_obj, calbindin_concs, ach_times, ach_concs=ach_levels, run_time=sim_dur)
        buff_figs = plot_buffer_comparison(buff_dict, t_ignore=200)
        all_figs += buff_figs

        b_names = ['buffer_comp_soma_v', 'buffer_comp_ca_cyt', 'buffer_comp_ca_er', 'buffer_comp_time2peak',
                   'buffer_comp_integrated_ca_cyt_2sec']

        for ach_conc in ach_levels:
            buff_names = [b_str+'_ACh_{0}'.format(ach_conc) for b_str in b_names]
            all_names += buff_names

    elif simulation_choice == 'calcium wave':
        # tl_dict = get_line_segs_calcium(config_obj)

        time_steps = [1500, 1600, 1700, 1800, 1900, 2000, 3000, 4000]
        ach_times = [1510, 1560]
        ach_conc = 100
        t_ig = 500
        line_dict, res_dict = run_calcium_wave(config_obj, time_steps, ach_times, ach_conc, run_time=4000)

        cal_figs = plot_calcium_wave(res_dict, line_dict, t_ignore=t_ig)

        all_figs += cal_figs
        cw_figs = ['calc_wave_time_{}'.format(time_step-t_ig) for time_step in time_steps]
        all_names += ['calc_wave_calc_cyt', 'calc_wave_calc_er'] + cw_figs

    #####################
    ## Plot all Figures
    #####################

    for f_name, fig in zip(all_names, all_figs):
        if 'vs_location' in f_name:
            set_size = False
        elif 'accel_test_v' in f_name:
            set_size = False
        elif 'calc_wave_time' in f_name:
            set_size = False
        else:
            set_size = True
        plt_res.configure_plot(fig, set_size=set_size)

    ########################
    # Plot Cell Morphology
    ########################

    plot_morph = False
    color_code = 'diameter'

    if plot_morph:

        t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

        cell = ca1p.MyCell(t_path, config_obj)

        my_line_dict = cell.get_line_segs()

        my_color_dict = {}

        if color_code == 'diameter':  # Color section based upon section diameter

            for sec_name in my_line_dict:
                if my_line_dict[sec_name]['diam'] > 5.0:
                    my_color_dict[sec_name] = 'purple'
                elif my_line_dict[sec_name]['diam'] > 1.0:
                    my_color_dict[sec_name] = 'red'
                elif my_line_dict[sec_name]['diam'] > 0.5:
                    my_color_dict[sec_name] = 'orange'
                else:
                    my_color_dict[sec_name] = 'black'

        elif color_code == 'rxd':
            for sec_name in my_line_dict:
                if sec_name in cell.cal_names:
                    my_color_dict[sec_name] = 'red'
                else:
                    my_color_dict[sec_name] = 'black'

        all_figs += [plot_cell_bokeh(my_line_dict, my_color_dict, title='Sections Using RXD')]
        all_names += ['cell_morphology']

    if save_figs:

        driver_path = r'/usr/local/bin/geckodriver'
        my_opts = webdriver.firefox.options.Options()
        my_opts.add_argument('start-maximized')
        my_opts.add_argument('disable-infobars')
        my_opts.add_argument('--disable-extensions')
        my_opts.add_argument('window-size=1200,1000')
        my_opts.headless = True

        my_driver = webdriver.Firefox(options=my_opts, executable_path=driver_path)

        for fig, name in zip(all_figs, all_names):
            fig_path = os.path.join(res_dir, '{}.png'.format(name))
            print('Saving {}'.format(fig_path))

            bkio.export_png(fig, filename=fig_path) # , webdriver=my_driver)

            # fig_path = os.path.join(res_dir, '{}.svg'.format(name))
            # fig.output_backend = "svg"
            # bkio.export_svgs(fig, filename=fig_path)

    if all_figs and not save_figs:
        bkio.show(blay.column(all_figs))