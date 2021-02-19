import numpy as np
import os
import sys

from neuron import h
from os.path import join

import matplotlib.pyplot as plt

import bokeh.io as bkio
import bokeh.layouts as blay
import bokeh.models as bmod
import bokeh.plotting as bplt
import bokeh.palettes as bpal
import frequency_analysis as fan
from bokeh.palettes import Category20 as palette
from bokeh.palettes import Category20b as paletteb
import colorcet as cc
from selenium import webdriver

colrs = palette[20] + paletteb[20]

import plot_results as plt_res

module_path = os.getcwd()  # os.path.abspath(os.path.join('../..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import ca1_pyramidal as ca1p

mu = u'\u03BC'
delta = u'\u0394'
ohm = u'\u03A9'


def plot_atype_test(result_d, t_ignore=0, plot_inds=[]):

    inhib_vals = result_d['inhibition_values']

    ret_figs = []

    s_v_fig = bplt.figure(title='Somatic Membrane Potential')
    ret_figs.append(s_v_fig)

    a_v_fig = bplt.figure(title='Apical Membrane Potential')
    ret_figs.append(a_v_fig)

    s_ika_fig = bplt.figure(title='Somatic A-Type Current')
    ret_figs.append(s_ika_fig)
    a_ika_fig = bplt.figure(title='Apical A-Type Current')
    ret_figs.append(a_ika_fig)

    s_v_items = []
    a_v_items = []
    s_i_items = []
    a_i_items = []

    i_ig = np.squeeze(np.where(result_d['t'] > t_ignore))
    t_arr = np.array(result_d['t'][i_ig]) - t_ignore

    for i_i, i_val in enumerate(inhib_vals):
        leg_str = '{0} %'.format(i_val*100)
        s_v_line = s_v_fig.line(t_arr, result_d['soma_v_dict'][i_i][i_ig], color=colrs[i_i], line_width=3)
        s_v_items.append((leg_str, [s_v_line,] ))
        a_v_line = a_v_fig.line(t_arr, result_d['apical_v_dict'][i_i][i_ig], color=colrs[i_i], line_width=3)
        a_v_items.append((leg_str, [a_v_line,] ))

        s_i_line = s_ika_fig.line(t_arr, result_d['soma_ika_dict'][i_i][i_ig], color=colrs[i_i], line_width=3)
        s_i_items.append((leg_str, [s_i_line,]))
        a_i_line = a_ika_fig.line(t_arr, result_d['apical_ika_dict'][i_i][i_ig], color=colrs[i_i], line_width=3)
        a_i_items.append((leg_str, [a_i_line,]))

    s_v_leg = bmod.Legend(items=s_v_items, location='center')
    s_v_fig.add_layout(s_v_leg, 'right')

    a_v_leg = bmod.Legend(items=a_v_items, location='center')
    a_v_fig.add_layout(a_v_leg, 'right')

    s_i_leg = bmod.Legend(items=s_i_items, location='center')
    s_ika_fig.add_layout(s_i_leg, 'right')

    a_i_leg = bmod.Legend(items=a_i_items, location='center')
    a_ika_fig.add_layout(a_i_leg, 'right')

    return ret_figs


def run_atype_test(param_dict={}, inhib_levels=[0.0, 1.0], run_time=300, current_dict={'amp': 0.0, 'start': 0.0, 'dur': 0.0}):

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    my_cells = []
    my_clamps = []
    soma_v_vecs = []
    s_ika_vecs = []
    apic_v_vecs = []
    a_ika_vecs = []

    for i_i, i_val in enumerate(inhib_levels):
        my_cells.append(ca1p.MyCell(t_path, True, param_dict))

        for sec in my_cells[-1].apical:
            for seg in sec:
                pre_val = seg.gkabar_kad
                seg.gkabar_kad = i_val*pre_val

        my_clamps.append(h.IClamp(my_cells[-1].somatic[0](0.5)))

        my_clamps[-1].amp = current_dict['amp']
        my_clamps[-1].delay = current_dict['start']
        my_clamps[-1].dur = current_dict['dur']

        soma_v_vecs.append(h.Vector().record(my_cells[-1].somatic[0](0.5)._ref_v))
        apic_v_vecs.append(h.Vector().record(my_cells[-1].apical[35](0.5)._ref_v))
        s_ika_vecs.append(h.Vector().record(my_cells[-1].somatic[0](0.5)._ref_ik_kap))
        a_ika_vecs.append(h.Vector().record(my_cells[-1].apical[35](0.5)._ref_ik_kad))

    t_vec = h.Vector().record(h._ref_t)
    print('Cells and Records created')

    cv_act = 1
    h.cvode.active(cv_act)

    h.v_init = -69.4
    h.celsius = 34.0

    print('cvode active: {0}'.format(cv_act))

    s_v_dict = {}
    a_v_dict = {}
    s_ika_dict = {}
    a_ika_dict = {}

    h.stdinit()
    print('Initialized')

    h.continuerun(run_time)
    print('Finished Running')

    for i_i, i_val in enumerate(inhib_levels):
        s_v_dict[i_i] = np.array(soma_v_vecs[i_i])
        a_v_dict[i_i] = np.array(apic_v_vecs[i_i])

        s_ika_dict[i_i] = np.array(s_ika_vecs[i_i])
        a_ika_dict[i_i] = np.array(a_ika_vecs[i_i])

    res_dict = {
                'soma_v_dict': s_v_dict,
                'apical_v_dict': a_v_dict,
                'soma_ika_dict': s_ika_dict,
                'apical_ika_dict': a_ika_dict,
                'inhibition_values': np.array(inhib_levels),
                't': np.array(t_vec),
                }

    return res_dict



if __name__ == "__main__":
    my_dict_ach = {'frac_cyt': 0.9,
                   'frac_er': 0.1,
                   'ca ext val': 2.0,
                   'ca_cyt_val': 0.1 * 0.001,
                   'ca_er_val': 175.0 * 0.001,
                   'apical calcium mult': 60.0,
                   'soma calcium mult': 120.0,
                   'soma kca mult': 20.0,
                   'apical kca mult': 20.0,
                   'soma im mult': 1.0,
                   'axon im mult': 1.0,
                   'ca_diff': 0.03,
                   'ip3_init': 0.000003,
                   'ip3k_total': 0.007,
                   'ip5p_total': 0.0002,
                   'cbd_total': 0.045,
                   'car_total': 86.0,
                   'dye identity': 'ogb1',
                   'ogb1_total': 0.05,
                   'ogb5_total': 0.0,
                   'g_ip3r': 1.5e6,
                   'g_er_leak': 2100.0,  # 0.75, # 1300.0,  # 1.8,
                   'g_serca': 63525.0,
                   'g ext leak': 3.2e-5,
                   'total pmca': 2.56e-5,  # 0.0015,
                   'total ncx': 1.6e-6,  # 0.0005,
                   'ip3 lock': False,
                   'tau_leak_soma': 5.0,  # 2.0,
                   'tau_leak_apical': 4.0,  # 0.95,
                   'm1 init': 15.87 * 5.0,
                   'g init': 40.0,
                   'plc init': 3.1198 * 5.0,
                   'pip2 init': 3232,
                   'pi4p init': 4540,
                   'pi init': 226975.0,
                   }

    sim_dur = 500

    # 0.525 = ~10 Hz

    c_dict = {'amp': 0.6, 'start': 0.0, 'dur': sim_dur}

    r_dict = run_atype_test(my_dict_ach, run_time=sim_dur, inhib_levels=[1.0], current_dict=c_dict)

    all_figs = []

    all_figs += plot_atype_test(r_dict, t_ignore=50)

    for fig in all_figs:
        plt_res.configure_plot_notebook(fig)

    bkio.show(blay.column(all_figs))