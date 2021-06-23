import config_utils as conf
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from mpl_toolkits.mplot3d import Axes3D
from itertools import cycle
from neuron import h, crxd as rxd
try:
    from neuron.crxd import rxdmath as rxdm
except:
    print('Could not import rxdmath')
from os.path import join

import bokeh.plotting as bplt
from bokeh.palettes import Category20 as palette
from bokeh.palettes import Category20b as paletteb

colrs = palette[20] + paletteb[20]

mu = u'\u03BC'

h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")
h.load_file("importCell_asc.hoc")


class MyCell(object):
    def __init__(self,
                 morphfilename,
                 config_obj,
                 locn=[[0, 0, 0], [0, 1, 0]],
                 seed=1,
                 cell_num=0,
                 ):

        self.config_obj = config_obj
        self.morph_file = morphfilename
        self.rxd_sim = True
        self.len_dict = {}

        self.c_type = ''

        self.all = []
        self.somatic = []
        self.axonal = []
        self.apical = []
        self.basal = []

        self.cell_orgn = np.array(locn[0])
        self.cell_vect = np.array(locn[1])
        self.cell_num = cell_num
        # self.gid_count = (cell_num + 1)*50000
        self.my_seed = seed
        self.seed_count = 0
        self.loc_mult = 1000000

        l_styles = ["-", "--", "-.", ":"]
        self.l_cycler = cycle(l_styles)

        self.records = {}
        self.syns = {}

        if morphfilename == 'none' or 'default' in morphfilename:
            pass
        elif morphfilename.endswith('.xml'):
            self.load_morphml()
        elif morphfilename.endswith('.swc'):
            self.load_swc()
        elif morphfilename.endswith('.asc'):
            self.load_asc()
        else:
            print('Invalid Morphology File')
            exit()

        self.axon_replacement()

        self.all = self.somatic + self.apical + self.basal + self.axonal

        self.divide_sections_into_segments(10.0)

        apical_trunk_inds = self.config_obj['RXD']['trunk_sections']

        self.cal_list = [self.axonal[0], self.axonal[1], self.somatic[0]]
        self.trunk_list = []
        self.trunk_names = []

        for a_i in apical_trunk_inds:
            a_i = int(a_i)
            self.cal_list.append(self.apical[a_i])
            self.trunk_list.append(self.apical[a_i])
            self.trunk_names.append(self.apical[a_i].name())

        self.cal_names = [sec.name() for sec in self.cal_list]

        self.insert_channels()

        # if self.rxd_sim:
        #     self.insert_rxd()

        self.layer_dict = self.sort_sections()

    def axon_replacement(self):

        len_target = 60
        nseg0 = 5

        nseg_total = 2 * nseg0
        chunk_size = len_target / nseg_total

        count_seg = 0

        seg_diams = []
        seg_lens = []

        for sec in self.c.axonal:

            sec.nseg = 1 + 2 * int(sec.L / chunk_size / 2.0)

            for seg in sec:
                if seg.x > 0 and seg.x < 1.0:

                    seg_diams.append(seg.diam)
                    seg_lens.append(sec.L / sec.nseg)

                    count_seg += 1

                    if count_seg == nseg_total:
                        break
            if count_seg == nseg_total:
                break

        for sec in self.c.axonal:
            h.delete_section(sec=sec)

        self.axonal = [h.Section(name='axon0'), h.Section(name='axon1')]
        self.axonal[0].connect(self.somatic[0](1.0), 0.0)
        self.axonal[1].connect(self.axonal[0](1.0), 0.0)

        count_seg = 0
        len_real = 0

        for sec in self.axonal:

            sec.L = len_target / 2
            sec.nseg = nseg0

            for seg in sec:
                if seg.x > 0 and seg.x < 1.0:
                    seg.diam = seg_diams[count_seg]
                    len_real += seg_lens[count_seg]

                    count_seg += 1

        # print('Target stub axon length: {0} um, equivalent length: {1} um'.format(len_target, len_real))

    def delete_sections(self):

        for sec in self.axonal:
            sec = None
            # h.delete_section(sec=sec)

        self.c.somatic = h.SectionList()
        self.c.apical = h.SectionList()
        self.c.basal = h.SectionList()
        self.c.axonal = h.SectionList()

        self.somatic = []
        self.apical = []
        self.basal = []
        self.axonal = []

        h('forall delete_section()')

        print('All sections deleted')

    def divide_sections_into_segments(self, chunk_size):

        for sec in self.all:
            s_len = sec.L

            nseg = int(round(s_len / chunk_size))
            if not nseg:
                nseg = 1

            sec.nseg = nseg

            # print('{0}: L={1}, {2} segments'.format(sec.name(), sec.L, sec.nseg))

    def insert_channels(self):

        h.distance(0, self.somatic[0](0.5))

        channel_dict = self.config_obj['Ion_Channel']['Apical']

        for sec in self.c.apical:
            sec.cm = channel_dict['cm']
            sec.Ra = channel_dict['Ra']
            sec.insert('pas')
            sec.g_pas = channel_dict['g_pas']

            sec.insert('kdr')
            sec.gkdrbar_kdr = channel_dict['gkdrbar_kdr']
            sec.insert('nax')
            sec.gbar_nax = channel_dict['gbar_nax']
            sec.insert('cal')
            sec.gcalbar_cal = channel_dict['gcalbar_cal']
            sec.insert('can')
            sec.gcanbar_can = channel_dict['gcanbar_can']
            sec.insert('cat')
            sec.gcatbar_cat = channel_dict['gcatbar_cat']
            sec.insert('kca')
            sec.gbar_kca = channel_dict['gbar_kca']
            sec.insert('cagk')
            sec.gbar_cagk = channel_dict['gbar_cagk']

            if not self.rxd_sim or sec.name() not in self.trunk_names:
                sec.insert('cacum')
            elif self.rxd_sim and sec.name() in self.trunk_names:
                sec.gcalbar_cal = channel_dict['calcium_mult'] * channel_dict['gcalbar_cal']
                sec.gcatbar_cat = channel_dict['calcium_mult'] * channel_dict['gcatbar_cat']
                sec.gcanbar_can = channel_dict['calcium_mult'] * channel_dict['gcanbar_can']

                sec.gbar_kca = channel_dict['kca_mult'] * channel_dict['gbar_kca']
            sec.insert('kad')
            sec.insert('hd')

            sec.ena = 50
            sec.ek = -90

            # print(sec.name())
            for seg in sec:
                seg_dist = h.distance(sec(seg.x))
                seg.ghdbar_hd = (1.0 + 3.0 / 100.0 * seg_dist) * channel_dict['ghdbar_hd_min']
                seg.e_pas = channel_dict['e_pas_mult'] - (5.0 * seg_dist / 150.0)
                seg.gkabar_kad = 15.0 / (1.0 + np.exp((300 - seg_dist) / 50.0)) * channel_dict['gkabar_kad_mult']

        channel_dict = self.config_obj['Ion_Channel']['Axonal']

        for sec in self.axonal:
            sec.cm = channel_dict['cm']
            sec.Ra = channel_dict['Ra']
            sec.insert('nax')
            sec.gbar_nax = channel_dict['gbar_nax']

            sec.insert('kdr')
            sec.gkdrbar_kdr = channel_dict['gkdrbar_kdr']

            sec.insert('kmb_inh')
            sec.gbar_kmb_inh = channel_dict['gbar_kmb_inh']

            sec.insert('kap')
            sec.gkabar_kap = channel_dict['gkabar_kap']

            sec.insert('pas')
            sec.g_pas = channel_dict['g_pas']
            sec.e_pas = channel_dict['e_pas']

            sec.ena = channel_dict['ena']
            sec.ek = channel_dict['ek']

        channel_dict = self.config_obj['Ion_Channel']['Basal']

        for sec in self.c.basal:
            sec.cm = channel_dict['cm']
            sec.Ra = channel_dict['Ra']
            sec.insert('kdr')
            sec.gkdrbar_kdr = channel_dict['gkdrbar_kdr']
            sec.insert('nax')
            sec.gbar_nax = channel_dict['gbar_nax']
            sec.insert('can')
            sec.gcanbar_can = channel_dict['gcanbar_can']
            sec.insert('cal')
            sec.gcalbar_cal = channel_dict['gcalbar_cal']
            sec.insert('cat')
            sec.gcatbar_cat = channel_dict['gcatbar_cat']
            sec.insert('kca')
            sec.gbar_kca = channel_dict['gbar_kca']
            sec.insert('cagk')
            sec.gbar_cagk = channel_dict['gbar_cagk']
            sec.insert('cacum')

            sec.insert('pas')
            sec.g_pas = channel_dict['g_pas']

            sec.insert('kad')
            sec.insert('hd')

            sec.ena = channel_dict['ena']
            sec.ek = channel_dict['ek']

            # print(sec.name())
            for seg in sec:
                seg_dist = h.distance(sec(seg.x))
                seg.ghdbar_hd = (1.0 + 3.0 / 100.0 * seg_dist) * channel_dict['ghdbar_hd_min']
                seg.e_pas = channel_dict['e_pas_mult'] - (5.0 * seg_dist / 150.0)
                seg.gkabar_kad = 15.0 / (1.0 + np.exp((300 - seg_dist) / 50.0)) * channel_dict['gkabar_kad_mult']
                # print('dist: {0:.4}, g_hd: {1:.4}, e_pas: {2:.4}, g_kad: {3:.4}'.format(seg_dist,seg.ghdbar_hd,seg.e_pas,seg.gkabar_kad))

        channel_dict = self.config_obj['Ion_Channel']['Somatic']

        for sec in self.c.somatic:
            sec.cm = channel_dict['cm']
            sec.Ra = channel_dict['Ra']

            sec.insert('pas')
            sec.g_pas = channel_dict['g_pas']

            sec.insert('kdr')
            sec.gkdrbar_kdr = channel_dict['gkdrbar_kdr']
            sec.insert('nax')
            sec.gbar_nax = channel_dict['gbar_nax']

            sec.insert('kmb_inh')
            sec.gbar_kmb_inh = channel_dict['gbar_kmb_inh']

            sec.insert('kap')
            sec.gkabar_kap = channel_dict['gkabar_kap']
            sec.insert('can')
            sec.gcanbar_can = channel_dict['gcanbar_can']
            sec.insert('cal')
            sec.gcalbar_cal = channel_dict['gcalbar_cal']
            sec.insert('cat')
            sec.gcatbar_cat = channel_dict['gcatbar_cat']
            sec.insert('kca')
            sec.gbar_kca = channel_dict['gbar_kca']
            sec.insert('cagk')
            sec.gbar_cagk = channel_dict['gbar_cagk']

            if self.rxd_sim:
                sec.gcalbar_cal = channel_dict['calcium_mult'] * channel_dict['gcalbar_cal']
                sec.gcanbar_can = channel_dict['calcium_mult'] * channel_dict['gcanbar_can']
                sec.gcatbar_cat = channel_dict['calcium_mult'] * channel_dict['gcatbar_cat']

                sec.gbar_kca = channel_dict['kca_mult'] * channel_dict['gbar_kca']
            else:
                sec.insert('cacum')

            sec.ena = channel_dict['ena']
            sec.ek = channel_dict['ek']

            sec.insert('hd')

            # print(sec.name())
            for seg in sec:
                seg_dist = h.distance(sec(seg.x))
                seg.ghdbar_hd = (1.0 + 3.0 / 100.0 * seg_dist) * channel_dict['ghdbar_hd_min']
                seg.e_pas = channel_dict['e_pas_mult'] - (5.0 * seg_dist / 150.0)

    def check_seg_is_in_list(self, sec_list, in_node):
        for sec in sec_list:
            if in_node.satisfies(sec):
                return 1.0

        return 0.0

    def insert_rxd(self):
        rxd_dict = self.config_obj['RXD']

        ################
        # Create Regions
        ################

        frac_cyt = rxd_dict['fraction_cyt']
        frac_er = rxd_dict['fraction_er']

        cal_list = self.cal_list

        # spine_psd_list = []
        spine_head_list = []

        for sec_name, loc_dict in self.syns.items():
            for num, sec_dict in loc_dict.items():
                # spine_psd_list.append(sec_dict['psd'])
                spine_head_list.append(sec_dict['head'])

        comb_spine_list = spine_head_list  # + spine_psd_list

        self.cyt = rxd.Region(cal_list, nrn_region='i', geometry=rxd.FractionalVolume(frac_cyt, surface_fraction=1))
        self.er = rxd.Region(cal_list, geometry=rxd.FractionalVolume(frac_er))
        self.cyt_er_membrane = rxd.Region(cal_list, geometry=rxd.ScalableBorder(0.1, on_cell_surface=False))

        ###################
        # Parameter Values
        ###################

        ca_cyt_val = rxd_dict['ca_cyt_val']
        ca_er_val = rxd_dict['ca_er_val']
        ca_diff = rxd_dict['ca_diff']

        ####################
        # M1 mAChR Receptor
        ####################

        kf_L1 = 2.78  # (mM^-1 msec^-1) 0.00278
        kb_L1 = 2.15e-3  # (msec^-1) Altered from original to match ACh binding
        kf_L2 = 2.78  # (mM^-1 msec^-1)
        kb_L2 = 2.15e-5  # (msec^-1) Altered from original to match ACh binding
        kf_G1 = 0.00000027  # (um^2 ms^-1)
        kb_G1 = 0.0068  # (msec^-1)
        kf_G2 = 0.0000027  # (um^2 msec^-1)
        kb_G2 = 0.00068  # (msec^-1)
        k_NX_RLG = 0.00065  # *10.0  # (msec^-1)
        k_NX_G = 0.000000015  # (msec^-1)
        k_NX_P = 0.0047  # (msec^-1)
        k_GTPase1 = 0.000026  # (msec^-1)
        k_GTPase2 = 0.015  # (msec^-1)
        k_PLCassoc = 0.001 * 10.0  # 10.0  # (um^2 ms^-1) Alterred to increase rate of PLC
        k_PLCdiss = 0.00071 * 10.0  # 10.0  # 0.00071  # (ms^-1)
        k_reconst = 0.001  # (um^2 msec^-1)
        k_PLC = 0.0003 * 100.0  # (um^2 msec^-1)

        #################################
        # Phosphatidylinositol Synthesis
        #################################

        k_4k = 0.0000008 * 24.9  # 0.0000228  # 0.0000008*28.5  # (msec^-1)
        k_4p = 0.00012 * 8.3  # 0.00114  # 0.00012*9.5  # (msec^-1)
        k_5k = 0.00002 * 8.3  # 0.00019  # 0.00002*9.5 # (msec^-1)
        k_5p = 0.000028 * 8.3  # 0.000266  # 0.000028*9.5  # (msec^-1)

        kf_pip2 = 0.001  # *10.0  # (msec^-1) [4]
        fold_PIP2 = 3.0  # [4,5]
        k_DAGase = 0.0002  # 0.0002 (ms^-1)

        KD_KCNQ_PIP2 = 2000
        kf_kcnq = 0.00005
        kb_kcnq = KD_KCNQ_PIP2 * kf_kcnq

        pip2_soma_init = rxd_dict['pip2_init']
        pip2_axon_init = rxd_dict['pip2_init']
        b_kcnq_soma = 1724.1 * 0.005449  # (channels*um^-2) (based on 5.8 pS/channel [1] with 0.005449 S/cm^2)
        b_kcnq_axon = 1724.1 * 0.02647  # (channels*um^-2) (based on 5.8 pS S/channel with 0.02647 S/cm^2
        tot_kcnq_soma = b_kcnq_soma * (KD_KCNQ_PIP2 + pip2_soma_init) / pip2_soma_init
        tot_kcnq_axon = b_kcnq_axon * (KD_KCNQ_PIP2 + pip2_axon_init) / pip2_axon_init
        ub_kcnq_soma = tot_kcnq_soma - b_kcnq_soma
        ub_kcnq_axon = tot_kcnq_axon - b_kcnq_axon

        perc_kcnq_base_soma = (b_kcnq_soma / tot_kcnq_soma) ** 2
        perc_kcnq_base_axon = (b_kcnq_axon / tot_kcnq_axon) ** 2

        ######
        # IP3
        ######

        ip3_init = rxd_dict['ip3_init']
        ip3_diff = 0.3

        ################
        # IP3 Breakdown
        ################

        k_ip5p_f = 59.0  # (mM^-1*msec^-1)
        k_ip5p_r = 0.072  # (msec^-1)
        KD_ip5p = k_ip5p_r / k_ip5p_f
        k_ip2 = 0.018  # (msec^-1)

        k_ip3k_ca_f = 1111.1  # (mM^-1*msec^-1)
        k_ip3k_ca_r = 0.1  # (msec^-1)
        KD_ip3k_ca = np.sqrt(k_ip3k_ca_r / k_ip3k_ca_f)

        k_ip3k_ip3_f = 500  # (mM^-1*msec^-1)
        k_ip3k_ip3_r = 0.08  # (msec^-1)
        KD_ip3k_ip3 = k_ip3k_ip3_r / k_ip3k_ip3_f

        k_ip4 = 0.02  # (msec^-1)

        tot_ip5p = rxd_dict['ip5p_total']
        b_ip5p = tot_ip5p * ip3_init / (KD_ip5p + ip3_init)
        ub_ip5p = tot_ip5p - b_ip5p

        #################
        # Calbindin-D28k
        #################

        total_cbd = rxd_dict['cbd_total']

        KD_cbdh = 0.000237
        kf_cbdh = 11.0
        kb_cbdh = KD_cbdh * kf_cbdh

        KD_cbdl = 0.000411
        kf_cbdl = 87.0
        kb_cbdl = KD_cbdl * kf_cbdl

        h0l0_init = total_cbd * 0.319
        h0l1_init = total_cbd * 0.157
        h0l2_init = total_cbd * 0.019
        h1l0_init = total_cbd * 0.236
        h1l1_init = total_cbd * 0.164
        h1l2_init = total_cbd * 0.0197
        h2l0_init = total_cbd * 0.0459
        h2l1_init = total_cbd * 0.0344
        h2l2_init = total_cbd * 0.00414

        ##############
        # Calmodulin
        ##############

        total_cam = rxd_dict['cam_total']
        kf_cam = 10e7
        kb_cam = 10.0
        KD_cam = kb_cam / kf_cam
        camca_init = total_cam * ca_cyt_val ** 3.0 / (KD_cam + ca_cyt_val ** 3.0)
        cam_init = total_cam - camca_init

        ###########################
        # Calcium Binding Proteins
        ###########################

        total_cbp = rxd_dict['cbp_total']
        kf_cbp = 247.0
        kb_cbp = 4.0
        KD_cbp = kb_cbp / kf_cbp
        cbpca_init = total_cbp * ca_cyt_val / (KD_cbp + ca_cyt_val)
        cbp_init = total_cbp - cbpca_init

        ###############
        # Calreticulin
        ###############

        total_car = rxd_dict['car_total']
        KD_car = 2.0
        kf_car = 0.01
        kb_car = KD_car * kf_car

        carca_init = total_car * ca_er_val / (KD_car + ca_er_val)
        car_init = total_car - carca_init

        ########
        # OGB-1
        ########

        total_ogb1 = rxd_dict['ogb1_total']
        KD_ogb1 = 0.000430
        kf_ogb1 = 10.0  # 10.0
        kb_ogb1 = KD_ogb1 * kf_ogb1

        ogb1ca_init = total_ogb1 * ca_cyt_val / (KD_ogb1 + ca_cyt_val)
        ogb1_init = total_ogb1 - ogb1ca_init

        #########
        # OGB-5N
        #########
        total_ogb5 = rxd_dict['ogb5_total']
        KD_ogb5 = 0.01
        kf_ogb5 = 10.0
        kb_ogb5 = KD_ogb5 * kf_ogb5

        ogb5ca_init = total_ogb5 * ca_cyt_val / (KD_ogb5 + ca_cyt_val)
        ogb5_init = total_ogb5 - ogb5ca_init

        #######
        # IP3R
        #######

        # Bicknell and Goodhill Model
        # kf1_ip3r_bg = 50.0  # (mM^-1*msec^-1)
        # kb1_ip3r_bg = 2.5*0.001  # (msec^-1)
        # kf2_ip3r_bg = 0.035  # (mM^-1*msec^-1)
        # kb2_ip3r_bg = 1.25 * 0.001  # (msec^-1)
        # kf3_ip3r_bg = 2.5 * (12.5/3.5) / (2.5/50.0) / (1.25/0.035) * 0.001  # (mM^-1*msec^-1)
        # kb3_ip3r_bg = 2.5 * 0.001  # (msec^-1)
        # kf4_ip3r_bg = 3.5  # (mM^-1*msec^-1)
        # kb4_ip3r_bg = 12.5 * 0.001  # (msec^-1)
        # kf5_ip3r_bg = 65.0  # (mM^-1*msec^-1)
        # kb5_ip3r_bg = 10.0 * 0.001  # (msec^-1)
        # kf6_ip3r_bg = 25.0  # (mM^-1*msec^-1)
        # kb6_ip3r_bg = 25.0 * (10.0/65.0) * (0.25/10) / (2.5/50) * 0.001 # (msec^-1)
        # kf7_ip3r_bg = 10.0  # (mM^-1*msec^-1)
        # kb7_ip3r_bg = 0.25 * 0.001  # (msec^-1)
        # kf8_ip3r_bg = 0.035  # (mM^-1*msec^-1)
        # kb8_ip3r_bg = 0.035 * (0.25/10.0) * (2.5/1.25) / (0.2/0.15) * 0.001  # (msec^-1)
        # kf9_ip3r_bg = 0.15*0.001  # (msec^-1)
        # kb9_ip3r_bg = 0.2 * 0.001  # (msec^-1)
        # kf10_ip3r_bg = 1.25 * 0.001  # (mM^-1*msec^-1)
        # kb10_ip3r_bg = 2.5 * 0.001  # (msec^-1)
        # kf11_ip3r_bg = 110.0*0.001  # (msec^-1)
        # kb11_ip3r_bg = 20.0*0.001  # (msec^-1)

        # Doi et al., 2005
        kf1_ip3r = 8000.0  # (mM^-1*msec^-1)
        kb1_ip3r = 2.0  # (msec^-1)
        kf2_ip3r = 1000.0  # (mM^-1*msec^-1)
        kb2_ip3r = 25.8  # (msec^-1)
        e_ip3r = 3
        kf_ip3r = 2.22
        kb_ip3r = 2.25 * kf_ip3r
        kf3_ip3r = 4 * kf_ip3r  # (mM^-1*msec^-1)
        kb3_ip3r = kb_ip3r  # 0.005  # (msec^-1)
        kf4_ip3r = e_ip3r * 3 * kf_ip3r  # (mM^-1*msec^-1)
        kb4_ip3r = 2 * kb_ip3r  # 0.010  # (msec^-1)
        kf5_ip3r = e_ip3r ** 2 * 2 * kf_ip3r  # (mM^-1*msec^-1)
        kb5_ip3r = 3 * kb_ip3r  # 0.015  # (msec^-1)
        kf6_ip3r = e_ip3r ** 3 * kf_ip3r  # (mM^-1*msec^-1)
        kb6_ip3r = 4 * kb_ip3r  # 0.02  # (msec^-1)

        ####################
        # Calcium Extrusion
        ####################

        tot_pmca = rxd_dict['pmca_total']
        tot_ncx = rxd_dict['ncx_total']
        kf_pmca_bind = 25000.0  # (msec^-1 mM^-1)
        kb_pmca_bind = 2.0  # (msec^-1)
        kf_pmca_release = 500.0  # (msec^-1)

        kf_ncx_bind = 93.827  # (msec^-1 mM^-1)
        kb_ncx_bind = 4.0  # (msec^-1)
        kf_ncx_release = 1000.0  # (msec^-1 mM^-1)

        pmca_ca_init = tot_pmca * ca_cyt_val / ((kb_pmca_bind / kf_pmca_bind) + ca_cyt_val)
        pmca_init = tot_pmca - pmca_ca_init

        ncx_2ca_init = tot_ncx * ca_cyt_val ** 2 / ((kb_ncx_bind / kf_ncx_bind) + ca_cyt_val ** 2)
        ncx_init = tot_ncx - ncx_2ca_init

        #############################
        # Channel Conductance Values
        #############################

        g_ip3r = rxd_dict['g_ip3r']
        g_leak = rxd_dict['g_er_leak']
        g_serca = rxd_dict['g_serca']

        ##################
        #################
        # Create Species
        #################
        ##################

        self.ca = rxd.Species([self.cyt, self.er], d=ca_diff, charge=2, name='ca',
                              initial=lambda nd: ca_cyt_val if nd.region == self.cyt else ca_er_val, atolscale=1e-6)

        self.cbd_h0l0 = rxd.Species(self.cyt, initial=h0l0_init)
        self.cbd_h1l0 = rxd.Species(self.cyt, initial=h1l0_init)
        self.cbd_h2l0 = rxd.Species(self.cyt, initial=h2l0_init)
        self.cbd_h0l1 = rxd.Species(self.cyt, initial=h0l1_init)
        self.cbd_h1l1 = rxd.Species(self.cyt, initial=h1l1_init)
        self.cbd_h2l1 = rxd.Species(self.cyt, initial=h2l1_init)
        self.cbd_h0l2 = rxd.Species(self.cyt, initial=h0l2_init)
        self.cbd_h1l2 = rxd.Species(self.cyt, initial=h1l2_init)
        self.cbd_h2l2 = rxd.Species(self.cyt, initial=h2l2_init)

        cam_init_func = lambda nd: cam_init * self.check_seg_is_in_list(comb_spine_list, nd)
        self.cam = rxd.Species(self.cyt, initial=cam_init_func)
        camca_init_func = lambda nd: camca_init * self.check_seg_is_in_list(comb_spine_list, nd)
        self.camca = rxd.Species(self.cyt, initial=camca_init_func)

        cbp_init_func = lambda nd: cbp_init * self.check_seg_is_in_list(comb_spine_list, nd)
        self.cbp = rxd.Species(self.cyt, initial=cbp_init_func)
        cbpca_init_func = lambda nd: cbpca_init * self.check_seg_is_in_list(comb_spine_list, nd)
        self.cbpca = rxd.Species(self.cyt, initial=cbpca_init_func)

        if rxd_dict['ogb1_total'] > 0.0:
            self.dye = rxd.Species(self.cyt, initial=ogb1_init)
            self.dyeca = rxd.Species(self.cyt, initial=ogb1ca_init)
        elif rxd_dict['ogb5_total'] > 0.0:
            self.dye = rxd.Species(self.cyt, initial=ogb5_init)
            self.dyeca = rxd.Species(self.cyt, initial=ogb5ca_init)

        self.car = rxd.Species(self.er, initial=car_init)
        self.carca = rxd.Species(self.er, initial=carca_init)

        self.ip3 = rxd.Species(self.cyt, d=ip3_diff, charge=0,
                               initial=ip3_init)

        self.ach = rxd.Species(self.cyt, charge=0, initial=0.0)

        ##########
        # M1 AChR
        ##########

        self.r_m1 = rxd.Species(self.cyt, initial=rxd_dict['m1_init'])
        self.g_m1 = rxd.Species(self.cyt, initial=rxd_dict['g_init'])
        self.gbg_m1 = rxd.Species(self.cyt, initial=0.0)
        self.rl_m1 = rxd.Species(self.cyt, initial=0.0)
        self.rg_m1 = rxd.Species(self.cyt, initial=0.0)
        self.rgbg_m1 = rxd.Species(self.cyt, initial=0.0)
        self.rlg_m1 = rxd.Species(self.cyt, initial=0.0)
        self.rlgbg_m1 = rxd.Species(self.cyt, initial=0.0)

        self.ga_gdp_m1 = rxd.Species(self.cyt, initial=0.0)
        self.ga_gtp_m1 = rxd.Species(self.cyt, initial=3.84e-6)  # 9.64e-6)
        self.plc_m1 = rxd.Species(self.cyt, initial=rxd_dict['plc_init'])
        self.ga_gdp_plc_m1 = rxd.Species(self.cyt, initial=0.0)
        self.ga_gtp_plc_m1 = rxd.Species(self.cyt, initial=3.097e-5)  # 1.32e-5)

        #################################
        # Phosphatidylinositol Synthesis
        #################################

        self.pi = rxd.Species(self.cyt, initial=rxd_dict['pi_init'])
        self.pi4p = rxd.Species(self.cyt, initial=rxd_dict['pi4p_init'])
        self.pip2 = rxd.Species(self.cyt, charge=0,
                                initial=rxd_dict['pip2_init'])

        self.pip2_bound = rxd.Species(self.cyt,
                                      initial=rxd_dict['pip2_init'] * (fold_PIP2 - 1.0))
        self.dag = rxd.Species(self.cyt, initial=13.0)

        ########################
        # PIP2 and KCNQ Binding
        ########################

        for sec in self.somatic:
            for seg in sec:
                sec(seg.x).base_perc_kmb_inh = perc_kcnq_base_soma

        for sec in self.axonal:
            for seg in sec:
                self.axonal[0](seg.x).base_perc_kmb_inh = perc_kcnq_base_axon

        self.kcnq = rxd.Species(self.cyt, charge=0, name='kcnq',
                                initial=lambda nd: ub_kcnq_axon if (nd.satisfies(self.axonal[0]) or nd.satisfies(
                                    self.axonal[1])) else ub_kcnq_soma)
        self.pip2_kcnq = rxd.Species(self.cyt, charge=0, name='pip2_kcnq',
                                     initial=lambda nd: b_kcnq_axon if (nd.satisfies(self.axonal[0]) or nd.satisfies(
                                         self.axonal[1])) else b_kcnq_soma)

        ################
        # IP3 Breakdown
        ################
        total_ip3k = rxd_dict['ip3k_total']
        b_ip3k_ca = total_ip3k * ca_cyt_val / (KD_ip3k_ca + ca_cyt_val)
        ub_ip3k_ca = total_ip3k - b_ip3k_ca

        b_ip3k_ip3 = b_ip3k_ca * ip3_init / (KD_ip3k_ip3 + ip3_init)

        self.ip5p = rxd.Species(self.cyt, initial=lambda nd: ub_ip5p * nd.surface_area / nd.volume)  # 9.97e-4)
        self.ip5p_ip3 = rxd.Species(self.cyt, initial=lambda nd: b_ip5p * nd.surface_area / nd.volume)  # 2.32e-6)

        self.ip3k = rxd.Species(self.cyt, initial=lambda nd: ub_ip3k_ca * nd.surface_area / nd.volume)  # 0.0128)
        self.ip3k_2ca = rxd.Species(self.cyt, initial=lambda nd: b_ip3k_ca * nd.surface_area / nd.volume)  # 1.458e-7)
        self.ip3k_2ca_ip3 = rxd.Species(self.cyt,
                                        initial=lambda nd: b_ip3k_ip3 * nd.surface_area / nd.volume)  # 2.582e-9)

        # #######
        # # IP3R
        # #######
        #
        # # Doi et al., 2005
        self.r_ip3r = rxd.State(self.cyt_er_membrane, initial=0.814)
        self.ri_ip3r = rxd.State(self.cyt_er_membrane, initial=8.913e-5)
        self.ro_ip3r = rxd.State(self.cyt_er_membrane, initial=3.594e-5)
        self.rc_ip3r = rxd.State(self.cyt_er_membrane, initial=0.146)
        self.rc2_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0294)
        self.rc3_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0079)
        self.rc4_ip3r = rxd.State(self.cyt_er_membrane, initial=0.00239)

        # Bicknell and Goodhill, 2016
        # self.x1_ip3r = rxd.State(self.cyt_er_membrane, initial=0.5)
        # self.x2_ip3r = rxd.State(self.cyt_er_membrane, initial=0.3)
        # self.x3_ip3r = rxd.State(self.cyt_er_membrane, initial=0.2)
        # self.x4_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)
        # self.x5_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)
        # self.x6_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)
        # self.x7_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)
        # self.x8_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)
        # self.x9_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)
        # self.x10_ip3r = rxd.State(self.cyt_er_membrane, initial=0.0)  # Active Open State

        ####################
        # Calcium Extrusion
        ####################

        self.pmca = rxd.Species(self.cyt, initial=pmca_init)
        self.pmca_ca = rxd.Species(self.cyt, initial=pmca_ca_init)

        self.ncx = rxd.Species(self.cyt, initial=ncx_init)
        self.ncx_2ca = rxd.Species(self.cyt, initial=ncx_2ca_init)

        #############
        ############
        # Reactions
        ############
        #############

        ###########
        # M1 mAChR
        ###########

        self.L1 = rxd.Reaction(self.r_m1[self.cyt], self.rl_m1[self.cyt],
                               kf_L1 * self.ach[self.cyt], kb_L1, membrane=self.cyt)
        self.L2 = rxd.Reaction(self.rg_m1[self.cyt], self.rlg_m1[self.cyt],
                               kf_L2 * self.ach[self.cyt], kb_L2, membrane=self.cyt)
        self.L2b = rxd.Reaction(self.rgbg_m1[self.cyt], self.rlgbg_m1[self.cyt],
                                kf_L2 * self.ach[self.cyt], kb_L2, membrane=self.cyt)
        self.G1 = rxd.Reaction(self.g_m1 + self.r_m1, self.rg_m1, kf_G1, kb_G1, region=[self.cyt])
        self.G2 = rxd.Reaction(self.g_m1 + self.rl_m1, self.rlg_m1, kf_G2, kb_G2, region=[self.cyt])
        self.G1b = rxd.Reaction(self.gbg_m1 + self.r_m1, self.rgbg_m1, kf_G1, kb_G1, region=[self.cyt])
        self.G2b = rxd.Reaction(self.gbg_m1 + self.rl_m1, self.rlgbg_m1, kf_G2, kb_G2, region=[self.cyt])

        self.nx_g = rxd.Reaction(self.g_m1 > self.gbg_m1 + self.ga_gtp_m1, k_NX_G, region=[self.cyt])
        self.nx_rg = rxd.Reaction(self.rg_m1 > self.rgbg_m1 + self.ga_gtp_m1, k_NX_G, region=[self.cyt])
        self.nx_rlg = rxd.Reaction(self.rlg_m1 > self.rlgbg_m1 + self.ga_gtp_m1, k_NX_RLG, region=[self.cyt])
        self.nx_p = rxd.Reaction(self.ga_gdp_plc_m1 > self.ga_gtp_plc_m1, k_NX_P, region=[self.cyt])

        self.gtpase1 = rxd.Reaction(self.ga_gtp_m1 > self.ga_gdp_m1, k_GTPase1, region=[self.cyt])
        self.gtpase2 = rxd.Reaction(self.ga_gtp_plc_m1 > self.ga_gdp_plc_m1, k_GTPase2, region=[self.cyt])

        self.plc_assoc = rxd.Reaction(self.plc_m1 + self.ga_gtp_m1 > self.ga_gtp_plc_m1, k_PLCassoc, region=[self.cyt])
        self.plc_diss = rxd.Reaction(self.ga_gdp_plc_m1 > self.ga_gdp_m1 + self.plc_m1, k_PLCdiss, region=[self.cyt])
        self.g_reconst = rxd.Reaction(self.gbg_m1 + self.ga_gdp_m1 > self.g_m1, k_reconst, region=[self.cyt])

        #################################
        # Phosphatidylinositol Synthesis
        #################################

        self.k4 = rxd.Rate(self.pi4p, k_4k * self.pi, regions=[self.cyt])
        self.p4 = rxd.Rate(self.pi4p, -k_4p * self.pi4p, regions=[self.cyt])
        self.k5_p5 = rxd.Reaction(self.pi4p, self.pip2, k_5k, k_5p, region=[self.cyt])

        ###############################
        # IP3 Production and Breakdown
        ###############################

        self.pip2_hydrolysis_dag = rxd.Reaction(self.pip2 > self.dag,
                                                k_PLC * fold_PIP2 * self.ga_gtp_plc_m1[self.cyt],
                                                region=[self.cyt])

        scale = 602214.129
        self.ip3_param = rxd.Parameter(self.cyt, value=lambda nd: nd.surface_area / scale / nd.volume)

        self.ip3_flux = rxd.Rate(self.ip3, self.pip2 * k_PLC * fold_PIP2 * self.ga_gtp_plc_m1 * self.ip3_param,
                                 regions=[self.cyt])

        self.ip5p_binding = rxd.Reaction(self.ip3 + self.ip5p, self.ip5p_ip3, k_ip5p_f, k_ip5p_r, region=self.cyt)
        self.ip2_form = rxd.Reaction(self.ip5p_ip3 > self.ip5p, k_ip2, region=self.cyt)

        self.ip3k_cal_bind = rxd.Reaction(self.ip3k + 2 * self.ca, self.ip3k_2ca, k_ip3k_ca_f, k_ip3k_ca_r,
                                          region=self.cyt)
        self.ip3k_ip3_bind = rxd.Reaction(self.ip3k_2ca + self.ip3, self.ip3k_2ca_ip3, k_ip3k_ip3_f, k_ip3k_ip3_r,
                                          region=self.cyt)
        self.ip4_form = rxd.Reaction(self.ip3k_2ca_ip3 > self.ip3k_2ca, k_ip4, region=self.cyt)

        self.dagase = rxd.Rate(self.dag, -k_DAGase * self.dag, regions=[self.cyt])

        #####################################
        # PIP2 Buffering and Binding to KCNQ
        #####################################

        self.pip2_buffer = rxd.Reaction(self.pip2, self.pip2_bound, (fold_PIP2 - 1) * kf_pip2, kf_pip2,
                                        regions=[self.cyt])
        self.kcnq_bind_pip2 = rxd.Reaction(self.pip2 + self.kcnq, self.pip2_kcnq, kf_kcnq,
                                           kb_kcnq, regions=[self.cyt])

        ############
        # Calbindin
        ############

        self.cbdh_h0l0bindh = rxd.Reaction(self.cbd_h0l0 + self.ca, self.cbd_h1l0, 2.0 * kf_cbdh, kb_cbdh,
                                           regions=[self.cyt])
        self.cbdh_h0l1bindh = rxd.Reaction(self.cbd_h0l1 + self.ca, self.cbd_h1l1, 2.0 * kf_cbdh, kb_cbdh,
                                           regions=[self.cyt])
        self.cbdh_h0l2bindh = rxd.Reaction(self.cbd_h0l2 + self.ca, self.cbd_h1l2, 2.0 * kf_cbdh, kb_cbdh,
                                           regions=[self.cyt])
        self.cbdh_h1l0bindh = rxd.Reaction(self.cbd_h1l0 + self.ca, self.cbd_h2l0, kf_cbdh, 2.0 * kb_cbdh,
                                           regions=[self.cyt])
        self.cbdh_h1l1bindh = rxd.Reaction(self.cbd_h1l1 + self.ca, self.cbd_h2l1, kf_cbdh, 2.0 * kb_cbdh,
                                           regions=[self.cyt])
        self.cbdh_h1l2bindh = rxd.Reaction(self.cbd_h1l2 + self.ca, self.cbd_h2l2, kf_cbdh, 2.0 * kb_cbdh,
                                           regions=[self.cyt])

        self.cbdh_h0l0bindl = rxd.Reaction(self.cbd_h0l0 + self.ca, self.cbd_h0l1, 2.0 * kf_cbdl, kb_cbdl,
                                           regions=[self.cyt])
        self.cbdh_h1l0bindl = rxd.Reaction(self.cbd_h1l0 + self.ca, self.cbd_h1l1, 2.0 * kf_cbdl, kb_cbdl,
                                           regions=[self.cyt])
        self.cbdh_h2l0bindl = rxd.Reaction(self.cbd_h2l0 + self.ca, self.cbd_h2l1, 2.0 * kf_cbdl, kb_cbdl,
                                           regions=[self.cyt])
        self.cbdh_h0l1bindl = rxd.Reaction(self.cbd_h0l1 + self.ca, self.cbd_h0l2, kf_cbdl, 2.0 * kb_cbdl,
                                           regions=[self.cyt])
        self.cbdh_h1l1bindl = rxd.Reaction(self.cbd_h1l1 + self.ca, self.cbd_h1l2, kf_cbdl, 2.0 * kb_cbdl,
                                           regions=[self.cyt])
        self.cbdh_h2l1bindl = rxd.Reaction(self.cbd_h2l1 + self.ca, self.cbd_h2l2, kf_cbdl, 2.0 * kb_cbdl,
                                           regions=[self.cyt])

        ################
        # Calreticulin
        ################

        self.car_bind = rxd.Reaction(self.car + self.ca, self.carca, kf_car, kb_car, regions=[self.er])

        ##############
        # Calmodulin
        ##############

        self.cam_bind = rxd.Reaction(3 * self.ca + self.cam, self.camca, kf_cam, kb_cam, regions=[self.cyt])

        ###########################
        # Calcium Binding Protein
        ###########################

        self.cbp_bind = rxd.Reaction(self.ca + self.cbp, self.cbpca, kf_cbp, kb_cbp, regions=[self.cyt])

        ####################
        # Calcium Indicator
        ####################

        if rxd_dict['ogb1_total'] > 0:
            self.dye_bind = rxd.Reaction(self.dye + self.ca, self.dyeca, kf_ogb1, kb_ogb1, regions=[self.cyt])
        elif rxd_dict['ogb5_total'] > 0:
            self.dye_bind = rxd.Reaction(self.dye + self.ca, self.dyeca, kf_ogb5, kb_ogb5, regions=[self.cyt])

        ###################################################
        # Calcium Extrusion and Extracellular Calcium Leak
        ###################################################

        self.pmca_bind = rxd.Reaction(self.pmca + self.ca, self.pmca_ca, kf_pmca_bind, kb_pmca_bind, region=[self.cyt])
        self.pmca_release = rxd.Reaction(self.pmca_ca > self.pmca, kf_pmca_release, region=[self.cyt])

        self.ncx_bind = rxd.Reaction(self.ncx + 2 * self.ca, self.ncx_2ca, kf_ncx_bind, kb_ncx_bind, region=[self.cyt])
        self.ncx_release = rxd.Reaction(self.ncx_2ca > self.ncx, kf_ncx_release, region=[self.cyt])

        g_ext_leak = rxd_dict['g_ext_leak']
        self.leak_cyt = rxd.Rate(self.ca, g_ext_leak * (rxd_dict['ca_ext_val'] - self.ca[self.cyt]),
                                 regions=[self.cyt])

        #############
        # SERCA Flux
        #############
        not_spine = lambda nd: not self.check_seg_is_in_list(comb_spine_list, nd)
        self.scale_serca = rxd.Parameter(self.cyt_er_membrane,
                                         value=lambda nd: g_serca * nd.segment.diam / 4.0 * not_spine(nd))
        self.serca_flux = rxd.MultiCompartmentReaction(self.ca[self.cyt] > self.ca[self.er],
                                                       self.scale_serca * self.ca[self.cyt] ** 2 / (
                                                                   self.ca[self.cyt] ** 2 + 0.0013 ** 2),
                                                       custom_dynamics=True, membrane=self.cyt_er_membrane)

        ############################################
        # SOCE refill of ER during depolarizations
        ############################################
        v_init = -65.0
        soce_g = 0.00001
        er_cent_val = 0.1
        k = 0.015
        soce_rate = soce_g * rxdm.log(1 + rxdm.exp(rxd.v - v_init)) * (
            rxdm.exp(-(self.ca[self.er] - er_cent_val) / k)) / (
                            1 + rxdm.exp(-(self.ca[self.er] - er_cent_val) / k))
        self.soce = rxd.Rate(self.ca, soce_rate, regions=[self.er])

        ##########################
        # Leak From ER to Cytosol
        ##########################

        scale_leak = rxd.Parameter(self.cyt_er_membrane,
                                   value=lambda nd: not_spine(nd) * g_leak * nd.segment.diam / 4.0)
        self.leak_er = rxd.MultiCompartmentReaction(self.ca[self.er] > self.ca[self.cyt], scale_leak,
                                                    membrane=self.cyt_er_membrane)

        ##############
        # IP3R Gating
        ##############
        # Bicknell and Goodhill, 2016
        # react1 = self.ip3[self.cyt] * kf1_ip3r_bg * self.x2_ip3r - kb1_ip3r_bg * self.x5_ip3r
        # react2 = self.ca[self.cyt] * kf2_ip3r_bg * self.x5_ip3r - kb2_ip3r_bg * self.x6_ip3r
        # react3 = self.ip3[self.cyt] * kf3_ip3r_bg * self.x3_ip3r - kb3_ip3r_bg * self.x6_ip3r
        # react4 = self.ca[self.cyt] * kf4_ip3r_bg * self.x2_ip3r - kb4_ip3r_bg * self.x3_ip3r
        # react5 = self.ca[self.cyt] * kf5_ip3r_bg * self.x4_ip3r - kb5_ip3r_bg * self.x5_ip3r
        # react5b = self.ca[self.cyt] * kf5_ip3r_bg * self.x7_ip3r - kb5_ip3r_bg * self.x8_ip3r
        # react6 = self.ca[self.cyt] * kf6_ip3r_bg * self.x1_ip3r - kb6_ip3r_bg * self.x2_ip3r
        # react7 = self.ip3[self.cyt] * kf7_ip3r_bg * self.x1_ip3r - kb7_ip3r_bg * self.x4_ip3r
        # react8 = self.ca[self.cyt] * kf8_ip3r_bg * self.x8_ip3r - kb8_ip3r_bg * self.x9_ip3r
        # react9 = kf9_ip3r_bg * self.x4_ip3r - kb9_ip3r_bg * self.x7_ip3r
        # react9b = kf9_ip3r_bg * self.x5_ip3r - kb9_ip3r_bg * self.x8_ip3r
        # react10 = kf10_ip3r_bg * self.x6_ip3r - kb10_ip3r_bg * self.x9_ip3r
        # react11 = kf11_ip3r_bg * self.x5_ip3r - kb11_ip3r_bg * self.x10_ip3r
        #
        # self.dx1_ip3r = rxd.Rate(self.x1_ip3r, -react6 - react7, regions=[self.cyt_er_membrane])
        # self.dx2_ip3r = rxd.Rate(self.x2_ip3r, react6 - react1 - react4, regions=[self.cyt_er_membrane])
        # self.dx3_ip3r = rxd.Rate(self.x3_ip3r, react4 - react3, regions=[self.cyt_er_membrane])
        # self.dx4_ip3r = rxd.Rate(self.x4_ip3r, react7 - react5 - react9, regions=[self.cyt_er_membrane])
        # self.dx5_ip3r = rxd.Rate(self.x5_ip3r, react5 + react1 - react9b - react2 - react11, regions=[self.cyt_er_membrane])
        # self.dx6_ip3r = rxd.Rate(self.x6_ip3r, react2 + react3 - react10, regions=[self.cyt_er_membrane])
        # self.dx7_ip3r = rxd.Rate(self.x7_ip3r, react9 - react5b, regions=[self.cyt_er_membrane])
        # self.dx8_ip3r = rxd.Rate(self.x8_ip3r, react5b + react9b - react8, regions=[self.cyt_er_membrane])
        # self.dx9_ip3r = rxd.Rate(self.x9_ip3r, react8 + react10, regions=[self.cyt_er_membrane])
        # self.dx10_ip3r = rxd.Rate(self.x10_ip3r, react11, regions=[self.cyt_er_membrane])

        # Doi et al., 2005
        react1 = self.ri_ip3r * self.ca[self.cyt] * kf1_ip3r - kb1_ip3r * self.ro_ip3r
        react2 = self.r_ip3r * self.ip3[self.cyt] * kf2_ip3r - kb2_ip3r * self.ri_ip3r
        react3 = self.r_ip3r * self.ca[self.cyt] * kf3_ip3r - self.rc_ip3r * kb3_ip3r
        react4 = self.rc_ip3r * self.ca[self.cyt] * kf4_ip3r - self.rc2_ip3r * kb4_ip3r
        react5 = self.rc2_ip3r * self.ca[self.cyt] * kf5_ip3r - self.rc3_ip3r * kb5_ip3r
        react6 = self.rc3_ip3r * self.ca[self.cyt] * kf6_ip3r - self.rc4_ip3r * kb6_ip3r

        self.dr_ip3r = rxd.Rate(self.r_ip3r, - react2 - react3, regions=[self.cyt_er_membrane])
        self.dri_ip3r = rxd.Rate(self.ri_ip3r, -react1 + react2, regions=[self.cyt_er_membrane])
        self.dro_ip3r = rxd.Rate(self.ro_ip3r, react1, regions=[self.cyt_er_membrane])
        self.drc_ip3r = rxd.Rate(self.rc_ip3r, react3 - react4, regions=[self.cyt_er_membrane])
        self.drc2_ip3r = rxd.Rate(self.rc2_ip3r, react4 - react5, regions=[self.cyt_er_membrane])
        self.drc3_ip3r = rxd.Rate(self.rc3_ip3r, react5 - react6, regions=[self.cyt_er_membrane])
        self.drc4_ip3r = rxd.Rate(self.rc4_ip3r, react6, regions=[self.cyt_er_membrane])

        self.ip3r_dense = rxd.Parameter(self.cyt_er_membrane,
                                        value=lambda nd: 1.2 * nd.segment.diam / 4.0 if nd.satisfies(
                                            self.somatic[0]) else not_spine(nd) * nd.segment.diam / 4.0)
        # k_ip3r = g_ip3r * (4*(1.0-self.x10_ip3r)*self.x10_ip3r**3 + self.x10_ip3r**4) * self.ip3r_dense  # For Bicknell and Goodhill 2016
        k_ip3r = g_ip3r * self.ip3r_dense * self.ro_ip3r  # For Doi et al., 2005
        self.ip3r_flux = rxd.MultiCompartmentReaction(self.ca[self.er], self.ca[self.cyt], k_ip3r,
                                                      membrane=self.cyt_er_membrane)

        print('RXD Setup Complete')

    def get_line_segs(self):

        seg_dict = {}

        min_diam = 1000
        max_diam = 0

        for sec in self.c.all:

            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': 'k', 'pts': np.array(pts), 'diam': sec.diam}

        return seg_dict

    def get_lines_for_mayavi(self):

        seg_dict = {}

        color_dict = {'somatic': 'k',
                      'basal': {'prox': 'g',
                                'dist': '#355e3b'},
                      'apical': {'sr': {'thin': '#000080',
                                        'thick_prox': '#89cff0',
                                        'thick_dist': '#4F97A3',
                                        'thick_med': 'c'},
                                 'slm': {'thin': '#CD2626',
                                         'med': '#FF3300',
                                         'thick': '#FF8000'},
                                 'trunk': {'all': 'k'},
                                 'obl_base': {'all': 'c'}
                                 }
                      }

        x_list = []
        y_list = []
        z_list = []
        d_list = []
        line_ind = []

        pt_ind = 0

        for sec in self.somatic:

            color = color_dict['somatic']

            n_pts = h.n3d(sec=sec)
            for i in range(int(n_pts)):
                x_list.append(h.x3d(i, sec=sec))
                y_list.append(h.y3d(i, sec=sec))
                z_list.append(h.z3d(i, sec=sec))
                d_list.append(h.diam3d(i, sec=sec))

            line_ind.append(np.vstack([np.arange(pt_ind, pt_ind + n_pts - 1),
                                       np.arange(pt_ind + 1, pt_ind + n_pts)]).T)

            pt_ind += n_pts

        for sec in self.apical:
            color = color_dict['somatic']

            n_pts = h.n3d(sec=sec)
            for i in range(int(n_pts)):
                x_list.append(h.x3d(i, sec=sec))
                y_list.append(h.y3d(i, sec=sec))
                z_list.append(h.z3d(i, sec=sec))
                d_list.append(h.diam3d(i, sec=sec))

            line_ind.append(np.vstack([np.arange(pt_ind, pt_ind + n_pts - 1),
                                       np.arange(pt_ind + 1, pt_ind + n_pts)]).T)

            pt_ind += n_pts

        for sec in self.basal:
            color = color_dict['somatic']

            n_pts = h.n3d(sec=sec)
            for i in range(int(n_pts)):
                x_list.append(h.x3d(i, sec=sec))
                y_list.append(h.y3d(i, sec=sec))
                z_list.append(h.z3d(i, sec=sec))
                d_list.append(h.diam3d(i, sec=sec))

            line_ind.append(np.vstack([np.arange(pt_ind, pt_ind + n_pts - 1),
                                       np.arange(pt_ind + 1, pt_ind + n_pts)]).T)

            pt_ind += n_pts

        x = np.hstack(x_list)
        y = np.hstack(y_list)
        z = np.hstack(z_list)
        s = np.hstack(d_list)
        connections = np.vstack(line_ind)

        return x, y, z, s, connections

    def get_line_segs_by_type(self):

        seg_dict = {}

        color_dict = {'somatic': 'b',
                      'basal': 'r',
                      'apical': 'orange',
                      'axon': 'g'
                      }

        min_diam = 1000
        max_diam = 0

        for sec in self.c.somatic:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            color = color_dict['somatic']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        for sec in self.c.apical:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            color = color_dict['apical']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        for sec in self.c.basal:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            color = color_dict['basal']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        # for sec in self.axonal:
        #     pts = []
        #     sname = sec.name()
        #     m_list = re.findall(r'\d*', sname)
        #     m_filt = [n_str for n_str in m_list if n_str]
        #     snum = m_filt[-1]
        #
        #     color = color_dict['axon']
        #
        #     for i in range(int(h.n3d(sec=sec))):
        #         pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])
        #
        #     seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        return seg_dict

    def get_line_segs_by_type(self):

        seg_dict = {}

        color_dict = {'somatic': 'b',
                      'basal': 'r',
                      'apical': 'orange',
                      'axon': 'g'
                      }

        min_diam = 1000
        max_diam = 0

        for sec in self.c.somatic:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            color = color_dict['somatic']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        for sec in self.c.apical:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            color = color_dict['apical']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        for sec in self.c.basal:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            color = color_dict['basal']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        # for sec in self.axonal:
        #     pts = []
        #     sname = sec.name()
        #     m_list = re.findall(r'\d*', sname)
        #     m_filt = [n_str for n_str in m_list if n_str]
        #     snum = m_filt[-1]
        #
        #     color = color_dict['axon']
        #
        #     for i in range(int(h.n3d(sec=sec))):
        #         pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])
        #
        #     seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        return seg_dict

    def get_line_segs_by_rxd(self):

        seg_dict = {}

        color_dict = {'rxd': 'r',
                      'not': 'k',
                      }

        min_diam = 1000
        max_diam = 0

        for sec in self.c.somatic:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            if sname in self.cal_names:
                color = color_dict['rxd']
            else:
                color = color_dict['not']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        for sec in self.c.apical:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            if sname in self.cal_names:
                color = color_dict['rxd']
            else:
                color = color_dict['not']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        for sec in self.c.basal:
            pts = []
            sname = sec.name()
            m_list = re.findall(r'\d*', sname)
            m_filt = [n_str for n_str in m_list if n_str]
            snum = m_filt[-1]

            if sname in self.cal_names:
                color = color_dict['rxd']
            else:
                color = color_dict['not']

            for i in range(int(h.n3d(sec=sec))):
                pts.append([h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)])

            seg_dict[sname] = {'num': snum, 'color': color, 'pts': np.array(pts), 'diam': sec.diam}

        return seg_dict

    def get_sec_by_layer(self, target_layer):
        """
        Returns a neuron section that is within the target layer
        :param target_layer: layer of ca1 which will be used in
        :return: randomly chosen neuron section that resides within target layer, target_layer
        """

        if target_layer in self.layer_dict:
            rand_num = np.random.uniform()

            sec_i = np.digitize(rand_num, self.sp_prob_dict[target_layer])

            return self.get_sec_by_name(self.layer_dict[target_layer]['sec_names'][sec_i])
        else:
            print
            "Invalid layer specified: {0}".format(target_layer)
            exit(1)

    def get_sec_by_name(self, sec_name):
        """
        Returns a section if it has the given name
        :param sec_name: name of section you desire with or without the cell name
        :return: section
        """

        l_dict = self.layer_dict

        if sec_name in l_dict['axon']['sec_names']:
            return self.axon[l_dict['axon']['sec_names'].index(sec_name)]
        elif sec_name in l_dict['so']['sec_names']:
            return self.basal[l_dict['so']['sec_names'].index(sec_name)]
        elif sec_name in l_dict['apical']['sec_names']:
            return self.apical[l_dict['apical']['sec_names'].index(sec_name)]
        elif sec_name in l_dict['soma']['sec_names']:
            return self.somatic[l_dict['soma']['sec_names'].index(sec_name)]

        sec_split = sec_name.split('_')

        if not hasattr(self, 'soma_base'):
            soma_name = self.somatic[0].name()
            # FLAG
            apical_name = self.apical[0].name()
            basal_name = self.basal[0].name()

            self.bracket_surround = soma_name[-1] == ']'
            self.soma_base = soma_name.rstrip('[]0123456789')
            # FLAG
            self.apical_base = apical_name.rstrip('[]0123456789')
            self.basal_base = basal_name.rstrip('[]0123456789')

        if self.bracket_surround:
            num = '[' + sec_split[1] + ']'
        else:
            num = sec_split[1]

        if sec_split[0] in 'somatic':
            rev_name = self.soma_base + num
        elif sec_split[0] in 'apical':
            rev_name = self.apical_base + num
        elif sec_split[0] in 'basal':
            rev_name = self.basal_base + num

        if rev_name in l_dict['axon']['sec_names']:
            return self.axon[l_dict['axon']['sec_names'].index(rev_name)]
        elif rev_name in l_dict['so']['sec_names']:
            return self.basal[l_dict['so']['sec_names'].index(rev_name)]
        elif rev_name in l_dict['apical']['sec_names']:
            return self.apical[l_dict['apical']['sec_names'].index(rev_name)]
        elif rev_name in l_dict['soma']['sec_names']:
            return self.somatic[l_dict['soma']['sec_names'].index(rev_name)]

        print(rev_name)
        print('{0} unknown section'.format(sec_name))
        return 0

    def load_asc(self):

        self.c = h.mkcell_asc(self.morph_file)

        for sec in self.c.somatic:
            self.somatic.append(sec)
        for sec in self.c.apical:
            self.apical.append(sec)
        for sec in self.c.basal:
            self.basal.append(sec)

    def load_swc(self):

        self.c = h.mkcell_swc(self.morph_file)

        self.c.somatic = list(self.c.somatic)
        self.c.apical = list(self.c.apical)
        self.c.basal = list(self.c.apical)

    def make_synapse_sections(self, par_section, input_var, loc, syn_type):
        syn_params = self.config_obj[syn_type]

        sec_name = par_section.name()

        if sec_name not in self.syns:
            self.syns[sec_name] = {}

        i_syn = len(self.syns[sec_name])

        self.syns[sec_name][i_syn] = {}

        # Make Post Synaptic Density (PSD)
        # self.syns[sec_name][i_syn]['psd'] = h.Section(name='{0}_psd_{1}'.format(sec_name, i_syn))
        # diam and L chosen to achieve Thin Spine Surface Area and Volume as measured by Stewart et al., 2005
        # PSD: Surface Area = 0.132 um^2, Volume = 0.0032 um^3

        # self.syns[sec_name][i_syn]['psd'].diam = 0.0032
        # self.syns[sec_name][i_syn]['psd'].L =  0.132
        # self.syns[sec_name][i_syn]['psd'].cm = syn_params['cm']
        # self.syns[sec_name][i_syn]['psd'].Ra = syn_params['Ra']
        #
        # self.syns[sec_name][i_syn]['psd'].insert('pas')
        # self.syns[sec_name][i_syn]['psd'].g_pas = syn_params['g_pas']
        # self.syns[sec_name][i_syn]['psd'].e_pas = syn_params['e_pas']

        # Make Spine Head
        self.syns[sec_name][i_syn]['head'] = h.Section(name='{0}_head_{1}'.format(sec_name, i_syn))
        # diam and L chosen to achieve Mushroom Spine Surface Area and Volume as measured by Stewart et al., 2005
        # Head: Surface Area = 0.671 um^2, Volume = 0.04 um^3
        self.syns[sec_name][i_syn]['head'].diam = 0.297
        self.syns[sec_name][i_syn]['head'].L = 0.720
        self.syns[sec_name][i_syn]['head'].cm = syn_params['cm']
        self.syns[sec_name][i_syn]['head'].Ra = syn_params['Ra']

        self.syns[sec_name][i_syn]['head'].insert('pas')
        self.syns[sec_name][i_syn]['head'].g_pas = syn_params['g_pas']
        self.syns[sec_name][i_syn]['head'].e_pas = syn_params['e_pas']

        self.syns[sec_name][i_syn]['head'].insert('cat')
        self.syns[sec_name][i_syn]['head'].gcatbar_cat = syn_params['gcatbar_cat']

        # Make Spine Neck
        self.syns[sec_name][i_syn]['neck'] = h.Section(name='{0}_neck_{1}'.format(sec_name, i_syn))
        self.syns[sec_name][i_syn]['neck'].diam = 0.1
        self.syns[sec_name][i_syn]['neck'].L = 0.5
        self.syns[sec_name][i_syn]['neck'].cm = syn_params['cm']
        self.syns[sec_name][i_syn]['neck'].Ra = syn_params['Ra']

        self.syns[sec_name][i_syn]['neck'].insert('pas')
        self.syns[sec_name][i_syn]['neck'].g_pas = syn_params['g_pas']
        self.syns[sec_name][i_syn]['neck'].e_pas = syn_params['e_pas']

        # self.syns[sec_name][i_syn]['psd'].connect(self.syns[sec_name][i_syn]['head'], 1.0, 0.0)
        self.syns[sec_name][i_syn]['head'].connect(self.syns[sec_name][i_syn]['neck'], 1.0, 0.0)
        self.syns[sec_name][i_syn]['neck'].connect(par_section, loc, 0.0)

        # self.cal_list.append(self.syns[sec_name][i_syn]['psd'])
        self.cal_list.append(self.syns[sec_name][i_syn]['head'])

        curr_sec = self.syns[sec_name][i_syn]['neck']
        self.cal_list.append(curr_sec)
        self.cal_names.append(curr_sec.name())

        while h.SectionRef(sec=curr_sec).parent not in self.cal_list:
            curr_sec = h.SectionRef(sec=curr_sec).parent
            self.cal_list.append(curr_sec)
            self.cal_names.append(curr_sec.name())

        # Create Mechanisms for AMPA and NMDA
        # Due to AMPA and NMDA mechanisms being structured as a Point Process, they are saved separately
        # and not inserted
        # self.syns[sec_name][i_syn]['ampa'] = h.AMPAGluR2flip(0.5, sec=self.syns[sec_name][i_syn]['head'])
        # self.syns[sec_name][i_syn]['ampa'].nbAMPAR = syn_params['nbampa']
        self.syns[sec_name][i_syn]['nmdanr2a'] = h.NMDANR2A(0.5, sec=self.syns[sec_name][i_syn]['head'])
        self.syns[sec_name][i_syn]['nmdanr2a'].nbNMDAR = syn_params['nbnmdanr2a']
        # self.syns[sec_name][i_syn]['nmdanr2b'] = h.NMDANR2B(0.5, sec=self.syns[sec_name][i_syn]['psd'])
        # self.syns[sec_name][i_syn]['nmdanr2b'].nbNMDAR = syn_params['nbnmdanr2b']

        # Connect mechanisms to allow activation using Netcon
        # diffusion mechanism
        self.syns[sec_name][i_syn]['diff'] = h.Diffusion(loc, sec=par_section)
        # This Netcon stimulates neurotransmitter calculation in the diffusion mechanism
        self.syns[sec_name][i_syn]['diff_netcon'] = h.NetCon(input_var,
                                                             self.syns[sec_name][i_syn]['diff'],
                                                             syn_params['diff_threshold'],
                                                             syn_params['delay'],
                                                             1.0
                                                             )

        self.syns[sec_name][i_syn]['diff'].distanceAMPA = syn_params['dist_ampa']
        self.syns[sec_name][i_syn]['diff'].distanceNMDA_A = syn_params['dist_nmdanr2a']
        # self.syns[sec_name][i_syn]['diff'].distanceNMDA_B = syn_params['dist_nmdanr2b']

        # Connect Glu to AMPA mechanism
        # h.setpointer(self.syns[sec_name][i_syn]['diff']._ref_NTconcAMPA,
        #              'Glu',
        #              self.syns[sec_name][i_syn]['ampa'])
        # Connect Glu to NMDA NR2A mechanism
        h.setpointer(self.syns[sec_name][i_syn]['diff']._ref_NTconcNMDA_A,
                     'Glu',
                     self.syns[sec_name][i_syn]['nmdanr2a'])
        # # Connect Glu to NMDA NR2B mechanism
        # h.setpointer(self.syns[sec_name][i_syn]['diff']._ref_NTconcNMDA_B,
        #              'Glu',
        #              self.syns[sec_name][i_syn]['nmdanr2b'])

    def print_number_of_sections(self):
        print('{0} apical sections'.format(len(self.apical)))
        print('{0} basal sections'.format(len(self.basal)))
        print('{0} somatic sections'.format(len(self.somatic)))
        print('{0} total sections'.format(len(self.all)))

    def plot_cell2d(self, pretty=False, with_sec_names=False, simplify=False, bullet=False):
        """

        :param pretty:
        :param with_sec_names:
        :param simplify:
        :param bullet:
        :return:
        """

        fig = plt.figure()
        mid_pt_shift = np.array([1, 0, 0])

        if with_sec_names:
            text_dict = {}

        if not pretty:
            try:
                line_dict = self.get_line_segs()
            except Exception:
                print("Error when getting line segments")

            for sec_n, sec_d in line_dict.items():
                sec_pts = sec_d['pts']

                if simplify:
                    plt.plot(sec_pts[[0, -1], 0], sec_pts[[0, -1], 1], color='k')
                else:
                    plt.plot(sec_pts[:, 0], sec_pts[:, 1], color='k')

                if with_sec_names:
                    if simplify:
                        mid_pt = sec_pts[0] + 0.5 * (sec_pts[-1] - sec_pts[0]) + mid_pt_shift
                    else:
                        mid_i = int(np.floor(sec_pts.shape[0] / 2.0))
                        mid_pt = sec_pts[mid_i, :] + mid_pt_shift

                    text_dict[sec_n] = plt.text(mid_pt[0], mid_pt[1], sec_d['num'], fontsize=12)
        else:

            try:
                line_dict = self.get_line_segs_pretty()
            except Exception:
                print("Error when getting pretty line segments")

            for sec_n, sec_d in line_dict.items():
                sec_pts = sec_d['pts']

                if simplify:
                    plt.plot(sec_pts[[0, -1], 0], sec_pts[[0, -1], 1], color=sec_d['color'])
                else:
                    plt.plot(sec_pts[:, 0], sec_pts[:, 1], color=sec_d['color'])

                if with_sec_names:
                    if simplify:
                        mid_pt = sec_pts[0] + 0.5 * (sec_pts[-1] - sec_pts[0]) + mid_pt_shift
                    else:
                        mid_i = int(np.floor(sec_pts.shape[0] / 2.0))
                        mid_pt = sec_pts[mid_i, :] + mid_pt_shift

                    text_dict[sec_n] = plt.text(mid_pt[0], mid_pt[1], sec_d['num'], color=sec_d['color'], fontsize=12)

            plt.axis('scaled')
            plt.draw()

        return fig

    def plot_cell(self, color_code='one color'):

        fig = plt.figure()
        ax = Axes3D(fig)

        if color_code == 'one color':
            try:
                lines = self.get_line_segs()
            except Exception:
                print("Error when getting 2D line segments")

        elif color_code == 'type':
            try:
                lines = self.get_line_segs_by_type()
            except Exception:
                print("Error when getting 3D line segments")

        elif color_code == 'layer':
            lines = self.get_line_segs_by_layer()
        elif color_code == 'rxd':
            lines = self.get_line_segs_by_rxd()

        for sec_name, sec_dict in lines.items():
            ax.plot(sec_dict['pts'][:, 0], sec_dict['pts'][:, 1], zs=sec_dict['pts'][:, 2], color=sec_dict['color'],
                    linewidth=sec_dict['diam'])

    def sort_sections(self, dist=0):
        """

        :param dist: distance in micrometers from the soma along cell's orientation vector (self.cell_vect) that will\n
        be the dividing line for sorting sections into stratum radium or stratum lacunosum moleculare
        :return: N/A
        """

        lay_dict = {'sr': {'sec_names': [], 'subclass': []},
                    'cross': {'sec_names': [], 'subclass': []},
                    'slm': {'sec_names': [], 'subclass': []},
                    'soma': {'sec_names': [], 'subclass': []},
                    'apical': {'sec_names': []},
                    'so': {'sec_names': [], 'subclass': []},
                    'trunk': {'sec_names': []},
                    'oblique_base': {'sec_names': []},
                    'base_dict': {},
                    'trunk_dict': {},
                    'axon': {}
                    }

        # tree_vec = np.array([(self.cell_vect[0]-self.cell_orgn[0]),
        #                      (self.cell_vect[1]-self.cell_orgn[1]),
        #                      (self.cell_vect[2]-self.cell_orgn[2])])

        tree_vec = self.cell_vect

        if not dist:
            max_proj = 0
            for sec in self.apical:
                n = h.n3d(sec=sec) - 1

                pt_n = np.subtract(np.array([h.x3d(n, sec=sec), h.y3d(n, sec=sec), h.z3d(n, sec=sec)]), self.cell_orgn)
                amp_n = np.linalg.norm(pt_n)
                proj_n = np.dot(pt_n, tree_vec)  # * amp_n

                if proj_n > max_proj:
                    max_proj = proj_n
            dist = 0.666 * max_proj

        lay_dict['axon']['sec_names'] = [sec.name() for sec in self.axonal]

        for sec in self.somatic:
            lay_dict['soma']['sec_names'].append(sec.name())

        for sec in self.basal:
            lay_dict['so']['sec_names'].append(sec.name())

            parent = h.SectionRef(sec=sec).parent()

            if 'soma' in parent.sec.name():
                lay_dict['so']['subclass'].append('prox')
            else:
                lay_dict['so']['subclass'].append('dist')

        for sec in self.apical:
            lay_dict['apical']['sec_names'].append(sec.name())

            n = h.n3d(sec=sec) - 1

            pt_0 = np.subtract(np.array([h.x3d(0, sec=sec), h.y3d(0, sec=sec), h.z3d(0, sec=sec)]), self.cell_orgn)
            pt_n = np.subtract(np.array([h.x3d(n, sec=sec), h.y3d(n, sec=sec), h.z3d(n, sec=sec)]), self.cell_orgn)

            amp_0 = np.linalg.norm(pt_0)
            amp_n = np.linalg.norm(pt_n)

            proj_0 = np.dot(pt_0, tree_vec)
            proj_n = np.dot(pt_n, tree_vec)

            if proj_0 < dist and proj_n < dist:
                lay_dict['sr']['sec_names'].append(sec.name())

                sec_diam = h.diam3d(0, sec=sec)

                nchild = h.SectionRef(sec=sec).nchild()

                if sec_diam < 1 or nchild == 0:
                    lay_dict['sr']['subclass'].append('thin')

                else:
                    if proj_n < 100:
                        lay_dict['sr']['subclass'].append('thick_prox')
                    elif proj_n < 250:
                        lay_dict['sr']['subclass'].append('thick_med')
                    else:
                        lay_dict['sr']['subclass'].append('thick_dist')

            elif proj_0 < dist and proj_n > dist:
                lay_dict['cross']['sec_names'].append(sec.name())
            elif proj_0 > dist and proj_n > dist:

                lay_dict['slm']['sec_names'].append(sec.name())

                sec_diam = h.diam3d(0, sec=sec)

                nchild = h.SectionRef(sec=sec).nchild()

                if nchild == 0:
                    lay_dict['slm']['subclass'].append('thin')

                elif sec_diam < 1:
                    lay_dict['slm']['subclass'].append('med')
                else:
                    lay_dict['slm']['subclass'].append('thick')

        if lay_dict['slm']['sec_names']:

            t_slm_secn = lay_dict['slm']['sec_names'][0]

            t_slm_sec = self.apical[lay_dict['apical']['sec_names'].index(t_slm_secn)]

            t_sec = t_slm_sec

            while h.SectionRef(sec=t_sec).has_parent() and t_sec.name() not in lay_dict['soma']['sec_names']:
                lay_dict['trunk']['sec_names'].append(t_sec.name())

                # if t_sec.name() in lay_dict['sr']['sec_names'] or t_sec.name() in lay_dict['cross']['sec_names']:
                #     lay_dict['trunk']['sec_names'].append(t_sec.name())

                t_sec = h.SectionRef(sec=t_sec).parent

        else:
            print('No SLM sections found')

        for sec in self.apical:

            s_name = sec.name()
            if s_name not in lay_dict['trunk']['sec_names']:

                p_sec = h.SectionRef(sec=sec).parent

                if p_sec.name() in lay_dict['trunk']['sec_names']:
                    lay_dict['oblique_base']['sec_names'].append(s_name)

                else:
                    while p_sec.name() not in lay_dict['trunk']['sec_names'] and p_sec.name() not in lay_dict['soma'][
                        'sec_names']:
                        prev_sec = p_sec
                        p_sec = h.SectionRef(sec=p_sec).parent

                    lay_dict['base_dict'][s_name] = prev_sec.name()

                lay_dict['trunk_dict'][s_name] = p_sec.name()

        return lay_dict


def run_current_injection(targ_sec, rxd_sim, param_dict={}, sim_dur=500, c_int=[50, 100], c_amp=1):
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

    cell = MyCell(t_path, rxd_sim, param_dict)

    print('Cell')

    count = 0

    for sec in h.allsec():
        count += 1
    print('Number of Sections: {0}'.format(count))

    t_sec = cell.get_sec_by_name(targ_sec)

    t_curr = h.IClamp(t_sec(0.5))

    t_curr.delay = c_int[0]
    t_curr.amp = c_amp
    t_curr.dur = c_int[1] - c_int[0]

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values
    i_rec = h.Vector().record(t_curr._ref_i)
    t = h.Vector().record(h._ref_t)
    soma_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    apic_v = h.Vector().record(cell.apical[9](0.5)._ref_v)

    if not rxd_sim:
        s_ca_cyt = h.Vector().record(cell.somatic[0](0.5)._ref_cai)
        a_ca_cyt = h.Vector().record(cell.apical[9](0.5)._ref_cai)
    else:
        s_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.somatic[0])[0]._ref_concentration)
        a_ca_cyt = h.Vector().record(cell.ca[cell.cyt].nodes(cell.apical[9])[0]._ref_concentration)
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

    h.v_init = -68.0

    h.tstop = sim_dur
    h.celsius = 34.0

    # h.load_file('negative_init.hoc')
    # h.init()

    h.stdinit()

    print('Running current injection, amplitude = {0} nA'.format(c_amp))

    h.continuerun(sim_dur)

    res_dict = {'t': np.array(t.as_numpy()),
                'i_stim': np.array(i_rec.as_numpy()),
                'soma_v': np.array(soma_v),
                'soma_cyt': np.array(s_ca_cyt),
                'apic9_v': np.array(apic_v),
                'apic9_cyt': np.array(a_ca_cyt),

                }

    cell.delete_sections()

    if rxd_sim:
        res_dict['soma_er'] = np.array(s_ca_er)
        res_dict['apic9_er'] = np.array(a_ca_er)
        res_dict['soma_dyeca'] = np.array(s_dyeca)
        res_dict['apic9_dyeca'] = np.array(s_dyeca)
        res_dict['soma_cbdhca'] = np.array(s_cbdhca)
        res_dict['apic9_cbdhca'] = np.array(s_cbdhca)
        res_dict['soma_cbdlca'] = np.array(s_cbdlca)
        res_dict['apic9_cbdlca'] = np.array(s_cbdlca)
        res_dict['soma_carca'] = np.array(s_carca)
        res_dict['apic9_carca'] = np.array(a_carca)
    return res_dict


def test_input_resistance(target_secs, current_values, rxd_sim, param_dict={}, sim_dur=500, c_int=[50, 100]):
    r_dict = {}

    peak_dp = np.zeros((len(current_values), len(target_secs)))

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = MyCell(t_path, rxd_sim, param_dict)

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values

    t = h.Vector().record(h._ref_t)
    soma_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    apic_v = h.Vector().record(cell.apical[9](0.5)._ref_v)

    h.cvode.active(1)
    h.v_init = -68.0

    h.tstop = sim_dur
    h.celsius = 34.0

    r_dict = {}

    for t_i, t_sec in enumerate(target_secs):

        r_dict[t_i] = {}
        t_sec = cell.get_sec_by_name(t_sec)

        t_curr = h.IClamp(t_sec(0.5))

        t_curr.delay = c_int[0]

        t_curr.dur = c_int[1] - c_int[0]

        i_rec = h.Vector().record(t_curr._ref_i)

        for c_i, c_val in enumerate(current_values):
            t_curr.amp = c_val

            h.stdinit()

            print('Running current injection, amplitude = {0} nA'.format(c_val))

            h.continuerun(sim_dur)

            r_dict[t_i][c_i] = {'t': np.array(t.as_numpy()),
                                'i_stim': np.array(i_rec.as_numpy()),
                                'soma_v': np.array(soma_v),
                                'apic9_v': np.array(apic_v),
                                }

            i_stim = np.squeeze(
                np.intersect1d(np.where(r_dict[t_i][c_i]['t'] > c_int[0]), np.where(r_dict[t_i][c_i]['t'] < c_int[1])))

            peak_dp[c_i, t_i] = np.max(r_dict[t_i][c_i]['soma_v'][i_stim]) - r_dict[t_i][c_i]['soma_v'][i_stim[0] - 1]

    dep_vs_i_fig = bplt.figure(title='Depolarization vs Current Amplitude, RXD={0}'.format(rxd_sim))
    dep_vs_i_fig.xaxis.axis_label = 'Current Amplitude (nA)'
    dep_vs_i_fig.yaxis.axis_label = 'Membrane Depolarization (mV)'

    for t_i, t_sec in enumerate(target_secs):
        dep_vs_i_fig.line(current_values, peak_dp[:, t_i], color=colrs[t_i], line_width=3,
                          legend='{}[0.5]'.format(t_sec))
        dep_vs_i_fig.circle(current_values, peak_dp[:, t_i], color=colrs[t_i], line_width=3)

    return dep_vs_i_fig


def test_spike_frequency_vs_current(target_secs, current_values, rxd_sim, param_dict={}, sim_dur=500, c_int=[50, 100]):
    r_dict = {}

    spike_c = np.zeros((len(current_values), len(target_secs)))

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')

    cell = MyCell(t_path, rxd_sim, param_dict)

    # apical[9][0.5] distance to soma[0][1.0] = 105.53 um

    # Record Values

    t = h.Vector().record(h._ref_t)
    soma_v = h.Vector().record(cell.somatic[0](0.5)._ref_v)
    apic_v = h.Vector().record(cell.apical[9](0.5)._ref_v)

    sp_vec = h.Vector()
    sp_count = h.NetCon(cell.somatic[0](0.5)._ref_v, None, sec=cell.somatic[0])
    sp_count.threshold = 0.0
    sp_count.record(sp_vec)

    h.cvode.active(1)
    h.v_init = -68.0

    h.tstop = sim_dur
    h.celsius = 34.0

    r_dict = {}

    for t_i, t_sec in enumerate(target_secs):

        r_dict[t_i] = {}
        t_sec = cell.get_sec_by_name(t_sec)

        t_curr = h.IClamp(t_sec(0.5))
        t_curr.delay = c_int[0]

        t_curr.dur = c_int[1] - c_int[0]

        i_rec = h.Vector().record(t_curr._ref_i)

        for c_i, c_val in enumerate(current_values):
            t_curr.amp = c_val

            h.stdinit()

            print('Running current injection, amplitude = {0} nA'.format(c_val))

            h.continuerun(sim_dur)

            r_dict[t_i][c_i] = {'spike_times': np.array(sp_vec),
                                't': np.array(t),
                                'soma_v': np.array(soma_v),
                                'apic9_v': np.array(apic_v)
                                }
            print(r_dict[t_i][c_i]['spike_times'].shape[0])

            spike_c[c_i, t_i] = r_dict[t_i][c_i]['spike_times'].shape[0]

    dep_vs_i_fig = bplt.figure(title='Number of Action Potentials vs Current Amplitude, RXD={0}'.format(rxd_sim))
    dep_vs_i_fig.xaxis.axis_label = 'Current Amplitude (nA)'
    dep_vs_i_fig.yaxis.axis_label = 'Membrane Depolarization (mV)'

    for t_i, t_sec in enumerate(target_secs):
        dep_vs_i_fig.line(current_values, spike_c[:, t_i], color=colrs[t_i], line_width=3,
                          legend='{}[0.5]'.format(t_sec))
        dep_vs_i_fig.circle(current_values, spike_c[:, t_i], color=colrs[t_i], line_width=3)

    soma_v_fig = bplt.figure(title='Soma Membrane Potential vs Time, RXD = {0}'.format(rxd_sim))

    count = 0
    for t_i, t_sec in enumerate(target_secs):
        for c_i, c_val in enumerate(current_values):
            soma_v_fig.line(r_dict[t_i][c_i]['t'], r_dict[t_i][c_i]['soma_v'], line_width=3, color=colrs[count],
                            legend='{0}, {1} nA'.format(t_sec, c_val))
            count += 1
    return dep_vs_i_fig, soma_v_fig


def print_hoc_distance():
    h.load_file('cell_seed3_0-pyr-08.hoc')

    t_cell2 = h.CA1_PC_cAC_sig6()

    h.distance(0, t_cell2.soma[0](0.5))

    for sec in t_cell2.somatic:
        print(sec.name())
        for seg in sec:
            seg_dist = h.distance(sec(seg.x))
            print('dist: {0:.4}, g_hd: {1:.4}, e_pas: {2:.4}'.format(seg_dist, seg.ghdbar_hd, seg.e_pas))

    for sec in t_cell2.apical:
        print(sec.name())
        for seg in sec:
            seg_dist = h.distance(sec(seg.x))
            print('dist: {0:.4}, g_hd: {1:.4}, e_pas: {2:.4}, g_kad: {3:.4}'.format(seg_dist, seg.ghdbar_hd,
                                                                                    seg.e_pas, seg.gkabar_kad))
    for sec in t_cell2.basal:
        print(sec.name())
        for seg in sec:
            seg_dist = h.distance(sec(seg.x))
            print('dist: {0:.4}, g_hd: {1:.4}, e_pas: {2:.4}, g_kad: {3:.4}'.format(seg_dist, seg.ghdbar_hd,
                                                                                    seg.e_pas, seg.gkabar_kad))


if __name__ == "__main__":
    print(os.getcwd())

    t_path = join(os.getcwd(), 'morphologies/mpg141209_A_idA.asc')
    # t_path = join(os.getcwd(), 'morphologies/mpg141208_B_idA.asc')

    path_parts = os.getcwd().split('/')
    home_i = path_parts.index('ca1_dendritic_experiments')
    home_path = '/'.join(path_parts[:(home_i + 1)])

    config_name = 'config_mpg141209_A_idA'
    # config_name = 'config_mpg141208_B_idA'
    spec_name = config_name + '_spec'
    conf_path = os.path.join(home_path, 'config_files', config_name)
    spec_path = conf_path + '_spec'
    conf_path += '.txt'
    spec_path += '.txt'

    config_obj = conf.load_config(conf_path, spec_path)

    t_cell = MyCell(t_path, config_obj)
    # t_cell.make_synapse_sections(t_cell.apical[22], 0.5, 'CA3toCA1')

    # t_cell.insert_rxd()

    for sec in t_cell.cal_list:
        print(sec.name())

    # t_cell.plot_cell(color_code='one color')
    t_cell.plot_cell2d(pretty=False)

    plt.show()
