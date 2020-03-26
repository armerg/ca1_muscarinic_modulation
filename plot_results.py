import os
import h5py
import bokeh.plotting as bplt
import bokeh.io as bkio
import bokeh.layouts as blay
import bokeh.models as bmod
import colorcet as cc
import frequency_analysis as fa
import numpy as np


mu = u'\u03BC'
sigma = u'\u03C3'
line_styles = ['', 'dot', 'dash', 'dashdot']
# ply_colors = ['rgb(31, 119, 180, 1.0)',
#               'rgb(44, 160, 44, 1.0)',
#               'rgb(214, 39, 40, 1.0)',
#               'rgb(148, 103, 189, 1.0)',
#               'rgb(140, 86, 75, 1.0)',
#               'rgb(255, 127, 14, 1.0)',
#               'rgb(23, 190, 207, 1.0)',
#               'rgb(227, 119, 194, 1.0)']
full_name = {'CA1P': 'CA1 Pyramidal Cell',
             'BCCCK': 'CCK Basket Cell',
             'BCPV': 'PV Basket Cell',
             'SC': 'CA3 SR',
             'RO': 'CA3 SO'}

# color_dict = {'CA1P': {'line': 'rgba(31, 119, 180, 1.0)', 'fill': 'rgba(31, 119, 180, 0.5)',
#                        'variations': ['rgba(27,79,114, 1.0)',
#                                       'rgba(40,116,166,1.0)',
#                                       'rgba(133,193,233,1.0)',
#                                       'rgba(52,152,219,1.0)',
#                                       'rgba(174,214,241,1.0)',
#                                       ]
#                        },
#               'BCCCK': {'line': 'rgba(214, 39, 40, 1.0)', 'fill': 'rgba(213, 39, 40, 0.5)',
#                         'variations': ['rgba(120,40,31, 1.0)',
#                                        'rgba(176,58,46,1.0)',
#                                        'rgba(231,76,60,1.0)',
#                                        'rgba(241,148,138,1.0)',
#                                        'rgba(250,219,216,1.0)',
#                                        ]
#                         },
#               'BCPV': {'line': 'rgba(44, 160, 44, 1.0)', 'fill': 'rgba(44, 160, 44, 0.5)',
#                        'variations': ['rgba(20,90,50,1.0)',
#                                       'rgba(30,132,73,1.0)',
#                                       'rgba(125,206,160,1.0)',
#                                       'rgba(212,239,223,1.0)',
#                                       ]
#                        }
#               }

color_dict = {'CA1P': {'line': 'blue', 'fill': 'blue',
                       'variations': ['powderblue',
                                      'cyan',
                                      'deepskyblue',
                                      'cornflowerblue',
                                      'slateblue',
                                      'mediumblue',
                                      'midnightblue',
                                      ]
                       },
              'BCCCK': {'line': 'red', 'fill': 'red',
                        'variations': ['lightsalmon',
                                       'orange',
                                       'lightcoral',
                                       'pink',
                                       'crimson',
                                       'firebrick',
                                       'darkred',
                                       ]
                        },
              'BCPV': {'line': 'green', 'fill': 'green',
                       'variations': ['lawngreen',
                                      'lime',
                                      'lightgreen',
                                      'mediumseagreen',
                                      'olive',
                                      'forestgreen',
                                      'darkgreen',
                                      ]
                       },
              'SC': {'line': 'purple', 'fill': 'purple',
                       'variations': ['#ccffff',
                                      '99ffff',
                                      '#66ffff',
                                      '#33ffff',
                                      '#00e6e6',
                                      '#00b3b3',
                                      '#006666',
                                      ]
                       },
              'RO': {'line': 'purple', 'fill': 'purple',
                       'variations': ['#ccffff',
                                      '99ffff',
                                      '#66ffff',
                                      '#33ffff',
                                      '#00e6e6',
                                      '#00b3b3',
                                      '#006666',
                                      ]
                       },

              }


def configure_plot(in_plot, set_size=True):
    """

    :param in_plot:
    :return:
    """
    if set_size:
        in_plot.plot_width = 950
        in_plot.plot_height = 600

    in_plot.xaxis.axis_label_text_font = 'arial'
    in_plot.xaxis.axis_label_text_font_style = 'bold'
    in_plot.xaxis.axis_label_text_font_size = '22pt'

    in_plot.yaxis.axis_label_text_font = 'arial'
    in_plot.yaxis.axis_label_text_font_style = 'bold'
    in_plot.yaxis.axis_label_text_font_size = '22pt'

    in_plot.title.text_font = 'arial'
    in_plot.title.text_font_size = '24pt'
    in_plot.title.text_font_style = 'bold'
    in_plot.legend.label_text_font_size = '18pt'

    in_plot.xaxis.major_label_text_font = 'arial'
    in_plot.xaxis.major_label_text_font_size = "18pt"
    in_plot.xaxis.major_label_text_font_style = 'bold'
    in_plot.yaxis.major_label_text_font = 'arial'
    in_plot.yaxis.major_label_text_font_size = "18pt"
    in_plot.yaxis.major_label_text_font_style = 'bold'

    in_plot.min_border_right = 50

    in_plot.xgrid.grid_line_color = None
    in_plot.ygrid.grid_line_color = None

    in_plot.toolbar.logo = None
    in_plot.toolbar_location = None

    in_plot.background_fill_color = None
    in_plot.border_fill_color = None
    in_plot.outline_line_color = None

    return in_plot