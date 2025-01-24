import pandas as pd                 # for manipulating data
import pcdl                         # physicell data loader library
import math, os, sys, re
import os.path
import glob
import multiprocessing
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
# from matplotlib.colors import LogNorm
import xml.etree.ElementTree as ET  # for accessing PhysiCell_settings.xml
from mpl_toolkits.axes_grid1 import make_axes_locatable


def make_still (mcds: object, intervention_name: str, output_path: str):

    fig, ax = plt.subplots(figsize = (8, 5.5))

    cell_df = mcds.get_cell_df()

    number_of_live_cells = len(cell_df[cell_df['dead'] == False])
    time = mcds.get_time()
    time = time/60 # (convert to hours)

    # iN theory these oculd be set once - don't worr yabout it. 

    plot_extend = mcds.get_mesh_mnp_range() # returns in a list of tuples the lowest and highest x-axis, y-axis, and z-axis mesh center value.

    # This assumes dx = dy = 20
    neg_plot_x_extend = plot_extend[0][0] - 10
    pos_plot_x_extend = plot_extend[0][1] + 10
    neg_plot_y_extend = plot_extend[1][0] - 10
    pos_plot_y_extend = plot_extend[1][1] + 10

    # substrate
    conc_df = mcds.get_conc_df() # also brilliant - first substrate willl always be ... 7th column (0-indexed)
    substrate = conc_df.columns[6] # gives me the substrate name

    substrate_conc_grid = mcds.get_concentration(substrate, 0.0)
    substrate_im = ax.imshow(substrate_conc_grid, cmap='coolwarm', interpolation='bicubic', origin='lower', vmin=0, vmax=1, 
            extent=[neg_plot_x_extend, pos_plot_x_extend, neg_plot_y_extend, pos_plot_y_extend])

    # Cells
    for cell in cell_df.itertuples(): 
        circ_nucleus = Circle((cell.position_x, cell.position_y), radius=cell.nuclear_radius, color='purple', alpha=0.5)
        circ_cytoplasm = Circle((cell.position_x, cell.position_y), radius=cell.radius, color='purple', alpha=0.1, edgecolor=None)
        ax.add_artist(circ_cytoplasm)
        ax.add_artist(circ_nucleus)
        # ax.scatter(cell.position_x, cell.position_y, s=cell.radius, color='blue', alpha=0.7)
        ax.set_aspect('equal')


    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.10)
    cb = fig.colorbar(substrate_im, cax=cax, format='%.3f')
    cb.set_label('Concentration [AU]')
    ax.set_xlabel('x position [$\mu$m]')
    ax.set_ylabel('y position [$\mu$m]')
    ax.set_title(intervention_name + '\n' + str(time) + ' hours: ' + str(number_of_live_cells) + ' live cells', fontsize=20)
    ax.set_xlim(neg_plot_x_extend, pos_plot_x_extend)
    ax.set_ylim(neg_plot_y_extend, pos_plot_y_extend)
    plt.tight_layout()

    name = int(mcds.get_time())
    name = f"{name:0>8}"
    name = "still_" + str(name)
    name = os.path.join(output_path, name)
    plt.savefig(name + '.png', dpi=300)
    plt.close()

def make_movie (mcds_ts: object, input_path: str, intervention_name: str, output_path: str):
    mcds_ts.make_movie(input_path, interface='png')
    # mv movie to output directory and rename
    os.system('mv output_png12.mp4' + ' ' + output_path + intervention_name + '.mp4')
    # os.system(f"ffmpeg -y -r 12 -i {input_path}/%08d.png -vcodec libx264 -pix_fmt yuv420p -vf scale=1920:-2 {output_path}/movie.mp4") # this is probably close but fix later
        #   os.system(
        #     'ffmpeg -start_number ' + str(
        #         start_file_index) + ' -y -framerate 24 -i ' + save_path + 'output%08d.png' + ' -frames:v ' + str(
        #         number_frames) + ' -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')


def main(input_dir, output_dir, intervention_name, plot_title_name):
    # Assuming the below directory adjustments and path listings are similar to your initial setup
    rel_input_dir = input_dir
    # rel_output_dir = output_dir

    full_path = os.getcwd()
    print("Current Working Directory:", full_path)

    input_dir = os.path.join(full_path, rel_input_dir)
    print("Input Directory:", input_dir)

    # output_dir = os.path.join(full_path, rel_output_dir)
    # print("Output Directory:", output_dir)

    # input("Press Enter to continue...\n Press Ctrl + C to exit...")

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load the data
    mcds_ts = pcdl.TimeSeries(input_dir, microenv=True, settingxml="PhysiCell_settings.xml", graph=False, verbose=False)
    
    # Get the data
    mcds_list = mcds_ts.get_mcds_list()

    args_list = [(mcds, input_dir, plot_title_name) for mcds in mcds_list]

    # Use multiprocessing to process each intervention in parallel
    multiprocessing.Queue(1000)
    with multiprocessing.Pool() as pool:
        pool.map(make_still, args_list)

    make_movie(mcds_ts, input_dir, intervention_name, output_dir)

if __name__ == '__main__':
    input_dir = sys.argv[1] # where files are coming from (its also whre the images go back to)
    output_dir = sys.argv[2] # where the movie will be saved
    intervention_name = sys.argv[3] # movie name
    plot_title_name = sys.argv[4] # BE SURE TO escape any slashes!!! 'GAP$\\rightarrow$RAS$=$1'
    main(input_dir, output_dir, intervention_name, plot_title_name)