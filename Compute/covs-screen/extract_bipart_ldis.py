import os
import sys
sys.path.append("/home/ma629/git")

from rouls.struct_utils import read_ct_file, write_ct_file, plot_arc


fse_structs_dir = "long_fse_structs"
fse_ldi_structs_dir = "ldi_fse_structs"

accs = ["NC_045512.2", "NC_004718.3", "NC_014470.1"]

bipart_ldi_pairs = {(188, 1306), (196, 1302)}


def extract_ldis():
    if not os.path.isdir(fse_ldi_structs_dir):
        os.mkdir(fse_ldi_structs_dir)
    for acc in accs:
        ct_file = os.path.join(fse_structs_dir, f"{acc}.ct")
        pairs, paired, seq = read_ct_file(ct_file, title_mode="number")
        n_ldis = 0
        for name, i_pairs in pairs.items():
            if i_pairs & bipart_ldi_pairs == bipart_ldi_pairs:
                n_ldis += 1
                struct_name = f"{acc}_{n_ldis}"
                ct_out = os.path.join(fse_ldi_structs_dir, f"{struct_name}.ct")
                write_ct_file(ct_out, seq, {struct_name: i_pairs}, overwrite=True)
                arc_plot_file = os.path.join(fse_ldi_structs_dir, f"{struct_name}_arc.png")
                plot_arc(arc_plot_file, seq, i_pairs, title=struct_name)

extract_ldis()

