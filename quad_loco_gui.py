#!/usr/bin/env python3

import os
import json
import tkinter as tk
from tkinter import ttk
import subprocess

def args_to_string(args):
    result = ""
    for key, value in args.items():
        if value is True:
            result += f" -{key}"
        elif value is not False:
            result += f" -{key} {value}"
    return result

def update_command_display(cmd):
    command_display.delete(1.0, tk.END)
    command_display.insert(tk.END, cmd)

def run_model(args):
    cmd = f"g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o{args_to_string(args)} >dat 2>ppp"
    update_command_display(cmd)
    os.system(cmd)

def run_gnuplot(plot_file):
    cmd = f"gnuplot -p -e 'load \"{plot_file}\"'"
    update_command_display(cmd)
    subprocess.Popen(cmd, shell=True)

def save_and_run():
    config = {}
    for key, entry in entries.items():
        value = entry.get()
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass
        config[key] = value

    config["ipsi_sw_inhib"] = ipsi_check_var.get()
    config["contra_sw_inhib"] = contra_check_var.get()
    config["diff_kv"] = diff_check_var.get()

    preset_var = preset_combobox.get()
    if preset_var == "Walk":
        config["presetWALK"] = True
    elif preset_var == "Trot":
        config["presetTROT"] = True
    elif preset_var == "Random":
        config["presetRANDOM"] = True

    with open(config_file, "w") as f:
        json.dump(config, f)

    run_model(config)

config_file = "quad_loco_config.json"

if os.path.exists(config_file):
    with open(config_file, "r") as f:
        config = json.load(f)
else:
    config = {}

root = tk.Tk()
root.title("Quadruped Locomotion Model Interface")

# Set the window icon
# root.iconbitmap('model_icon.ico')  # For Windows
root.iconphoto(True, tk.PhotoImage(file='model_icon.png'))  # For other platforms

parameter_labels = [
    ("T", "Locomotion time"),
    ("fps", "FPS"),
    ("v0", "Initial velocity"),
    ("v1", "Final velocity"),
    ("Gu0", "Initial unloading parameter"),
    ("Gu1", "Final unloading parameter"),
    ("kv0", "Initial balance parameter"),
    ("kv1", "Final balance parameter"),
    ("tsw0", "Initial swing duration"),
    ("tsw1", "Final swing duration"),
    ("kr", "Rotation resistance"),
]

label_font = ("Helvetica", 12)
entry_font = ("Helvetica", 12)
# botton_font = ("Helvetica", 12)
title_font = ("Helvetica", 12)

entries = {}
for idx, (arg_key, label) in enumerate(parameter_labels):
    row, col = divmod(idx, 2)  # Divide the index by 2 to create two columns
    frame = ttk.Frame(root)
    frame.grid(row=row, column=col, padx=10, pady=10)

    ttk.Label(frame, text=label, font=label_font).grid(row=0, column=0, padx=5, pady=5)
    entry = ttk.Entry(frame, width=10, font=entry_font)
    entry.grid(row=0, column=1, padx=5, pady=5)
    if arg_key in config:
        entry.insert(0, str(config[arg_key]))
    entries[arg_key] = entry


checkbox_frame = ttk.Frame(root)
checkbox_frame.grid(row=row+1, column=0, columnspan=2, padx=5, pady=5)

ipsi_check_var = tk.BooleanVar()
ipsi_check_var.set(config.get("ipsi_sw_inhib", False))
ipsi_check = ttk.Checkbutton(checkbox_frame, text="Ipsilateral swing inhibition", variable=ipsi_check_var)
ipsi_check.grid(row=0, column=0, padx=5, pady=5)

contra_check_var = tk.BooleanVar()
contra_check_var.set(config.get("contra_sw_inhib", False))
contra_check = ttk.Checkbutton(checkbox_frame, text="Contralateral swing inhibition", variable=contra_check_var)
contra_check.grid(row=0, column=1, padx=5, pady=5)

diff_check_var = tk.BooleanVar()
diff_check_var.set(config.get("diff_kv", False))
diff_check = ttk.Checkbutton(checkbox_frame, text="Diagonal swing excitation", variable=diff_check_var)
diff_check.grid(row=0, column=2, padx=5, pady=5)

preset_frame = ttk.Frame(root)
preset_frame.grid(row=row+2, column=0, columnspan=2, padx=5, pady=5)

preset_label = ttk.Label(preset_frame, text="Preset step cycle:", font=label_font)
preset_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E)

preset_combobox = ttk.Combobox(preset_frame, values=["Walk", "Trot", "Random"])
preset_combobox.set("Walk" if config.get("presetWALK") else "Trot" if config.get("presetTROT") else "Random")
preset_combobox.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W)


button_frame = ttk.Frame(root)
button_frame.grid(row=row+3, column=0, columnspan=2, padx=5, pady=5)

run_button = ttk.Button(button_frame, text="Run Model", command=save_and_run, width=50)
run_button.grid(row=0, column=0, columnspan=3, padx=5, pady=5)

plot_buttons = [
    ("Plot total load", "totload.plot"),
    ("Plot step cycles", "stepcycle.plot"),
    ("Plot stride frequency", "stridefreqv.plot"),
    ("Plot stride length", "stridelengthv.plot"),
    ("Plot duty factors", "dutyfacv.plot"),
    ("Plot swing/stance time", "swstfreq.plot"),
]

for idx, (label, plot_file) in enumerate(plot_buttons):
    row, col = divmod(idx, 3)  # Divide the index by 3 to create two rows
    plot_button = ttk.Button(button_frame, text=label, command=lambda plot_file=plot_file: run_gnuplot(plot_file))
    plot_button.grid(row=row+1, column=col, padx=5, pady=5)


command_display_frame = ttk.Frame(root)
command_display_frame.grid(row=row+8, column=0, columnspan=2, padx=5, pady=5)

command_display_label = ttk.Label(command_display_frame, text="Running Command:", font=title_font)
command_display_label.grid(row=0, column=0, padx=5, pady=5)

command_display = tk.Text(command_display_frame, width=60, height=4, wrap=tk.WORD)
command_display.grid(row=1, column=0, padx=5, pady=5)

root.mainloop()

