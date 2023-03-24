#!/usr/bin/env python3

import os
import json
import tkinter as tk
from tkinter import ttk

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
    cmd = f"g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o {args_to_string(args)} >dat 2>ppp"
    update_command_display(cmd)
    os.system(cmd)

def run_gnuplot(plot_file):
    cmd = f"gnuplot load \"{plot_file}\""
    update_command_display(cmd)
    os.system(cmd)

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

    inhib_var = inhib_check_var.get()
    config["inhib"] = inhib_var

    preset_var = preset_combobox.get()
    if preset_var == "Walk":
        config["presetWALK"] = True
    elif preset_var == "Trot":
        config["presetTROT"] = True

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

entries = {}
for idx, (arg_key, label) in enumerate(parameter_labels):
    row, col = divmod(idx, 2)
    frame = ttk.Frame(root)
    frame.grid(row=row, column=col, padx=5, pady=5)

    ttk.Label(frame, text=label).grid(row=0, column=0, padx=5, pady=5)
    entry = ttk.Entry(frame, width=10)
    entry.grid(row=0, column=1, padx=5, pady=5)
    if arg_key in config:
        entry.insert(0, str(config[arg_key]))
    entries[arg_key] = entry

inhib_check_var = tk.BooleanVar()
inhib_check_var.set(config.get("inhib", False))
inhib_check = ttk.Checkbutton(root, text="Inhib", variable=inhib_check_var)
inhib_check.grid(row=row+1, column=0, padx=5, pady=5)

preset_combobox = ttk.Combobox(root, values=["Walk", "Trot", "Random"])
preset_combobox.set("Walk" if config.get("presetWALK") else "Trot" if config.get("presetTROT") else "Random")
preset_combobox.grid(row=row+1, column=1, padx=5, pady=5)

button_frame = ttk.Frame(root)
button_frame.grid(row=row+2, column=0, columnspan=2, padx=5, pady=5)

command_display_frame = ttk.Frame(root)
command_display_frame.grid(row=row+3, column=0, columnspan=2, padx=5, pady=5)

command_display_label = ttk.Label(command_display_frame, text="Generated Command:")
command_display_label.grid(row=0, column=0, padx=5, pady=5)

command_display = tk.Text(command_display_frame, width=60, height=2, wrap=tk.WORD)
command_display.grid(row=1, column=0, padx=5, pady=5)

run_button = ttk.Button(button_frame, text="Run Model", command=save_and_run)
run_button.grid(row=0, column=0, padx=5, pady=5)

plot_total_load_button = ttk.Button(button_frame, text="Plot total load", command=lambda: run_gnuplot("totload.plot"))
plot_total_load_button.grid(row=0, column=1, padx=5, pady=5)

plot_step_cycles_button = ttk.Button(button_frame, text="Plot step cycles", command=lambda: run_gnuplot("stepcycle.plot"))
plot_step_cycles_button.grid(row=0, column=2, padx=5, pady=5)

plot_duty_factors_button = ttk.Button(button_frame, text="Plot duty factors", command=lambda: run_gnuplot("dutyfacv.plot"))
plot_duty_factors_button.grid(row=0, column=3, padx=5, pady=5)


root.mainloop()

