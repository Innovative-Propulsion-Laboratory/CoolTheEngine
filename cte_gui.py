"""
-*- coding: utf-8 -*-
"""

import tkinter as tk
from tkinter import ttk
from tkinter import simpledialog
from tkinter import filedialog

from PIL import Image, ImageTk


import sys
import time
import os
import traceback


class InputsWin:
    def __init__(self):
        # -----------------------------------------------------------------------------------------------------------------

        self.list_inputs_features = [
            "Engine model", "Channel global", "Channel dimensions", "Coolant properties", "Advanced"]
        inputs_features_function = [
            self.engine_model, self.channel_global, self.channel_dimension, self.coolant_properties, self.advanced]

        self.list_settings_features = ["Plots", "Other"]
        settings_features_function = [self.fplot, self.fother]

        self.material_list = ["pure copper", "cucrzr", "inconel"]

        self.coolant_list = ["Methane", "Water"]

        self.plot_detail_list = ["0", "1", "2", "3"]

        self.engine_name_list = ["Arrax","Viserion"]

        self.mesh_size_list = ["0.25", "0.5", "1"]

        self.cea_list = ["Viserion_2023.txt","Arrax_2023.txt"]

        self.boolean = ["False", "True"]

        # -----------------------------------------------------------------------------------------------------------------

        self.inputs_function_list = [
            inputs_features_function, settings_features_function]

        self.entry_name_dict = {}
        self.entry_dict = {}

    def engine_model(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Channel global")

        self.create_separator(frame_name)

        self.create_combobox(frame_name, text="Engine",
                             var_name="engine_name", var_list=self.engine_name_list)

        self.create_combobox(frame_name, text="Mesh size",
                             var_name="mesh_size", var_list=self.mesh_size_list)

        self.create_combobox(frame_name, text="CEA parameters",
                             var_name="input_CEA_data", var_list=self.cea_list)

    def channel_global(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Channel global")

        self.create_separator(frame_name)

        self.create_entry(
            frame_name, text="Number", var_name="nbc")
        self.create_entry(
            frame_name, text="Manifold position", var_name="manifold_pos")

        self.create_entry(frame_name, text="Roughness", var_name="roughness")

        self.create_combobox(frame_name, text="Material",
                             var_name="material_name", var_list=self.material_list)

    def channel_dimension(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Channel dimension")

        self.create_separator(frame_name)

        self.create_label(frame_name, text="At the injection plate :")
        self.create_entry(frame_name, text="Width",
                          var_name="lrg_inj")
        self.create_entry(frame_name, text="Height",
                          var_name="ht_inj")

        self.create_separator(frame_name)

        self.create_label(frame_name, text="At the end of the chamber :")
        self.create_entry(frame_name, text="Width",
                          var_name="lrg_conv")
        self.create_entry(frame_name, text="Height",
                          var_name="ht_conv")
        self.create_entry(frame_name, text="Thickness",
                          var_name="e_conv")

        self.create_separator(frame_name)

        self.create_label(frame_name, text="In the throat :")
        self.create_entry(frame_name, text="Width",
                          var_name="lrg_col")
        self.create_entry(frame_name, text="Height",
                          var_name="ht_col")
        self.create_entry(frame_name, text="Thickness",
                          var_name="e_col")

        self.create_separator(frame_name)

        self.create_label(frame_name, text="At the manifold :")
        self.create_entry(frame_name, text="Width",
                          var_name="lrg_tore")
        self.create_entry(frame_name, text="Height",
                          var_name="ht_tore")
        self.create_entry(frame_name, text="Thickness",
                          var_name="e_tore")

    def coolant_properties(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Coolant properties")

        self.create_separator(frame_name)

        self.create_combobox(frame_name, text="Fluid",
                             var_name="fluid", var_list=self.coolant_list)
        self.create_entry(frame_name, text="Density",
                          var_name="density_cool_init")

        self.create_separator(frame_name)

        self.create_entry(frame_name, text="Initial temperature",
                          var_name="Temp_cool_init")
        self.create_entry(frame_name, "Initial pressure", "Pressure_cool_init")

    def advanced(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Geometry coefficients")

        self.create_separator(frame_name)

        self.create_entry(frame_name, text="Width convergent", var_name="n1")
        self.create_entry(frame_name, text="Width divergent", var_name="n2")

        self.create_separator(frame_name)

        self.create_entry(frame_name, text="Height convergent", var_name="n3")
        self.create_entry(frame_name, text="Height divergent", var_name="n4")

        self.create_separator(frame_name)

        self.create_entry(
            frame_name, text="Thickness convergent", var_name="n5")
        self.create_entry(
            frame_name, text="Thickness divergent", var_name="n6")

    def fplot(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Plot settings")

        self.create_separator(frame_name)

        self.create_combobox(frame_name, text="Plot detail",
                             var_name="plot_detail", var_list=self.plot_detail_list)

        self.create_combobox(frame_name, text="3D plots",
                             var_name="show_3d_plots", var_list=self.boolean)

        self.create_combobox(frame_name, text="2D temperature",
                             var_name="show_2D_temperature", var_list=self.boolean)

        self.create_combobox(frame_name, text="Final 3D plot",
                             var_name="do_final_3d_plot", var_list=self.boolean)

    def fother(self, frame_name):
        self.init_feature()
        self.create_label(frame_name, text="Other settings")

        self.create_separator(frame_name)

        self.create_combobox(frame_name, text="Write in CSV",
                             var_name="write_in_csv", var_list=self.boolean)

    def get_entry(self):
        try:
            for entry in self.entry_name_dict.items():
                entry_params = entry[1].get()
                if entry_params == "False":
                    entry_params = False
                self.entry_dict[entry[0]] = entry_params
            print(
                "█ Changes has been saved                                                   █")

        except:
            print("Impossible to get entry")
        # print(self.entry_dict)

    def init_feature(self):
        self.line = 0
        self.entry_name_dict = {}

    def save_params(self):
        save_name = simpledialog.askstring(
            "Input", "Save name :")
        if save_name != None and save_name != "":
            file = open("save/" + str(save_name) + ".scte", "w")
            for key, value in self.entry_dict.items():
                file.write(str(key) + "$$" + str(value) + "\n")
            file.close()
            print(
                "█ Settings has been saved                                                  █")

    def open_params(self):
        filename = filedialog.askopenfilename(
            initialdir="save", title="Select File", filetypes=(("scte files", "*.scte"), ("all files", "*.*")))
        file = open(str(filename), "r")
        for line in file.readlines():
            line = line.split("$$")
            self.entry_dict[line[0]] = line[1][:-1]
            if self.entry_dict[line[0]] == "False":
                self.entry_dict[line[0]] = False
        file.close()
        print(
            "█ Data has been imported                                                   █")

    def create_label(self, frame_name, text):
        tk.Label(frame_name, text=text).grid(
            row=self.line, columnspan=2, sticky="nw")
        self.line += 1

    def create_entry(self, frame_name, text, var_name):
        var = tk.IntVar()
        try:
            var.set(self.entry_dict[var_name])
        except:
            var.set("")
        tk.Label(frame_name, text=text).grid(
            row=self.line, column=0, sticky="nw")
        var_entry = tk.Entry(frame_name, textvariable=var, width=10)
        var_entry.grid(row=self.line, column=1, sticky="nw")
        self.entry_name_dict[var_name] = var_entry
        self.line += 1

    def create_large_entry(self, frame_name, var_name):
        var = tk.IntVar()
        try:
            var.set(self.entry_dict[var_name])
        except:
            var.set("")
        var_entry = tk.Entry(frame_name, textvariable=var)
        var_entry.grid(row=self.line, columnspan=2, sticky="nw")
        self.entry_name_dict[var_name] = var_entry
        self.line += 1

    def create_combobox(self, frame_name, text, var_name, var_list):
        tk.Label(frame_name, text=text).grid(
            row=self.line, column=0, sticky="nw")

        """def selected(event):
            select = combo_list.get()
            self.entry_dict[var_name] = select"""
        var_combo = ttk.Combobox(
            frame_name, values=var_list, state="readonly", width=10)
        try:
            index = var_list.index(self.entry_dict[var_name])
        except:
            index = 0
        var_combo.current(index)
        var_combo.grid(row=self.line, column=1, sticky="nw")
        self.entry_name_dict[var_name] = var_combo
        self.line += 1

    def create_separator(self, frame_name):
        ttk.Separator(frame_name, orient="horizontal").grid(
            row=self.line, columnspan=2, padx=10, pady=10, sticky="nesw")
        self.line += 1


class Output:
    def __init__(self, plot_dir="plot_cache"):
        self.plot_dir = plot_dir

        self.list_output_features = os.listdir(self.plot_dir)

    def list_plot(self, frame_name, index, display_frame, display_label):
        self.display_frame = display_frame
        self.display_label = display_label
        self.listbox = tk.Listbox(frame_name)
        self.listbox.grid(row=0, column=0, sticky="nesw")

        self.target_dir = os.path.join(
            self.plot_dir, self.list_output_features[index])

        self.list_plot_file = os.listdir(self.target_dir)
        self.list_plot_name = [name[:-4] for name in self.list_plot_file]

        for i in range(len(self.list_plot_name)):
            self.listbox.insert(i, self.list_plot_name[i])

        self.listbox.bind("<<ListboxSelect>>", self.selection)

    def selection(self, event):
        selected = event.widget.curselection()

        if selected:
            index = selected[0]
            file_path = os.path.join(
                self.target_dir, self.list_plot_file[index])
            self.display_plot(file_path, self.display_frame,
                              self.display_label)

    def display_plot(self, file_path, display_frame, display_label):
        def resize(event):
            width = event.width
            height = event.height

            aspect_ratio = plot.width / plot.height

            new_width = width
            new_height = int(new_width / aspect_ratio)

            if new_height > height:
                new_height = height
                new_width = int(new_height * aspect_ratio)

            resized_plot = plot.resize(
                (new_width, new_height))
            resized_plot_tk = ImageTk.PhotoImage(resized_plot)

            display_label.config(image=resized_plot_tk)
            display_label.image = resized_plot_tk

        def initial_resize():
            width = display_frame.winfo_width()
            height = display_frame.winfo_height()

            aspect_ratio = plot.width / plot.height

            new_width = width
            new_height = int(new_width / aspect_ratio)

            if new_height > height:
                new_height = height
                new_width = int(new_height * aspect_ratio)

            resized_plot = plot.resize(
                (new_width, new_height))
            resized_plot_tk = ImageTk.PhotoImage(resized_plot)

            display_label.config(image=resized_plot_tk)
            display_label.image = resized_plot_tk

        plot = Image.open(file_path)

        initial_resize()

        display_frame.bind("<Configure>", resize)


class OutRedirection:
    def __init__(self, out_frame):
        self.out_frame = out_frame

    def write(self, str):
        self.out_frame.insert("end", str)
        self.out_frame.see("end")
        self.out_frame.update()

    def flush(self):
        pass


class Run:
    def __init__(self, process_class, inputs_class):
        self.process_class = process_class()
        self.entry_dict = inputs_class.entry_dict

    def run_process(self):
        self.process_class.process(self.entry_dict)


class MainGUI(tk.Tk):
    def __init__(self, process_class):
        tk.Tk.__init__(self)

        self.index_temp = None
        self.set_selected = None
        self.entry_params = {}

        """
        w = self.winfo_screenwidth()
        h = self.winfo_screenheight()
        self.geometry("%ix%i" % (w/2, h/2))
        """
        self.geometry("1200x600")

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.inputs_class = InputsWin()
        self.output_class = Output()

        self.inputs_function_list = self.inputs_class.inputs_function_list

        self.menu_bar()

        try:
            self.runprocess = Run(process_class, self.inputs_class)
        except:
            print("Failed to load main process class in MainGUI")

        self.primary()
        self.info_frame()

        sys.stdout = OutRedirection(self.secondary_r1B_frame)

    def menu_bar(self):
        main_menu = tk.Menu(self, tearoff=0)
        menu_file = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="File", menu=menu_file)
        menu_file.add_command(label="Quit", command=self.quit)
        menu_file.add_separator()
        menu_file.add_command(
            label="Save", command=self.inputs_class.save_params)
        menu_file.add_command(
            label="Open", command=self.inputs_class.open_params)
        self.config(menu=main_menu)

    def primary(self):
        self.primary_frame = tk.Frame(self)
        self.primary_frame.grid(row=0, column=0, sticky="nesw")

        self.primary_frame.grid_rowconfigure(0, weight=1)
        self.primary_frame.grid_columnconfigure(0, weight=15)
        self.primary_frame.grid_columnconfigure(1, weight=0)
        self.primary_frame.grid_columnconfigure(2, weight=85)

        self.secondary_l()
        self.secondary_r()

    def secondary_l(self):
        self.secondary_l_frame = tk.Frame(
            self.primary_frame, bd=1, relief="flat")
        self.secondary_l_frame.grid(row=0, column=0, sticky="nesw")

        self.secondary_l_frame.grid_rowconfigure(1, weight=1)
        self.secondary_l_frame.grid_rowconfigure(3, weight=1)
        self.secondary_l_frame.grid_rowconfigure(5, weight=1)
        self.secondary_l_frame.grid_columnconfigure(0, weight=1)

        tk.Label(self.secondary_l_frame, text="Inputs").grid(row=0)
        self.tertiary_l1()
        tk.Label(self.secondary_l_frame, text="Settings").grid(row=2)
        self.tertiary_l2()
        tk.Label(self.secondary_l_frame, text="Outputs").grid(row=4)
        self.tertiary_l3()

    def tertiary_l1(self):
        self.inputs = tk.Listbox(self.secondary_l_frame)
        self.inputs.grid(row=1, sticky="nesw")

        list_inputs_features = self.inputs_class.list_inputs_features
        for i in range(len(list_inputs_features)):
            self.inputs.insert(i, list_inputs_features[i])

        self.inputs.bind("<<ListboxSelect>>", self.l1_selected)

    def tertiary_l2(self):
        self.settings = tk.Listbox(self.secondary_l_frame)
        self.settings.grid(row=3, sticky="nesw")

        list_settings_features = self.inputs_class.list_settings_features
        for i in range(len(list_settings_features)):
            self.settings.insert(i, list_settings_features[i])

        self.settings.bind("<<ListboxSelect>>", self.l2_selected)

    def tertiary_l3(self):
        self.outputs = tk.Listbox(self.secondary_l_frame)
        self.outputs.grid(row=5, sticky="nesw")

        list_output_features = self.output_class.list_output_features
        for i in range(len(list_output_features)):
            self.outputs.insert(i, list_output_features[i])

        self.outputs.bind("<<ListboxSelect>>", self.l3_selected)

    def l1_selected(self, event):
        selected = event.widget.curselection()
        if selected:
            self.set_selected = 0
            self.listbox_selection(event)

    def l2_selected(self, event):
        selected = event.widget.curselection()
        if selected:
            self.set_selected = 1
            self.listbox_selection(event)

    def l3_selected(self, event):
        selected = event.widget.curselection()
        if selected:
            self.set_selected = 2
            self.output_selection(event)

    def listbox_selection(self, event):
        selection = event.widget.curselection()
        index = selection[0]

        if selection and self.index_temp != index:
            self.destroy_secondary_m()
            self.index_temp = index

            self.create_input_win(self.inputs_function_list,
                                  self.set_selected, index)
        else:
            self.destroy_secondary_m()
            self.inputs.selection_clear(0, "end")
            self.settings.selection_clear(0, "end")
            self.index_temp = None

    def output_selection(self, event):
        selection = event.widget.curselection()
        index = selection[0]

        if selection and self.index_temp != index:
            self.destroy_secondary_m()
            self.index_temp = index

            self.create_output_win(index=index)
        else:
            self.destroy_secondary_m()
            self.inputs.selection_clear(0, "end")
            self.settings.selection_clear(0, "end")
            self.index_temp = None

    def destroy_secondary_m(self):
        try:
            for objects in self.secondary_m_frame.winfo_children():
                objects.destroy()
            self.secondary_m_frame.destroy()

        except:
            pass
        self.primary_frame.grid_columnconfigure(0, weight=15)
        self.primary_frame.grid_columnconfigure(1, weight=0)
        self.primary_frame.grid_columnconfigure(2, weight=85)

    def create_input_win(self, function_list, set_selected, index):
        self.destroy_secondary_m()
        self.secondary_m()
        function = function_list[set_selected][index]
        function(self.subframe_m)

    def create_output_win(self, index):
        self.destroy_secondary_m()
        self.secondary_m(button=False)

        self.output_class.list_plot(
            self.secondary_m_frame, index=index, display_frame=self.secondary_r1A_frame, display_label=self.display_label)

    def secondary_r(self):
        self.secondary_r_frame = tk.Frame(
            self.primary_frame, bd=1, relief="flat")
        self.secondary_r_frame.grid(row=0, column=2, sticky="nesw")

        self.secondary_r_frame.grid_rowconfigure(1, weight=1)
        # self.secondary_r_frame.grid_rowconfigure(2, weight=1)
        self.secondary_r_frame.grid_columnconfigure(0, weight=1)

        tk.Label(self.secondary_r_frame, text="Results").grid(row=0, column=0)
        self.secondary_r1()

    def secondary_r1(self):
        self.secondary_r1_frame = tk.PanedWindow(
            self.secondary_r_frame, bd=0, relief="flat", orient="vertical")
        self.secondary_r1_frame.grid(row=1, column=0, sticky="nesw")

        self.secondary_r1_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1_frame.grid_rowconfigure(1, weight=1)

        self.secondary_r1A()
        self.secondary_r1B()

        # self.secondary_r1_frame.grid_rowconfigure(0, weight=8)
        # self.secondary_r1_frame.grid_rowconfigure(1, weight=2)
        # self.secondary_r1_frame.grid_columnconfigure(0, weight=1)

    def secondary_r1A(self):
        self.secondary_r1A_frame = tk.Frame(
            self.secondary_r1_frame, bd=2, relief="groove")
        self.secondary_r1A_frame.grid(row=0, column=0, sticky="nesw")

        self.secondary_r1A_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1A_frame.grid_columnconfigure(0, weight=1)

        self.secondary_r1_frame.add(self.secondary_r1A_frame, height=400)

        self.display_label = tk.Label(
            self.secondary_r1A_frame, background="white")
        self.display_label.grid(row=0, column=0, sticky="nesw")

        self.output_class.display_plot(
            "logo_ipl.png", self.secondary_r1A_frame, self.display_label)

    def secondary_r1B(self):
        self.secondary_r1B_frame = tk.Text(
            self.secondary_r1_frame, bd=2, relief="groove", bg="white")
        # self.secondary_r1B_frame.grid(row=1, column=0, sticky="nesw")

        self.secondary_r1B_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1B_frame.grid_columnconfigure(0, weight=1)

        self.secondary_r1_frame.add(self.secondary_r1B_frame, height=100)

    def secondary_m(self, button=True):
        if button:
            self.secondary_m_frame = tk.Frame(
                self.primary_frame, bd=2, relief="groove", width=200)
        else:
            self.secondary_m_frame = tk.Frame(
                self.primary_frame, bd=2, relief="flat", width=200)
        self.secondary_m_frame.grid(
            row=0, column=1, sticky="nesw", pady=(22, 0))

        self.secondary_m_frame.grid_propagate(False)
        self.secondary_m_frame.pack_propagate(False)

        self.primary_frame.grid_columnconfigure(0, weight=20)
        self.primary_frame.grid_columnconfigure(1, weight=20)
        self.primary_frame.grid_columnconfigure(2, weight=60)

        self.secondary_m_frame.grid_rowconfigure(0, weight=1)
        self.secondary_m_frame.grid_columnconfigure(0, weight=1)

        self.subframe_m = tk.Frame(self.secondary_m_frame)
        self.subframe_m.grid(row=0, column=0, sticky="nesw")

        if button:
            tk.Button(self.secondary_m_frame, text="Save changes",
                      command=self.inputs_class.get_entry).grid(row=1, sticky="sw")

    def info_frame(self):
        self.info = tk.Frame(self, borderwidth=2, relief="groove")
        self.info.grid(row=1, column=0, sticky="nesw")

        self.info.grid_rowconfigure(0, weight=1)
        self.info.grid_columnconfigure(0, weight=1)
        self.info.grid_columnconfigure(1, weight=1)

        tk.Label(self.info, text="CTE - V2.0.0   © IPL").grid(row=0,
                                                              column=0, sticky="w")

        tk.Button(self.info, text="Run", command=self.run_button_pressed).grid(
            row=0, column=1, sticky="e")

    def run_button_pressed(self):
        try:
            self.runprocess.run_process()
            self.output_class = Output()
            self.tertiary_l3()

        except Exception as e:
            print("Failed to run process : ")
            print(e)
            """error = traceback.format_exc()
            print("##################\n\n", error, "\n\n##################")"""
            self.output_class = Output()
        # self.runprocess.run_process()


if __name__ == "__main__":
    try:
        runprocess = Run(None)
    except:
        print("Failed to load process class in Run")
    gui = MainGUI(None)
    gui.title("CTE")

    while True:
        try:
            gui.update()
            gui.mainloop()
            exit()
        except tk.TclError:
            break
