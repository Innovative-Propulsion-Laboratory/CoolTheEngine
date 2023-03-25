"""
-*- coding: utf-8 -*-
"""

import tkinter as tk
from tkinter import ttk
import sys
import time


class InputsWin:
    def __init__(self):
        self.list_inputs_features = [
            "Engine model", "Channel global", "Channel dimensions", "Coolant properties"]
        self.list_settings_features = ["Plots", "Other"]

        self.inputs_function_list = [[self.engine_model, self.channel_global,
                                      self.channel_dimension, self.coolant_properties], [self.fplot, self.fother]]

        self.entry_name_dict = {}
        self.entry_dict = {}

    def get_entry(self):
        for entry in self.entry_name_dict.items():
            entry_params = entry[1].get()
            self.entry_dict[entry[0]] = entry_params
        print("Data saved")

    def persistent_entry(self, variable_name):
        var = tk.IntVar()
        try:
            var.set(self.entry_dict[str(variable_name)])
        except:
            var.set("")
        return var

    def engine_model(self, frame_name):
        tk.Label(frame_name, text="Engine model :").pack(anchor="w")
        tk.Label(frame_name, text="Entry test :").pack(anchor="w")

        entry1 = self.persistent_entry("entry1")
        self.test_entry1 = tk.Entry(
            frame_name, textvariable=entry1)
        self.test_entry1.pack(anchor="w")

        entry2 = self.persistent_entry("entry2")
        self.test_entry2 = tk.Entry(
            frame_name, textvariable=entry2)
        self.test_entry2.pack(anchor="w")

        self.entry_name_dict = {
            "entry1": self.test_entry1, "entry2": self.test_entry2}

    def channel_global(self, frame_name):
        tk.Label(frame_name, text="Channel global").pack(anchor="w")

        ttk.Separator(frame_name, orient="horizontal").pack(
            fill="x", padx=10, pady=10)

        tk.Label(frame_name, text="Number of channels :").pack(anchor="w")

        nbc = self.persistent_entry("nbc")
        self.nbc = tk.Entry(frame_name, textvariable=nbc)
        self.nbc.pack(anchor="w")

        tk.Label(frame_name, text="Position of the manifold from the throat (in m) :").pack(
            anchor="w")

        manifold_pos = self.persistent_entry("manifold_pos")
        self.manifold_pos = tk.Entry(frame_name, textvariable=manifold_pos)
        self.manifold_pos.pack(anchor="w")

        self.entry_name_dict = {
            "nbc": self.nbc, "manifold_pos": self.manifold_pos}

    def channel_dimension(self, frame_name):
        tk.Label(frame_name, text="Channel dimension").pack(anchor="w")

        ttk.Separator(frame_name, orient="horizontal").pack(
            fill="x", padx=10, pady=10)

        tk.Label(frame_name, text="Widths").pack(anchor="w")

        tk.Label(frame_name, text="Injection plate").pack(
            anchor="w")

        lrg_inj = self.persistent_entry("lrg_inj")
        self.lrg_inj = tk.Entry(frame_name, textvariable=lrg_inj)
        self.lrg_inj.pack(anchor="w")

        self.entry_name_dict = {
            "lrg_inj": self.lrg_inj}

    def coolant_properties(self, frame_name):
        tk.Label(frame_name, text="Coolant properties :" +
                 "\n"+"n/a").pack(anchor="w")

    def fplot(self, frame_name):
        tk.Label(frame_name, text="Plot settings :" +
                 "\n"+"n/a").pack(anchor="w")

    def fother(self, frame_name):
        tk.Label(frame_name, text="Other settings :" +
                 "\n"+"n/a").pack(anchor="w")


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

        main_menu = tk.Menu(self, tearoff=0)
        menu_file = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="File", menu=menu_file)
        menu_file.add_command(label="Quit", command=self.quit)
        self.config(menu=main_menu)

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.inputs_class = InputsWin()

        self.inputs_function_list = self.inputs_class.inputs_function_list

        try:
            self.runprocess = Run(process_class, self.inputs_class)
        except:
            print("Failed to load main process class in MainGUI")

        self.primary()
        self.info_frame()

        sys.stdout = OutRedirection(self.secondary_r1B_frame)

    def primary(self):
        self.primary_frame = tk.Frame(self)
        self.primary_frame.grid(row=0, column=0, sticky="nesw")

        self.primary_frame.grid_rowconfigure(0, weight=1)
        self.primary_frame.grid_columnconfigure(0, weight=15)
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

    def destroy_secondary_m(self):
        try:
            self.secondary_m_frame.destroy()
            for objects in self.secondary_m_frame.winfo_children():
                objects.destroy()
        except:
            pass
        self.primary_frame.grid_columnconfigure(0, weight=15)
        self.primary_frame.grid_columnconfigure(1, weight=0)
        self.primary_frame.grid_columnconfigure(2, weight=85)

    def create_input_win(self, function_list, set_selected, index):
        self.destroy_secondary_m()
        self.secondary_m()
        function = function_list[set_selected][index]
        function(self.secondary_m_frame)

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
        # self.secondary_r1A_frame.grid(row=0, column=0, sticky="nesw")

        self.secondary_r1A_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1A_frame.grid_columnconfigure(0, weight=1)

        self.secondary_r1_frame.add(self.secondary_r1A_frame, height=100)

        tk.Label(self.secondary_r1A_frame,
                 text="Results plot").grid(row=0, column=0)

    def secondary_r1B(self):
        self.secondary_r1B_frame = tk.Text(
            self.secondary_r1_frame, bd=2, relief="groove", bg="white")
        # self.secondary_r1B_frame.grid(row=1, column=0, sticky="nesw")

        self.secondary_r1B_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1B_frame.grid_columnconfigure(0, weight=1)

        self.secondary_r1_frame.add(self.secondary_r1B_frame, height=500)

    def secondary_m(self):
        self.secondary_m_frame = tk.Frame(
            self.primary_frame, bd=2, relief="groove", width=200)
        self.secondary_m_frame.grid(
            row=0, column=1, sticky="nesw", pady=(22, 0))
        self.secondary_m_frame.grid_propagate(False)
        self.secondary_m_frame.pack_propagate(False)

        self.primary_frame.grid_columnconfigure(0, weight=20)
        self.primary_frame.grid_columnconfigure(1, weight=20)
        self.primary_frame.grid_columnconfigure(2, weight=60)

        self.secondary_m_frame.grid_rowconfigure(0, weight=1)
        self.secondary_m_frame.grid_columnconfigure(0, weight=1)

        tk.Button(self.secondary_m_frame, text="Save data",
                  command=self.inputs_class.get_entry).grid(row=1, sticky="sw")

    def info_frame(self):
        self.info = tk.Frame(self, borderwidth=2, relief="groove")
        self.info.grid(row=1, column=0, sticky="nesw")

        self.info.grid_rowconfigure(0, weight=1)
        self.info.grid_columnconfigure(0, weight=1)
        self.info.grid_columnconfigure(1, weight=1)

        tk.Label(self.info, text="CTE GUI - In development   © IPL").grid(row=0,
                                                                          column=0, sticky="w")

        tk.Button(self.info, text="Run", command=self.run_button_pressed).grid(
            row=0, column=1, sticky="e")

    def run_button_pressed(self):
        try:
            self.runprocess.run_process()
        except:
            print("Failed to run process")


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
