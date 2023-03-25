"""
-*- coding: utf-8 -*-
"""

import tkinter as tk
from tkinter import ttk
import sys
import time


class InputsWin:
    def __init__(self):
        # -----------------------------------------------------------------------------------------------------------------

        self.list_inputs_features = [
            "Engine model", "Channel global", "Channel dimensions", "Coolant properties"]
        inputs_features_function = [
            self.engine_model, self.channel_global, self.channel_dimension, self.coolant_properties]

        self.list_settings_features = ["Plots", "Other"]
        settings_features_function = [self.fplot, self.fother]

        # -----------------------------------------------------------------------------------------------------------------

        self.inputs_function_list = [
            inputs_features_function, settings_features_function]

        self.entry_name_dict = {}
        self.entry_dict = {}

    def engine_model(self, frame_name):
        self.init_features()
        self.create_label(frame_name, "Test")
        self.create_separator(frame_name)
        self.create_entry(frame_name, text="entry test", var_name="test1")
        self.create_entry(frame_name, text="entry test2", var_name="test2")

    def channel_global(self, frame_name):
        self.init_features()
        self.create_label(frame_name, text="Channel global")
        self.create_separator(frame_name)
        self.create_entry(
            frame_name, text="Number of channels :", var_name="nbc")
        self.create_entry(
            frame_name, text="Position of the manifold", var_name="manifold_pos")

    def channel_dimension(self, frame_name):
        self.init_features()
        self.create_label(frame_name, text="Channel dimension")
        self.create_separator(frame_name)
        self.create_label(frame_name, text="Widths")
        self.create_entry(frame_name, text="Injection plate",
                          var_name="lrg_inj")

    def coolant_properties(self, frame_name):
        tk.Label(frame_name, text="Coolant properties :" +
                 "\n"+"n/a").pack(anchor="w")

    def fplot(self, frame_name):
        tk.Label(frame_name, text="Plot settings :" +
                 "\n"+"n/a").pack(anchor="w")

    def fother(self, frame_name):
        tk.Label(frame_name, text="Other settings :" +
                 "\n"+"n/a").pack(anchor="w")

    def get_entry(self):
        try:
            for entry in self.entry_name_dict.items():
                entry_params = entry[1].get()
                self.entry_dict[entry[0]] = entry_params
            print("Data saved")
        except:
            print("Impossible to get entry")

    def init_features(self):
        self.line = 0
        self.entry_name_dict = {}

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
        var_entry = tk.Entry(frame_name, textvariable=var)
        var_entry.grid(row=self.line, column=1, sticky="nw")
        self.entry_name_dict[var_name] = var_entry
        self.line += 1

    def create_separator(self, frame_name):
        ttk.Separator(frame_name, orient="horizontal").grid(
            row=self.line, columnspan=2, padx=10, pady=10, sticky="nesw")
        self.line += 1


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
        function(self.subframe_m)

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

        self.subframe_m = tk.Frame(self.secondary_m_frame)
        self.subframe_m.grid(row=0, column=0, sticky="nesw")

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
