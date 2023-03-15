"""
-*- coding: utf-8 -*-
"""

import tkinter as tk
import sys


class Inputs_win(tk.Tk):
    def __init__(self):
        self.list_inputs_features = [
            "Engine model", "Channel global", "Channel dimensions", "Coolant properties"]
        self.list_settings_features = ["Plots", "Other"]

    def engine_model(self):
        tk.Label(self, text="Engine model :"+"\n" +
                 "n/a").grid(row=0, sticky="nw")

    def channel_global(self):
        tk.Label(self, text="Channel global :"+"\n" +
                 "n/a").grid(row=0, sticky="nw")

    def channel_dimension(self):
        tk.Label(self, text="Channel dimension :" +
                 "\n"+"n/a").grid(row=0, sticky="nw")

    def coolant_properties(self):
        tk.Label(self, text="Coolant properties :" +
                 "\n"+"n/a").grid(row=0, sticky="nw")

    def fplot(self):
        tk.Label(self, text="Plot settings :" +
                 "\n"+"n/a").grid(row=0, sticky="nw")

    def fother(self):
        tk.Label(self, text="Other settings :" +
                 "\n"+"n/a").grid(row=0, sticky="nw")


class Out_redirection:
    def __init__(self, out_frame):
        self.out_frame = out_frame

    def write(self, str):
        self.out_frame.insert("end", str)
        self.out_frame.see("end")


class Main_GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)

        self.index_temp = None
        self.set_selected = None

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

        self.inputs_class = Inputs_win()

        self.primary()
        self.info_frame()

        sys.stdout = Out_redirection(self.secondary_r1B_frame)

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

        def destroy_secondary_m():
            try:
                self.secondary_m_frame.destroy()
                for objects in self.secondary_m_frame.winfo_children():
                    objects.destroy()
            except:
                pass
            self.primary_frame.grid_columnconfigure(0, weight=15)
            self.primary_frame.grid_columnconfigure(1, weight=0)
            self.primary_frame.grid_columnconfigure(2, weight=85)

        if selection and self.index_temp != index:
            self.index_temp = index
            print(index)
            print(self.set_selected)

            if self.set_selected == 0 and index == 0:
                destroy_secondary_m()
                self.secondary_m()
                Inputs_win.engine_model(self.secondary_m_frame)
            elif self.set_selected == 0 and index == 1:
                destroy_secondary_m()
                self.secondary_m()
                Inputs_win.channel_global(self.secondary_m_frame)
            elif self.set_selected == 0 and index == 2:
                destroy_secondary_m()
                self.secondary_m()
                Inputs_win.channel_dimension(self.secondary_m_frame)
            elif self.set_selected == 0 and index == 3:
                destroy_secondary_m()
                self.secondary_m()
                Inputs_win.coolant_properties(self.secondary_m_frame)
            elif self.set_selected == 1 and index == 0:
                destroy_secondary_m()
                self.secondary_m()
                Inputs_win.fplot(self.secondary_m_frame)
            elif self.set_selected == 1 and index == 1:
                destroy_secondary_m()
                self.secondary_m()
                Inputs_win.fother(self.secondary_m_frame)
            else:
                destroy_secondary_m()

        else:
            self.inputs.selection_clear(0, "end")
            self.settings.selection_clear(0, "end")
            destroy_secondary_m()
            self.index_temp = None

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

        self.secondary_r1_frame.add(self.secondary_r1A_frame)

        tk.Label(self.secondary_r1A_frame,
                 text="Results plot").grid(row=0, column=0)

    def secondary_r1B(self):
        self.secondary_r1B_frame = tk.Text(
            self.secondary_r1_frame, bd=2, relief="groove", bg="white")
        # self.secondary_r1B_frame.grid(row=1, column=0, sticky="nesw")

        self.secondary_r1B_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1B_frame.grid_columnconfigure(0, weight=1)

        self.secondary_r1_frame.add(self.secondary_r1B_frame)

        """
        tk.Label(self.secondary_r1B_frame,
                 text="Results text").grid(row=0, column=0)
        """

    def secondary_m(self):
        self.secondary_m_frame = tk.Frame(
            self.primary_frame, bd=2, relief="groove", width=200)
        self.secondary_m_frame.grid(
            row=0, column=1, sticky="nesw", pady=(23, 0))
        self.secondary_m_frame.grid_propagate(False)

        self.primary_frame.grid_columnconfigure(0, weight=20)
        self.primary_frame.grid_columnconfigure(1, weight=20)
        self.primary_frame.grid_columnconfigure(2, weight=60)

        self.secondary_m_frame.grid_rowconfigure(0, weight=1)
        self.secondary_m_frame.grid_columnconfigure(0, weight=1)

    def info_frame(self):
        self.info = tk.Frame(self, borderwidth=2, relief="groove")
        self.info.grid(row=1, column=0, sticky="nesw")

        self.info.grid_rowconfigure(0, weight=1)
        self.info.grid_columnconfigure(0, weight=1)

        self.text_info = tk.Label(
            self.info, text="CTE GUI - In development   © IPL")
        self.text_info.grid(row=0, column=0, sticky="sw")


if __name__ == "__main__":
    gui = Main_GUI()
    gui.title("CTE")
    gui.mainloop()
    exit()
