import tkinter as tk


class Main_GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)

        w = self.winfo_screenwidth()
        h = self.winfo_screenheight()
        self.geometry("%ix%i" % (w/2, h/2))

        main_menu = tk.Menu(self, bg="blue", tearoff=0)
        menu_file = tk.Menu(main_menu, tearoff=0)
        main_menu.add_cascade(label="File", menu=menu_file)
        menu_file.add_command(label="Quit", command=self.quit)
        self.config(menu=main_menu)

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.primary()
        self.info_frame()

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

        list_inputs_features = [
            "Engine model", "Channel global", "Channel dimensions", "Coolant properties"]
        for i in range(len(list_inputs_features)):
            self.inputs.insert(i, list_inputs_features[i])

        self.index_temp = None

        def callback(event):
            selection = event.widget.curselection()
            index = selection[0]
            if selection and self.index_temp != index:
                self.secondary_m()
                self.index_temp = index
                print(index)
            else:
                self.inputs.selection_clear(0, "end")
                self.secondary_m_frame.destroy()
                self.primary_frame.grid_columnconfigure(0, weight=15)
                self.primary_frame.grid_columnconfigure(1, weight=0)
                self.primary_frame.grid_columnconfigure(2, weight=85)
                self.index_temp = None

        self.inputs.bind("<<ListboxSelect>>", callback)

    def tertiary_l2(self):
        self.inputs = tk.Listbox(self.secondary_l_frame)
        self.inputs.grid(row=3, sticky="nesw")

        list_settings_features = ["Plots", "Other"]
        for i in range(len(list_settings_features)):
            self.inputs.insert(i, list_settings_features[i])

    def tertiary_l3(self):
        self.inputs = tk.Listbox(self.secondary_l_frame)
        self.inputs.grid(row=5, sticky="nesw")

    def secondary_r(self):
        self.secondary_r_frame = tk.Frame(
            self.primary_frame, bd=2, relief="flat")
        self.secondary_r_frame.grid(row=0, column=2, sticky="nesw")

        self.secondary_r_frame.grid_rowconfigure(1, weight=9)
        self.secondary_r_frame.grid_rowconfigure(2, weight=1)
        self.secondary_r_frame.grid_columnconfigure(0, weight=1)

        tk.Label(self.secondary_r_frame, text="Results").grid(row=0, column=0)
        self.secondary_r1()
        self.secondary_r2()

    def secondary_r1(self):
        self.secondary_r1_frame = tk.Frame(
            self.secondary_r_frame, bd=2, relief="groove")
        self.secondary_r1_frame.grid(row=1, column=0, sticky="nesw")

        self.secondary_r1_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r1_frame.grid_columnconfigure(0, weight=1)

        tk.Label(self.secondary_r1_frame,
                 text="Results png").grid(row=0, column=0)

    def secondary_r2(self):
        self.secondary_r2_frame = tk.Frame(
            self.secondary_r_frame, bd=2, relief="groove", bg="white")
        self.secondary_r2_frame.grid(row=2, column=0, sticky="nesw")

        self.secondary_r2_frame.grid_rowconfigure(0, weight=1)
        self.secondary_r2_frame.grid_columnconfigure(0, weight=1)

        tk.Label(self.secondary_r2_frame,
                 text="Results text").grid(row=0, column=0)

    def secondary_m(self):
        self.secondary_m_frame = tk.Frame(
            self.primary_frame, bd=2, relief="groove")
        self.secondary_m_frame.grid(row=0, column=1, sticky="nesw")

        self.primary_frame.grid_columnconfigure(0, weight=15)
        self.primary_frame.grid_columnconfigure(1, weight=25)
        self.primary_frame.grid_columnconfigure(2, weight=60)

    def info_frame(self):
        self.info = tk.Frame(self, borderwidth=2, relief="groove")
        self.info.grid(row=1, column=0, sticky="nesw")

        self.info.grid_rowconfigure(0, weight=1)
        self.info.grid_columnconfigure(0, weight=1)

        self.text_info = tk.Label(
            self.info, text="In development")
        self.text_info.grid(row=0, column=0, sticky="sw")


if __name__ == "__main__":
    gui = Main_GUI()
    gui.title("CTE")
    gui.mainloop()
    exit()
