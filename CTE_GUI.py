import tkinter as tk

class Main_GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.label = tk.Label(self, text="In development")
        self.label.pack()


if __name__ == "__main__":
    gui = Main_GUI()
    gui.title("CTE")
    gui.mainloop()
