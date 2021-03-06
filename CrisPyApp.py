import matplotlib.pyplot as plt
from tkinter import *
from tkinter import filedialog, simpledialog
import csv
import CrisPy

class App(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()
        self.mut_freq_dict = {}
        # user can set a default target sequence here
        #self.target_sequence = 'CCGGCAAGCTGCCCGTGCCC'
        self.target_sequence = 'GGGCACGGGCAGCTTGCCGG'
        #self.target_range = [13,14,15,16]
        self.target_range = [4,5,6]

    def createWidgets(self):
        # Text box
        self.title = LabelFrame(self, text="CrisPy")
        self.title.pack(fill="both", expand=1)
        self.left = Label(self.title, text=textBoxText)
        self.left["justify"] = "left"
        self.left["width"] = 50
        self.left.pack(fill="both", expand=1)

        # File select button
        self.load = Button(self)
        self.load["text"] = "Load Data"
        self.load["command"] = self.selectFiles
        self.load["width"] = 10
        self.load["height"] = 3
        self.load["border"] = 3
        self.load["font"] = ["calibri", 24, "bold"]
        self.load["anchor"] = "center"
        self.load.pack(side="left", padx=30, pady=30, fill="both")

        # analyze select button
        self.analyze = Button(self)
        self.analyze["text"] = "Analyze \nSequences"
        self.analyze["command"] = self.analyzeData
        self.analyze["width"] = 10
        self.analyze["height"] = 3
        self.analyze["border"] = 3
        self.analyze["font"] = ["calibri", 24, "bold"]
        self.analyze["anchor"] = "center"
        self.analyze.pack(side = "right", padx=30, pady=30, fill="both")

    def selectFiles(self):
        # Delete the selection buttons and update text box
        self.load.destroy()
        self.analyze.destroy()
        self.left["text"] = "Select Ref and Experiment files, and a target sequence.\n"

        self.ref_path = filedialog.askopenfilename(title = "Select reference file", filetypes = [("ab1 files","*.ab1")])
        self.test_path = filedialog.askopenfilename(title = "Select experiment file", filetypes = [("ab1 files","*.ab1")])
        input_sequence = simpledialog.askstring("Target Sequence","Input your nucleotide sequence:")
        if (len(input_sequence)>0):
              self.target_sequence = input_sequence
        for letter in self.target_sequence:
            if letter not in ('A','C','G','T','a','c','g','t'):
                raise "not a valid nucleotide sequence"
        input_range = simpledialog.askstring("Target Range","Start and end index(ex: 13,15)")
        if (len(input_range)>0):
            self.target_range = input_range
        try:
            if not (isinstance(self.target_range, list)): # by default sets this range to the CAMERA experiment's
                self.target_range = list(map(int, self.target_range.split(',')))
                if len(self.target_range) != 2:
                    raise "not a valid range"
                elif self.target_range[0] > self.target_range[1]:
                    raise "not a valid range"
                elif self.target_range[0] not in range(0,len(self.target_sequence)):
                    raise "not a valid range"
                elif self.target_range[1] not in range(0,len(self.target_sequence)):
                    raise "not a valid range"
                else:
                    self.target_range = list(range(self.target_range[0], self.target_range[1]+1))
        except ValueError:
            raise "not a valid range"

        self.next = Button(self, text="Continue",
                                 command=self.doAgain, width=10, height=3,
                                 border=3, font=["calibri",24,"bold"])
        self.next.pack(side="left", padx=30, pady=30, fill="both", expand=1)

    def analyzeData(self):
        self.load.destroy()
        self.analyze.destroy()

        # Calls on CrisPy
        # Initializes SeqDoc object with reference and test file
        self.seqdoc = CrisPy.SeqDoc(self.ref_path, self.test_path)
        align_length, self.diffs = self.seqdoc.get_all_data()

        # Initializes OfftargetFinder
        offtargeter = CrisPy.OfftargetFinder(self.seqdoc.ref_trace, self.target_sequence)
        match_dict = offtargeter.get_targets()

        # Initializes Sequalizer object
        sequalizer = CrisPy.Sequalizer(self.seqdoc.ref_trace, self.seqdoc.test_trace, self.diffs, self.target_sequence, self.target_range)
        self.mut_freq_dict[0] = sequalizer.get_mutation_freq()
        print('Target Sequence:')
        print('Score is ', 1)
        print(self.mut_freq_dict[0], '\n')
        for k in sorted(match_dict):
            self.mut_freq_dict[k] = sequalizer.get_mutation_freq(match_override=match_dict[k])
            print('Off-Target Sequence:')
            print('Score is ' , 1-k)
            print(self.mut_freq_dict[k],'\n')


        self.left["text"] = "Mutation frequency calcultations succesful.\n"
        self.save = Button(self, text="Save Output",
                                 command=self.saveOutput, width=10, height=3,
                                 border=3, font=["calibri",24,"bold"])
        self.save.pack(side="left", padx=30, pady=30, fill="both", expand=1)
        self.display = Button(self, text="Display",
                                    command=self.displayOutput, width=10, height=3,
                                    border=3, font=["calibri",24,"bold"])
        self.display.pack(side="right", padx=30, pady=30, fill="both", expand=1)

    def displayOutput(self):
        # display difference data on a plot
        plt.xticks(self.seqdoc.ref_trace['base_pos'])
        plt.title('Difference Between Traces')
        plt.xlabel('Base Index')
        plt.ylabel('Relative Frequency')
        line1 = plt.plot(self.diffs['A'], label = 'Adenine')
        line2 = plt.plot(self.diffs['G'], label = 'Guanine')
        line3 = plt.plot(self.diffs['C'], label = 'Cytosine')
        line4 = plt.plot(self.diffs['T'], label = 'Thymine')
        plt.legend()
        plt.show()
        pass

    def saveOutput(self):
        self.save.destroy()
        self.display.destroy()

        with open('CrisPy_Output.csv', 'w') as f:  # Just use 'w' mode in 3.x
            w = csv.DictWriter(f, self.mut_freq_dict.keys())
            w.writeheader()
            w.writerow(self.mut_freq_dict)

        self.left["text"] += "Would you like to  analyze another file, or are you done?"

        self.again = Button(self, text="Analyze\nNew Sequences",
                            command=self.doAgain, width=16, height=3,
                            border=3, font=["calibri",24,"bold"])
        self.again.pack(side="left", padx=30, pady=30, fill="both", expand=1)
        
        self.quit = Button(self, text="Quit",
                           command=root.destroy, width=10, height=3,
                           border=3, font=["calibri",24,"bold"])
        self.quit.pack(side="right", padx=30, pady=30, fill="both", expand=1)

    def doAgain(self):
        try:
            self.again.destroy()
        except:
            pass
        try:
            self.quit.destroy()
        except:
            pass
        try:
            self.title.destroy()
        except:
            pass
        try:
            self.left.destroy()
        except:
            pass
        try:
            self.next.destroy()
        except:
            pass

        self.createWidgets()



# The main loop of the program
root = Tk()
origType = None
incorrectAttempts = 0
textBoxText = "Welcome to CrisPy!\nPlease select your sanger files and input target sequence"

# create the application
myapp = App(master=root)
myapp.master.title("Pitt iGEM 2018")
myapp.master.minsize(640, 320)
myapp.mainloop()   
