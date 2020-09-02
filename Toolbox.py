
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 18:32:51 2020

@author: malou
"""

#!/usr/bin/python3


import matplotlib.pyplot as plt
import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from numpy import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
import os
from itertools import permutations
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import threading
import time
windowbg = 'pink'
Buttoncolor = '#e9c4ee'
Activebuttoncolor = '#faceff'

root=Tk()
root.title('Toolbox')
root.geometry('160x290')
root.configure(background = windowbg)


Aminoaciddictionary = {'AAA':'K','AAT':'N','AAC':'N','AAG':'K',\
                       'ATA':'I','ATT':'I','ATC':'I','ATG':'M',\
                       'ACA':'T','ACT':'T','ACC':'T','ACG':'T',\
                       'AGA':'R','AGT':'S','AGC':'S','AGG':'R',\
                       'TAA':'STOP','TAT':'Y','TAC':'Y','TAG':'STOP',\
                       'TTA':'L','TTT':'F','TTC':'F','TTG':'L',\
                       'TCA':'S','TCT':'S','TCC':'S','TCG':'S',\
                       'TGA':'STOP','TGT':'C','TGC':'C','TGG':'W',\
                       'CAA':'Q','CAT':'H','CAC':'H','CAG':'Q',\
                       'CTA':'L','CTT':'L','CTC':'L','CTG':'L',\
                       'CCA':'P','CCT':'P','CCC':'P','CCG':'P',\
                       'CGA':'R','CGT':'R','CGC':'R','CGG':'R',\
                       'GAA':'Z','GAT':'D','GAC':'D','GAG':'Z',\
                       'GTA':'V','GTT':'V','GTC':'V','GTG':'V',\
                       'GCA':'A','GCT':'A','GCC':'A','GCG':'A',\
                       'GGA':'G','GGT':'G','GGC':'G','GGG':'G'}
def destroy():
    root8_1.destroy()

def factorial(x):
    perms = 1
    for i in range(1,x+1):
        perms *= i
    return(perms)

def Permutations():

    Sequence = T2.get(1.0, END).rstrip().upper()
    Seqset = set(Sequence)
    Lettercount = []
    for i in Seqset:
        Lettercount.append(Sequence.count(i))

    Totperms = factorial(len(Sequence))
    fraction = 1
    for j in Lettercount:
        fraction *= factorial(j)
    Uniqueperms = int(Totperms/fraction)
    def oligobox():
        global listbox8_1
        try:
            root8_1.destroy()
        except:
            pass
        root8 = Toplevel(bg = windowbg)
        root8.title('Sequences')
        root8.geometry('150x320')
        L8_1 = Label(root8, text = '{} permutations\n found'.format(Uniqueperms), bg = windowbg)
        L8_1.place(x = 15, y = 10)
        listbox8_1 = Listbox(root8)
        listbox8_1.place(x = 10, y = 60, width = 130, height = 200)
        scrollbar8_1 = Scrollbar(root8, orient = "vertical")
        scrollbar8_1.config(command = listbox8_1.yview, bg = windowbg, activebackground = Activebuttoncolor, elementborderwidth = 2, troughcolor = Buttoncolor)
        scrollbar8_1.place(x = 115, y = 60, width = 25, height = 200)
        listbox8_1.config(yscrollcommand = scrollbar8_1.set)
        perm = [''.join(p) for p in permutations(Sequence)]
        perm = set(perm)
        perm = list(perm)
        perm = sort(perm)
        for i in perm:
            listbox8_1.insert(END, i)

        B8_1 = Button(root8, text = 'Use as primer', bg = Buttoncolor, activebackground = Activebuttoncolor, borderwidth = 2, command = Primer)
        B8_1.place(x = 15, y = 270)

    def warningbox():
        global root8_1
        root8_1 = Toplevel(bg = windowbg)
        root8_1.geometry('300x200')
        root8_1.title('Warning')
        B8_1 = Button(root8_1, text = 'Cancel', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command = destroy )
        B8_1.place(x = 70, y = 150, width = 80, anchor = 'center')
        B8_2 = Button(root8_1, text = 'Continue anyway', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command=lambda:[destroy(),oligobox()])
        B8_2.place(x = 200, y = 150, width = 150, anchor = 'center')
        L8_1 = Label(root8_1, text = 'The choosen sequence will generate \n{} unique permutations out of \n{} possible.'.format(Uniqueperms,Totperms), bg = windowbg)
        L8_1.place(x = 150, y = 40, anchor = 'center')
        L8_2 = Label(root8_1, text = 'This might take some time to generate.', bg = windowbg)
        L8_2.place(x = 150, y = 110, anchor = 'center')
    if len(Sequence) > 10:
        warningbox()

    else:
        oligobox()

def Oligo():
    global T2
    root6 = Toplevel(bg = windowbg)
    root6.title('Oligo')
    root6.geometry('160x130')
    L3 = Label(root6, text = 'Enter oligo:', bg = windowbg)
    L3.place(x = 10, y = 10)
    T2 = Text(root6)
    T2.place(x = 10, y = 40, width = 140, height = 30)
    B3 = Button(root6, text = 'Get permutations', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command = Permutations)
    B3.place(x = 10, y = 75, width = 140, height = 35)
def Baselist():
    root7_1 = Toplevel(bg = windowbg)
    root7_1.title('Baselist')
    root7_1.geometry('150x310')
    listframe = Frame(root7_1, bg = Activebuttoncolor)
    listframe.place(x = 10, y = 10 ,width = 130, height = 290)
    L7_3 = Label(listframe, text = 'R = A or G', bg = Activebuttoncolor)
    L7_3.place(x = 10, y = 10)
    L7_4 = Label(listframe, text = 'M = A or C', bg = Activebuttoncolor)
    L7_4.place(x = 10, y = 35)
    L7_5 = Label(listframe, text = 'W = A or T', bg = Activebuttoncolor)
    L7_5.place(x = 10, y = 60)
    L7_6 = Label(listframe, text = 'Y = C or T', bg = Activebuttoncolor)
    L7_6.place(x = 10, y = 85)
    L7_7 = Label(listframe, text = 'S = C or G', bg = Activebuttoncolor)
    L7_7.place(x = 10, y = 110)
    L7_8 = Label(listframe, text = 'K = G or T', bg = Activebuttoncolor)
    L7_8.place(x = 10, y = 135)
    L7_9 = Label(listframe, text = 'B = C, G or T', bg = Activebuttoncolor)
    L7_9.place(x = 10, y = 160)
    L7_10 = Label(listframe, text = 'H = A, C or T', bg = Activebuttoncolor)
    L7_10.place(x = 10, y = 185)
    L7_11 = Label(listframe, text = 'V = A, C or G', bg = Activebuttoncolor)
    L7_11.place(x = 10, y = 210)
    L7_12 = Label(listframe, text = 'D = A, G or T', bg = Activebuttoncolor)
    L7_12.place(x = 10, y = 235)
    L7_13 = Label(listframe, text = 'N = A, G, C or T', bg = Activebuttoncolor)
    L7_13.place(x = 10, y = 260)

def Primeseq():
    global sequencelist
    Basedict = {'R': ['A','G'], 'M': ['A','C'],'W': ['A','T'],\
     'Y':['C', 'T'],'S': ['C', 'G'], 'K': ['G','T'],'B':['C','G','T'],\
     'H': ['A','C','T'],'V':['A','C','G'], 'D':['A','G','T'],'N': ['A','C','G','T']}
    primerseq = T7_1.get(1.0,END).upper()
    sequencelist = []
    for base in primerseq:
        if base in Basedict:
            baselist = Basedict[base]
            baselen = len(baselist)
            if len(sequencelist)!=0:
                for i in range(len(sequencelist)):
                    addseq = sequencelist[i]
                    sequencelist[i] = addseq + baselist[0]
                    for j in range(1,baselen):
                        sequencelist.append(addseq+baselist[j])
            else:
                for nucleotide in baselist:
                    sequencelist.append('{}'.format(nucleotide))

        else:
            if len(sequencelist)!=0:

                for i in range(len(sequencelist)):
                    sequencelist[i] += base
            else:
                sequencelist.append('{}'.format(base))
    def Retrieve_as_primer():
        selected_primer = listbox7_1.get(listbox7_1.curselection()[0])
        T7_1.delete(1.0,END)
        T7_1.insert(1.0, selected_primer)
        root7_1.destroy()

    root7_1 = Toplevel(bg = windowbg)
    root7_1.title('Primerlist')
    root7_1.geometry('150x320')
    L7_11 = Label(root7_1, text = '{} possible \n sequences '.format(len(sequencelist)), bg = windowbg)
    L7_11.place(x = 75, y = 20, anchor = 'center')
    B7_11 = Button(root7_1, text = 'Use as primer', bg = Buttoncolor, activebackground = Activebuttoncolor, borderwidth = 2, command = Retrieve_as_primer)
    B7_11.place(x = 75, y = 290, anchor = 'center')
    listbox7_1 = Listbox(root7_1)
    listbox7_1.place(x = 10, y = 60, width = 130, height = 200)
    scrollbar7_1 = Scrollbar(root7_1, orient = "vertical")
    scrollbar7_1.config(command = listbox7_1.yview, bg = windowbg, activebackground = Activebuttoncolor, elementborderwidth = 2, troughcolor = Buttoncolor)
    scrollbar7_1.place(x = 115, y = 60, width = 25, height = 200)
    listbox7_1.config(yscrollcommand = scrollbar7_1.set)
    for item in sequencelist:
        item = item.rstrip()
        listbox7_1.insert(END, item)

def Primersearch():

    try:
        x = len(sequencelist)
    except:
        Basedict = {'R': ['A','G'], 'M': ['A','C'],'W': ['A','T'],\
         'Y':['C', 'T'],'S': ['C', 'G'], 'K': ['G','T'],'B':['C','G','T'],\
         'H': ['A','C','T'],'V':['A','C','G'], 'D':['A','G','T'],'N': ['A','C','G','T']}
        primerseq = T7_1.get(1.0,END).upper()
        sequencelist = []
        for base in primerseq:
            if base in Basedict:
                baselist = Basedict[base]
                baselen = len(baselist)
                if len(sequencelist)!=0:
                    for i in range(len(sequencelist)):
                        addseq = sequencelist[i]
                        sequencelist[i] = addseq + baselist[0]
                        for j in range(1,baselen):
                            sequencelist.append(addseq+baselist[j])
                else:
                    for nucleotide in baselist:
                        sequencelist.append('{}'.format(nucleotide))

            else:
                if len(sequencelist)!=0:

                    for i in range(len(sequencelist)):
                        sequencelist[i] += base
                else:
                    sequencelist.append('{}'.format(base))

    Querysequence = T1.get(1.0,END).rstrip().upper()
    end = -1
    Primerfile = open('/home/malou/Pythonprojekt/Primer.txt', 'w')
    for sequence in sequencelist:
        start = 0
        positionlist = []
        for i in range(len(Querysequence)):
            position = Querysequence.find(sequence.rstrip(), start, end)
            if position == -1:
                break
            start = position +1
            positionlist.append(position)
        if len(positionlist) == 0:
            Primerfile.write('No matches found for sequence {}\n'.format(sequence))
        else:
            Primerfile.write('Sequence {} found at position/s:\n'.format(sequence.rstrip()))
            for position in positionlist:
                Primerfile.write('{}\n'.format(position))
        Primerfile.write('\n')
    Primerfile.close()
    want_run = ('"gedit /home/malou/Pythonprojekt/Primer.txt"')
    os.system('bash -c '  + want_run)



def Primer():
    global T7_1, T1, L4
    try:
        selected_oligo = listbox8_1.get(listbox8_1.curselection()[0])
        T7_1.delete(1.0,END)
        T7_1.insert(1.0, selected_oligo)
    except:
        root7 = Toplevel(bg = windowbg)
        root7.title('Primer')
        root7.geometry('300x380')
        L7_1 = Label(root7, text = 'Enter primer:', bg = windowbg)
        L7_1.place(x = 10, y = 10)
        T7_1 = Text(root7)
        T7_1.place(x = 10, y = 40, width = 280, height = 30)
        L7_2 = Label(root7, text = 'Choose sequence from file', bg = windowbg)
        L7_2.place(x = 10, y = 140)
        L4 = Label(root7, text = '', bg = windowbg)
        L4.place(x = 130, y = 300)
        F_1 = Frame(root7, bg = 'Grey')
        F_1.place(x = 0, y = 120, width = 410, height = 1)
        B7_2 = Button(root7, text = 'Open file', borderwidth = 2, bg = Buttoncolor, command = Openfile, activebackground = Activebuttoncolor)
        B7_2.place(x = 200, y = 135)
        T1 = Text(root7)
        T1.place(x = 10, y = 195, width = 280, height = 100)
        L7_3 = Label(root7, text = 'Or enter a sequence:', bg = windowbg)
        L7_3.place(x = 10, y = 170)
        B7_1 = Button(root7, text = 'Retrieve sequences', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command = Primeseq)
        B7_1.place(x = 90, y = 95, width = 160, anchor = 'center')
        B8_1 = Button(root7, text = 'List of bases', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command = Baselist)
        B8_1.place(x = 235, y = 95, width = 110, anchor = 'center')
        B7_2 = Button(root7, text = 'Search for primer', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command = Primersearch)
        B7_2.place(x = 150, y = 345, anchor = 'center')

        try:
            selected_oligo = listbox8_1.get(listbox8_1.curselection()[0])
            T7_1.delete(1.0,END)
            T7_1.insert(1.0, selected_oligo)
        except:
            pass


def Openfile():
    style = ttk.Style()
    style.configure('TFrame', background = windowbg)
    style.configure('TLabel', background = windowbg)
    style.map('TLabel', background = [('active',windowbg), ('disabled', windowbg)])
    style.configure('TMenubutton', background = windowbg)
    style.map('TMenubutton', background = [('active',windowbg), ('disabled', windowbg)])
    style.configure('TButton',  background = Buttoncolor, activebackground = Activebuttoncolor)
    style.configure('Horizontal.TScrollbar', troughcolor=Activebuttoncolor, bordercolor=Buttoncolor, background='#ba9cbe', lightcolor='white', darkcolor=Buttoncolor)
    file = filedialog.askopenfilename(initialdir = "/home/malou/Pythonprojekt")
    root3 = Toplevel(bg = windowbg)
    root3.title('ID-list')
    root3.geometry('260x420')

    def ChooseSeq():
        listitem = listbox1.get(listbox1.curselection()[0])
        listitem = listitem.split(' ')[1]
        T1.insert(END, IDdict[listitem])
        try:
            L4.configure(text = 'Sequence lenght: {}'.format((len(IDdict[listitem]))))
        except:
            pass
        root3.destroy()

    label3_1 = Label(root3, text = "Several sequence-ID's found,\n please choose one", bg = windowbg)
    label3_1.place(x = 20, y =10)
    listbox1 = Listbox(root3)
    listbox1.place(x = 20, y = 50, height = 300, width = 200)
    scrollbar1 = Scrollbar(root3, orient="vertical")
    scrollbar1.config(command=listbox1.yview, bg = windowbg, activebackground = Activebuttoncolor, elementborderwidth = 2, troughcolor = Buttoncolor)
    scrollbar1.place(x=220,y=50,height=300, width = 25)
    listbox1.config(yscrollcommand=scrollbar1.set)
    B31 = Button(root3, text = 'Choose sequence', command = ChooseSeq, bg = Buttoncolor, activebackground = Activebuttoncolor, borderwidth = 3)
    B31.place(x = 50, y = 360, width = 170, height = 50)
    with open('{}'.format(file),'r') as F1:
        T1.delete(1.0,END)
        IDdict = {}
        for line in F1:
            if line.startswith('>'):
                line = line.split('\t')
                IDline = line[0]
                IDdict[IDline] = ''
                listbox1.insert(END, 'IDname:' + ' ' +'{}'.format(IDline))
            else:
                if bool(IDdict) == False:
                    T1.insert(END, line.rstrip())
                    try:
                        L4.configure(text = 'Sequence lenght: {}'.format(len(line.rstrip())))
                    except:
                        continue
                else:
                    IDdict[IDline] += line.rstrip()
        if bool(IDdict) == False:
            root3.destroy()

def GCcontent():
    global T1, L4
#============================= Window set-up ===================================
    root2 = Toplevel(bg = windowbg)
    root2.title('GC-content')
    root2.geometry('320x310')
#===============================================================================
    def GCplot():
        Sequence = T1.get(1.0,END).upper()
        Windowsize = int(E1.get())
        Steps = int(E2.get())
        Iterationlength = len(Sequence)-Windowsize
        Procentlist = []
        for i in range(0, Iterationlength, Steps):
            GCcontent = (Sequence[i:(i+Windowsize)].count('G') + Sequence[i:(i+Windowsize)].count('C'))/Windowsize
            Procentlist.append(GCcontent)
        Positionlist = linspace(Windowsize/2, Iterationlength + (Windowsize/2), len(Procentlist))
        Procentlist = [n*100 for n in Procentlist]
        root22 = Toplevel(bg = windowbg)
        root22.title('Plot')
        root22.geometry('510x410')
        figure = Figure(figsize = (5,4), dpi = 100)
        plot_1 = figure.add_subplot(1,1,1)
        plot_1.set_ylabel('Procent')
        plot_1.set_xlabel('Position')
        plot_1.set_title('Windowsize {} and stepsize {}'.format(Windowsize,Steps))
        plot_1.plot(Positionlist,Procentlist)
        figure.patch.set_facecolor(windowbg)
        canvas = FigureCanvasTkAgg(figure, root22)
        canvas.get_tk_widget().place(x = 2, y = 2)
#================================== Widgets ====================================
    L2 = Label(root2, text = 'Choose sequence from file', bg = windowbg)
    L2.place(x = 10, y = 15)
    B2 = Button(root2, text = 'Open file', borderwidth = 2, bg = Buttoncolor, command = Openfile, activebackground = Activebuttoncolor)
    B2.place(x = 200, y = 10)
    L1 = Label(root2, text = 'Or enter a sequence to calculate GC-content:', bg = windowbg)
    L1.place(x = 10, y = 50)
    T1 = Text(root2)
    T1.place(x = 10, y = 80, width = 300, height = 100)
    L4 = Label(root2, text = '', bg = windowbg)
    L4.place(x = 150, y = 190)
    E1 = Entry(root2)
    E1.place(x = 110, y = 230, width = 30)
    E2 = Entry(root2)
    E2.place(x = 110, y = 260, width = 30)
    L2 = Label(root2, text = 'Windowsize', bg = windowbg)
    L2.place(x = 20, y = 232)
    L3 = Label(root2, text = 'Steps', bg = windowbg)
    L3.place(x = 20, y = 262)
    B1 = Button(root2, text = 'Calculate', borderwidth = 5, bg = Buttoncolor, activebackground = Activebuttoncolor, command = GCplot)
    B1.place(x = 160, y = 225, width = 140, height = 60)
#===============================================================================

def start_submit_thread(event):
    global submit_thread, root9_1, progressbar
    root9_1 = Toplevel(bg = windowbg)
    root9_1.title('Loading')
    # Gets the requested values of the height and widht.
    windowWidth = root9_1.winfo_reqwidth()
    windowHeight = root9_1.winfo_reqheight()

    # Gets both half the screen width/height and window width/height
    positionRight = int(root9_1.winfo_screenwidth()/2 - windowWidth/2)
    positionDown = int(root9_1.winfo_screenheight()/2 - windowHeight/2)

    # Positions the window in the center of the page.
    root9_1.geometry("300x100+{}+{}".format(positionRight, positionDown))

    # root9_1.geometry('300x100')
    L9_1 = Label(root9_1,text = 'Loading, please wait', bg = windowbg)
    L9_1.place(x = 150, y = 30, anchor = 'center')
    s = ttk.Style()
    s.theme_use('clam')
    s.configure("bar.Horizontal.TProgressbar", troughcolor=Activebuttoncolor, bordercolor=Buttoncolor, background='#ba9cbe', lightcolor='white', darkcolor=Buttoncolor)
    progressbar = ttk.Progressbar(root9_1, style = "bar.Horizontal.TProgressbar", mode='indeterminate')
    progressbar.place(x = 10, y = 60, width = 280)


    submit_thread = threading.Thread(target=Blastbutton)
    submit_thread.daemon = True
    progressbar.start()
    submit_thread.start()
    root9.after(20, check_submit_thread)

def root9destroy():
    root9_2.destroy()
    want_run = ('"gedit {}"'.format(filename1))
    os.system('bash -c '  + want_run)

def check_submit_thread():
    if submit_thread.is_alive():
        root9.after(20, check_submit_thread)
    else:
        return
def Blastbutton():
    global root9_2, filename1
    try:
        E_VALUE_THRESH = float(E1.get())
    except:
        E_VALUE_THRESH = 150
    sequence = T1.get(1.0, END).upper()
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    with open('my_blast.xml', 'w') as out_handle:
        out_handle.write(result_handle.read())
    result_handle = open('my_blast.xml')
    blast_record = NCBIXML.read(result_handle)
    progressbar.stop()

    style = ttk.Style()
    style.configure('TFrame', background = windowbg)
    style.configure('TLabel', background = windowbg)
    style.map('TLabel', background = [('active',windowbg), ('disabled', windowbg)])
    style.configure('TMenubutton', background = windowbg)
    style.map('TMenubutton', background = [('active',windowbg), ('disabled', windowbg)])
    style.configure('TButton',  background = Buttoncolor, activebackground = Activebuttoncolor)
    style.configure('Horizontal.TScrollbar', troughcolor=Activebuttoncolor, bordercolor=Buttoncolor, background='#ba9cbe', lightcolor='white', darkcolor=Buttoncolor)


    filename = open('/home/malou/Pythonprojekt/Blastfile.txt','w')
    filename1 = '/home/malou/Pythonprojekt/Blastfile.txt'
    for alignment in blast_record.alignments:
         for hsp in alignment.hsps:
             if hsp.expect < E_VALUE_THRESH:
                 filename.write("****Alignment****\n")
                 filename.write('{}\n'.format(alignment.title))
                 filename.write('{}\n'.format(alignment.length))
                 filename.write('{}\n'.format(hsp.expect))
                 filename.write('{}\n'.format(hsp.query))
                 filename.write('{}\n'.format(hsp.match))
                 filename.write('{}\n'.format(hsp.sbjct))
    result_handle.close()
    filename.close()
    root9_1.destroy()
    root9_2 = Toplevel(bg = windowbg)
    root9_2.title('Finished')
    # Gets the requested values of the height and widht.
    windowWidth = root9_2.winfo_reqwidth()
    windowHeight = root9_2.winfo_reqheight()

    # Gets both half the screen width/height and window width/height
    positionRight = int(root9_2.winfo_screenwidth()/2 - windowWidth/2)
    positionDown = int(root9_2.winfo_screenheight()/2 - windowHeight/2)

    # Positions the window in the center of the page.
    root9_2.geometry("300x100+{}+{}".format(positionRight, positionDown))
    L9_1 = Label(root9_2,text = 'Search successful', bg = windowbg)
    L9_1.place(x = 150, y = 30, anchor = 'center')
    B9_1 = Button(root9_2, text = 'Ok', borderwidth = 2, bg = Buttoncolor, activebackground = Activebuttoncolor, command = root9destroy)
    B9_1.place(x = 150, y = 70, anchor = 'center')


def Blast():
    global T1, E1, root9
    root9 = Toplevel(bg = windowbg)
    root9.title('Blast')
    root9.geometry('320x300')
    L2 = Label(root9, text = 'Choose sequence from file', bg = windowbg)
    L2.place(x = 10, y = 15)
    B2 = Button(root9, text = 'Open file', borderwidth = 2, bg = Buttoncolor, command = Openfile, activebackground = Activebuttoncolor)
    B2.place(x = 200, y = 10)
    L1 = Label(root9, text = 'Or enter a sequence to do a BLAST search:', bg = windowbg)
    L1.place(x = 10, y = 50)
    T1 = Text(root9)
    T1.place(x = 10, y = 80, width = 300, height = 80)
    E1 = Entry(root9)
    E1.place(x = 170, y = 180, width = 50)
    # E2 = Entry(root9)
    # E2.place(x = 110, y = 230, width = 30)
    L2 = Label(root9, text = 'set E-value treshhold', bg = windowbg)
    L2.place(x = 20, y = 182)
    # L3 = Label(root9, text = 'Steps', bg = windowbg)
    # L3.place(x = 20, y = 232)

    B1 = Button(root9, text = 'Do a blastn search', borderwidth = 5, bg = Buttoncolor, activebackground = Activebuttoncolor, command=lambda:start_submit_thread(None))
    B1.place(x = 90, y = 220, width = 150, height = 60)




def Complement(x):                      # Function that takes a string as input and translates
    string1 = x.lower().rstrip()                # Makes all chars lower case for the replace method to work
    string2 = re.sub('a','T',string1)
    string2 = re.sub('t','A',string2)
    string2 = re.sub('c','G',string2)
    string2 = re.sub('g','C',string2)
    return(string2[::-1])               # Inverts the string.


def Translate():
    global T1
    root5 = Toplevel(bg = windowbg)
    root5.title('Translate')
    root5.geometry('300x270')

    def TranslateSeq():
        style = ttk.Style()
        style.configure('TFrame', background = windowbg)
        style.configure('TLabel', background = windowbg)
        style.map('TLabel', background = [('active',windowbg), ('disabled', windowbg)])
        style.configure('TMenubutton', background = windowbg)
        style.map('TMenubutton', background = [('active',windowbg), ('disabled', windowbg)])
        style.configure('TButton',  background = Buttoncolor, activebackground = Activebuttoncolor)
        style.configure('Horizontal.TScrollbar', troughcolor=Activebuttoncolor, bordercolor=Buttoncolor, background='#ba9cbe', lightcolor='white', darkcolor=Buttoncolor)

        saveas = filedialog.asksaveasfilename(initialdir = '/home/malou/Pythonprojekt')
        Seq = T1.get(1.0,END).upper().rstrip()
        Seqcomplement = Complement(Seq)
        Seqreading = {0:'',1:'',2:''}
        complreading = {0:'',1:'',2:''}
        for i in range(3):
            for j in range(int(len(Seq)/3)):
                Seqcodon = Seq[i+j*3:(j+1)*3+i]
                Compcodon = Seqcomplement[i+j*3:(j+1)*3+i]
                if len(Seqcodon) == 3:
                    if Aminoaciddictionary[Seqcodon] == 'STOP':
                        Seqreading[i] += '-'
                    else:
                        Seqreading[i] += Aminoaciddictionary[Seqcodon]
                    if Aminoaciddictionary[Compcodon] == 'STOP':
                        complreading[i] += '-'
                    else:
                        complreading[i] += Aminoaciddictionary[Compcodon]
        with open(saveas, 'w') as O1:
            O1.write("5'3 Frame 1:\n{} \n \n".format(Seqreading[0]))
            O1.write("5'3 Frame 2:\n{} \n \n".format(Seqreading[1]))
            O1.write("5'3 Frame 3:\n{} \n \n".format(Seqreading[2]))
            O1.write("3'5 Frame 1:\n{} \n \n".format(complreading[0]))
            O1.write("3'5 Frame 2:\n{} \n \n".format(complreading[1]))
            O1.write("3'5 Frame 3:\n{} \n \n".format(complreading[2]))

        want_run = ('"gedit {}"'.format(saveas))
        os.system('bash -c '  + want_run)

    L2 = Label(root5, text = 'Choose sequence from file', bg = windowbg)
    L2.place(x = 10, y = 15)
    B2 = Button(root5, text = 'Open file', borderwidth = 2, bg = Buttoncolor, command = Openfile, activebackground = Activebuttoncolor)
    B2.place(x = 200, y = 10)
    L1 = Label(root5, text = 'Or enter a sequence to translate', bg = windowbg)
    L1.place(x = 10, y = 50)
    T1 = Text(root5)
    T1.place(x = 10, y = 80, width = 280, height = 100)
    B1 = Button(root5, text = 'Translate', borderwidth = 5, bg = Buttoncolor, activebackground = Activebuttoncolor, command = TranslateSeq)
    B1.place(x = 85, y = 195, width = 140, height = 60)


#========================= Buttonsection root window ===========================
B = {}
ButtonNames = ['Oligo', 'Primer', 'GC-content', 'Blast', 'Translate']
CommandList = [Oligo, Primer, GCcontent, Blast, Translate]
for i in range(len(ButtonNames)):
    B[i] = Button(root, text = ButtonNames[i], borderwidth = 4, command = CommandList[i], background  = Buttoncolor, activebackground = Activebuttoncolor)
    B[i].place(x = 30, y = 10 + 55*i, height = 50, width = 100)
#===============================================================================



root.mainloop()
